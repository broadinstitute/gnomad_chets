import hail as hl
import logging
import timeit
import argparse

from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_LEAST_CONSEQUENCE,
    DEFAULT_MAX_FREQ,
    DEFAULT_TMP_DIR,
    SITES_FIELD_CLINVAR,
    SITES_FIELD_PANGOLIN,
    SITES_FIELD_SPLICEAI,
    get_phasing_resources,
)
from gnomad_chets.v4.utils import filter_for_testing

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("phase_gnomad")
logger.setLevel(logging.INFO)


_EM_PSEUDO = [0, 0, 0, 0, 1, 0, 0, 0, 0]
"""One pseudo-count in the double-het (AaBb) cell for the ``em_plus_one`` variant."""


def get_em_expr(gt_counts):
    gt_counts_int32 = gt_counts.map(lambda x: hl.int32(x))
    hap_counts = hl.experimental.haplotype_freq_em(gt_counts_int32)
    return hl.bind(
        lambda x: hl.struct(
            hap_counts=x,
            # p_chet collapses to 0/0 = NaN when a stratum has no double-carrier
            # pairs (common for small per-pop strata / empty groups). Surface it
            # as missing so downstream is-defined / threshold logic treats it as
            # "not phased" rather than propagating NaN.
            p_chet=hl.bind(
                lambda p: hl.or_missing(~hl.is_nan(p), p),
                (x[1] * x[2]) / (x[0] * x[3] + x[1] * x[2]),
            ),
        ),
        hap_counts
    )


def _em_by_pop_expr(gt_counts_by_pop):
    """Map :func:`get_em_expr` (raw + adj) over a per-pop gt-counts dict.

    ``gt_counts_by_pop`` is ``dict<pop, struct{raw, adj}>`` (from the count
    steps' ``--stratify-by-pop`` output). Returns a parallel
    ``dict<pop, struct{raw, adj}>`` of EM results, mirroring v2's
    ``phasing.get_phased_gnomad_ht`` ``map_values`` approach.
    """
    return (
        gt_counts_by_pop.map_values(
            lambda c: hl.struct(raw=get_em_expr(c.raw), adj=get_em_expr(c.adj))
        ),
        gt_counts_by_pop.map_values(
            lambda c: hl.struct(
                raw=get_em_expr(c.raw + _EM_PSEUDO),
                adj=get_em_expr(c.adj + _EM_PSEUDO),
            )
        ),
    )


def get_phased_gnomad_ht(
        ht: hl.Table
) -> hl.Table:

    phased = dict(
        em=hl.struct(
            raw=get_em_expr(ht.gt_counts_raw),
            adj=get_em_expr(ht.gt_counts_adj),
        ),
        em_plus_one=hl.struct(
            raw=get_em_expr(ht.gt_counts_raw + _EM_PSEUDO),
            adj=get_em_expr(ht.gt_counts_adj + _EM_PSEUDO),
        )
    )
    # When the counts carry a per-population breakdown (--stratify-by-pop),
    # run EM per pop too (keyed by the same pop labels, incl. "all").
    if "gt_counts_by_pop" in ht.row:
        em_by_pop, em_plus_one_by_pop = _em_by_pop_expr(ht.gt_counts_by_pop)
        phased["em_by_pop"] = em_by_pop
        phased["em_plus_one_by_pop"] = em_plus_one_by_pop
    return phased


def _get_variant_ann_expr(
    locus,
    alleles,
    sites_ht: hl.Table,
    variant_filter_ht: hl.Table,
    include_full_vep: bool = False,
    include_full_freq: bool = False,
    include_transcripts: bool = False,
):
    """Build the per-variant annotation struct for one endpoint of a pair.

    Always attaches a compact ``canonical_transcript`` summary
    (gene_symbol, transcript_id, consequence_terms, lof) — picked from the
    VEP canonical flag if available, else the first transcript. Set
    ``include_transcripts=True`` to also attach the full per-transcript
    array. Missing scores / ClinVar rows come through as ``NA`` naturally
    via the joins.
    """
    sites_row = sites_ht[locus, alleles]
    filter_row = variant_filter_ht[locus, alleles]
    tcs = sites_row.vep.transcript_consequences
    # Representative transcript for the compact ``canonical_transcript``
    # summary: prefer protein-coding Ensembl (matches the pipeline's
    # variant-filter definition), then canonical=1, then first. This
    # avoids picking a non-coding ``LNC`` or pseudogene transcript
    # (which lack a ``gene_symbol``) as the representative — a bug that
    # blanked out MUC16 / FBN3 in a6/c10/e3 symbol columns.
    pc_tcs = tcs.filter(
        lambda tc: (tc.biotype == "protein_coding") & (tc.source == "Ensembl")
    )
    canonical_tc = pc_tcs.filter(lambda tc: tc.canonical == 1)
    canonical_or_first = (
        hl.case()
        .when(hl.len(canonical_tc) > 0, canonical_tc[0])
        .when(hl.len(pc_tcs) > 0, pc_tcs[0])
        .when(hl.len(tcs) > 0, tcs[0])
        .or_missing()
    )

    # Coalesce source/gene_id to empty containers so downstream group_by
    # keys (b3, c7 priority-source; c11 v2_match mask) don't drop pairs
    # whose variant is present in sites_ht but absent from variant_filter_ht
    # (a PASS variant that's AC0 in release makes it into the sites HT but
    # never enters the release-scoped variant filter).
    ann = {
        "source": hl.or_else(filter_row.source, hl.empty_set(hl.tstr)),
        "gene_id": hl.or_else(filter_row.gene_id, hl.empty_array(hl.tstr)),
        "an_pct": filter_row.an_pct,
        "ac": sites_row.ac,
        "af": sites_row.af,
        "an": sites_row.an,
        "most_severe_consequence": sites_row.vep.most_severe_consequence,
        "gene_symbols": hl.array(hl.set(tcs.map(lambda tc: tc.gene_symbol))),
        "canonical_transcript": hl.or_missing(
            hl.is_defined(canonical_or_first),
            hl.struct(
                gene_id=canonical_or_first.gene_id,
                gene_symbol=canonical_or_first.gene_symbol,
                transcript_id=canonical_or_first.transcript_id,
                consequence_terms=canonical_or_first.consequence_terms,
                lof=canonical_or_first.lof,
            ),
        ),
        SITES_FIELD_SPLICEAI: sites_row[SITES_FIELD_SPLICEAI],
        SITES_FIELD_PANGOLIN: sites_row[SITES_FIELD_PANGOLIN],
        SITES_FIELD_CLINVAR: sites_row[SITES_FIELD_CLINVAR],
    }
    if include_transcripts:
        # Full per-transcript compact info (~15 transcripts per variant
        # mean; kept in VEP emission order). Bloats the HT ~10x, gated
        # behind --include-transcripts for consumers that need it.
        ann["transcripts"] = tcs.map(lambda tc: hl.struct(
            gene_id=tc.gene_id,
            gene_symbol=tc.gene_symbol,
            transcript_id=tc.transcript_id,
            canonical=tc.canonical == 1,
            mane_select=hl.is_defined(tc.mane_select) & (tc.mane_select != ""),
            consequence_terms=tc.consequence_terms,
            lof=tc.lof,
        ))
    if include_full_vep:
        ann["vep"] = sites_row.vep
    if include_full_freq:
        # freq is on the raw freq HT, not the assembled sites HT — sites
        # only projects freq[0] scalars. Full freq requires a separate join
        # (see main()); here we only carry through if the user pre-joined it.
        ann["freq"] = sites_row.freq
    return hl.struct(**ann)


def annotate_phased_ht_with_sites(
    ht: hl.Table,
    sites_ht: hl.Table,
    variant_filter_ht: hl.Table,
    include_full_vep: bool = False,
    include_full_freq: bool = False,
    include_transcripts: bool = False,
) -> hl.Table:
    """Attach per-variant annotation structs (``v1_ann`` / ``v2_ann``) to a phased pair HT.

    Also adds a pair-level ``shared_gene_ids`` field with the intersection
    of ``v1_ann.gene_id`` and ``v2_ann.gene_id`` — the "which gene do both
    endpoints hit" answer downstream consumers almost always want.

    Pure transformation: takes and returns HTs, no I/O.
    """
    v1_ann = _get_variant_ann_expr(
        ht.locus1, ht.alleles1, sites_ht, variant_filter_ht,
        include_full_vep=include_full_vep,
        include_full_freq=include_full_freq,
        include_transcripts=include_transcripts,
    )
    v2_ann = _get_variant_ann_expr(
        ht.locus2, ht.alleles2, sites_ht, variant_filter_ht,
        include_full_vep=include_full_vep,
        include_full_freq=include_full_freq,
        include_transcripts=include_transcripts,
    )
    ht = ht.annotate(v1_ann=v1_ann, v2_ann=v2_ann)
    ht = ht.annotate(
        shared_gene_ids=hl.array(
            hl.set(ht.v1_ann.gene_id).intersection(hl.set(ht.v2_ann.gene_id))
        ),
    )
    return ht

def main(args):
    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix or ""
    data_type = args.data_type
    test = args.test
    
    hl.init(
        log="/create_vp_matrix.log",
        tmp_dir=tmp_dir,
    )
    
    logger.info(
        f"""
        Running script with the following parameters:

            Data type: {data_type}
            Test: {test}
            Output postfix: {output_postfix}
            Overwrite: {overwrite}
            Tmp dir: {tmp_dir}
        """
    )
    
    resources = get_phasing_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir,
        output_postfix=output_postfix,
        overwrite=overwrite,
    )
    
    if args.phase:
        logger.info("Phasing variant pairs...")
        res=resources.phase

        #phase variant pairs
        ht=hl.read_table(args.file_to_phase)
        print(ht.describe())
        # Repartition before EM only when --em-partitions > 0. With the
        # Tier 3-corrected counts, the EM's earlier stalling behavior
        # (pathological negative gt_counts driving haplotype_freq_em into
        # an infinite oscillation) shouldn't reoccur, so the input's
        # native partitioning is often fine. Passing --em-partitions 0
        # skips the repartition shuffle entirely.
        if args.em_partitions and args.em_partitions > 0:
            logger.info(
                "Repartitioning to %d partitions before EM.", args.em_partitions,
            )
            ht = ht.repartition(args.em_partitions).checkpoint(
                hl.utils.new_temp_file("phase_gnomad_repartitioned", "ht")
            )
        else:
            logger.info(
                "Skipping repartition (--em-partitions=%s); using input's "
                "native partitioning (%d partitions).",
                args.em_partitions, ht.n_partitions(),
            )
        phased_dict=get_phased_gnomad_ht(ht)
        logger.info("Phasing complete. Now annotating phased data...")

        ht = ht.annotate(**dict(phased_dict)).checkpoint(
            hl.utils.new_temp_file("get_phased_gnomad", "ht")
        )
        logger.info("Annotating complete. Now writing phased data...")
        #write phased data
        ht=ht.write(res.phase.path,overwrite=overwrite)

    if args.annotate_phased_ht:
        logger.info("Annotating phased HT with per-variant sites annotations...")
        res = resources.annotate_phased_ht

        phased_ht = hl.read_table(args.phased_ht_path or res.phased.path)
        sites_ht = hl.read_table(args.sites_ht_path or res.sites_ht.path)
        variant_filter_ht = hl.read_table(
            args.variant_filter_ht_path or res.variant_filter_ht.path
        )

        if args.include_full_freq:
            # sites_ht only projects freq[0] scalars; pull the full freq
            # array from the raw freq HT. Semi-join to sites keys so test
            # / postfix filtering is inherited from the sites HT itself.
            from gnomad_qc.v4.resources.annotations import get_freq
            freq_ht = get_freq(data_type=data_type).ht()
            sites_ht = sites_ht.annotate(
                freq=freq_ht[sites_ht.locus, sites_ht.alleles].freq
            )

        annotated = annotate_phased_ht_with_sites(
            phased_ht,
            sites_ht,
            variant_filter_ht,
            include_full_vep=args.include_full_vep,
            include_full_freq=args.include_full_freq,
            include_transcripts=args.include_transcripts,
        )
        annotated.write(res.annotated_phase.path, overwrite=overwrite)
        logger.info("Wrote annotated phased HT to %s", res.annotated_phase.path)

    stop = timeit.default_timer()
    logger.info(f"Time taken to run the script is {stop - start} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tmp-dir",
        default=DEFAULT_TMP_DIR,
        help="Temporary directory for intermediate files.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to PCNT gene (chr21:46324141-46445769) for testing purposes.",
    )
    parser.add_argument(
        "--output-postfix",
        help=(
            'Postfix to append to output file names (e.g., "pcnt_test" for files like '
            "exomes.vp_list.pcnt_test.ht)."
        ),
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Whether to overwrite existing files."
    )
    parser.add_argument(
        "--data-type",
        default=DEFAULT_DATA_TYPE,
        choices=DATA_TYPE_CHOICES,
        help=(
            f'Data type to use. Must be one of {", ".join(DATA_TYPE_CHOICES)}. Default '
            f"is {DEFAULT_DATA_TYPE}.",
        ),
    )
    parser.add_argument(
        "--phase",
        action="store_true",
        help="Whether to phase variant pairs.",
    )
    
    parser.add_argument(
        "--file-to-phase",
        help="input file for phasing",
    )
    parser.add_argument(
        "--em-partitions",
        type=int,
        default=10000,
        help=(
            "Number of partitions to repartition the gt_counts HT to before "
            "running haplotype_freq_em. EM is row-local with uncapped "
            "iterations; finer partitions bound per-task wall time on "
            "slow-converging (near-degenerate) rows."
        ),
    )
    parser.add_argument(
        "--annotate-phased-ht",
        action="store_true",
        help=(
            "Join per-variant annotations (source tags, gene_id, an_pct, "
            "AC/AN/AF, most_severe_consequence, gene_symbols, SpliceAI, "
            "Pangolin, ClinVar) onto the phased pair HT. Writes to "
            "phased.annotated.{postfix}.ht."
        ),
    )
    parser.add_argument(
        "--include-full-vep",
        action="store_true",
        help=(
            "When --annotate-phased-ht is set, also carry the full VEP "
            "struct (transcript_consequences) per endpoint. Off by default "
            "— bloats the HT ~10x for what most consumers don't need."
        ),
    )
    parser.add_argument(
        "--include-full-freq",
        action="store_true",
        help=(
            "When --annotate-phased-ht is set, also carry the full freq "
            "array (per-population breakdown) per endpoint. Off by default."
        ),
    )
    parser.add_argument(
        "--include-transcripts",
        action="store_true",
        help=(
            "When --annotate-phased-ht is set, also carry the full "
            "per-transcript array (gene_symbol, transcript_id, "
            "consequence_terms, lof, canonical/mane flags) per endpoint. "
            "Off by default — bloats the HT ~10x. canonical_transcript "
            "is always included regardless of this flag."
        ),
    )
    parser.add_argument(
        "--phased-ht-path",
        help=(
            "Override path to the phased HT (input to --annotate-phased-ht). "
            "Defaults to the postfix-scoped path from get_phasing_resources."
        ),
    )
    parser.add_argument(
        "--sites-ht-path",
        help=(
            "Override path to the sites HT (input to --annotate-phased-ht). "
            "Use e.g. the production gs://gnomad/v4.1/variant_cooccurrence/"
            "exomes.sites.ht when the postfix-scoped one isn't materialised."
        ),
    )
    parser.add_argument(
        "--variant-filter-ht-path",
        help=(
            "Override path to the variant_filter HT (input to "
            "--annotate-phased-ht). Defaults to the postfix-scoped path."
        ),
    )



    args = parser.parse_args()
    main(args)

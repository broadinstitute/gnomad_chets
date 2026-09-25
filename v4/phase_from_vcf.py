"""Derive per-pair cis/trans calls from an externally phased VCF.

Statistical phasers (Eagle, Beagle, SHAPEIT) emit a VCF of phased genotypes.
This script turns one into the same per-pair representation the co-occurrence
pipeline already uses, so an external phaser can be compared against trio
phase-by-transmission and against our EM estimates on equal terms.

No new phasing logic is needed. :func:`~gnomad_chets.v4.compute_vp_counts.\
_count_phase_from_sets` already converts per-sample physical phase into per-pair
``n_phased_cis`` / ``n_phased_trans``, keyed off the ``PGT``/``PID`` entry
fields that ``encode_genotypes`` feature-detects. A statistically phased VCF is
one phase set per chromosome, so setting ``PGT = GT`` and a constant ``PID``
makes the existing kernel do the work: the ``PID`` equality test becomes
vacuous and the haplotype comparison carries the call.

Two differences from the read-backed path are worth stating. The resulting
call means "statistically phased" rather than "read-backed", and because every
heterozygote is phased the counts cover the whole double-het cell rather than
refining part of it.

Steps::

    --import-vcf      phased VCF -> MatrixTable (cached; the VCFs carry no index)
    --build-vp-list   MatrixTable -> within-gene pair list from its OWN variants
    --call-vcf-chet   MatrixTable + pair list -> per-pair cis/trans Table

Use ``--build-vp-list`` when the production exome pair list overlaps the phased
data poorly. On the 2020 chr20 genome pilot only 2,365 of the pair-list
endpoints survived, dropping 94.6% of pairs, because that run removed
singletons while the exome pair list is dominated by rare coding variants.
Building the list from the phaser's own variants keeps the pairs it actually
phased, at the cost of no longer matching the production list.

Usage (Dataproc)::

    hailctl dataproc submit <CLUSTER> --project broad-mpg-gnomad --region us-central1 \\
      --pyfiles=/abs/path/to/gnomad_chets \\
      /abs/path/to/gnomad_chets/v4/phase_from_vcf.py -- \\
      --import-vcf --call-vcf-chet \\
      --vcf-path gs://... --vp-list-path gs://... \\
      --output-postfix chr21_eagle_v2 --overwrite
"""

import argparse
import logging
import timeit
from uuid import uuid4
from typing import Optional

import hail as hl
from gnomad.resources.grch38.reference_data import gencode
from gnomad.utils.annotations import get_adj_expr

from gnomad_chets.v4.compute_vp_counts import (
    count_all_pairs_auto_split,
    encode_genotypes,
)
from gnomad_chets.v4.create_vp_list import create_variant_pair_ht
from gnomad_chets.v4.resources import DEFAULT_TMP_DIR

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("phase_from_vcf")
logger.setLevel(logging.INFO)


def import_phased_vcf(
    vcf_path: str,
    mt_path: str,
    reference_genome: str,
    min_partitions: int,
    overwrite: bool,
) -> hl.MatrixTable:
    """Import a phased VCF and persist it as a MatrixTable.

    The phased VCFs we have carry no tabix index, so each read re-parses the
    whole file. Every downstream step needs the same MatrixTable, so the parse
    is paid once. Only ``GT`` is retained; in the Eagle outputs every other
    FORMAT field is ``.`` throughout.

    :param vcf_path: Source phased VCF (BGZF).
    :param mt_path: Where to write the MatrixTable.
    :param reference_genome: Reference build of the VCF.
    :param min_partitions: Minimum partitions for the import.
    :param overwrite: Whether to overwrite an existing MatrixTable.
    :return: The phased MatrixTable.
    """
    logger.info("Importing %s", vcf_path)
    mt = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        reference_genome=reference_genome,
        min_partitions=min_partitions,
    )
    mt = mt.select_entries("GT")
    mt.write(mt_path, overwrite=overwrite)
    logger.info("Wrote %s", mt_path)
    return hl.read_matrix_table(mt_path)


def prepare_phased_mt(
    mt: hl.MatrixTable, phase_set: str, adj_mask_mt: hl.MatrixTable = None
) -> hl.MatrixTable:
    """Present externally phased GTs to the pipeline as read-backed phase.

    ``encode_genotypes`` detects physical phase through ``PGT``/``PID`` and
    ``_count_phase_from_sets`` requires both sides of a pair to share a ``PID``
    before comparing haplotypes. A statistically phased VCF is a single phase
    set for the whole chromosome, so a constant ``PID`` is the faithful
    representation and the ``PID`` check passes for every pair.

    ``adj`` defaults to "any defined genotype", because a phased VCF carries no
    per-genotype quality fields to re-derive it from. That is not free: a
    phaser fills in every genotype it was not given, so if the input was
    adj-masked, the masked calls come back imputed and would be scored as
    observations -- inflating the phaser's coverage on exactly the genotypes
    the trio and EM arms discard. Pass ``adj_mask_mt`` (the pre-phasing dense
    MatrixTable, which still has GQ/DP/AD) to require the original genotype to
    have passed adj, which is the like-for-like comparison.

    :param mt: Phased MatrixTable with a ``GT`` entry field.
    :param phase_set: Value to use as the ``PID``, e.g. the contig name.
    :param adj_mask_mt: Optional pre-phasing MatrixTable supplying GQ/DP/AD.
    :return: MatrixTable with ``PGT``, ``PID`` and ``adj`` entry fields.
    """
    adj_expr = hl.is_defined(mt.GT)
    if adj_mask_mt is not None:
        src = adj_mask_mt.annotate_entries(
            _adj=get_adj_expr(
                adj_mask_mt.GT, adj_mask_mt.GQ, adj_mask_mt.DP, adj_mask_mt.AD
            )
        )
        adj_expr = adj_expr & hl.coalesce(
            src[mt.row_key, mt.col_key]._adj, False
        )
    return mt.annotate_entries(
        PGT=mt.GT,
        PID=hl.if_else(mt.GT.phased, hl.literal(phase_set), hl.missing(hl.tstr)),
        adj=adj_expr,
    )


def annotate_vcf_chet(ht: hl.Table) -> hl.Table:
    """Collapse the phase counts into a tri-state call per pair.

    Matches the convention of ``trio_phasing._chet_call``: ``True`` is trans,
    ``False`` is cis, and missing means the pair carries no evidence either way
    or contradictory evidence across samples.

    :param ht: Table with ``n_phased_cis`` and ``n_phased_trans``.
    :return: ``ht`` with a ``vcf_chet`` boolean and ``n_vcf_phased`` total.
    """
    return ht.annotate(
        n_vcf_phased=ht.n_phased_cis + ht.n_phased_trans,
        vcf_chet=hl.case()
        .when((ht.n_phased_trans > 0) & (ht.n_phased_cis == 0), True)
        .when((ht.n_phased_cis > 0) & (ht.n_phased_trans == 0), False)
        .or_missing(),
    )


def call_vcf_chet(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    tmp_dir: str,
    phase_set: str,
    adj_mask_mt: hl.MatrixTable = None,
    count_partitions: int = 1000,
) -> hl.Table:
    """Produce per-pair cis/trans counts from a phased MatrixTable.

    :param mt: Phased MatrixTable.
    :param vp_ht: Variant pair list, keyed however; must carry ``locus1``,
        ``alleles1``, ``locus2``, ``alleles2``.
    :param tmp_dir: Directory for the encoder's intermediates.
    :param phase_set: ``PID`` value, see :func:`prepare_phased_mt`.
    :param count_partitions: Partitions for the all-light counting plan.
        ``count_all_pairs_auto_split`` repartitions the pair list to this before
        the indexed count, so the default of 1000 -- sized for a gene like PCNT
        with millions of pairs -- shatters a small locus into partitions of a
        few dozen pairs and the following shuffle fails non-deterministically
        with ``zip: length mismatch``. Scale it to the pair count.
    :param adj_mask_mt: Optional pre-phasing MatrixTable; see
        :func:`prepare_phased_mt`.
    :return: Table keyed by the pair locus tuple with ``n_phased_cis`` /
        ``n_phased_trans`` / ``vcf_chet``.
    """
    mt = prepare_phased_mt(mt, phase_set, adj_mask_mt)

    # The encoder requires the MatrixTable to hold exactly the pair variants.
    endpoints = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1).select()
    endpoints = endpoints.union(
        vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2).select()
    ).distinct()
    mt = mt.semi_join_rows(endpoints)

    # Unique per invocation. This used to be a fixed path, which made two
    # concurrent extractions overwrite each other's var_idx.ht mid-read -- a
    # 404 on a parts file if you were lucky, silently wrong counts if not.
    encode_dir = f"{tmp_dir}/phase_from_vcf_encoded_{uuid4().hex[:12]}"
    logger.info("Encoder scratch: %s", encode_dir)
    encode_genotypes(mt, vp_ht, encode_dir, use_precomputed_adj=True)
    var_idx_ht = hl.read_table(f"{encode_dir}/var_idx.ht")
    encoded_gt_ht = hl.read_table(f"{encode_dir}/encoded_gt_sets_by_var_idx.ht")

    ht = count_all_pairs_auto_split(
        vp_ht,
        var_idx_ht,
        encoded_gt_ht,
        fallback_partitions=count_partitions,
        emit_phase=True,
    )
    return annotate_vcf_chet(ht)


def build_vp_list_from_mt(
    mt: hl.MatrixTable,
    max_freq: float,
    gene_names: Optional[list] = None,
) -> hl.Table:
    """Build a within-gene variant pair list from a MatrixTable's own variants.

    Assigns genes by overlap with GENCODE protein-coding gene bodies rather than
    by consequence, so intronic variants are kept. That widens the comparison
    beyond coding pairs, which is appropriate when scoring a phaser on genome
    data and is the regime the grant's noncoding aim cares about.

    ``an_pct`` is set to 100 for every variant: the phaser already applied its
    own call-rate filter upstream, and its output carries no missing genotypes
    to re-derive a per-variant AN from.

    :param mt: Phased MatrixTable keyed by ``(locus, alleles)``.
    :param max_freq: Keep variants at or below this alternate allele frequency.
    :param gene_names: Optional GENCODE gene-name allowlist.
    :return: Variant pair Table from :func:`create_variant_pair_ht`.
    """
    genes = gencode.ht()
    genes = genes.filter(
        (genes.feature == "gene") & (genes.gene_type == "protein_coding")
    )
    if gene_names:
        genes = genes.filter(hl.literal(set(gene_names)).contains(genes.gene_name))

    # Allele frequency, not carrier frequency: a hom-alt sample contributes two
    # alternate alleles, so fraction(is_non_ref) would understate common variants.
    mt = mt.annotate_rows(
        _af=hl.agg.sum(mt.GT.n_alt_alleles())
        / (2 * hl.agg.count_where(hl.is_defined(mt.GT)))
    )
    rows = mt.rows()
    idx = genes.index(rows.locus, all_matches=True)
    filter_ht = rows.select(
        gene_id=hl.set(hl.or_else(idx.gene_id, hl.empty_array(hl.tstr))),
        an_pct=hl.int32(100),
        _af=rows._af,
    )
    filter_ht = filter_ht.filter(
        (hl.len(filter_ht.gene_id) > 0) & (filter_ht._af <= max_freq)
    )
    filter_ht = filter_ht.drop("_af")
    logger.info("Variants eligible for pairing: %d", filter_ht.count())

    mt = mt.semi_join_rows(filter_ht)
    return create_variant_pair_ht(mt, filter_ht)


def main(args):
    """Run the requested steps."""
    start = timeit.default_timer()
    hl.init(log="/tmp/phase_from_vcf.log", tmp_dir=args.tmp_dir)

    postfix = f".{args.output_postfix}" if args.output_postfix else ""
    mt_path = args.mt_path or f"{args.tmp_dir}/phased_vcf{postfix}.mt"
    out_path = args.output_path or f"{args.tmp_dir}/vcf_phase{postfix}.ht"

    logger.info(
        "phase_from_vcf\n"
        "  vcf_path      : %s\n"
        "  mt_path       : %s\n"
        "  vp_list_path  : %s\n"
        "  output_path   : %s\n"
        "  reference     : %s\n"
        "  interval      : %s\n"
        "  phase_set     : %s",
        args.vcf_path,
        mt_path,
        args.vp_list_path,
        out_path,
        args.reference_genome,
        args.interval,
        args.phase_set,
    )

    if args.import_vcf:
        import_phased_vcf(
            args.vcf_path,
            mt_path,
            args.reference_genome,
            args.min_partitions,
            args.overwrite,
        )

    vp_list_path = args.vp_list_path or f"{args.tmp_dir}/vcf_vp_list{postfix}.ht"

    if args.build_vp_list:
        mt = hl.read_matrix_table(mt_path)
        if args.interval:
            mt = hl.filter_intervals(
                mt,
                [
                    hl.parse_locus_interval(
                        args.interval, reference_genome=args.reference_genome
                    )
                ],
            )
        vp = build_vp_list_from_mt(mt, args.max_freq, args.genes)
        vp = vp.checkpoint(vp_list_path, overwrite=args.overwrite)
        logger.info("Built %d pairs -> %s", vp.count(), vp_list_path)

    if args.call_vcf_chet:
        mt = hl.read_matrix_table(mt_path)
        vp_ht = hl.read_table(vp_list_path)
        if args.interval:
            interval = hl.parse_locus_interval(
                args.interval, reference_genome=args.reference_genome
            )
            # Both endpoints must fall inside the requested region.
            vp_ht = vp_ht.filter(
                interval.contains(vp_ht.locus1) & interval.contains(vp_ht.locus2)
            )
            mt = hl.filter_intervals(mt, [interval])
        n_pairs = vp_ht.count()
        logger.info("Calling phase for %d pairs.", n_pairs)

        if args.exclude_samples_fam:
            fam = hl.import_fam(args.exclude_samples_fam, delimiter="\t")
            drop = fam.filter(hl.is_defined(fam.id))
            drop = drop.key_by(s=drop.id).select()
            before = mt.count_cols()
            mt = mt.anti_join_cols(drop)
            logger.info(
                "Excluded probands from the phaser vote: %d -> %d samples",
                before, mt.count_cols(),
            )

        adj_mask_mt = (
            hl.read_matrix_table(args.adj_mask_mt) if args.adj_mask_mt else None
        )
        if adj_mask_mt is not None:
            logger.info("Requiring pre-phasing adj from %s", args.adj_mask_mt)
        ht = call_vcf_chet(
            mt,
            vp_ht,
            args.tmp_dir,
            args.phase_set,
            adj_mask_mt,
            count_partitions=args.count_partitions,
        )
        ht = ht.checkpoint(out_path, overwrite=args.overwrite)

        summary = ht.aggregate(
            hl.struct(
                n=hl.agg.count(),
                n_called=hl.agg.count_where(hl.is_defined(ht.vcf_chet)),
                n_trans=hl.agg.count_where(ht.vcf_chet),
                n_cis=hl.agg.count_where(~ht.vcf_chet),
                n_no_double_het=hl.agg.count_where(ht.n_vcf_phased == 0),
            )
        )
        logger.info(
            "Per-pair phase: %d pairs, %d called (%d trans, %d cis), "
            "%d with no double-het carrier.",
            summary.n,
            summary.n_called,
            summary.n_trans,
            summary.n_cis,
            summary.n_no_double_het,
        )
        logger.info("Wrote %s", out_path)

    logger.info("Finished in %.1f s", timeit.default_timer() - start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--import-vcf", action="store_true", help="Import the phased VCF to a MT."
    )
    parser.add_argument(
        "--build-vp-list",
        action="store_true",
        help="Build a within-gene pair list from the MatrixTable's own variants.",
    )
    parser.add_argument(
        "--max-freq",
        type=float,
        default=0.05,
        help="Allele-frequency ceiling for --build-vp-list.",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        help="Restrict --build-vp-list to these GENCODE gene names.",
    )
    parser.add_argument(
        "--call-vcf-chet",
        action="store_true",
        help="Derive per-pair cis/trans calls from the phased MT.",
    )
    parser.add_argument("--vcf-path", help="Phased VCF to import.")
    parser.add_argument("--mt-path", help="MatrixTable path (derived if unset).")
    parser.add_argument("--vp-list-path", help="Variant pair list Table.")
    parser.add_argument("--output-path", help="Output Table path (derived if unset).")
    parser.add_argument(
        "--reference-genome",
        default="GRCh38",
        choices=["GRCh37", "GRCh38"],
        help="Reference build of the VCF. The 2018 v2 chr21 run is GRCh37.",
    )
    parser.add_argument(
        "--interval", help="Restrict to a locus interval, e.g. chr20:5000000-10000000."
    )
    parser.add_argument(
        "--phase-set",
        default="statistical",
        help="Constant PID for the whole VCF; see prepare_phased_mt.",
    )
    parser.add_argument(
        "--min-partitions", type=int, default=500, help="Minimum import partitions."
    )
    parser.add_argument(
        "--count-partitions",
        type=int,
        default=1000,
        help=(
            "Partitions for the all-light counting plan. The default suits a "
            "gene with millions of pairs; a small locus needs far fewer or the "
            "counting shuffle fails with 'zip: length mismatch'."
        ),
    )
    parser.add_argument(
        "--exclude-samples-fam",
        help=(
            "Pedigree .fam whose probands are dropped before the per-pair "
            "vote. The phasers keep probands in the phasing cohort while the "
            "EM arm removes every trio member, so by default the phaser is "
            "scored in-sample and EM out-of-sample. Dropping probands here "
            "makes the phaser answer the same question EM does: does the "
            "population's phase for this pair hold in someone who did not "
            "contribute to it. Pairs carried only by probands go uncalled, "
            "which is the honest deployment answer."
        ),
    )
    parser.add_argument(
        "--adj-mask-mt",
        help=(
            "Pre-phasing dense MatrixTable (with GQ/DP/AD). A phaser fills in "
            "every genotype it was not given, so without this a genotype "
            "masked out of the phasing input returns imputed and is scored as "
            "an observation. Supplying it requires the original genotype to "
            "have passed adj, matching the trio and EM arms."
        ),
    )
    parser.add_argument("--output-postfix", help="Postfix for derived output paths.")
    parser.add_argument(
        "--tmp-dir", default=DEFAULT_TMP_DIR, help="Temporary directory."
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing outputs."
    )
    args = parser.parse_args()

    if not (args.import_vcf or args.build_vp_list or args.call_vcf_chet):
        parser.error(
            "Specify at least one step (--import-vcf, --build-vp-list, "
            "--call-vcf-chet)."
        )
    if args.import_vcf and not args.vcf_path:
        parser.error("--import-vcf requires --vcf-path.")
    if args.call_vcf_chet and not (args.vp_list_path or args.build_vp_list):
        parser.error("--call-vcf-chet requires --vp-list-path or --build-vp-list.")

    main(args)

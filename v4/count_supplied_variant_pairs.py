"""
Genotype counts for a user-supplied list of variant pairs.

Unlike ``create_vp_list.py`` + ``compute_vp_counts.py``, which *discover* co-occurring
variant pairs genome-wide and write into the pipeline's resource paths, this script
takes a pre-built list of specific pairs you already care about, counts them, and hands
back your own table with the counts annotated on.

Written for Rachel Ungar's ``vp_genotype_matrix.py`` use case (branch ``add_phasing``):
same inputs, same output columns, same one-command invocation -- but the counting is
done by importing ``compute_vp_counts`` rather than reimplementing it, so this
automatically gets the light/heavy split, the size-balanced partitioning, the v4 high-AB
het correction and the sex-ploidy adjustment. The standalone script has none of those,
and its per-pair indexed join is what makes a single gene take hours: profiling FKRP
(10,190 pairs, 2,390 variants) showed 25 GiB of encoded-set payload dragged across that
join, ~98% of the count stage's wall time.

Input
-----
A Hail Table with ``locus1``, ``alleles1``, ``locus2``, ``alleles2`` row fields. The
table's key can be anything -- e.g. keyed by ``locus, alleles, gene, gene_id`` with
``locus1``/``alleles1`` duplicating ``locus``/``alleles``. Every other column is
preserved and returned as-is.

Outputs
-------
``--output`` gives the input table with 18 columns appended, one per two-variant
genotype class, for both raw and adj genotypes::

    raw_AABB, raw_AABb, raw_AAbb, raw_AaBB, raw_AaBb, raw_Aabb, raw_aaBB, raw_aaBb, raw_aabb
    adj_AABB, adj_AABb, adj_AAbb, adj_AaBB, adj_AaBb, adj_Aabb, adj_aaBB, adj_aaBb, adj_aabb

"A/a" is variant 1 (hom-ref/het/hom-var), "B/b" is variant 2, so ``AaBb`` is the
double-het compound-het candidate cell. Pairs whose variants aren't in the callset keep
their row but get missing counts -- distinct from a zero count, which means the variant
was found and nobody carried it.

``--emit-em-phase`` adds EM haplotype counts and ``p_chet`` (the probability the pair
is in trans) to the ``--output`` table.

``--gt-counts-output`` optionally writes the same counts a second time in the
pipeline's array-shaped schema (``gt_counts_raw`` / ``gt_counts_adj``), which is what
``run_in_trans_oe.py --gt-counts-ht-path`` consumes.

Pair ordering is left exactly as supplied. Counts are correct for whatever orientation
you give, but note the published tables order each pair by locus position (tie-broken on
alt allele); a pair supplied the other way round won't join against them, and matching
unordered leaves the cells transposed (AABb<->AaBB, AAbb<->aaBB, Aabb<->aaBb).

Example
-------
.. code-block:: bash

    hailctl dataproc submit CLUSTER \\
      --pyfiles=/abs/path/to/gnomad_chets \\
      /abs/path/to/gnomad_chets/v4/count_supplied_variant_pairs.py -- \\
      --variant-pair-list-ht gs://my-bucket/FKRP_variant_pairs.ht \\
      --output gs://my-bucket/FKRP_variant_pairs.counts.ht \\
      --tmp-dir gs://my-tmp-bucket/fkrp
"""
import argparse
import logging
from typing import List, Optional

import hail as hl
from gnomad.utils.file_utils import file_exists
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds

from gnomad_chets.v4.compute_vp_counts import (
    TARGET_HEAVY_PARTITION_BYTES,
    _size_info_to_heavy_variants,
    build_variant_size_info_ht,
    compute_counts_heavy,
    compute_counts_light,
    count_all_pairs_via_index,
    densify_encode_input_mt,
    encode_genotypes,
)
from gnomad_chets.v4.phase_gnomad import get_em_expr
from gnomad_chets.v4.resources import DATA_TYPE_CHOICES, DEFAULT_DATA_TYPE

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("count_supplied_variant_pairs")
logger.setLevel(logging.INFO)

GENOTYPE_CLASSES = [
    "AABB", "AABb", "AAbb",
    "AaBB", "AaBb", "Aabb",
    "aaBB", "aaBb", "aabb",
]
"""The 9 two-variant genotype classes, in gt_counts array order."""

PAIR_FIELDS = ("locus1", "alleles1", "locus2", "alleles2")


def get_variants_ht(vp_ht: hl.Table) -> hl.Table:
    """
    Get the distinct variants referenced by a variant pair Table.

    :param vp_ht: Variant pair Table with locus1/alleles1/locus2/alleles2.
    :return: Distinct Table keyed by (locus, alleles).
    """
    v1 = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1).select().distinct()
    v2 = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2).select().distinct()

    return v1.union(v2).distinct()


def get_covering_intervals(
    variants_ht: hl.Table, padding: int = 100
) -> List[hl.utils.Interval]:
    """
    Compute one small locus interval per exact variant.

    One tight interval per variant rather than one bounding box per contig: a bounding
    box still scans every partition between two variants that happen to be far apart,
    which for a scattered pair list is most of the callset. The exact variant filter
    downstream does the real precision filtering; this only bounds what gets scanned.

    Overlapping and adjacent intervals are merged before returning, so a compact gene
    collapses to a handful rather than one interval per variant. Merging is loss-free:
    the union of covered loci is identical either way, so a scattered pair list keeps
    exactly the tight per-variant pruning it needs.

    Do not expect this to be faster. Measured on FKRP, collapsing 2,009 per-variant
    intervals to 10 left the densify unchanged (~104 min merged against ~106 min
    unmerged, on matching 2-worker clusters). The stage is bound by reading and
    densifying the VDS across ~730k samples, not by interval bookkeeping. The merge is
    kept because it is loss-free and bounds the interval count for pathological inputs
    -- it is not an optimisation, and interval count is not the lever worth pulling.

    (A sweep over sorted start/end events, the same shape as ``split_interval_ht`` in
    tgg_methods' ``gnomad_small_variant_list_query`` -- that one splits overlapping
    intervals into disjoint segments and maps each back to its variants, which is the
    opposite accumulation from the merge wanted here.)

    :param variants_ht: Table keyed by (locus, alleles).
    :param padding: Bp padding either side of each variant.
    :return: Minimal list of merged intervals covering every variant.
    """
    # Collect the expression rather than re-keying the Table: `ht.key_by()` returns a
    # new source object, so referring to `variants_ht.locus` afterwards raises
    # "Cannot combine expressions from different source objects".
    loci = variants_ht.locus.collect()
    rg = variants_ht.locus.dtype.reference_genome
    intervals = []
    for locus in loci:
        contig = locus.contig
        start = max(1, locus.position - padding)
        end = min(rg.contig_length(contig), locus.position + padding)
        intervals.append(
            hl.utils.Interval(
                hl.Locus(contig, start, reference_genome=rg),
                hl.Locus(contig, end, reference_genome=rg),
                includes_end=True,
            )
        )

    contig_order = rg.contigs
    intervals.sort(key=lambda i: (contig_order.index(i.start.contig), i.start.position))
    merged = [intervals[0]]
    for interval in intervals[1:]:
        last = merged[-1]
        # includes_end=True, so touching intervals (start == end + 1) also merge.
        if (
            interval.start.contig == last.end.contig
            and interval.start.position <= last.end.position + 1
        ):
            if interval.end.position > last.end.position:
                merged[-1] = hl.utils.Interval(
                    last.start, interval.end, includes_end=True
                )
        else:
            merged.append(interval)

    logger.info(
        "Merged %d per-variant intervals into %d covering interval(s).",
        len(intervals), len(merged),
    )

    return merged


def _reuse_or_write(build_fn, path: str, overwrite: bool) -> hl.Table:
    """Write ``build_fn()`` to ``path``, or read it back if it's already there.

    Lets a rerun pick up from the last completed step instead of redoing the densify,
    the size-info build or a finished count. Cutoff-dependent artifacts are written to
    cutoff-tagged paths by the caller, so retuning the cutoff picks up fresh ones
    without ``--overwrite-intermediates`` and without another densify.
    """
    if not overwrite and file_exists(f"{path}/_SUCCESS"):
        logger.info("Reusing existing %s", path)
        return hl.read_table(path)

    return build_fn().checkpoint(path, overwrite=True)


def count_supplied_pairs(
    vp_ht: hl.Table,
    tmp_dir: str,
    data_type: str = DEFAULT_DATA_TYPE,
    *,
    interval_padding: int = 100,
    heavy_contribution_cutoff: Optional[int] = None,
    release_only: bool = True,
    overwrite: bool = False,
) -> hl.Table:
    """
    Count genotypes for every pair in ``vp_ht``.

    Densifies just this pair list's variants out of the VDS, encodes them once into
    per-variant sample-index sets, then counts with the light/heavy split so a hub
    variant's payload isn't replicated into every pair it participates in.

    :param vp_ht: Pair Table carrying locus1/alleles1/locus2/alleles2 (any key).
    :param tmp_dir: Directory for the encode intermediates, which are reused on rerun
        unless ``overwrite`` is set.
    :param data_type: 'exomes' or 'genomes'.
    :param interval_padding: Bp padding per variant for partition pruning; negative
        disables interval pruning entirely.
    :param heavy_contribution_cutoff: Bytes; a variant is heavy at or above this
        ``degree x payload`` contribution. Defaults to TARGET_HEAVY_PARTITION_BYTES.
    :param release_only: Restrict to release samples. Matches the published counts.
    :param overwrite: Recompute the encode even if intermediates already exist.
    :return: Pair Table keyed by the 4 pair fields with gt_counts_raw / gt_counts_adj.
    """
    if heavy_contribution_cutoff is None:
        heavy_contribution_cutoff = TARGET_HEAVY_PARTITION_BYTES

    # Deliberately UNkeyed. compute_counts_heavy does
    # ``select("v1_idx", "v2_idx", "locus1", ...)``, and Hail's check_keys rejects a
    # select() naming the table's own key fields -- so handing it a table keyed on the
    # pair fields fails there. The production pair lists aren't keyed on them either.
    pair_key_ht = vp_ht.key_by().select(*PAIR_FIELDS)

    encode_dir = f"{tmp_dir}/genotype_count_intermediates"
    var_idx_path = f"{encode_dir}/var_idx.ht"
    encoded_path = f"{encode_dir}/encoded_gt_sets_by_var_idx.ht"

    if not overwrite and file_exists(f"{encoded_path}/_SUCCESS"):
        logger.info("Reusing existing encode intermediates in %s", encode_dir)
    else:
        variants_ht = get_variants_ht(vp_ht).checkpoint(
            f"{tmp_dir}/variants.ht", overwrite=True
        )
        n_variants = variants_ht.count()
        intervals = (
            get_covering_intervals(variants_ht, padding=interval_padding)
            if interval_padding >= 0
            else None
        )
        logger.info(
            "Densifying %d variants out of the %s VDS (%d covering intervals).",
            n_variants, data_type, len(intervals) if intervals else 0,
        )
        get_vds_func = (
            get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
        )
        mt = densify_encode_input_mt(
            get_vds_func,
            variants_ht,
            data_type,
            release_only=release_only,
            filter_intervals=intervals,
        )
        # min_an_pct=-1 keeps every supplied pair: this script's contract is a row per
        # pair the caller asked for, so the AN floor the production pipeline uses to
        # bound set sizes is deliberately not applied here.
        encode_genotypes(mt, pair_key_ht, encode_dir, min_an_pct=-1)

    var_idx_ht = hl.read_table(var_idx_path)
    encoded_gt_ht = hl.read_table(encoded_path)

    # The size-info HT carries no heavy/light decision, so it is cutoff-independent
    # and reused across cutoffs. Everything downstream of the cutoff is written to a
    # cutoff-tagged path instead, so retuning --heavy-contribution-cutoff costs a
    # recount but never another densify.
    cutoff_tag = f"cutoff{heavy_contribution_cutoff}"
    size_info_ht = _reuse_or_write(
        lambda: build_variant_size_info_ht(encoded_gt_ht, pair_key_ht, var_idx_ht),
        f"{tmp_dir}/variant_size_info.ht", overwrite,
    )
    heavy_variants = _reuse_or_write(
        lambda: _size_info_to_heavy_variants(size_info_ht, heavy_contribution_cutoff),
        f"{tmp_dir}/heavy_variants.{cutoff_tag}.ht", overwrite,
    )

    n_heavy = heavy_variants.count()
    total_contribution = size_info_ht.aggregate(
        hl.agg.sum(size_info_ht._contribution)
    )
    logger.info(
        "Total contribution %.2f GiB across %d variants; %d heavy at the %.0f MB "
        "cutoff.",
        total_contribution / 1024**3, size_info_ht.count(), n_heavy,
        heavy_contribution_cutoff / 1024**2,
    )

    if n_heavy == 0:
        # compute_counts_light assumes a non-empty heavy set; with nothing heavy there
        # is no split to plan and the indexed lookup is the intended path.
        logger.info("No heavy variants; counting every pair via indexed lookup.")
        return count_all_pairs_via_index(pair_key_ht, var_idx_ht, encoded_gt_ht)

    light_ht = _reuse_or_write(
        lambda: compute_counts_light(
            pair_key_ht, var_idx_ht, encoded_gt_ht, heavy_variants,
            size_info_ht=size_info_ht,
        ),
        f"{tmp_dir}/counts_light.{cutoff_tag}.ht", overwrite,
    )
    heavy_ht = _reuse_or_write(
        lambda: compute_counts_heavy(
            pair_key_ht, var_idx_ht, encoded_gt_ht, heavy_variants,
        ),
        f"{tmp_dir}/counts_heavy.{cutoff_tag}.ht", overwrite,
    )
    logger.info(
        "Counted %d light pairs and %d heavy pairs.",
        light_ht.count(), heavy_ht.count(),
    )

    # Project both sides to the same row type before unioning: the two count paths
    # don't guarantee identical field order, and Table.union requires an exact match.
    # Both are keyed by the 4 pair fields, so selecting the non-key count columns
    # leaves the key untouched.
    count_cols = ["gt_counts_raw", "gt_counts_adj"]

    return light_ht.select(*count_cols).union(heavy_ht.select(*count_cols))


def main(args):
    """Annotate a supplied variant pair list with genotype counts."""
    hl.init(
        log="/tmp/count_supplied_variant_pairs.log",
        tmp_dir=args.tmp_dir,
        default_reference="GRCh38",
    )

    vp_ht = hl.read_table(args.variant_pair_list_ht)
    missing = [f for f in PAIR_FIELDS if f not in vp_ht.row]
    if missing:
        raise ValueError(
            f"--variant-pair-list-ht is missing required field(s): {missing}. "
            f"Needs {list(PAIR_FIELDS)}."
        )
    orig_key = list(vp_ht.key)
    logger.info(
        "Read %d variant pairs from %s", vp_ht.count(), args.variant_pair_list_ht
    )

    counts_ht = count_supplied_pairs(
        vp_ht,
        tmp_dir=args.tmp_dir,
        data_type=args.data_type,
        interval_padding=args.interval_padding,
        heavy_contribution_cutoff=args.heavy_contribution_cutoff,
        release_only=not args.all_high_quality_samples,
        overwrite=args.overwrite_intermediates,
    )

    if args.gt_counts_output:
        # Write the array-shaped counts too. run_in_trans_oe.py's --gt-counts-ht-path
        # wants gt_counts_raw / gt_counts_adj, not the flat per-class columns, and this
        # is exactly that shape before it gets flattened -- so checkpointing here saves
        # the caller a conversion pass and feeds the flattening below from disk.
        counts_ht = counts_ht.checkpoint(
            args.gt_counts_output, overwrite=args.overwrite
        )
        logger.info("Wrote array-shaped genotype counts to %s", args.gt_counts_output)

    if args.emit_em_phase:
        # Reuses the pipeline's EM (phase_gnomad.get_em_expr) rather than a second
        # implementation, so p_chet here means exactly what it means everywhere else.
        counts_ht = counts_ht.annotate(
            em_raw=get_em_expr(counts_ht.gt_counts_raw),
            em_adj=get_em_expr(counts_ht.gt_counts_adj),
        )
        counts_ht = counts_ht.transmute(
            hap_counts_raw=counts_ht.em_raw.hap_counts,
            p_chet_raw=counts_ht.em_raw.p_chet,
            hap_counts_adj=counts_ht.em_adj.hap_counts,
            p_chet_adj=counts_ht.em_adj.p_chet,
        )
        logger.info("Annotated EM haplotype counts and p_chet (raw + adj).")

    em_cols = (
        ["hap_counts_raw", "p_chet_raw", "hap_counts_adj", "p_chet_adj"]
        if args.emit_em_phase
        else []
    )
    counts_ht = counts_ht.select(
        *em_cols,
        **{
            f"raw_{name}": counts_ht.gt_counts_raw[i]
            for i, name in enumerate(GENOTYPE_CLASSES)
        },
        **{
            f"adj_{name}": counts_ht.gt_counts_adj[i]
            for i, name in enumerate(GENOTYPE_CLASSES)
        },
    )
    out_ht = vp_ht.annotate(
        **counts_ht[vp_ht.locus1, vp_ht.alleles1, vp_ht.locus2, vp_ht.alleles2]
    ).key_by(*orig_key)

    if args.output_format == "ht":
        out_ht.write(args.output, overwrite=args.overwrite)
    else:
        # Unkey first: the pair fields are usually part of the key, and Hail refuses to
        # annotate over a key field. Arrays are delimited rather than left as Hail's
        # bracketed repr so the result opens cleanly in a spreadsheet.
        out_ht = out_ht.key_by()
        str_exprs = {
            f: hl.str(out_ht[f]) for f in ("locus1", "locus2") if f in out_ht.row
        }
        str_exprs.update(
            {
                f: hl.delimit(out_ht[f], ",")
                for f in ("alleles1", "alleles2")
                if f in out_ht.row
            }
        )
        if args.emit_em_phase:
            str_exprs.update(
                {
                    f: hl.delimit(out_ht[f].map(hl.str), ",")
                    for f in ("hap_counts_raw", "hap_counts_adj")
                }
            )
        out_ht = out_ht.annotate(**str_exprs)
        out_ht.flatten().export(args.output)
    logger.info("Wrote genotype counts to %s", args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--variant-pair-list-ht",
        required=True,
        help=(
            "Path to a Hail Table of variant pairs with locus1/alleles1/locus2/"
            "alleles2 row fields. Any other columns are preserved on the output."
        ),
    )
    parser.add_argument("--output", required=True, help="Path to write the result to.")
    parser.add_argument(
        "--emit-em-phase",
        action="store_true",
        help=(
            "Also emit EM haplotype counts and p_chet (raw and adj) per pair, via the "
            "same phase_gnomad.get_em_expr the rest of the pipeline uses. p_chet is "
            "the probability the two variants are in trans; it is missing where the EM "
            "denominator collapses, which is every pair with no double-het carrier."
        ),
    )
    parser.add_argument(
        "--gt-counts-output",
        help=(
            "Optional second output path. Writes the same counts in the pipeline's "
            "array-shaped schema (gt_counts_raw / gt_counts_adj, keyed by the four "
            "pair fields), which is what run_in_trans_oe.py --gt-counts-ht-path "
            "expects. The --output table keeps the flat per-class columns and the "
            "caller's own columns."
        ),
    )
    parser.add_argument(
        "--output-format",
        default="ht",
        choices=("ht", "tsv"),
        help="Write a Hail Table, or a flattened text export. Default is ht.",
    )
    parser.add_argument(
        "--tmp-dir",
        required=True,
        help=(
            "Directory for the encode intermediates and per-step checkpoints. Reused "
            "on rerun unless --overwrite-intermediates is passed, so a killed run "
            "picks up from the last completed step instead of re-densifying."
        ),
    )
    parser.add_argument(
        "--data-type",
        default=DEFAULT_DATA_TYPE,
        choices=DATA_TYPE_CHOICES,
        help=f"Data type to count against. Default is {DEFAULT_DATA_TYPE}.",
    )
    parser.add_argument(
        "--interval-padding",
        type=int,
        default=100,
        help=(
            "Bp padding around each variant for the partition-pruning intervals. "
            "Negative disables interval pruning. Default 100."
        ),
    )
    parser.add_argument(
        "--heavy-contribution-cutoff",
        type=int,
        default=None,
        help=(
            "Bytes. A variant is heavy at or above this degree x payload "
            "contribution, and its pairs are counted on the split-aware heavy path "
            f"instead of the light one. Defaults to {TARGET_HEAVY_PARTITION_BYTES} "
            "(500 MB)."
        ),
    )
    parser.add_argument(
        "--all-high-quality-samples",
        action="store_true",
        help=(
            "Count over all high-quality samples rather than release samples only. "
            "Off by default: the published co-occurrence counts are release-only, and "
            "high-quality additionally includes related and non-releasable samples."
        ),
    )
    parser.add_argument(
        "--overwrite-intermediates",
        action="store_true",
        help=(
            "Recompute every intermediate in --tmp-dir (encode, size info, heavy "
            "set, both count tables) even if it already exists. Not needed to change "
            "--heavy-contribution-cutoff: cutoff-dependent artifacts are written to "
            "cutoff-tagged paths, so a new cutoff recounts without re-densifying."
        ),
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite the output if it exists."
    )

    main(parser.parse_args())

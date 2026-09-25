"""Densify a region out of the gnomAD VDS and export a VCF for an external phaser.

Produces the input for running Eagle (or Beagle/SHAPEIT) ourselves, so a phaser
can be scored on genome data where it has a common-variant scaffold to work
with. The motivating question: our chr21 exome result put Eagle 18 points below
gnomAD EM, but exome capture gives a phaser almost no common non-coding
variants to anchor rare ones to, so that result may be about the data type
rather than the tool.

Scope choices worth knowing:

* **Singletons are dropped.** They are ~48% of polymorphic SNPs in the target
  region and a phaser provably cannot place them -- there is no second carrier
  whose haplotype background could be matched. Including them would roughly
  double runtime for variants whose answer is already known to be near-random
  (SHAPEIT5 reports 35.2% switch error on exome singletons). The resulting
  comparison therefore covers AC>=2 only.
* **The cohort is the ``high_quality`` set, not the release.** gnomAD release
  sets are pruned to unrelated individuals, so a complete trio never survives
  into one -- trio work has to run on the pre-release ``high_quality`` callset
  including related and unreleasable samples. For v4 genomes that is 140,508
  samples against 76,215 in the release, so it is also the better phasing
  cohort. All 6,743 complete genome trios are inside it.
* **Trio parents are excluded from that cohort**, not the probands. That leaves
  every proband leakage-free for the trio-truth comparison while keeping the
  cohort large, which matters because phasing accuracy scales with the number
  of haplotypes available to copy from. The 2020 pilot phased probands only and
  gave up an order of magnitude of cohort size to get the same property.
* **Variants are selected on release frequency** even though the cohort is
  larger, because the release table is the standard reference and needs no
  densify to read. A variant that is a singleton in the release but polymorphic
  in the fuller cohort is therefore excluded; that is conservative, not neutral.
* **Only the scaffold region needs this job.** The downstream EM, trio-truth and
  read-backed arms need the gene alone, which is a far smaller densify; do not
  persist this region's dense MatrixTable for them.

Usage (Dataproc)::

    hailctl dataproc submit <CLUSTER> --project broad-mpg-gnomad --region us-central1 \\
      --pyfiles=/abs/path/to/gnomad_chets \\
      /abs/path/to/gnomad_chets/v4/export_for_phasing.py -- \\
      --region chr21:44324141-46709983 --data-type genomes \\
      --exclude-trio-parents --pedigree-path gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_trios.fam \\
      --out-vcf gs://gnomad-tmp-30day/pcnt_region.genomes.vcf.bgz
"""

import argparse
import logging
import timeit
from typing import Optional

import hail as hl

from gnomad.utils.annotations import get_adj_expr

from gnomad_chets.v4.compute_vp_counts import _split_variant_data_keeping_phase
from gnomad_chets.v4.resources import DEFAULT_TMP_DIR

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("export_for_phasing")
logger.setLevel(logging.INFO)


def build_phasing_variant_ht(
    region: str,
    data_type: str,
    min_ac: int,
    min_an_frac: float,
    snps_only: bool,
    pass_only: bool = False,
) -> hl.Table:
    """Select the variants an external phaser should be given.

    Built from the release frequency table rather than from the densified data,
    so the variant set is known before any densify happens and the job size is
    bounded up front.

    :param region: Locus interval string.
    :param data_type: ``exomes`` or ``genomes``.
    :param min_ac: Minimum alternate allele count; 2 drops singletons.
    :param min_an_frac: Call-rate floor as a fraction of the maximum AN seen.
    :param snps_only: Restrict to SNPs, matching how Eagle is usually run.
    :param pass_only: Drop variants failing variant QC, so every arm is scored
        on the same variants. Follows the pipeline convention: PASS, or
        filtered for ``AC0`` alone.
    :return: Table keyed by ``(locus, alleles)``.
    """
    from gnomad_qc.v4.resources.annotations import get_freq

    interval = hl.parse_locus_interval(region, reference_genome="GRCh38")
    ht = get_freq(data_type=data_type).ht()
    ht = hl.filter_intervals(ht, [interval])
    ht = ht.annotate(_ac=ht.freq[0].AC, _an=ht.freq[0].AN)

    max_an = ht.aggregate(hl.agg.max(ht._an))
    logger.info("Max AN in region: %d (~%d samples)", max_an, max_an // 2)

    keep = (ht._ac >= min_ac) & (ht._an >= min_an_frac * max_an)
    if snps_only:
        keep = keep & hl.is_snp(ht.alleles[0], ht.alleles[1])
    ht = ht.filter(keep)
    ht = ht.select()

    if pass_only:
        from gnomad.resources.grch38.gnomad import public_release

        sites = public_release(data_type).ht()
        sites = hl.filter_intervals(sites, [interval])
        sites = sites.filter(
            sites.filters.difference(hl.set(["AC0"])).length() == 0
        ).select()
        before = ht.count()
        ht = ht.semi_join(sites)
        logger.info(
            "PASS filter: %d of %d variants kept", ht.count(), before
        )

    n = ht.count()
    logger.info("Variants selected for phasing: %d", n)
    if n == 0:
        raise RuntimeError("No variants passed the filters; check --region.")
    return ht


def trio_members_ht(pedigree_path: str) -> hl.Table:
    """Every individual named in a pedigree: probands and both parents."""
    fam = hl.import_fam(pedigree_path, delimiter="\t")
    parts = []
    for field in ("id", "pat_id", "mat_id"):
        sub = fam.filter(hl.is_defined(fam[field]))
        parts.append(sub.key_by(s=sub[field]).select())
    return parts[0].union(*parts[1:]).distinct()


def get_trio_parents_ht(pedigree_path: str) -> hl.Table:
    """Collect every parent in a pedigree, as samples to exclude.

    Excluding parents rather than probands keeps the cohort large while making
    each proband's phase independent of its own family, so trio transmission
    stays a clean reference.

    :param pedigree_path: ``.fam`` file path.
    :return: Table keyed by ``s`` holding the parent sample IDs.
    """
    fam = hl.import_fam(pedigree_path, delimiter="\t")
    # Rebind after each filter: chaining filter().key_by() referencing the
    # pre-filter handle breaks Hail's source tracking (CLAUDE.md).
    pats = fam.filter(hl.is_defined(fam.pat_id))
    pats = pats.key_by(s=pats.pat_id).select()
    mats = fam.filter(hl.is_defined(fam.mat_id))
    mats = mats.key_by(s=mats.mat_id).select()
    parents = pats.union(mats).distinct()
    logger.info("Trio parents to exclude: %d", parents.count())
    return parents


def get_high_quality_samples_ht(data_type: str) -> hl.Table:
    """The ``high_quality`` sample set, related and unreleasable included.

    ``high_quality`` is the gnomAD sample-QC flag -- not hard-filtered and no
    sample-QC metric outlier flags -- which is what the release frequencies are
    computed on. The release set additionally prunes to unrelated individuals
    and drops unreleasable samples, so it contains no complete trios;
    phase-by-transmission therefore needs this wider set, but no wider.

    NOT ``~sample_filters.hard_filtered``: that is 149,623 v4 genomes against
    140,508 here, admitting 9,115 samples that fail ``high_quality`` and would
    put the phasing cohort on a different sample basis than the release
    frequencies the EM arm is compared against.

    :param data_type: ``exomes`` or ``genomes``.
    :return: Table keyed by ``s``.
    """
    from gnomad_qc.v4.resources.meta import meta

    ht = meta(data_type=data_type).ht()
    ht = ht.filter(ht.high_quality)
    ht = ht.select()
    logger.info("high_quality %s samples: %d", data_type, ht.count())
    return ht


def build_phasing_cohort_ht(
    hq_ht: hl.Table, pedigree_path: str, target_n: int, seed: int
) -> hl.Table:
    """Cap the cohort at ``target_n`` while keeping every trio proband.

    Phasing runtime is roughly linear in samples, so a large cohort can be
    unaffordable. Probands are what gets evaluated against trio truth, so they
    are kept unconditionally; parents are dropped so probands stay leakage-free;
    the remainder is filled with unrelated samples chosen reproducibly.

    Cutting the cohort does cost phasing accuracy -- fewer haplotypes are
    available to copy from -- so the resulting number is a floor for what this
    phaser would achieve at full cohort size. Say so alongside any figure.

    :param hq_ht: High-quality samples, keyed by ``s``.
    :param pedigree_path: ``.fam`` supplying probands and parents.
    :param target_n: Target cohort size.
    :param seed: Seed for the fill sample, for reproducibility.
    :return: Table keyed by ``s``.
    """
    fam = hl.import_fam(pedigree_path, delimiter="\t")
    probands = fam.filter(hl.is_defined(fam.id))
    probands = probands.key_by(s=probands.id).select()
    probands = hq_ht.semi_join(probands)
    n_probands = probands.count()

    members = trio_members_ht(pedigree_path)
    fillers = hq_ht.anti_join(members)
    n_fill = max(target_n - n_probands, 0)
    n_avail = fillers.count()
    if n_fill > n_avail:
        logger.warning("Only %d fill samples available; using all.", n_avail)
        n_fill = n_avail
    # Deterministic random subset: seeded uniform, then take the smallest.
    fillers = fillers.annotate(_r=hl.rand_unif(0, 1, seed=seed))
    fillers = fillers.order_by("_r").head(n_fill)
    fillers = fillers.key_by("s").select()

    cohort = probands.union(fillers).distinct()
    logger.info(
        "Phasing cohort: %d probands + %d unrelated fill = %d",
        n_probands,
        n_fill,
        cohort.count(),
    )
    return cohort


def _export_vcf_from_mt(mt: hl.MatrixTable, args) -> None:
    """Drop to GT, optionally adj-masking first, and write the VCF."""
    if args.adj_only:
        # The phasers emit only GT, so a phased VCF carries no quality
        # fields and `phase_from_vcf` has to treat every defined genotype
        # as adj. Masking here is what puts the phaser arms on the same
        # genotype-quality bar as the trio and EM arms, which derive adj
        # from GQ/DP/AD themselves off the checkpointed MT.
        mt = mt.annotate_entries(
            GT=hl.or_missing(get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD), mt.GT)
        )
        logger.info("adj mask applied to the exported genotypes")
    vcf_mt = mt.select_entries("GT").select_rows()
    hl.export_vcf(vcf_mt, args.out_vcf, tabix=args.tabix)
    logger.info("Wrote %s", args.out_vcf)


def main(args):
    """Densify the region and write a phaser-ready VCF."""
    start = timeit.default_timer()
    hl.init(log="/tmp/export_for_phasing.log", tmp_dir=args.tmp_dir)

    if args.from_checkpoint:
        # Export-only half of a split run. The MT on disk already carries the
        # cohort, the variant set and the filters, so none of the selection
        # work above is repeated -- that is the point, since it is what makes
        # the densify and the export separately costable.
        logger.info("Reading dense MT: %s", args.from_checkpoint)
        mt = hl.read_matrix_table(args.from_checkpoint)
        logger.info("Dense MT: %d samples", mt.count_cols())
        _export_vcf_from_mt(mt, args)
        logger.info("Finished in %.1f s", timeit.default_timer() - start)
        return

    logger.info(
        "export_for_phasing\n"
        "  region      : %s\n"
        "  data_type   : %s\n"
        "  min_ac      : %d (2 drops singletons)\n"
        "  min_an_frac : %s\n"
        "  out_vcf     : %s",
        args.region,
        args.data_type,
        args.min_ac,
        args.min_an_frac,
        args.out_vcf,
    )

    variant_ht = build_phasing_variant_ht(
        args.region,
        args.data_type,
        args.min_ac,
        args.min_an_frac,
        not args.keep_indels,
        args.pass_only,
    )

    exclude_ht: Optional[hl.Table] = None
    if args.exclude_trio_parents:
        if not args.pedigree_path:
            raise RuntimeError("--exclude-trio-parents requires --pedigree-path.")
        exclude_ht = get_trio_parents_ht(args.pedigree_path)
    if args.exclude_trio_probands:
        if not args.pedigree_path:
            raise RuntimeError("--exclude-trio-probands requires --pedigree-path.")
        fam = hl.import_fam(args.pedigree_path, delimiter="\t")
        probands = fam.filter(hl.is_defined(fam.id))
        probands = probands.key_by(s=probands.id).select()
        exclude_ht = probands if exclude_ht is None else exclude_ht.union(probands)
        exclude_ht = exclude_ht.distinct()
        logger.info("Excluding probands too; exclusion set: %d", exclude_ht.count())

    from gnomad_qc.v4.resources.basics import (
        get_gnomad_v4_genomes_vds,
        get_gnomad_v4_vds,
    )

    intervals = [hl.parse_locus_interval(args.region, reference_genome="GRCh38")]

    samples_ht = None
    if args.cohort == "high_quality":
        samples_ht = get_high_quality_samples_ht(args.data_type)
        if args.target_n:
            samples_ht = build_phasing_cohort_ht(
                samples_ht, args.pedigree_path, args.target_n, args.seed
            )

    # densify_encode_input_mt calls its loader with `high_quality_only`, which
    # only the exomes loader accepts -- its documented genomes support has never
    # been exercised. Read the VDS here and hand it over through the `vds=`
    # escape hatch (the same one --vds-subset uses), which skips the loader call
    # and assumes release/interval restriction is already applied.
    if args.data_type == "exomes":
        vds = get_gnomad_v4_vds(
            release_only=args.cohort == "release",
            high_quality_only=args.cohort == "high_quality",
            split=False,
            filter_intervals=intervals,
            split_reference_blocks=False,
        )
    else:
        # The genomes loader has no high_quality_only, and its
        # remove_hard_filtered_samples default reads
        # gs://gnomad-autoclass/sample_qc/ht/genomes_v3.1/gnomad_v3.1_hard_filtered_samples.ht,
        # which does not exist. Supply the sample list from the metadata instead.
        vds = get_gnomad_v4_genomes_vds(
            release_only=args.cohort == "release",
            remove_hard_filtered_samples=False,
            filter_samples_ht=samples_ht,
            split=False,
            filter_intervals=intervals,
            split_reference_blocks=False,
        )

    if exclude_ht is not None:
        # Safe on the unsplit VDS and without remove_dead_alleles; note this
        # prunes rows left with no carriers among the kept samples, which is
        # what we want for a phaser's input.
        vds = hl.vds.filter_samples(vds, exclude_ht, keep=False)

    # Deliberately NOT densify_encode_input_mt: that function exists to apply
    # the v4 *exome* corrections and joins project_meta.fixed_homalt_model,
    # which genomes metadata does not carry. Its other annotations (release af,
    # sex karyotype) are not needed to export genotypes either, and the
    # sex-ploidy adjustment is a no-op on an autosome. What we do need from it
    # is the phase-preserving split, so use that directly.
    variant_mt = _split_variant_data_keeping_phase(
        vds.variant_data, variant_ht, ["GT", "GQ", "DP", "AD", "PGT", "PID"]
    )
    vds = hl.vds.VariantDataset(vds.reference_data, variant_mt)
    mt = hl.vds.to_dense_mt(vds)
    logger.info("Dense MT: %d samples", mt.count_cols())

    if args.n_partitions:
        mt = mt.naive_coalesce(args.n_partitions)

    # The checkpoint keeps every entry field, including PGT/PID: the trio-truth,
    # EM and read-backed arms all read it and need read-backed phase and the
    # quality fields that adj is derived from. Only the VCF drops down to GT.
    if args.checkpoint_mt:
        mt = mt.checkpoint(args.checkpoint_mt, overwrite=args.overwrite)
        logger.info("Checkpointed dense MT: %s", args.checkpoint_mt)

    if args.out_vcf:
        _export_vcf_from_mt(mt, args)
    logger.info("Finished in %.1f s", timeit.default_timer() - start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", required=True, help="e.g. chr21:44324141-46709983.")
    parser.add_argument("--data-type", default="genomes", choices=["exomes", "genomes"])
    parser.add_argument(
        "--out-vcf", help="Output VCF (.bgz); omit for a densify-only run."
    )
    parser.add_argument(
        "--cohort",
        default="high_quality",
        choices=["high_quality", "release"],
        help=(
            "Sample set. 'high_quality' is the gnomAD sample-QC flag and keeps "
            "related and unreleasable samples, which trio work requires because "
            "releases are pruned to unrelated individuals and contain no "
            "complete trios."
        ),
    )
    parser.add_argument(
        "--min-ac",
        type=int,
        default=2,
        help="Minimum AC; the default of 2 drops singletons a phaser cannot place.",
    )
    parser.add_argument(
        "--min-an-frac",
        type=float,
        default=0.9,
        help="Call-rate floor, as a fraction of the region's maximum AN.",
    )
    parser.add_argument(
        "--keep-indels", action="store_true", help="Keep indels as well as SNPs."
    )
    parser.add_argument(
        "--pass-only",
        action="store_true",
        help=(
            "Drop variants failing variant QC (PASS, or AC0-only). Without it "
            "the phasers are handed variants the scored pair list excludes."
        ),
    )
    parser.add_argument(
        "--adj-only",
        action="store_true",
        help=(
            "Mask non-adj genotypes to missing in the exported VCF, so the "
            "phasers see the same genotype-quality bar as the trio and EM arms."
        ),
    )
    parser.add_argument(
        "--exclude-trio-parents",
        action="store_true",
        help="Drop pedigree parents so probands stay leakage-free for trio truth.",
    )
    parser.add_argument(
        "--exclude-trio-probands",
        action="store_true",
        help=(
            "Drop probands as well, so they are absent from the phasing cohort "
            "entirely rather than merely excluded from the per-pair vote. This "
            "removes them from the reference panel other samples are phased "
            "against. With --exclude-trio-parents and no --target-n the cohort "
            "becomes exactly the EM cohort (120,279 v4 genomes); with "
            "--target-n it keeps that run's deterministic fill sample, minus "
            "the probands."
        ),
    )
    parser.add_argument("--pedigree-path", help=".fam supplying the pedigree.")
    parser.add_argument(
        "--target-n",
        type=int,
        help=(
            "Cap the cohort at this many samples, keeping every proband. "
            "Phasing runtime is ~linear in samples; a smaller cohort is a floor "
            "on accuracy, not a neutral choice."
        ),
    )
    parser.add_argument("--seed", type=int, default=42, help="Fill-sample seed.")
    parser.add_argument("--checkpoint-mt", help="Optional dense MatrixTable path.")
    parser.add_argument(
        "--from-checkpoint",
        help=(
            "Read this dense MatrixTable instead of rebuilding it, and go "
            "straight to the VCF export. Splits the run into two separately "
            "billable jobs -- filter+densify, then export -- so the cost of "
            "each can be attributed. Must name a path a previous "
            "--checkpoint-mt run wrote with the same filters."
        ),
    )
    parser.add_argument(
        "--n-partitions",
        type=int,
        help="Coalesce before export (fewer, larger shards).",
    )
    parser.add_argument(
        "--tabix", action="store_true", help="Also write a tabix index."
    )
    parser.add_argument("--tmp-dir", default=DEFAULT_TMP_DIR)
    parser.add_argument("--overwrite", action="store_true")
    parsed = parser.parse_args()
    if parsed.from_checkpoint:
        if not parsed.out_vcf:
            parser.error("--from-checkpoint is export-only; pass --out-vcf.")
        if parsed.checkpoint_mt:
            parser.error("--from-checkpoint reads an MT; --checkpoint-mt writes one.")
    elif not (parsed.out_vcf or parsed.checkpoint_mt):
        parser.error("Nothing to produce: pass --out-vcf and/or --checkpoint-mt.")
    main(parsed)

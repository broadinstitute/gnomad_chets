"""
Compute per-variant-pair genotype counts from the gnomAD v4 VariantDataset.

The upstream sites / variant-filter / filtered-VMT / variant-pair-list steps
(and their functions) now live in create_vp_list.py; this script consumes the
variant-pair list and produces per-pair genotype-count arrays.

Genotype-count steps:

1. Encoded genotypes (--encode-genotypes): densify only the pair-list variants
   out of the gnomAD v4 VDS and encode them into per-variant sample sets,
   reusable across the count steps. The dense MT is transient (checkpointed to
   scratch, not persisted). Read-backed phase (PGT + PID) is kept through the
   split via _split_variant_data_keeping_phase and stored as a per-variant
   phased-het sidecar, so downstream can recover the same physical-phase signal
   as the gnomAD MNV pipeline. The v4 high-AB het -> hom-alt correction (GATK
   <4.1.4.1 artifact) is applied to the adj call, matching gnomad_qc
   generate_freq (joins the release freq HT for per-variant AF and the meta HT
   for per-sample fixed_homalt_model).

2. Variant size-info HT (--build-variant-size-info): per-variant contribution
   info used to split pairs into light vs heavy count jobs.

3. Genotype counts (--compute-counts-light / --compute-counts-heavy /
   --combine-counts): per-pair genotype-count arrays (raw and adj), optionally
   stratified by genetic-ancestry group (--stratify-by-pop / --pops) and — for
   the full-cohort path (--emit-phase-counts, default on) — with n_phased_cis /
   n_phased_trans refining the double-het (AaBb) cell from physical phase.

Use --backend batch to run on Hail Query-on-Batch instead of Spark (local or
Dataproc). Requires hailctl auth login and hailctl config set batch/remote_tmpdir,
batch/billing_project, and query/backend batch (or pass --backend batch).
See https://hail.is/docs/0.2/cloud/query_on_batch.html.
"""

import argparse
import logging
import os
import tempfile
import timeit
from typing import Optional, Union

import hail as hl
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import get_adj_expr
from gnomad_qc.v4.resources.annotations import get_freq
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.sample_qc import pedigree, trios

from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_TMP_DIR,
    GLOBAL_POP,
    TEST_INTERVALS,
    _get_output_postfix,
    get_count_subset_vds_path,
    get_count_subsets,
    get_pops,
    get_sample_pop_ht,
    get_variant_filter_ht,
    get_variant_pair_genotype_counts_ht,
    get_variant_pair_resources,
)
from gnomad_chets.v4.size_info_report import build_report
from gnomad_chets.v4.utils import (
    calculate_partitions_by_size,
    complete_trio_samples,
    samples_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("compute_vp_counts")
logger.setLevel(logging.INFO)

# High-AB het -> hom-alt correction (GATK <4.1.4.1 artifact), mirroring
# gnomad_qc.v4.annotations.generate_freq: an adj het-ref call with allele
# balance above HIGH_AB_CUTOFF, that is NOT a true het-non-ref, on an
# unfixed-model sample, at a variant with adj AF above HIGH_AB_AF_THRESHOLD,
# is really a hom-alt. Applied to the adj call only (raw is never adjusted).
HIGH_AB_CUTOFF = 0.9
HIGH_AB_AF_THRESHOLD = 0.01


def filter_pairs_by_an_pct(ht: hl.Table, min_an_pct: int) -> hl.Table:
    """Drop pairs whose AN_percent is at or below ``min_an_pct`` on either side.

    With ~no callable samples at a locus (``an_pct == 0``) nearly every
    sample is no-entry, so the AABB (hom-ref/hom-ref) cell collapses to zero
    and the haplotype EM degenerates — those pairs carry no co-occurrence
    signal; higher floors additionally trade power for AN quality. Applied at
    pair-list
    consumption (``--min-an-pct``) rather than baked into the build, so the
    list stays the complete raw artifact and the floor can change without a
    rebuild. A negative ``min_an_pct`` is a no-op (keeps every pair).

    :param ht: Variant pair Table with ``an_pct1`` / ``an_pct2`` fields.
    :param min_an_pct: Exclusive AN_percent floor required on both sides.
    :return: ``ht`` with sub-floor pairs removed.
    """
    if min_an_pct < 0:
        return ht
    return ht.filter((ht.an_pct1 > min_an_pct) & (ht.an_pct2 > min_an_pct))


def _read_min_an_pct(t: Union[hl.Table, hl.MatrixTable]) -> int:
    """Return the ``min_an_pct`` floor stamped on ``t``'s globals, or -1.

    Artifacts written before AN_pct flooring existed carry no such global;
    -1 makes :func:`_assert_min_an_pct_not_lowered` permissive for them.
    """
    g = t.index_globals()
    return hl.eval(g.min_an_pct) if "min_an_pct" in g.dtype else -1


def _assert_min_an_pct_not_lowered(
    current: int, upstream: int, upstream_name: str
) -> None:
    """Fail if the ``current`` floor is below the ``upstream`` artifact's floor.

    Downstream steps may only raise the AN_percent floor: a lower floor would
    count pairs whose variants were never densified / encoded upstream,
    silently corrupting the genotype counts.
    """
    if current < upstream:
        raise ValueError(
            f"--min-an-pct ({current}) is below the floor baked into the "
            f"{upstream_name} ({upstream}). Downstream steps may only raise "
            f"the AN_percent floor, never lower it. Rebuild the {upstream_name} "
            f"with a floor <= {current}, or pass --min-an-pct >= {upstream}."
        )


def create_variant_pair_filter_ht(vp_ht: hl.Table) -> hl.Table:
    """
    Create a filter Table for variant pairs (unique variants appearing in any pair).

    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2, alleles2.
    :return: Filter Table keyed by locus, alleles.
    """
    v1_ht = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1).select().distinct()
    v2_ht = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2).select().distinct()
    n_partitions = vp_ht.n_partitions()
    ht = (
        v1_ht.union(v2_ht)
        .distinct()
        .repartition(n_partitions, shuffle=True)
        .checkpoint(hl.utils.new_temp_file("encode_genotypes.variants", "ht"))
    )
    return ht


# Local-entry names that must be remapped to their local (``L*``) form before
# a sparse split. gnomad_qc's loader only remaps ``{GT, AD, PL}`` — omitting
# ``PGT`` silently drops phase — so :func:`_split_variant_data_keeping_phase`
# adds ``PGT`` (matching the gnomad_mnv pipeline) to carry read-backed phase
# (``LPGT`` → ``PGT``) + the phase-set id (``PID``) through the split.
_LOCAL_ENTRY_REMAP = {"GT", "AD", "PL", "PGT"}


def _split_variant_data_keeping_phase(
    variant_mt: hl.MatrixTable,
    filter_variant_ht: Optional[hl.Table],
    entries_to_keep: list,
) -> hl.MatrixTable:
    """Split a VDS variant-data MT, preserving phased-GT (``PGT``) + ``PID``.

    Mirrors ``gnomad_qc``'s ``_split_and_filter_variant_data_for_loading`` but
    adds ``PGT`` to the local-entry remap so ``hl.experimental.sparse_split_multi``
    downcodes ``LPGT`` → ``PGT`` instead of dropping it (the gnomad_qc loader
    only remaps ``{GT, AD, PL}``). ``PID`` (phase-set id) is a shared FORMAT
    field and passes through unchanged. This is the read pattern the gnomAD MNV
    pipeline uses to recover physical phase.

    The variant restriction is applied here (locus pre-filter before the split,
    full ``(locus, alleles)`` semi-join after) because ``get_gnomad_v4_vds``
    rejects ``filter_variant_ht`` on unsplit reads — callers pass ``split=False``
    and hand the raw variant data to this helper.

    :param variant_mt: Unsplit VDS variant-data MT (``vds.variant_data``).
    :param filter_variant_ht: Optional ``(locus, alleles)``-keyed Table to
        restrict to. Applied as a locus pre-filter + post-split semi-join.
    :param entries_to_keep: Post-split global entry names to keep (e.g.
        ``["GT", "GQ", "DP", "AD", "PGT", "PID"]``).
    :return: Split variant-data MT keyed by ``(locus, alleles)``, with an added
        ``_het_non_ref`` entry flag (from the local GT before the split
        downcodes it) for the high-AB het correction.
    """
    split_entries = [
        "L" + e if e in _LOCAL_ENTRY_REMAP else e
        for e in (entries_to_keep + ["LA"])
    ]
    variant_mt = variant_mt.select_entries(*split_entries)
    # Capture het-non-ref (e.g. 1/2) from the LOCAL GT *before* the split
    # downcodes each alt to its own het-ref record; carried through
    # sparse_split_multi as a plain passthrough entry so the encoder can exempt
    # true het-non-ref calls from the high-AB het -> hom-alt correction.
    variant_mt = variant_mt.annotate_entries(
        _het_non_ref=variant_mt.LGT.is_het_non_ref()
    )
    if filter_variant_ht is not None:
        # Locus-only pre-filter before the split (cheaper than splitting the
        # whole interval). filter_variant_ht is (locus, alleles)-keyed hence
        # locus-sorted, so re-key to locus without a shuffle (gnomad_qc idiom).
        filter_locus_ht = hl.Table(
            hl.ir.TableKeyBy(filter_variant_ht._tir, ["locus"], is_sorted=True)
        )
        variant_mt = variant_mt.filter_rows(
            hl.is_defined(filter_locus_ht[variant_mt.locus])
        )
    variant_mt = hl.experimental.sparse_split_multi(
        variant_mt, filter_changed_loci=True
    )
    if filter_variant_ht is not None:
        variant_mt = variant_mt.semi_join_rows(filter_variant_ht)
    return variant_mt


def _intervals_span_sex_chromosomes(filter_intervals) -> bool:
    """True if ``filter_intervals`` is None (genome-wide) or references chrX/chrY.

    Used to skip the sex-ploidy adjustment on autosomal(-only) runs: it is a
    strict no-op there, but ``adjusted_sex_ploidy_expr``'s index optimisation
    broadcasts the full column table into the per-entry expression, bloating the
    dense-MT checkpoint write. Handles both interval strings (``--gene`` /
    ``--test-genes``, from ``TEST_INTERVALS``) and parsed interval expressions
    (``--test-chrom``).
    """
    if filter_intervals is None:
        return True
    sex_contigs = {"chrX", "chrY", "X", "Y"}
    for iv in filter_intervals:
        if isinstance(iv, str):
            if iv.split(":")[0] in sex_contigs:
                return True
        elif hl.eval(iv.start.contig) in sex_contigs:
            return True
    return False


def _read_count_subset_vds(
    subset: str,
    *,
    filter_intervals,
) -> hl.vds.VariantDataset:
    """Read a by-sample count-subset VDS by path, interval- and chr19-filtered.

    Returns the subset with **all** of its cohort's samples — the release /
    high-quality restriction is deliberately NOT applied here. ``hl.vds.filter_
    samples`` on the sparse VDS prunes any variant row that has zero entries
    among the kept samples, so restricting samples on the VDS would drop the
    row for every pair-list variant that has no carrier within this subset and
    lose its hom-ref baseline (that is exactly the bug that made a per-subset
    count undercount AABB). ``to_dense_mt`` instead PRESERVES those rows and
    fills them as hom-ref from the reference blocks, so the cohort restriction
    is applied AFTER densify as a dense-MT column filter (see
    :func:`densify_encode_input_mt`'s ``restrict_samples_ht``). The subsets are
    splits of the same v4.0 exomes VDS ``get_gnomad_v4_vds`` reads and preserve
    all variant rows (see ``analysis/ukb_vds_split_runs.md``), so every
    pair-list variant is a row here even when monomorphic within the subset.
    The excessively-multi-allelic chr19 site is dropped to match the loader.

    :param subset: Subset name (see :func:`resources.get_count_subsets`).
    :param filter_intervals: Interval restriction (``None`` for a full run).
    :return: Interval-restricted VDS (all cohort samples) ready for
        :func:`densify_encode_input_mt` (which applies the cohort filter post-
        densify via ``restrict_samples_ht``).
    """
    vds = hl.vds.read_vds(get_count_subset_vds_path(subset))
    if filter_intervals is not None:
        ivs = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            if isinstance(x, str)
            else x
            for x in filter_intervals
        ]
        vds = hl.vds.filter_intervals(vds, ivs, split_reference_blocks=False)
    return hl.vds.filter_intervals(
        vds,
        [hl.parse_locus_interval("chr19:5787204-5787205", reference_genome="GRCh38")],
        keep=False,
    )


def densify_encode_input_mt(
    get_vds_func,
    filter_variant_ht: hl.Table,
    data_type: str,
    *,
    release_only: bool,
    filter_intervals,
    exclude_samples_ht: Optional[hl.Table] = None,
    vds: Optional[hl.vds.VariantDataset] = None,
    restrict_samples_ht: Optional[hl.Table] = None,
) -> hl.MatrixTable:
    """Densify the (filtered) pair-list variants into the MT the encoder wants.

    The single source of truth for building the encode input, shared by
    ``--encode-genotypes`` and ``trio_phasing.py``'s gnomAD-counts path so the
    two can't drift (they did once: the trio path skipped the high-AB
    correction + phase). Reads the VDS unsplit (to keep local ``LPGT``), splits
    via :func:`_split_variant_data_keeping_phase` to preserve read-backed phase
    (``PGT`` + ``PID``), then annotates the fields the v4 corrections need:
    per-variant ``af`` (release ``get_freq().freq[0].AF``, the same source as
    create_vp_list) and per-sample ``fixed_homalt_model`` (``meta().project_meta``)
    for the high-AB het → hom-alt correction, plus per-sample ``sex_karyotype``
    (``meta().sex_imputation``) for the sex-ploidy adjustment. Finally applies
    ``adjusted_sex_ploidy_expr`` to ``GT`` (a no-op on autosomes), matching
    gnomad_qc generate_freq so per-chromosome sex-chr runs are handled the same
    way the release freq was computed. ``get_gnomad_v4_vds`` rejects
    ``filter_variant_ht`` on unsplit reads, so the variant restriction happens
    inside the split helper.

    :param get_vds_func: ``get_gnomad_v4_vds`` / ``get_gnomad_v4_genomes_vds``.
    :param filter_variant_ht: ``(locus, alleles)`` filter (pair-list variants).
    :param data_type: ``exomes`` / ``genomes`` (for the freq + meta joins).
    :param release_only: passed to the VDS loader (``high_quality_only`` is its
        negation, matching the count-cohort convention).
    :param filter_intervals: interval restriction (``None`` for the full run).
    :param exclude_samples_ht: samples to REMOVE before densifying (e.g. the
        trio path drops PBT members for the gnomAD-minus-PBT counts).
    :param vds: Optional pre-read, already interval-restricted unsplit VDS to
        use instead of calling ``get_vds_func`` (the ``--vds-subset`` path
        passes a :func:`_read_count_subset_vds` result, which holds ALL of the
        subset's cohort samples). When given, ``release_only`` /
        ``filter_intervals`` are assumed already applied to it.
    :param restrict_samples_ht: Optional ``s``-keyed Table to restrict the
        cohort to AFTER densify (the ``--vds-subset`` path passes the release /
        high-quality set). Applied as a dense-MT COLUMN filter, which — unlike
        ``hl.vds.filter_samples`` on the sparse VDS — never prunes variant rows,
        so a pair-list variant monomorphic within the subset keeps its dense
        (hom-ref) row. Must NOT be used to restrict the sparse VDS upstream.
    :return: dense MatrixTable with ``GT/GQ/DP/AD/PGT/PID/_het_non_ref`` entries
        (``GT`` sex-ploidy-adjusted), per-variant ``af`` row field, and per-sample
        ``fixed_homalt_model`` + ``sex_karyotype`` cols.
    """
    if vds is None:
        vds = get_vds_func(
            release_only=release_only,
            high_quality_only=not release_only,
            split=False,
            filter_intervals=filter_intervals,
            split_reference_blocks=False,
        )
    if exclude_samples_ht is not None:
        vds = hl.vds.filter_samples(vds, exclude_samples_ht, keep=False)
    variant_mt = _split_variant_data_keeping_phase(
        vds.variant_data, filter_variant_ht, ["GT", "GQ", "DP", "AD", "PGT", "PID"],
    )
    vds = hl.vds.VariantDataset(vds.reference_data, variant_mt)
    mt = hl.vds.to_dense_mt(vds)
    # Restrict to the release / high-quality cohort AFTER densify (a column
    # filter): to_dense_mt has already materialized every variant row (incl.
    # ones monomorphic within a subset, filled hom-ref from reference blocks),
    # and a dense-MT column filter keeps all rows — so the hom-ref baseline of a
    # subset-monomorphic pair-list variant survives. (Restricting the sparse VDS
    # upstream would prune those zero-entry rows.)
    if restrict_samples_ht is not None:
        mt = mt.filter_cols(hl.is_defined(restrict_samples_ht[mt.s]))
    # Project the freq + meta HTs to ONLY the fields we join in, BEFORE the join.
    # gnomad_qc's meta HT is huge (project_meta / sample_qc / population_inference
    # / sex_imputation / …); accessing nested fields via `meta_ht[mt.s].<struct>.<f>`
    # broadcasts the WHOLE row, and the full-meta broadcast lands in the dense-MT
    # checkpoint's write task — ~835 MB, over spark.rpc.message.maxSize. Selecting
    # the two needed fields first shrinks the broadcast to a couple of columns.
    # Same reasoning for the freq HT's large per-row `freq` array → keep only AF.
    freq_ht = get_freq(data_type=data_type).ht()
    freq_ht = freq_ht.select(_af=freq_ht.freq[0].AF)
    meta_ht = meta(data_type=data_type).ht()
    meta_ht = meta_ht.select(
        _fixed_homalt_model=meta_ht.project_meta.fixed_homalt_model,
        _sex_karyotype=meta_ht.sex_imputation.sex_karyotype,
    )
    mt = mt.annotate_rows(af=freq_ht[mt.locus, mt.alleles]._af)
    mt = mt.annotate_cols(
        fixed_homalt_model=meta_ht[mt.s]._fixed_homalt_model,
        sex_karyotype=meta_ht[mt.s]._sex_karyotype,
    )
    # Sex-ploidy adjustment, matching gnomad_qc generate_freq
    # densify_and_prep_vds_for_freq: on non-PAR X/Y, XY calls become haploid and
    # XY hets are dropped to missing; XX calls on Y become missing. Applied
    # BEFORE the encoder computes adj + the high-AB correction, so adj is
    # computed on the sex-adjusted GT and the correction can't fire on a
    # now-missing chrX XY het. On chrX/Y: a hemizygous alt is haploid hom-var →
    # the `aa`/`bb` hom cell; a dropped het is uncallable → excluded from AABB.
    #
    # SKIPPED on autosomal(-only) runs: it is a strict no-op there, but
    # adjusted_sex_ploidy_expr's index optimisation
    # (annotate_and_index_source_mt_for_sex_ploidy) broadcasts the source MT's
    # WHOLE column table into the per-entry expression — a ~835 MB task in the
    # dense-MT checkpoint write that blows past spark.rpc.message.maxSize even on
    # autosomes. So only pay that cost when the data actually spans chrX/chrY.
    if _intervals_span_sex_chromosomes(filter_intervals):
        mt = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype)
        )
    return mt


def _create_var_idx_ht(mt: hl.MatrixTable) -> hl.Table:
    """
    Assign a unique int64 ``var_idx`` to each variant in ``mt``.

    Integer keys are dramatically faster for joins than ``(locus, alleles)`` struct
    keys and let us co-partition the encoded genotype Table and variant pair Table
    on the same int64 space.

    :param mt: MatrixTable whose row key is ``(locus, alleles)``.
    :return: Table keyed by ``(locus, alleles)`` with a ``var_idx`` field.
    """
    return mt.rows().add_index("var_idx").select("var_idx")


def _materialize_pop_map(pop_ht: hl.Table) -> dict:
    """Driver-side ``{s: pop}`` dict from a ``(s → pop)`` table.

    ``pop_ht`` is keyed by ``s`` with a ``pop`` field (see
    :func:`resources.get_sample_pop_ht`). Collected once per count step so the
    population label can be attached to the encoded ``samples`` global (which
    already carries ``s`` per index) at count time — no need to bake ``pop``
    into the encoding, so per-pop counts run over an EXISTING encoded table.
    """
    return {r.s: r.pop for r in pop_ht.select("pop").collect()}


def _build_pop_stratification(samples, pop_map=None, requested_pops=None):
    """Per-pop sample-index sets + sizes from the driver-side ``samples`` global.

    ``samples`` is the evaluated value of the encoded / gt_info ``samples``
    global — a list of structs carrying ``s``, in index order matching the
    encoded sets. The pop label per sample comes from ``pop_map`` (``{s: pop}``
    from the meta HT) when given, else from a ``pop`` field on each sample
    struct. Groups the sample *indices* (0-based) by pop.

    ``requested_pops`` (list of group names, excluding ``GLOBAL_POP``) restricts
    the specific strata to that subset — only requested groups that are actually
    present in the cohort are emitted, in the requested order. When ``None``,
    every present group is used.

    Returns ``(pops, pop_index_sets, pop_sizes)`` where ``pops`` is
    ``[GLOBAL_POP] + specific groups``, ``pop_index_sets`` maps each specific pop
    to an ``hl.literal`` set of its sample indices, and ``pop_sizes`` maps each
    specific pop to its sample count. ``GLOBAL_POP`` ("all") is intentionally
    absent from the two dicts — it reuses the flat full-cohort counts rather
    than an index-set intersection.
    """
    if pop_map is not None:
        def _pop_of(smp):
            return pop_map.get(smp.s)
    elif len(samples) and "pop" in samples[0]:
        def _pop_of(smp):
            return smp.pop
    else:
        raise ValueError(
            "Per-pop stratification needs a pop_map (from the meta HT) or a "
            "`pop` field on the encoded `samples` global."
        )
    idx_by_pop = {}
    for i, smp in enumerate(samples):
        p = _pop_of(smp)
        if p is not None:
            idx_by_pop.setdefault(p, []).append(i)
    if requested_pops is not None:
        # Keep only requested groups that are present, in the requested order.
        specific = [p for p in requested_pops if p in idx_by_pop]
    else:
        specific = sorted(idx_by_pop)
    pops = [GLOBAL_POP] + specific
    pop_index_sets = {
        p: hl.literal(set(idx_by_pop[p]), hl.tset(hl.tint32)) for p in specific
    }
    pop_sizes = {p: len(idx_by_pop[p]) for p in specific}
    return pops, pop_index_sets, pop_sizes


def _pop_stratification_for(encoded_gt_ht: hl.Table, pop_ht: hl.Table, requested_pops=None):
    """Per-pop index sets for an encoded table, using the meta pop HT.

    Evaluates the encoded ``samples`` global (``s`` per index) and joins it to
    ``{s: pop}`` from ``pop_ht`` driver-side — so per-pop counts work over any
    existing encoding without re-encoding to bake in ``pop``. ``requested_pops``
    optionally restricts the specific strata (see :func:`_build_pop_stratification`).
    """
    samples = hl.eval(encoded_gt_ht.index_globals().samples)
    return _build_pop_stratification(
        samples, _materialize_pop_map(pop_ht), requested_pops
    )


def _subtract_pbt_counts(full_ht: hl.Table, pbt_ht: hl.Table) -> hl.Table:
    """Annotate ``full_ht`` with ``gt_counts_{raw,adj}_no_pbt = full − PBT∩cohort``.

    ``cohort \\ PBT = cohort − (PBT∩cohort)``, element-wise on the 9-cell arrays —
    exact for all 9 cells incl. AABB (disjoint-cohort additivity). ``full_ht`` is
    the full-cohort counts; ``pbt_ht`` is the counts over the PBT∩cohort sample
    restriction (from :func:`restrict_encoded_to_samples` + a normal count),
    keyed by the same pair key. Both come from the SAME encode (``pbt_ht`` is a
    re-indexed restriction of it), so a sample lands in the same genotype cell in
    both — which is what makes the subtraction valid.
    """
    p = pbt_ht[full_ht.key]
    return full_ht.annotate(
        gt_counts_raw_no_pbt=hl.zip(full_ht.gt_counts_raw, p.gt_counts_raw).map(
            lambda t: t[0] - t[1]
        ),
        gt_counts_adj_no_pbt=hl.zip(full_ht.gt_counts_adj, p.gt_counts_adj).map(
            lambda t: t[0] - t[1]
        ),
    )


def merge_subset_counts(hts: list) -> hl.Table:
    """Sum per-subset genotype-count HTs into full-cohort counts.

    The ``--vds-subset`` count HTs are computed over a disjoint by-sample
    partition of the cohort (``non_ukb`` + ``ukb.<group>``) that keeps every
    variant row, so every 9-cell count — AABB included — is additive across
    subsets (disjoint-cohort additivity, the same property the PBT subtraction
    relies on): merging is a plain element-wise sum on the shared
    ``(locus1, alleles1, locus2, alleles2)`` key. ``gt_counts_raw`` /
    ``gt_counts_adj`` are summed coordinate-wise; ``n_phased_cis`` /
    ``n_phased_trans`` are summed when every input carries them. Sums stay in
    ``int32`` to match the per-subset schema (full-cohort AABB ≈ n_samples fits).

    A full-outer union + group-by (rather than an inner join) means a pair
    present in only some subsets still sums correctly, but with the intended
    all-rows-kept subsets every pair appears in every subset.

    :param hts: Per-subset genotype-count Tables (same schema).
    :return: Full-cohort counts Table keyed by the pair 4-tuple.
    """
    common = [
        "locus1", "alleles1", "locus2", "alleles2",
        "gt_counts_raw", "gt_counts_adj",
    ]
    has_phase = all(
        "n_phased_cis" in t.row and "n_phased_trans" in t.row for t in hts
    )
    if has_phase:
        common += ["n_phased_cis", "n_phased_trans"]
    normed = [
        t.key_by().select(*common).key_by(
            "locus1", "alleles1", "locus2", "alleles2"
        )
        for t in hts
    ]
    unioned = normed[0]
    for t in normed[1:]:
        unioned = unioned.union(t)
    # Cache before the group_by shuffle (Spark-backend shuffle stability).
    unioned = unioned.cache()
    agg = dict(
        gt_counts_raw=hl.agg.array_sum(unioned.gt_counts_raw).map(hl.int32),
        gt_counts_adj=hl.agg.array_sum(unioned.gt_counts_adj).map(hl.int32),
    )
    if has_phase:
        agg["n_phased_cis"] = hl.int32(hl.agg.sum(unioned.n_phased_cis))
        agg["n_phased_trans"] = hl.int32(hl.agg.sum(unioned.n_phased_trans))
    return unioned.group_by(
        "locus1", "alleles1", "locus2", "alleles2"
    ).aggregate(**agg)


_ENCODED_SET_FIELDS = (
    "all_samples", "raw_hr_adj_missing",
    "raw_het", "raw_hv", "adj_het", "adj_hv",
)


def restrict_encoded_to_pops(
    encoded_gt_ht: hl.Table, pop_ht: hl.Table, requested_pops: list
) -> hl.Table:
    """Restrict an encoded GT table to the union of ``requested_pops``' samples.

    Applied ONCE per variant (not per pair): each per-variant sample-index set
    is intersected with the kept-sample set and **re-indexed** into a dense
    ``0..n_keep-1`` space, and the size fields + ``samples`` global are updated
    to the kept cohort. The result is a self-consistent encoding of only the
    requested-pop samples, so the light/heavy count shuffle moves the small
    per-pop sets instead of the full-cohort sets (which otherwise blows worker
    shuffle disk). All sets are positive form, so restriction is a plain
    intersect + re-index and the sizes are just the re-indexed lengths. The
    ``phased_het`` sidecar is re-indexed the same way (kept present + consistent
    so the downstream projection doesn't break; per-pop phase counts aren't
    emitted yet but would be correct if wired up).

    Pure transform (aside from the driver-side samples/meta materialization),
    intended to be checkpointed by the caller.
    """
    samples = hl.eval(encoded_gt_ht.index_globals().samples)
    pop_map = _materialize_pop_map(pop_ht)
    req = set(requested_pops)
    keep = [i for i, smp in enumerate(samples) if pop_map.get(smp.s) in req]
    if not keep:
        raise ValueError(
            f"No samples in the encoded cohort for requested pops {sorted(req)}."
        )
    logger.info(
        "restrict_encoded_to_pops: keeping %d of %d samples for pops %s",
        len(keep), len(samples), sorted(req),
    )
    return _restrict_encoded_to_indices(encoded_gt_ht, keep, samples)


def _restrict_encoded_to_indices(encoded_gt_ht: hl.Table, keep: list, samples) -> hl.Table:
    """Re-index an encoded GT table's sets to ``keep`` (orig sample indices) into
    a dense ``0..len(keep)-1`` space, updating the size fields + ``samples``
    global. Shared by :func:`restrict_encoded_to_pops` and
    :func:`restrict_encoded_to_samples`.

    The keep-set / remap **literals are used HERE, once per variant in this
    (caller-checkpointed) transform** — never in the per-pair count expression.
    That is the whole point: a large keep-set literal is fine as a one-shot
    broadcast over the encoded rows, but replicating it into every light/heavy
    count task blows past ``spark.rpc.message.maxSize`` (the no-PBT footgun).
    """
    remap = {orig: dense for dense, orig in enumerate(keep)}
    keep_samples = [samples[i] for i in keep]
    samples_dtype = encoded_gt_ht.index_globals().samples.dtype
    keep_set = hl.literal(set(keep), hl.tset(hl.tint32))
    remap_lit = hl.literal(remap, hl.tdict(hl.tint32, hl.tint32))

    def _reidx(s):
        return hl.set(
            s.filter(lambda i: keep_set.contains(i)).map(lambda i: remap_lit[i])
        )

    e = encoded_gt_ht
    reidx = {f: _reidx(e[f]) for f in _ENCODED_SET_FIELDS}

    # Positive within-kept sizes from the re-indexed stored sets.
    n_f = reidx["raw_hr_adj_missing"].length()
    annotations = dict(
        **reidx,
        n_with_data=hl.int32(reidx["all_samples"].length() + n_f),
        n_raw_hr_adj_missing=hl.int32(n_f),
    )
    # phased_het holds sample indices too, so re-index it into the same dense
    # 0..n_keep-1 space (drop kept-out keys, remap the rest) — keeping the
    # field present and consistent so the downstream projection / count path
    # doesn't break, and any future per-pop phase count stays correct.
    if "phased_het" in e.row:
        annotations["phased_het"] = hl.dict(
            e.phased_het.items()
            .filter(lambda kv: keep_set.contains(kv[0]))
            .map(lambda kv: (remap_lit[kv[0]], kv[1]))
        )
    return e.annotate(**annotations).annotate_globals(
        samples=hl.literal(keep_samples, samples_dtype)
    )


def restrict_encoded_to_samples(encoded_gt_ht: hl.Table, keep_samples_ht: hl.Table):
    """Restrict an encoded GT table to the samples in ``keep_samples_ht`` (keyed
    by ``s``) and re-index into a dense space — the no-PBT primitive (restrict to
    ``PBT∩cohort``, count it, subtract from the full-cohort counts).

    Uses the same re-index-once mechanism as :func:`restrict_encoded_to_pops`, so
    the keep-set literal stays in this checkpointed transform and OUT of the
    per-pair count expression (the ``hl.literal`` that blew up when the no-PBT
    restriction was applied inline). Returns ``(restricted_encoded, n_keep)``.
    """
    keep_s = set(keep_samples_ht.s.collect())
    samples = hl.eval(encoded_gt_ht.index_globals().samples)
    keep = [i for i, smp in enumerate(samples) if smp.s in keep_s]
    if not keep:
        raise ValueError(
            "restrict_encoded_to_samples: no encoded samples matched "
            "keep_samples_ht (PBT∩cohort is empty)."
        )
    logger.info(
        "restrict_encoded_to_samples: keeping %d of %d samples.",
        len(keep), len(samples),
    )
    return _restrict_encoded_to_indices(encoded_gt_ht, keep, samples), len(keep)


def _encode_genotype_sets_by_var_idx(
    mt: hl.MatrixTable,
    var_idx_ht: hl.Table,
    *,
    use_precomputed_adj: bool = False,
) -> hl.Table:
    """
    Encode genotypes as per-variant sample sets keyed by ``var_idx``.

    For each variant, produces sets of sample indices by genotype category for
    both the raw and adj call, plus an ``all_samples`` set used to derive
    hom-ref counts by inclusion/exclusion downstream.

    Raw/adj genotypes are encoded as:

        - missing (None) = hom-ref (space saving)
        - 0 = missing data (no GT call, or failed adj for ``adj_gt``)
        - 1 = het
        - 2 = hom-var

    Each sample falls in exactly one of these 7 disjoint per-variant
    categories (the implicit "adj-PASS-0/0" majority is left untracked):

        cat | GT                  | adj | raw_gt | adj_gt | stored in
        ----+---------------------+-----+--------+--------+----------------------
         1  | no entry            |  -  |   NA   |   NA   | all_samples (only)
         2  | 0/0                 | FAIL|   NA   |    0   | raw_hr_adj_missing
                                                          |   (only)
         3  | 0/1                 | FAIL|    1   |    0   | all_samples, raw_het
         4  | 0/1                 | PASS|    1   |    1   | all_samples, raw_het,
                                                          |   adj_het
         5  | 1/1                 | FAIL|    2   |    0   | all_samples, raw_hv
         6  | 1/1                 | PASS|    2   |    2   | all_samples, raw_hv,
                                                          |   adj_hv
         7  | 0/0                 | PASS|   NA   |   NA   | (implicit majority)

    ``all_samples`` stores cats 1, 3-6 directly (positive form), disjoint
    from ``raw_hr_adj_missing`` (cat 2). cat 7 (the adj-PASS-0/0 majority)
    is never stored — :func:`_count_from_sets` reconstructs it as
    ``N − (cats 1-6)`` from ``n_samples``.

    The count algebra in :func:`_count_from_sets` reconstructs
    ``D'_v = A_v ∪ F_v`` (= cats 1-6) by summing the four A/F cross terms.
    Storing positive sets keeps the decode a plain intersection — no
    complement bookkeeping, so no proper-complement invariant to violate.
    At low-coverage variants the positive ``all_samples`` (large cat-1
    no-entry) can approach ``n_samples``; ``--min-an-pct`` bounds this by
    dropping the lowest-AN endpoints.

    The v4 high-AB het correction (see module constants) can set ``adj_gt=2``
    for a call whose ``raw_gt=1`` — i.e. a corrected call is in ``raw_het``
    AND ``adj_hv`` (not the cat-4 ``raw==adj==1`` combo). ``raw`` sets and
    ``adj`` sets are consumed independently by :func:`_count_from_sets`, so
    this needs no special handling; the cat-1-7 table below describes the
    UNcorrected mapping.

    Output schema:

        ----------------------------------------
        Global fields:
            'samples': array<struct { s: str }>
        ----------------------------------------
        Row fields:
            'v_idx': int64
            'all_samples': set<int32>      # cats 1, 3-6 (positive form),
                                           # disjoint from raw_hr_adj_missing.
            'n_with_data': int32           # |cats 1-6| = samples NOT in
                                           # adj-PASS-0/0 majority
                                           # (= |all_samples| + |raw_hr_adj_missing|)
            'raw_hr_adj_missing': set<int32>  # cat 2 (positive form)
            'n_raw_hr_adj_missing': int32  # |cat 2|
            'raw_het': set<int32>          # all het (cats 3 ∪ 4)
            'raw_hv': set<int32>           # all hv  (cats 5 ∪ 6)
            'adj_het': set<int32>          # adj-pass het (cat 4 only)
            'adj_hv': set<int32>           # adj-pass hv  (cat 6 only)
            'phased_het': dict<int32,      # het sample_idx → (pid, gt0) for
                struct{pid:str,gt0:int32}> #   read-backed-phased het calls.
                                           #   Two het variants of one sample
                                           #   are cis iff same pid & same gt0.
                                           #   Empty when the MT had no PGT/PID.
        ----------------------------------------
        Key: ['v_idx']
        ----------------------------------------

    :param mt: MatrixTable with variant data. Row key must be ``(locus, alleles)``.
    :param var_idx_ht: Var-idx Table from :func:`_create_var_idx_ht` keyed by
        ``(locus, alleles)`` with a ``var_idx`` field.
    :param use_precomputed_adj: When ``True``, read the adj-PASS flag from the
        ``mt.adj`` entry field instead of recomputing it from ``GT/GQ/DP/AD``.
        Lets callers encode an MT that kept only ``GT`` + ``adj`` (e.g. the
        exploded PBT trio MatrixTable) without re-densifying for GQ/DP/AD.
    :return: Table keyed by ``v_idx`` with per-variant sample sets.
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )
    adj_pass_expr = (
        mt.adj if use_precomputed_adj else get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD)
    )
    # v4 high-AB het -> hom-alt correction (GATK <4.1.4.1 artifact): reclassify
    # an adj het-ref call as hom-var (2) when AB > cutoff, it isn't a true
    # het-non-ref, the sample's model isn't fixed, and the variant's adj AF is
    # above threshold — mirroring gnomad_qc generate_freq. Skipped unless the
    # dense MT carries the required fields (af / fixed_homalt_model /
    # _het_non_ref), e.g. the trio PBT MT doesn't, so it is left uncorrected.
    #
    # The v4 release corrects the RAW call stats too, not just adj: the released
    # freq is `ab_adjusted_freq`, whose correction adds the (adj-determined)
    # high-AB hom-alt count to every stratum including freq[1] = raw
    # (gnomad_qc generate_freq.correct_for_high_ab_hets, and the raw-group
    # aggregate of the adj-gated high_ab_het). So we reclassify an adj-passing
    # high-AB het to hom-var in BOTH the raw and adj genotypes. A non-adj
    # high-AB het is NOT in that correction set, so it stays het in raw — hence
    # the raw reclassification is gated on adj-pass.
    correct_high_ab = (
        "af" in mt.row
        and "fixed_homalt_model" in mt.col
        and "_het_non_ref" in mt.entry
    )
    if correct_high_ab:
        high_ab_homalt = (
            mt.GT.is_het_ref()
            & (mt.AD[1] / mt.DP > HIGH_AB_CUTOFF)
            & ~hl.coalesce(mt._het_non_ref, False)
            & ~hl.coalesce(mt.fixed_homalt_model, False)
            & hl.coalesce(mt.af > HIGH_AB_AF_THRESHOLD, False)
        )
        adj_gt_expr = hl.if_else(high_ab_homalt, 2, gt_count_expr)
        raw_gt_expr = hl.if_else(
            high_ab_homalt & adj_pass_expr, 2, gt_count_expr, missing_false=True
        )
    else:
        adj_gt_expr = gt_count_expr
        raw_gt_expr = gt_count_expr
    adj_gt_count_expr = hl.if_else(
        adj_pass_expr, adj_gt_expr, 0, missing_false=True
    )

    # Physical (read-backed) phase, when the dense MT carries it (PGT + PID
    # kept via _split_variant_data_keeping_phase). For a phased het-ref call
    # we record the phase-set id (PID) and which haplotype carries the alt
    # (PGT[0]); at count time two het variants of the same sample are cis iff
    # they share a PID and agree on PGT[0]. Reference-block (hom-ref) samples
    # have no PGT/PID, so phase is missing for them. Absent on MTs without
    # phase (e.g. the exploded PBT trio MT) → an empty phased_het dict.
    has_phase = "PGT" in mt.entry and "PID" in mt.entry
    select_entry_exprs = dict(raw_gt=raw_gt_expr, adj_gt=adj_gt_count_expr)
    if has_phase:
        select_entry_exprs["phase"] = hl.if_else(
            mt.GT.is_het()
            & hl.is_defined(mt.PGT)
            & mt.PGT.phased
            & hl.is_defined(mt.PID),
            hl.struct(pid=mt.PID, gt0=mt.PGT[0]),
            hl.missing(hl.tstruct(pid=hl.tstr, gt0=hl.tint32)),
        )
    mt = mt.select_entries(**select_entry_exprs)
    ht = mt.localize_entries("_entries", "samples")

    # Build per-variant sample-index sets. The category boundary uses:
    #   - is_missing(x[1])               => cat 1 (no entry)
    #   - is_missing raw_gt & adj_gt==0  => cat 2 (raw_hr_adj_missing)
    #   - raw_gt == 1, adj_gt == 0       => cat 3
    #   - raw_gt == 1, adj_gt == 1       => cat 4
    #   - raw_gt == 2, adj_gt == 0       => cat 5
    #   - raw_gt == 2, adj_gt == 2       => cat 6
    #   - is_missing raw_gt & is_missing adj_gt & defined entry => cat 7
    #
    # all_samples (stored) = cats 1, 3-6 — disjoint from raw_hr_adj_missing
    # (cat 2) so the same sample never appears in both. cat 7 (adj-PASS-0/0
    # majority) is left untracked and reconstructed by the count kernel.
    gt = hl.enumerate(ht._entries)
    # ``not_adj_hom_ref`` = cats 1-6 (used only to derive n_with_data;
    # not stored). The set stored as ``all_samples`` is the cats-1,3-6
    # subset (``not_adj_hom_ref_no_F``), disjoint from ``raw_hr_adj_missing``
    # (cat 2) so the same sample never appears in both stored sets.
    not_adj_hom_ref = gt.filter(
        lambda x: hl.is_missing(x[1])
        | hl.is_defined(x[1].raw_gt)
        | hl.is_defined(x[1].adj_gt)
    )
    # cats 1, 3-6: keep cat 1 (no entry) and cats 3-6 (raw_gt is defined
    # only for het/hom-var, i.e. cats 3-6). Drops cat 2 (hom-ref +
    # adj-FAIL) which is already in ``raw_hr_adj_missing``.
    not_adj_hom_ref_no_F = gt.filter(
        lambda x: hl.is_missing(x[1]) | hl.is_defined(x[1].raw_gt)
    )
    raw_hr_adj_missing = gt.filter(
        lambda x: hl.is_defined(x[1])
        & hl.is_missing(x[1].raw_gt)
        & (x[1].adj_gt == 0)
    )
    n_with = not_adj_hom_ref.length()
    n_raw_hr_adj_missing = raw_hr_adj_missing.length()
    # Per-variant phased-het sidecar: sample_idx → (pid, gt0) for het samples
    # with a valid phased call. Rides alongside the het sets; unphased hets
    # are simply absent. Empty when the MT carried no phase (has_phase=False).
    if has_phase:
        phased_het_expr = hl.dict(
            gt.filter(
                lambda x: hl.is_defined(x[1]) & hl.is_defined(x[1].phase)
            ).map(lambda x: (x[0], x[1].phase))
        )
    else:
        phased_het_expr = hl.empty_dict(
            hl.tint32, hl.tstruct(pid=hl.tstr, gt0=hl.tint32)
        )
    ht = ht.select(
        # all_samples = cats 1, 3-6 (positive form), disjoint from
        # raw_hr_adj_missing (cat 2). n_with_data = |cats 1-6|; the count
        # kernel reconstructs cat 7 (= N − cats 1-6) from n_samples.
        all_samples=hl.set(not_adj_hom_ref_no_F.map(lambda x: x[0])),
        n_with_data=hl.int32(n_with),
        raw_hr_adj_missing=hl.set(raw_hr_adj_missing.map(lambda x: x[0])),
        n_raw_hr_adj_missing=hl.int32(n_raw_hr_adj_missing),
        raw_het=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].raw_gt == 1).map(lambda x: x[0])
        ),
        raw_hv=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].raw_gt == 2).map(lambda x: x[0])
        ),
        adj_het=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].adj_gt == 1).map(lambda x: x[0])
        ),
        adj_hv=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].adj_gt == 2).map(lambda x: x[0])
        ),
        phased_het=phased_het_expr,
    )

    # Rekey by var_idx and drop the locus/alleles fields to shrink row size.
    ht = ht.annotate(v_idx=var_idx_ht[ht.locus, ht.alleles].var_idx)
    return ht.key_by("v_idx").drop("locus", "alleles")



def _count_from_sets(
    v1_het: hl.expr.SetExpression,
    v1_hv: hl.expr.SetExpression,
    v1_all: hl.expr.SetExpression,
    v1_n: hl.expr.Int32Expression,
    v1_F: hl.expr.SetExpression,
    v1_n_F: hl.expr.Int32Expression,
    v2_het: hl.expr.SetExpression,
    v2_hv: hl.expr.SetExpression,
    v2_all: hl.expr.SetExpression,
    v2_n: hl.expr.Int32Expression,
    v2_F: hl.expr.SetExpression,
    v2_n_F: hl.expr.Int32Expression,
    n_samples: hl.expr.Int32Expression,
    *,
    include_raw_hr_adj_missing: bool,
) -> hl.expr.ArrayExpression:
    """
    Compute 9-element genotype count array from per-variant sample sets.

    Per-variant storage (all sets are positive form — sample indices held
    directly, no complement encoding):

      - ``v_all`` = ``A_v`` = cats 1, 3-6 (no-entry ∪ het ∪ hom-var),
        disjoint from ``v_F``. ``v_n`` = ``|cats 1-6| = |A_v| + |F_v|`` =
        n_with_data, the "samples with any data" count.
      - ``v_F`` = ``raw_hr_adj_missing`` = cat 2 (raw-0/0, adj-fail).
        ``v_n_F`` = |cat 2|.

    The "real" not-adj-hom-ref set is ``D'_v = A_v ∪ F_v = cats 1-6``;
    intersections over D' decompose into the A/F cross terms.

    Two modes (Python-level branch via ``include_raw_hr_adj_missing``):

      - ``False`` (ADJ cells): hom-ref = adj-PASS-0/0 = N \\ D'. Pass
        ``v_het = adj_het``, ``v_hv = adj_hv``.
      - ``True`` (RAW cells): hom-ref = raw-0/0 = (adj-PASS-0/0) ∪ F.
        Pass ``v_het = raw_het`` (all het), ``v_hv = raw_hv`` (all hv).
        F adds the adj-fail-hom-ref samples to the hom-ref pool.

    Count array layout: ``[AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]``
    where A/a = v1 ref/alt, B/b = v2 ref/alt.

    :param v1_het, v1_hv: Sample sets for v1 het / hom-var.
    :param v1_all, v1_n: v1 ``all_samples`` (= A_v) set + ``n_with_data``.
    :param v1_F, v1_n_F: v1 ``raw_hr_adj_missing`` (= cat 2) set + size.
    :param v2_...: same for v2.
    :param n_samples: Total number of samples in the cohort.
    :param include_raw_hr_adj_missing: ``True`` for raw cells, ``False``
        for adj cells.
    :return: 9-element genotype count array.
    """
    # D'_v = A_v ∪ F_v (disjoint per variant). Decompose |D'_v1 ∩ D'_v2|
    # into the four cross terms; reuse the A-F and F-F pieces in the
    # raw-cells branch below.
    a1_isect_a2 = v1_all.intersection(v2_all).length()
    a1_isect_f2 = v1_all.intersection(v2_F).length()
    f1_isect_a2 = v1_F.intersection(v2_all).length()
    f1_isect_f2 = v1_F.intersection(v2_F).length()
    d1_isect_d2 = a1_isect_a2 + a1_isect_f2 + f1_isect_a2 + f1_isect_f2

    # |H_v1 ∩ H_v2| where H = N \ D'  (adj-PASS-0/0 at both).
    #   = N - |D'_v1 ∪ D'_v2| = N - |D'_v1| - |D'_v2| + |D'_v1 ∩ D'_v2|
    h1_isect_h2 = n_samples - v1_n - v2_n + d1_isect_d2

    # Carrier-carrier cells.
    het_het = v1_het.intersection(v2_het).length()
    het_hv = v1_het.intersection(v2_hv).length()
    hv_het = v1_hv.intersection(v2_het).length()
    hv_hv = v1_hv.intersection(v2_hv).length()

    # Edge cells (adj component): |carrier_v ∩ H_other|
    #   = |carrier_v| - |carrier_v ∩ D'_other|
    # where |carrier_v ∩ D'_other| = |carrier_v ∩ A_other|
    #                              + |carrier_v ∩ F_other|.
    v1_het_in_d2 = (
        v1_het.intersection(v2_all).length()
        + v1_het.intersection(v2_F).length()
    )
    v1_hv_in_d2 = (
        v1_hv.intersection(v2_all).length()
        + v1_hv.intersection(v2_F).length()
    )
    v2_het_in_d1 = (
        v2_het.intersection(v1_all).length()
        + v2_het.intersection(v1_F).length()
    )
    v2_hv_in_d1 = (
        v2_hv.intersection(v1_all).length()
        + v2_hv.intersection(v1_F).length()
    )
    v1_het_in_h2 = v1_het.length() - v1_het_in_d2
    v1_hv_in_h2 = v1_hv.length() - v1_hv_in_d2
    v2_het_in_h1 = v2_het.length() - v2_het_in_d1
    v2_hv_in_h1 = v2_hv.length() - v2_hv_in_d1

    if include_raw_hr_adj_missing:
        # RAW cells: hom-ref pool extended by F (cat 2) at each variant.
        # AABB_raw = |(H1 ∪ F1) ∩ (H2 ∪ F2)| (H_v and F_v disjoint per
        # variant, so the union sums) =
        #   |H1∩H2| + |H1∩F2| + |F1∩H2| + |F1∩F2|.
        # Reuse the cross terms already computed above:
        #   |D'_v1 ∩ F_v2| = |A_v1 ∩ F_v2| + |F_v1 ∩ F_v2|
        #   |F_v1 ∩ D'_v2| = |F_v1 ∩ A_v2| + |F_v1 ∩ F_v2|
        h1_isect_f2 = v2_n_F - (a1_isect_f2 + f1_isect_f2)
        f1_isect_h2 = v1_n_F - (f1_isect_a2 + f1_isect_f2)
        hom_ref_both = h1_isect_h2 + h1_isect_f2 + f1_isect_h2 + f1_isect_f2

        # Edge cells (raw): also extend hom-ref-at-other by F_other.
        v1_het_in_f2 = v1_het.intersection(v2_F).length()
        v1_hv_in_f2 = v1_hv.intersection(v2_F).length()
        v2_het_in_f1 = v2_het.intersection(v1_F).length()
        v2_hv_in_f1 = v2_hv.intersection(v1_F).length()

        v1_het_homref_v2 = v1_het_in_h2 + v1_het_in_f2
        v1_hv_homref_v2 = v1_hv_in_h2 + v1_hv_in_f2
        v2_het_homref_v1 = v2_het_in_h1 + v2_het_in_f1
        v2_hv_homref_v1 = v2_hv_in_h1 + v2_hv_in_f1
    else:
        # ADJ cells: hom-ref = H only.
        hom_ref_both = h1_isect_h2
        v1_het_homref_v2 = v1_het_in_h2
        v1_hv_homref_v2 = v1_hv_in_h2
        v2_het_homref_v1 = v2_het_in_h1
        v2_hv_homref_v1 = v2_hv_in_h1

    return hl.array([
        hom_ref_both,                       # AABB
        v2_het_homref_v1,                   # AABb
        v2_hv_homref_v1,                    # AAbb
        v1_het_homref_v2,                   # AaBB
        het_het,                            # AaBb
        het_hv,                             # Aabb
        v1_hv_homref_v2,                    # aaBB
        hv_het,                             # aaBb
        hv_hv,                              # aabb
    ])


def _count_phase_from_sets(
    v1_het: hl.expr.SetExpression,
    v1_phased_het: hl.expr.DictExpression,
    v2_het: hl.expr.SetExpression,
    v2_phased_het: hl.expr.DictExpression,
) -> hl.expr.StructExpression:
    """Physical-phase refinement of the double-het (``AaBb``) cell.

    Among samples het at both variants (``v1_het ∩ v2_het``), classify those
    read-backed-phased on both sides *within the same phase set* (matching
    ``pid``) as **cis** (alt on the same haplotype — ``gt0`` agrees) or
    **trans** (opposite). This is the gnomAD MNV pipeline's same-haplotype
    call, applied per pair from the per-variant phase sidecar. Samples without
    shared-``pid`` phase are neither; they are the unphased remainder
    (``AaBb − cis − trans``) the EM still resolves statistically.

    Pass adj carrier sets so the counts refine the adj-quality ``AaBb`` cell.
    Missing-key dict lookups return missing, so an unphased het at either side
    fails the ``pid`` equality and is excluded — the ``contains`` guards make
    that explicit.

    :param v1_het, v2_het: Het sample-index sets (adj) for the two variants.
    :param v1_phased_het, v2_phased_het: ``sample_idx → struct{pid, gt0}``
        phase sidecars from the encoder.
    :return: ``struct(n_phased_cis, n_phased_trans)`` (int32).
    """
    both_phased = v1_het.intersection(v2_het).filter(
        lambda s: v1_phased_het.contains(s)
        & v2_phased_het.contains(s)
        & (v1_phased_het[s].pid == v2_phased_het[s].pid)
    )
    n_cis = both_phased.filter(
        lambda s: v1_phased_het[s].gt0 == v2_phased_het[s].gt0
    ).length()
    return hl.struct(
        n_phased_cis=hl.int32(n_cis),
        n_phased_trans=hl.int32(both_phased.length() - n_cis),
    )


def _pop_restrict_variant(v, het, hv, pop_set):
    """Restrict one variant's :func:`_count_from_sets` inputs to ``pop_set``.

    ``het`` / ``hv`` are the raw *or* adj carrier sets (passed explicitly so the
    same helper serves both counts). Returns the 6-tuple of pop-restricted args
    in the order :func:`_count_from_sets` consumes per variant:
    ``(het, hv, all, n_with_data, F, n_F)`` — every positive-form set
    intersected with the pop's sample-index set, sizes recomputed within-pop.
    """
    all_p = v.all_samples.intersection(pop_set)
    f_p = v.raw_hr_adj_missing.intersection(pop_set)
    return (
        het.intersection(pop_set),
        hv.intersection(pop_set),
        all_p,
        hl.int32(all_p.length() + f_p.length()),
        f_p,
        hl.int32(f_p.length()),
    )


def _count_from_sets_by_pop(v1, v2, pops, pop_index_sets, pop_sizes, flat_raw, flat_adj):
    """Per-population 9-cell counts as ``dict<pop, struct{raw, adj}>``.

    Reuses :func:`_count_from_sets` for each specific (non-``GLOBAL_POP``) group
    by restricting both variants' sample sets to that group's sample-index set
    (:func:`_pop_restrict_variant`) and passing ``|pop|`` as ``n_samples``. The
    ``GLOBAL_POP`` ("all") entry reuses the already-computed flat full-cohort
    arrays, guaranteeing ``by_pop["all"]`` equals the flat ``gt_counts_*``.

    :param v1, v2: Encoded per-variant structs (indexed from the encoded GT
        table) carrying the ``_COUNT_FROM_SETS_FIELDS``.
    :param pops: Pop list incl. ``GLOBAL_POP`` (from :func:`_build_pop_stratification`).
    :param pop_index_sets: ``dict pop -> hl.literal(set<int32>)`` for specific pops.
    :param pop_sizes: ``dict pop -> int`` for specific pops.
    :param flat_raw, flat_adj: Flat full-cohort 9-cell arrays (the ``"all"`` value).
    :return: ``hl.dict`` keyed by pop → ``struct(raw, adj)`` 9-cell arrays.
    """
    entries = [(GLOBAL_POP, hl.struct(raw=flat_raw, adj=flat_adj))]
    for p in pops:
        if p == GLOBAL_POP:
            continue
        ps = pop_index_sets[p]
        sz = hl.int32(pop_sizes[p])
        raw = _count_from_sets(
            *_pop_restrict_variant(v1, v1.raw_het, v1.raw_hv, ps),
            *_pop_restrict_variant(v2, v2.raw_het, v2.raw_hv, ps),
            sz,
            include_raw_hr_adj_missing=True,
        )
        adj = _count_from_sets(
            *_pop_restrict_variant(v1, v1.adj_het, v1.adj_hv, ps),
            *_pop_restrict_variant(v2, v2.adj_het, v2.adj_hv, ps),
            sz,
            include_raw_hr_adj_missing=False,
        )
        entries.append((p, hl.struct(raw=raw, adj=adj)))
    return hl.dict(entries)


def _compute_counts_for_subset(
    vp_subset: hl.Table,
    encoded_gt_ht: hl.Table,
    n_samples: hl.expr.Int32Expression,
    label: str,
    n_partitions: int,
    heavy_variants_with_contribution: hl.Table,
    *,
    pops=None,
    pop_index_sets=None,
    pop_sizes=None,
    emit_phase: bool = False,
) -> hl.Table:
    """
    Compute genotype counts for a subset of variant pairs using split-aware
    ``(v_idx, _split_idx)`` co-partitioning.

    Why the ``_split_idx`` axis? ``calculate_partitions_by_size`` can only
    place boundaries between rows — it cannot split a single row whose
    weight already exceeds the target. For DYSF the worst single variant's
    predicted shuffle contribution is ~14.6 GB (one variant alone, 4,091
    partners × 3.6 MB payload), and the target is 500 MB. Without
    splitting, that whole variant's data would land on one partition and
    blow worker disk.

    Instead, each heavy variant ``v`` is duplicated ``split_count(v)`` times
    in the encoded view (one row per ``(v_idx, s)`` for s in
    ``[0, split_count)``). Each pair touching that variant gets routed to
    one of those duplicates via ``hash(opposite_v_idx) % split_count``, so
    its ~``degree/split_count`` partners spread evenly across the duplicate
    partitions. ``split_count`` is chosen so each per-split contribution is
    at most one target, and the standard size-balanced partition algorithm
    then produces partitions of roughly equal weight even in the presence
    of one ultra-heavy variant. Light variants stay at ``split_count = 1``
    and behave identically to single-key partitioning.

    Two shuffles total: one group_by ``(v1_idx, v1_split_idx)`` (the v1
    zip-join's co-partition), one rekey by ``(v2_idx, v2_split_idx)`` (the
    v2 zip-join's co-partition). Both joins are partition-local against the
    split-aware encoded view.

    :param vp_subset: Pair Table with ``v1_idx, v2_idx`` (Design C carries only
        the integer keys; the locus/alleles 4-tuple is rebuilt by the caller).
    :param encoded_gt_ht: Encoded genotype sets keyed by v_idx.
    :param n_samples: Total sample count expression.
    :param label: Label for temp-file naming (e.g., "heavy").
    :param n_partitions: Target number of join partitions.
    :param heavy_variants_with_contribution: Table keyed by ``v_idx`` with
        ``contribution`` (bytes) and ``split_count`` (int32). Light variants
        not in this table get ``contribution = 0`` and ``split_count = 1``
        via ``hl.or_else`` defaults.
    :return: Counts Table with gt_counts_raw and gt_counts_adj.
    """
    # --- 1. Split-aware encoded view, keyed by (v_idx, _split_idx) ---
    # Annotate each encoded row with its split_count, then explode to
    # split_count copies (light variants → one copy). Per-split contribution
    # is the variant's total contribution divided by its split_count, so
    # every row's predicted weight is at most ~TARGET_HEAVY_PARTITION_BYTES.
    enc = encoded_gt_ht.annotate(
        _split_count=hl.or_else(
            heavy_variants_with_contribution[encoded_gt_ht.v_idx].split_count,
            hl.int32(1),
        ),
        _contribution=hl.or_else(
            heavy_variants_with_contribution[encoded_gt_ht.v_idx].contribution,
            hl.int64(0),
        ),
    )
    enc = enc.annotate(
        _per_split_contribution=enc._contribution // hl.int64(enc._split_count),
        _split_idx=hl.range(0, enc._split_count),
    ).explode("_split_idx")
    # _split_idx came out of hl.range as int32; just promote it into the key.
    enc = enc.key_by("v_idx", "_split_idx")

    encoded_path = hl.utils.new_temp_file(f"encoded_split_{label}", "ht")
    enc.write(encoded_path, overwrite=True)
    enc = hl.read_table(encoded_path)

    partition_intervals = calculate_partitions_by_size(
        enc, n_partitions, size_field="_per_split_contribution",
    )
    encoded_v = hl.read_table(encoded_path, _intervals=partition_intervals)

    # --- 2-5: produce vp_collected = vp_with_v1 ---
    # One-shot resume hook: when _RESUME_HEAVY_VP_WITH_V1_PATH is set and
    # label == "heavy", skip the three-day v1-zip-join work and read the
    # surviving intermediate from the previous failed run.
    if _RESUME_HEAVY_VP_WITH_V1_PATH and label == "heavy":
        logger.warning(
            "[ONE-SHOT RESUME] Skipping steps 2-5 of heavy; reading "
            "vp_with_v1 from %s",
            _RESUME_HEAVY_VP_WITH_V1_PATH,
        )
        vp_collected = hl.read_table(_RESUME_HEAVY_VP_WITH_V1_PATH, _intervals=partition_intervals)
    else:
        # --- 2. Tag every pair with the split_idx it will land in on each side ---
        # v1_split_idx = hash(v2_idx) mod split_count(v1) — distributes a heavy
        # v1's 4,091 partners evenly across its split_count copies.
        # v2_split_idx = hash(v1_idx) mod split_count(v2) — symmetric.
        vp_subset = vp_subset.annotate(
            _v1_split_count=hl.or_else(
                heavy_variants_with_contribution[vp_subset.v1_idx].split_count,
                hl.int32(1),
            ),
            _v2_split_count=hl.or_else(
                heavy_variants_with_contribution[vp_subset.v2_idx].split_count,
                hl.int32(1),
            ),
        )
        # Scramble the opposite side's v_idx with Knuth's multiplicative hash
        # before the modulo. Hail 0.2.134 has no hl.hash, and plain
        # ``v2_idx % split_count`` could bias the assignment if partner v_idxs
        # are locally clustered (chromosome-local genes). Multiply-then-modulo
        # spreads them. v_idx values are non-negative int64 ≤ ~1e8, so the
        # product fits in int64 without overflow; the modulo against int32
        # split_count then safely casts back to int32.
        _KNUTH_HASH = hl.int64(2654435761)
        vp_subset = vp_subset.annotate(
            _v1_split_idx=hl.if_else(
                vp_subset._v1_split_count > 1,
                hl.int32(
                    (vp_subset.v2_idx * _KNUTH_HASH)
                    % hl.int64(vp_subset._v1_split_count)
                ),
                hl.int32(0),
            ),
            _v2_split_idx=hl.if_else(
                vp_subset._v2_split_count > 1,
                hl.int32(
                    (vp_subset.v1_idx * _KNUTH_HASH)
                    % hl.int64(vp_subset._v2_split_count)
                ),
                hl.int32(0),
            ),
        )

        # --- 3. Group by (v1_idx, _v1_split_idx) → collect partners ---
        # Partners carry their own _v2_split_idx forward so step 5 can rekey
        # without having to recompute it (and without joining heavy_variants
        # again post-v1-join).
        vp_collected = vp_subset.group_by(
            "v1_idx", "_v1_split_idx",
        ).aggregate(
            # Design C: carry ONLY the integer partner key (v2_idx) + its
            # routing split-idx through the shuffle — no locus/alleles. v1_idx
            # is the group key (recovered from v_idx below); the 4-tuple key is
            # rebuilt at the very end via _remap_vidx_pairs_to_loci.
            pairs=hl.agg.collect(
                hl.struct(
                    v2_idx=vp_subset.v2_idx,
                    _v2_split_idx=vp_subset._v2_split_idx,
                )
            ),
        )
        vp_collected = vp_collected.key_by(
            v_idx=vp_collected.v1_idx,
            _split_idx=vp_collected._v1_split_idx,
        ).drop("v1_idx", "_v1_split_idx")

        # --- 4. Batch each (v1, split) row's partner list into chunks ---
        # Keeps per-row size bounded even when degree(v1)/split_count is large.
        PAIR_BATCH_SIZE = 250
        vp_collected = vp_collected.annotate(
            pairs=hl.range(
                0, hl.len(vp_collected.pairs), PAIR_BATCH_SIZE,
            ).map(lambda i: vp_collected.pairs[i:i + PAIR_BATCH_SIZE]),
        ).explode("pairs")

        vp_collected_path = hl.utils.new_temp_file(f"vp_collected_v1_{label}", "ht")
        vp_collected.write(vp_collected_path, overwrite=True)
        vp_collected = hl.read_table(
            vp_collected_path, _intervals=partition_intervals,
        )

        # --- 5. v1 zip-join: annotate with the (v_idx, _split_idx) copy of v1 ---
        vp_collected = vp_collected.annotate(
            v1=encoded_v[vp_collected.v_idx, vp_collected._split_idx],
        )
        vp_with_v1_path = hl.utils.new_temp_file(f"vp_with_v1_sets_{label}", "ht")
        vp_collected.write(vp_with_v1_path, overwrite=True)
        vp_collected = hl.read_table(vp_with_v1_path)

    # --- 6. Explode partner batches, rekey by (v2_idx, _v2_split_idx) ---
    # We deliberately unkey + drop the v1-side (v_idx, _split_idx) before
    # creating the new v2-side key with the same names, otherwise the
    # key_by kwargs would collide with the existing v1-side key columns.
    vp_exploded = vp_collected.explode("pairs")
    vp_exploded = vp_exploded.transmute(
        v2_idx=vp_exploded.pairs.v2_idx,
        _v2_split_idx=vp_exploded.pairs._v2_split_idx,
    )
    # v_idx here is the v1-side idx (the group key) — capture it as v1_idx
    # before dropping the v1-side (v_idx, _split_idx) key, so the final result
    # can carry (v1_idx, v2_idx). v2_idx is kept as a value field alongside the
    # new v2-side key (v_idx == v2_idx).
    vp_exploded = vp_exploded.annotate(v1_idx=vp_exploded.v_idx)
    vp_exploded = vp_exploded.key_by().drop("v_idx", "_split_idx").cache()
    vp_exploded = vp_exploded.key_by(
        v_idx=vp_exploded.v2_idx,
        _split_idx=vp_exploded._v2_split_idx,
    ).drop("_v2_split_idx")

    vp_by_v2_path = hl.utils.new_temp_file(f"vp_by_v2_{label}", "ht")
    # The shuffle that feeds this write (rekey by v2_idx of the v1-joined,
    # exploded pair stream) has been the failure point on autoscaling
    # clusters — torn partition files from preemptible-secondary
    # decommission. Hail's "new shuffle" routes that shuffle through Hail's
    # own service instead of Spark's local-disk model, which sidesteps the
    # tornn-file failure mode. Turn it off immediately after so unrelated
    # downstream ops keep their default behavior.
    #hl._set_flags(use_new_shuffle="1")
    vp_exploded.write(vp_by_v2_path, overwrite=True)
    #hl._set_flags(use_new_shuffle=None)
    vp_exploded = hl.read_table(vp_by_v2_path, _intervals=partition_intervals).cache()
    encoded_v2 = hl.read_table(encoded_path, _intervals=partition_intervals).cache()

    # --- 7. v2 zip-join + per-pair counts ---
    v1 = vp_exploded.v1
    v2 = encoded_v2[vp_exploded.v_idx, vp_exploded._split_idx]
    gt_counts_raw = _count_from_sets(
        v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        n_samples,
        include_raw_hr_adj_missing=True,
    )
    gt_counts_adj = _count_from_sets(
        v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        n_samples,
        include_raw_hr_adj_missing=False,
    )
    count_fields = dict(gt_counts_raw=gt_counts_raw, gt_counts_adj=gt_counts_adj)
    if pops is not None:
        count_fields["gt_counts_by_pop"] = _count_from_sets_by_pop(
            v1, v2, pops, pop_index_sets, pop_sizes, gt_counts_raw, gt_counts_adj,
        )
    if emit_phase:
        phase = _count_phase_from_sets(
            v1.adj_het, v1.phased_het, v2.adj_het, v2.phased_het,
        )
        count_fields["n_phased_cis"] = phase.n_phased_cis
        count_fields["n_phased_trans"] = phase.n_phased_trans
    return vp_exploded.select(
        "v1_idx",
        "v2_idx",
        **count_fields,
    ).cache()


TARGET_HEAVY_PARTITION_BYTES = 500 * 1024 ** 2
"""Target per-partition shuffle-data budget for the heavy step (500 MB).

Used in two coupled places:

1. :func:`build_variant_size_info_ht` divides each heavy variant's
   contribution by this value to set ``split_count`` — the number of
   duplicate copies of that variant's encoded row to spread across
   partitions. A variant whose contribution exceeds the target gets
   ``split_count > 1`` so its data is spread across multiple partitions
   instead of piling up on one.

2. :func:`compute_counts_heavy` clamps the heavy step's partition count to
   ``ceil(total_heavy_contribution / TARGET_HEAVY_PARTITION_BYTES)``,
   floored at :data:`MIN_HEAVY_PARTITIONS` and capped at
   ``max_join_partitions``.
"""

MIN_HEAVY_PARTITIONS = 50
"""Floor on heavy-step partition count, regardless of total contribution."""

_RESUME_HEAVY_VP_WITH_V1_PATH: Optional[str] = None
"""ONE-SHOT resume hook for a heavy step that failed at step 6 with a
preemptible-secondary truncated shuffle. When set and ``label == "heavy"``,
:func:`_compute_counts_for_subset` skips steps 2-5 and reads ``vp_with_v1``
directly from this path. Leave ``None`` for normal runs — a stale path from
an unrelated postfix would inject the wrong intermediate into the current
heavy step.

Note: under Design C the collected ``pairs`` struct carries only
``(v2_idx, _v2_split_idx)`` (no locus/alleles), so a resume must point at an
intermediate produced by this Design-C code, not a baseline one."""

_RESUME_LIGHT_GT_PATH: Optional[str] = None
_RESUME_LIGHT_VP_PATH: Optional[str] = None
"""ONE-SHOT resume hooks for :func:`compute_counts_light`. Both intermediates
(``gt_light``, ``vp_light``) are written by the light step just before the
co-partition + zip-join. When these are set, the function skips the
v_idx-annotation, heavy-filter, ``semi_join``, and both writes; it reads
the two paths directly and jumps to co-partitioning. Reset to ``None`` after
the resume run completes — stale paths from a different postfix would inject
the wrong intermediates.

Note: under Design C ``vp_light`` is keyed by ``(v1_idx, v2_idx)`` and
carries no locus/alleles, so a resume must point at a Design-C-produced
``vp_light``."""


def encode_genotypes(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    output_dir: str,
    min_an_pct: int = -1,
    *,
    use_precomputed_adj: bool = False,
) -> None:
    """
    Step A: Encode genotypes and write the encoded GT table.

    Creates ``var_idx.ht`` and ``encoded_gt_sets_by_var_idx.ht`` in *output_dir*.
    Run once; reuse for light and heavy compute steps with any threshold.

    :param mt: Dense filtered MatrixTable.
    :param vp_ht: Variant pair list Table (unused directly, but mt must
        contain exactly the variants in vp_ht).
    :param output_dir: GCS directory for outputs.
    :param min_an_pct: AN_percent floor of the dense MT, carried onto the
        encoded table's globals so the count steps can refuse to lower it.
    :param use_precomputed_adj: Forwarded to
        :func:`_encode_genotype_sets_by_var_idx`; read adj from ``mt.adj``
        instead of recomputing from GQ/DP/AD.
    """
    var_idx_path = f"{output_dir}/var_idx.ht"
    logger.info("Writing var_idx to %s", var_idx_path)
    _create_var_idx_ht(mt).write(var_idx_path, overwrite=True)
    var_idx_ht = hl.read_table(var_idx_path)

    encoded_path = f"{output_dir}/encoded_gt_sets_by_var_idx.ht"
    logger.info("Encoding genotypes to %s", encoded_path)
    _encode_genotype_sets_by_var_idx(
        mt, var_idx_ht, use_precomputed_adj=use_precomputed_adj
    ).annotate_globals(min_an_pct=min_an_pct).write(encoded_path, overwrite=True)
    logger.info("Encoded genotypes written.")


_COUNT_FROM_SETS_FIELDS = (
    "raw_het", "raw_hv",
    "all_samples", "n_with_data",
    "raw_hr_adj_missing", "n_raw_hr_adj_missing",
    "adj_het", "adj_hv",
    "phased_het",
)


def _drop_pairs_missing_v_idx(vp_ht: hl.Table, caller: str) -> hl.Table:
    """Drop pairs whose v1_idx or v2_idx is missing, logging the count.

    Pairs land here when ``vp_ht``'s ``(locus, alleles)`` doesn't appear in
    ``var_idx_ht`` (the index over the encoded MT's rows). The most common
    cause is the pair list and the dense MT having a different
    multi-allelic split, or some pair variants being absent from the dense
    MT. Without this filter the silent failure mode is row loss downstream
    — ``hl.read_table(_intervals=partition_intervals)`` matches its
    intervals against the first key column and a NULL key falls outside
    every interval, so the row is dropped without warning.

    :param vp_ht: Pair Table annotated with ``v1_idx`` / ``v2_idx``.
    :param caller: Caller name for the log message (e.g. "compute_counts_light").
    :return: ``vp_ht`` filtered to rows with both v_idx fields defined.
    """
    n_missing = vp_ht.aggregate(
        hl.agg.count_where(
            hl.is_missing(vp_ht.v1_idx) | hl.is_missing(vp_ht.v2_idx)
        )
    )
    if n_missing > 0:
        logger.warning(
            "%s: %d pairs have v_idx missing for v1 and/or v2 "
            "(variant not in var_idx_ht); dropping.",
            caller, n_missing,
        )
        vp_ht = vp_ht.filter(
            hl.is_defined(vp_ht.v1_idx) & hl.is_defined(vp_ht.v2_idx)
        )
    return vp_ht


def _project_count_fields(encoded_ht: hl.Table) -> hl.Table:
    """Restrict an encoded table to the fields the count paths read.

    The encoder currently emits exactly ``_COUNT_FROM_SETS_FIELDS``, so this
    is effectively a no-op today — kept as a defensive projection so any
    future diagnostic column added to the encoder doesn't silently ride
    through every count shuffle.
    """
    return encoded_ht.select(*_COUNT_FROM_SETS_FIELDS)


def _empty_counts_ht(
    reference_genome: str = "GRCh38",
    *,
    stratify_by_pop: bool = False,
    emit_phase: bool = False,
) -> hl.Table:
    """Empty (locus1, alleles1, locus2, alleles2)-keyed counts HT.

    Used by :func:`compute_counts_heavy` to return a typed-empty result when
    no heavy variants exist, so its caller can ``.union(...)`` with the light
    side unconditionally. When ``stratify_by_pop`` is set, the schema also
    carries the ``gt_counts_by_pop`` dict so it unions with a per-pop light
    result; when ``emit_phase`` is set it carries ``n_phased_cis`` /
    ``n_phased_trans`` so it unions with a phase-annotated light result. (The
    ``gt_counts_*_no_pbt`` columns are added later, at ``--combine-counts``, by
    subtracting the PBT∩cohort counts — not carried on the light/heavy tables.)
    """
    # tint32 must match what _count_from_sets actually returns (its set
    # algebra produces array<int32>); a tint64 schema here makes
    # compute_counts_heavy's empty table fail to union with
    # compute_counts_light's result whenever the heavy filter is empty.
    fields = dict(
        locus1=hl.tlocus(reference_genome=reference_genome),
        alleles1=hl.tarray(hl.tstr),
        locus2=hl.tlocus(reference_genome=reference_genome),
        alleles2=hl.tarray(hl.tstr),
        gt_counts_raw=hl.tarray(hl.tint32),
        gt_counts_adj=hl.tarray(hl.tint32),
    )
    if stratify_by_pop:
        fields["gt_counts_by_pop"] = hl.tdict(
            hl.tstr,
            hl.tstruct(raw=hl.tarray(hl.tint32), adj=hl.tarray(hl.tint32)),
        )
    if emit_phase:
        fields["n_phased_cis"] = hl.tint32
        fields["n_phased_trans"] = hl.tint32
    return hl.Table.parallelize(
        [], schema=hl.tstruct(**fields),
        key=["locus1", "alleles1", "locus2", "alleles2"],
    )


def _remap_vidx_pairs_to_loci(counts_ht: hl.Table, var_idx_ht: hl.Table) -> hl.Table:
    """Rebuild the ``(locus1, alleles1, locus2, alleles2)`` key from v_idx pairs.

    Design C carries only the integer keys ``(v1_idx, v2_idx)`` through the
    count-step shuffles / intermediates (slimmer shuffled rows); the
    locus/alleles 4-tuple is reconstructed HERE, once, right before each count
    function returns. ``var_idx_ht`` (``(locus, alleles) → var_idx``) is joined
    twice — v1_idx → locus1/alleles1, v2_idx → locus2/alleles2 — so the WRITTEN
    count outputs keep the baseline ``(locus1, alleles1, locus2, alleles2)`` key
    and schema and ``--combine-counts`` / downstream are unchanged.

    :param counts_ht: Counts Table carrying ``v1_idx`` / ``v2_idx`` fields
        (key or value) plus the count columns.
    :param var_idx_ht: ``(locus, alleles) → var_idx`` lookup.
    :return: ``counts_ht`` keyed by ``(locus1, alleles1, locus2, alleles2)`` with
        ``v1_idx`` / ``v2_idx`` dropped.
    """
    var_idx_by_vidx = var_idx_ht.key_by("var_idx")
    v1 = var_idx_by_vidx[counts_ht.v1_idx]
    v2 = var_idx_by_vidx[counts_ht.v2_idx]
    counts_ht = counts_ht.annotate(
        locus1=v1.locus, alleles1=v1.alleles,
        locus2=v2.locus, alleles2=v2.alleles,
    )
    return counts_ht.key_by(
        "locus1", "alleles1", "locus2", "alleles2"
    ).drop("v1_idx", "v2_idx")


def _count_pairs_via_index(
    vp: hl.Table,
    encoded_gt_ht: hl.Table,
    n_samples,
    *,
    pops=None,
    pop_index_sets=None,
    pop_sizes=None,
    emit_phase: bool = False,
) -> hl.Table:
    """Per-pair _count_from_sets via direct table indexing on var_idx.

    No co-partition, no shuffle, no semi_join. Used whenever the encoded
    table is small enough that Hail's auto-join (hash / broadcast / sort-
    merge) is cheaper than an explicit shuffle setup. Equivalent to the
    benchmark's ``approach_b``.

    Design C: ``vp`` must carry ``v1_idx`` / ``v2_idx``; the result is keyed by
    ``(v1_idx, v2_idx)`` and carries ONLY the count columns (no locus/alleles).
    The caller reconstructs the 4-tuple key via :func:`_remap_vidx_pairs_to_loci`.

    When ``pops`` is provided (with ``pop_index_sets`` / ``pop_sizes`` from
    :func:`_build_pop_stratification`), additionally emits a
    ``gt_counts_by_pop`` dict<pop, struct{raw, adj}> alongside the flat
    ``gt_counts_raw`` / ``gt_counts_adj`` (which remain the full-cohort counts).
    """
    v1 = encoded_gt_ht[vp.v1_idx]
    v2 = encoded_gt_ht[vp.v2_idx]
    gt_counts_raw = _count_from_sets(
        v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        n_samples,
        include_raw_hr_adj_missing=True,
    )
    gt_counts_adj = _count_from_sets(
        v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        n_samples,
        include_raw_hr_adj_missing=False,
    )
    count_fields = dict(gt_counts_raw=gt_counts_raw, gt_counts_adj=gt_counts_adj)
    if pops is not None:
        count_fields["gt_counts_by_pop"] = _count_from_sets_by_pop(
            v1, v2, pops, pop_index_sets, pop_sizes, gt_counts_raw, gt_counts_adj,
        )
    if emit_phase:
        phase = _count_phase_from_sets(
            v1.adj_het, v1.phased_het, v2.adj_het, v2.phased_het,
        )
        count_fields["n_phased_cis"] = phase.n_phased_cis
        count_fields["n_phased_trans"] = phase.n_phased_trans
    # vp is keyed by (v1_idx, v2_idx); select only the count columns (the key is
    # auto-preserved — listing key fields positionally is rejected by Hail).
    vp = vp.select(**count_fields).cache()
    return vp.key_by("v1_idx", "v2_idx").cache()


def build_variant_size_info_ht(
    encoded_ht: hl.Table,
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    variant_filter_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """Build the per-variant size-info Table.

    Pure transform: HTs in → HT out. The caller is responsible for
    writing the result to a resource path.

    No heavy/light decision is made here. The decision lives entirely
    at consumption time (``--heavy-contribution-cutoff`` on the count
    steps) so the cutoff can be retuned by inspecting this HT without
    a rebuild.

    One row per variant in ``encoded_ht``, keyed by ``v_idx``:

    * ``locus``, ``alleles`` — joined from ``var_idx_ht``.
    * ``gene_id`` (array<str>) — joined from ``variant_filter_ht``
      (omitted when ``variant_filter_ht`` is ``None``).
    * ``_contribution`` — ``degree(v) × payload(v)`` bytes.
    * ``_cum_before`` — cumulative ``_contribution`` of variants ranked
      above v in the descending-by-contribution order. Useful for
      driving a greedy "minimum cutoff at which the remaining light
      contribution drops below X" calculation.
    * ``_bytes`` — partner-load metric: for each variant v, the maximum
      of ``Σ payload(partner)`` over (v1=v, partner=v2) pairs and
      (v1=partner, v2=v) pairs, restricted to partners whose own
      contribution exceeds ``TARGET_HEAVY_PARTITION_BYTES`` (i.e.
      partners that themselves would force ``split_count > 1``). This
      approximates the bucket size at the v1/v2 zip-join. The fixed
      ``TARGET`` gate keeps the metric cutoff-independent.
    * ``split_count`` — ``ceil(max(_contribution, _bytes) / TARGET)``
      where ``TARGET = TARGET_HEAVY_PARTITION_BYTES = 500 MB``.

    Globals: ``target_heavy_partition_bytes`` (the only "fixed"
    threshold used at build time).

    :param encoded_ht: Encoded GT table, post-``_project_count_fields``.
    :param vp_ht: Variant pair list Table.
    :param var_idx_ht: ``(locus, alleles) → v_idx`` lookup.
    :param variant_filter_ht: Variant filter HT (must have ``gene_id``).
        Optional; if omitted, ``gene_id`` is left off the output.
    :return: Materialized variant size-info Table keyed by ``v_idx``.
    """
    # === 1. Degree per variant ===
    vp_with_idx = vp_ht.annotate(
        _v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        _v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    side_rows = vp_with_idx.select(
        _e=hl.array([vp_with_idx._v1_idx, vp_with_idx._v2_idx])
    ).key_by().explode("_e")
    side_rows = side_rows.select(v_idx=side_rows._e)
    degree_ht = (
        side_rows.group_by("v_idx")
        .aggregate(degree=hl.agg.count())
        .cache()
    )

    # === 2. Payload + contribution per encoded row ===
    # Set elements are int32 (4 bytes). The phased_het sidecar rides through
    # the heavy shuffle too; each entry is an int32 key + struct{pid:str,
    # gt0:int32} — budget ~20 bytes/entry so phase-heavy variants aren't
    # under-sized (empty dict on phase-less encodings → adds nothing).
    payload_expr = hl.int64(
        hl.len(encoded_ht.all_samples)
        + hl.len(encoded_ht.raw_het) + hl.len(encoded_ht.raw_hv)
        + hl.len(encoded_ht.adj_het) + hl.len(encoded_ht.adj_hv)
        + hl.len(encoded_ht.raw_hr_adj_missing)
    ) * 4 + hl.int64(hl.len(encoded_ht.phased_het)) * 20
    enc = encoded_ht.select(
        _payload=payload_expr,
        _contribution=(
            hl.or_else(degree_ht[encoded_ht.v_idx].degree, hl.int64(0))
            * payload_expr
        ),
    ).cache()

    # === 3. Total + count (one fused pass); used as driver constants ===
    total, n_total = enc.aggregate((
        hl.agg.sum(enc._contribution),
        hl.agg.count(),
    ))

    # === 4. Sort desc by contribution; scan-cumsum ===
    enc_sorted = enc.order_by(hl.desc(enc._contribution))
    enc_sorted = enc_sorted.annotate(
        _cum_before=hl.scan.sum(enc_sorted._contribution),
    )
    enc_sorted = enc_sorted.key_by("v_idx").select(
        "_contribution", "_cum_before",
    )
    enc_sorted = enc_sorted.checkpoint(
        hl.utils.new_temp_file("variant_size_info_pre_bytes", "ht"),
    )

    # === 5. Cross-role partner load (_bytes) ===
    # For each variant v: sum payloads of partners whose OWN contribution
    # exceeds TARGET (i.e. partners that themselves force split_count>1
    # and so are guaranteed to carry their own encoded data through the
    # v1 or v2 zip-join). Doing this with a fixed TARGET gate keeps
    # _bytes cutoff-independent: any heavy cutoff ≥ TARGET produces the
    # same partner set; below TARGET only adds variants that don't
    # affect bucket size anyway.
    v1_enc = enc[vp_with_idx._v1_idx]
    v2_enc = enc[vp_with_idx._v2_idx]
    pairs = vp_with_idx.annotate(
        _v1_payload=hl.or_else(v1_enc._payload, hl.int64(0)),
        _v2_payload=hl.or_else(v2_enc._payload, hl.int64(0)),
        _v1_over_target=hl.or_else(
            v1_enc._contribution > TARGET_HEAVY_PARTITION_BYTES, False,
        ),
        _v2_over_target=hl.or_else(
            v2_enc._contribution > TARGET_HEAVY_PARTITION_BYTES, False,
        ),
    ).cache()

    pairs_v1heavy = pairs.filter(pairs._v1_over_target)
    v2_load = (
        pairs_v1heavy.group_by(v_idx=pairs_v1heavy._v2_idx)
        .aggregate(_v2_load=hl.agg.sum(pairs_v1heavy._v1_payload))
    )
    pairs_v2heavy = pairs.filter(pairs._v2_over_target)
    v1_load = (
        pairs_v2heavy.group_by(v_idx=pairs_v2heavy._v1_idx)
        .aggregate(_v1_load=hl.agg.sum(pairs_v2heavy._v2_payload))
    )

    load = v1_load.join(v2_load, how="outer")
    load = load.select(
        _bytes=hl.max(
            hl.or_else(load._v1_load, hl.int64(0)),
            hl.or_else(load._v2_load, hl.int64(0)),
        ),
    )

    # === 6. Assemble: contribution + cum + bytes + split ===
    size_info = enc_sorted.annotate(
        _bytes=hl.or_else(load[enc_sorted.v_idx]._bytes, hl.int64(0)),
    )
    size_info = size_info.annotate(
        split_count=hl.max(
            hl.int32(1),
            hl.int32(
                (hl.max(size_info._contribution, size_info._bytes)
                 + TARGET_HEAVY_PARTITION_BYTES - 1)
                // TARGET_HEAVY_PARTITION_BYTES
            ),
        ),
    )

    # === 7. Join locus/alleles/gene_id for human-readable diagnostics ===
    var_idx_by_v_idx = var_idx_ht.key_by("var_idx")
    size_info = size_info.annotate(
        locus=var_idx_by_v_idx[size_info.v_idx].locus,
        alleles=var_idx_by_v_idx[size_info.v_idx].alleles,
    )
    if variant_filter_ht is not None:
        size_info = size_info.annotate(
            gene_id=variant_filter_ht[size_info.locus, size_info.alleles].gene_id,
        )

    # === 8. Globals + diagnostic log ===
    size_info = size_info.annotate_globals(
        target_heavy_partition_bytes=hl.int64(TARGET_HEAVY_PARTITION_BYTES),
    )
    sum_bytes_all, max_split, n_over_target = size_info.aggregate((
        hl.agg.sum(size_info._bytes),
        hl.agg.max(size_info.split_count),
        hl.agg.count_where(
            size_info._contribution > TARGET_HEAVY_PARTITION_BYTES
        ),
    ))
    logger.info(
        "Variant size info: total=%.2f GB contrib across %d variants; "
        "%d variants with _contribution > TARGET (=split_count>1 by own "
        "contribution alone); Σ _bytes=%.2f GB; max split_count=%d "
        "(TARGET=%.0f MB).",
        total / 1024**3, n_total, n_over_target,
        sum_bytes_all / 1024**3, max_split or 0,
        TARGET_HEAVY_PARTITION_BYTES / 1024**2,
    )
    return size_info


def _size_info_to_heavy_variants(
    size_info_ht: hl.Table,
    heavy_contribution_cutoff: int,
    excluded_genes_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """Pure transform: filter the size-info HT to the heavy subset.

    Heavy iff ``_contribution >= heavy_contribution_cutoff``. When
    ``excluded_genes_ht`` is provided (keyed by ``gene_id``), variants
    whose ``gene_id`` array intersects the excluded set are dropped.

    Returns a Table keyed by ``v_idx`` with ``contribution`` and
    ``split_count`` — the schema consumed by
    :func:`_compute_counts_for_subset`.
    """
    if excluded_genes_ht is not None and "gene_id" in size_info_ht.row:
        excluded_set = excluded_genes_ht.aggregate(
            hl.agg.collect_as_set(excluded_genes_ht.gene_id)
        )
        size_info_ht = size_info_ht.filter(
            hl.len(
                hl.set(size_info_ht.gene_id).intersection(
                    hl.literal(excluded_set)
                )
            ) == 0
        )
    size_info_ht = size_info_ht.filter(
        size_info_ht._contribution >= heavy_contribution_cutoff
    )
    return size_info_ht.select(
        contribution=size_info_ht._contribution,
        split_count=size_info_ht.split_count,
    )


def _filter_pairs_by_excluded_genes(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    size_info_ht: hl.Table,
    excluded_genes_ht: hl.Table,
) -> hl.Table:
    """Drop pairs from ``vp_ht`` where EITHER side belongs to an excluded
    gene (gene_id intersection with the excluded set is non-empty).

    Uses the variant size-info HT for the variant → gene_id mapping
    (joined from variant_filter_ht at build time)."""
    excluded_set = excluded_genes_ht.aggregate(
        hl.agg.collect_as_set(excluded_genes_ht.gene_id)
    )
    excl_lit = hl.literal(excluded_set)
    gene_by_v_idx = size_info_ht.select("gene_id")
    vp_ht = vp_ht.annotate(
        _v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        _v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    vp_ht = vp_ht.annotate(
        _v1_genes=gene_by_v_idx[vp_ht._v1_idx].gene_id,
        _v2_genes=gene_by_v_idx[vp_ht._v2_idx].gene_id,
    )
    vp_ht = vp_ht.filter(
        (hl.len(hl.set(vp_ht._v1_genes).intersection(excl_lit)) == 0)
        & (hl.len(hl.set(vp_ht._v2_genes).intersection(excl_lit)) == 0)
    )
    return vp_ht.drop("_v1_idx", "_v2_idx", "_v1_genes", "_v2_genes")


def count_all_pairs_via_index(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    encoded_gt_ht: hl.Table,
    *,
    pop_ht: Optional[hl.Table] = None,
    pops: Optional[list] = None,
    emit_phase: bool = False,
) -> hl.Table:
    """Count every pair via the per-pair indexed-lookup plan.

    Use this entry point when no variants need splitting (i.e.
    ``heavy_variants`` from
    :func:`_size_info_to_heavy_variants` would be empty). All pairs
    are treated as light; Hail's table-indexing planner picks the
    join. Equivalent to the benchmark's ``approach_b``.

    Handles v_idx annotation, ``_drop_pairs_missing_v_idx``,
    :func:`_project_count_fields`, and the ``n_samples`` lookup —
    callers just pass in the raw HTs.

    :param vp_ht: Variant pair list Table.
    :param var_idx_ht: ``(locus, alleles) → v_idx`` lookup.
    :param encoded_gt_ht: Encoded GT table (pre-projection).
    :param pop_ht: When provided (``s`` → ``pop`` meta table), also emit a
        per-population ``gt_counts_by_pop`` dict. The pop label is joined to the
        encoded ``samples`` global at count time (no re-encode needed). The
        per-pop counts intersect each variant's stored sets with per-pop
        sample-index sets, so this multiplies the per-pair set-intersection
        work by the number of groups; cheap at gene/light scale, meaningful at
        full scale.
    :param pops: Optional subset of genetic-ancestry groups to stratify by
        (list of group names, e.g. ``["nfe", "afr"]``); ``None`` = all groups
        present. ``GLOBAL_POP`` ("all") is always included regardless.
    :param emit_phase: When ``True``, also emit ``n_phased_cis`` /
        ``n_phased_trans`` — the physically-phased refinement of the adj
        ``AaBb`` cell (see :func:`_count_phase_from_sets`). Requires the encoded
        table to carry ``phased_het``; contributes zeros when it is empty.
    :return: Counts Table with gt_counts_raw / gt_counts_adj (and
        gt_counts_by_pop when ``pop_ht`` is given, n_phased_cis / n_phased_trans
        when ``emit_phase``).
    """
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())
    pop_kwargs = {}
    if pop_ht is not None:
        strat_pops, pop_index_sets, pop_sizes = _pop_stratification_for(
            encoded_gt_ht, pop_ht, pops
        )
        pop_kwargs = dict(
            pops=strat_pops, pop_index_sets=pop_index_sets, pop_sizes=pop_sizes,
        )
    encoded_gt_ht = _project_count_fields(encoded_gt_ht)
    vp_ht = vp_ht.annotate(
        v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    vp_ht = _drop_pairs_missing_v_idx(vp_ht, "count_all_pairs_via_index")
    # Design C: carry only the integer keys through the count; the locus/alleles
    # 4-tuple is rebuilt at the end via _remap_vidx_pairs_to_loci.
    vp_ht = vp_ht.key_by("v1_idx", "v2_idx").select()
    result = _count_pairs_via_index(
        vp_ht, encoded_gt_ht, n_samples,
        emit_phase=emit_phase, **pop_kwargs,
    )
    return _remap_vidx_pairs_to_loci(result, var_idx_ht)


def compute_counts_light(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    encoded_gt_ht: hl.Table,
    heavy_variants: hl.Table,
    max_join_partitions: int = 10000,
    *,
    size_info_ht: Optional[hl.Table] = None,
    pop_ht: Optional[hl.Table] = None,
    pops: Optional[list] = None,
    emit_phase: bool = False,
) -> hl.Table:
    """
    Step B: Compute genotype counts for the light split.

    Pure transform: HTs in → counts HT out. **Assumes
    ``heavy_variants`` is non-empty** — i.e. that an actual light/heavy
    split is happening. When ``heavy_variants`` is empty, the caller
    should use :func:`count_all_pairs_via_index` directly (every pair
    is light, no shuffle plan needed).

    Filters pairs to light-only (neither side in ``heavy_variants``),
    builds the small ``gt_light`` table via ``semi_join``,
    co-partitions on v1_idx, then v1 zip-joins to count.

    The caller is responsible for: reading the encoded GT table,
    building/filtering ``heavy_variants`` (via
    :func:`_size_info_to_heavy_variants`), applying any gene exclusion
    to ``vp_ht`` (via :func:`_filter_pairs_by_excluded_genes`), and
    writing the result.

    :param vp_ht: Variant pair list Table.
    :param var_idx_ht: ``(locus, alleles) → v_idx`` lookup.
    :param encoded_gt_ht: Encoded GT table (pre-projection).
    :param heavy_variants: Heavy subset of the variant size-info HT,
        keyed by ``v_idx``. Must be non-empty.
    :param max_join_partitions: Upper bound on partition count.
    :param size_info_ht: Full variant size-info HT (keyed by ``v_idx`` with
        ``_contribution``). When provided, the light join's partition count and
        size-balance are both driven by ``_contribution`` (= degree × payload
        bytes, mirroring the heavy path) instead of the encoded HT's incidental
        partition count. When ``None`` (e.g. per-pop restricted counts), falls
        back to ``n_with_data``-weighted balancing at ``n_partitions() * 3``.
    :return: Counts Table with gt_counts_raw and gt_counts_adj. Internal
        ``hl.utils.new_temp_file`` writes are kept (algorithmic — needed
        for the partition-interval re-read pattern).
    """
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())
    pop_kwargs = {}
    if pop_ht is not None:
        strat_pops, pop_index_sets, pop_sizes = _pop_stratification_for(
            encoded_gt_ht, pop_ht, pops
        )
        pop_kwargs = dict(
            pops=strat_pops, pop_index_sets=pop_index_sets, pop_sizes=pop_sizes,
        )
    encoded_gt_ht = _project_count_fields(encoded_gt_ht)

    if _RESUME_LIGHT_GT_PATH and _RESUME_LIGHT_VP_PATH:
        logger.warning(
            "compute_counts_light: RESUMING from precomputed intermediates."
            " gt_light=%s, vp_light=%s."
            " Skipping v_idx-annotate + heavy-filter + semi_join + writes.",
            _RESUME_LIGHT_GT_PATH, _RESUME_LIGHT_VP_PATH,
        )
        gt_light_path = _RESUME_LIGHT_GT_PATH
        vp_light_path = _RESUME_LIGHT_VP_PATH
        gt_light = hl.read_table(gt_light_path)
        vp_light = hl.read_table(vp_light_path)
    else:
        # Add v_idx to pair table.
        vp_ht = vp_ht.annotate(
            v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
            v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
        )
        vp_ht = _drop_pairs_missing_v_idx(vp_ht, "compute_counts_light")

        logger.info("compute_counts_light: building light-pair shuffle plan.")
        vp_light = vp_ht.filter(
            ~hl.is_defined(heavy_variants[vp_ht.v1_idx])
            & ~hl.is_defined(heavy_variants[vp_ht.v2_idx])
        )
        # Design C: carry ONLY the integer keys through the shuffle; the
        # locus/alleles 4-tuple is rebuilt at the end via
        # _remap_vidx_pairs_to_loci.
        vp_light = vp_light.key_by("v1_idx", "v2_idx").select().cache()

        # Build small GT table for light variants only.
        light_v1 = vp_light.key_by(v_idx=vp_light.v1_idx).select().distinct()
        light_v2 = vp_light.key_by(v_idx=vp_light.v2_idx).select().distinct()
        light_variants = light_v1.union(light_v2).distinct().cache()
        gt_light = encoded_gt_ht.semi_join(light_variants)

        gt_light_path = hl.utils.new_temp_file("gt_light", "ht")
        gt_light.write(gt_light_path, overwrite=True)
        gt_light = hl.read_table(gt_light_path)

        # Write vp_light so we can re-read with partition intervals.
        vp_light_path = hl.utils.new_temp_file("vp_light", "ht")
        vp_light.write(vp_light_path, overwrite=True)
        vp_light = hl.read_table(vp_light_path)

    logger.info(
        "Light pairs: %d, light GT variants: %d",
        vp_light.count(), gt_light.count(),
    )

    # Co-partition on v1_idx for the v1 zip-join, size-balancing the encoded
    # rows so no single partition holds a disproportionate share of the join.
    if size_info_ht is not None:
        # Preferred: partition COUNT and size-balance both driven by
        # ``_contribution`` (= degree × payload bytes, the same metric the heavy
        # path uses). ``payload`` counts the stored ``all_samples`` set, so a
        # low-AN positive-form variant (set ≈ n_samples) is scored as expensive
        # even though its ``n_with_data`` is small — the case that otherwise
        # piles low-AN variants into one partition and stalls the join for
        # hours. n_parts = total light contribution / TARGET, mirroring
        # compute_counts_heavy — NOT the old ``n_partitions() * 3``, which
        # anchored the join width to the encoded HT's incidental partition count
        # (e.g. 4 → 12 partitions for a 54 GB no-floor light set).
        gt_light_for_parts = gt_light.annotate(
            _contribution=hl.or_else(
                size_info_ht[gt_light.v_idx]._contribution, hl.int64(0)
            )
        )
        total_light_contribution = gt_light_for_parts.aggregate(
            hl.agg.sum(gt_light_for_parts._contribution)
        )
        n_parts = min(
            max(
                gt_light.n_partitions(),
                int(
                    (total_light_contribution + TARGET_HEAVY_PARTITION_BYTES - 1)
                    // TARGET_HEAVY_PARTITION_BYTES
                ),
            ),
            max_join_partitions,
        )
        logger.info(
            "compute_counts_light: total light contribution %.2f GB → "
            "n_parts=%d (target %.0f MB/partition).",
            total_light_contribution / 1024**3, n_parts,
            TARGET_HEAVY_PARTITION_BYTES / 1024**2,
        )
        partition_intervals = calculate_partitions_by_size(
            gt_light_for_parts, n_parts, size_field="_contribution",
        )
    else:
        # Fallback (per-pop restricted counts / benchmarks, no size-info HT):
        # weight ``n_with_data`` by per-variant v1 pair degree so a genomic
        # cluster of high-degree variants gets subdivided across partitions
        # instead of collapsing into one hot spot. The restricted per-group sets
        # are small, so anchoring the width to ``n_partitions() * 3`` is fine.
        v1_deg_ht = vp_light.group_by("v1_idx").aggregate(
            pair_deg=hl.int64(hl.agg.count())
        ).cache()
        gt_light_for_parts = gt_light.annotate(
            _pair_deg=hl.or_else(v1_deg_ht[gt_light.v_idx].pair_deg, hl.int64(1))
        )
        n_parts = min(gt_light.n_partitions() * 3, max_join_partitions)
        partition_intervals = calculate_partitions_by_size(
            gt_light_for_parts, n_parts, size_field="n_with_data",
            weight_field="_pair_deg",
        )
    gt_light = hl.read_table(gt_light_path, _intervals=partition_intervals).cache()
    vp_light = hl.read_table(vp_light_path, _intervals=partition_intervals).cache()

    result = _count_pairs_via_index(
        vp_light, gt_light, n_samples,
        emit_phase=emit_phase, **pop_kwargs,
    )
    return _remap_vidx_pairs_to_loci(result, var_idx_ht)


def compute_counts_heavy(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    encoded_gt_ht: hl.Table,
    heavy_variants: hl.Table,
    max_join_partitions: int = 36000,
    *,
    pop_ht: Optional[hl.Table] = None,
    pops: Optional[list] = None,
    emit_phase: bool = False,
) -> hl.Table:
    """
    Step C: Compute genotype counts for the heavy split.

    Pure transform: HTs in → counts HT out. The caller is responsible
    for reading the encoded GT table, building/filtering
    ``heavy_variants`` (via :func:`_size_info_to_heavy_variants`),
    applying any gene exclusion to ``vp_ht``, and writing the result.

    Returns an empty (typed) counts Table when no heavy variants
    remain so callers can ``.union(...)`` with ``compute_counts_light``
    unconditionally.

    :param vp_ht: Variant pair list Table.
    :param var_idx_ht: ``(locus, alleles) → v_idx`` lookup.
    :param encoded_gt_ht: Encoded GT table (pre-projection).
    :param heavy_variants: Heavy subset of the variant size-info HT.
        **Must be non-empty.** When it's empty, the caller should use
        :func:`_empty_counts_ht` directly (no heavy work to do).
    :param max_join_partitions: Upper bound on partition count.
    :return: Counts Table with gt_counts_raw and gt_counts_adj.
    """
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())
    pop_kwargs = {}
    if pop_ht is not None:
        strat_pops, pop_index_sets, pop_sizes = _pop_stratification_for(
            encoded_gt_ht, pop_ht, pops
        )
        pop_kwargs = dict(
            pops=strat_pops, pop_index_sets=pop_index_sets, pop_sizes=pop_sizes,
        )
    encoded_gt_ht = _project_count_fields(encoded_gt_ht)

    n_heavy, total_heavy_contribution = heavy_variants.aggregate(
        (hl.agg.count(), hl.agg.sum(heavy_variants.contribution))
    )
    n_heavy = int(n_heavy)
    total_heavy_contribution = int(total_heavy_contribution or 0)

    # Drive n_partitions off the actual heavy contribution (not the
    # encoded HT's arbitrary input partition count).
    n_partitions = max(
        MIN_HEAVY_PARTITIONS,
        min(
            (total_heavy_contribution + TARGET_HEAVY_PARTITION_BYTES - 1)
            // TARGET_HEAVY_PARTITION_BYTES,
            max_join_partitions,
        ),
    )
    logger.info(
        "compute_counts_heavy: %d heavy variants, "
        "total contribution %.2f GB → n_partitions=%d (target %.0f MB/partition).",
        n_heavy, total_heavy_contribution / 1024**3, n_partitions,
        TARGET_HEAVY_PARTITION_BYTES / 1024**2,
    )

    # Filter pairs to the heavy set.
    vp_ht = vp_ht.annotate(
        v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    vp_ht = _drop_pairs_missing_v_idx(vp_ht, "compute_counts_heavy")
    vp_heavy = vp_ht.filter(
        hl.is_defined(heavy_variants[vp_ht.v1_idx])
        | hl.is_defined(heavy_variants[vp_ht.v2_idx])
    )
    # Design C: carry ONLY the integer keys through the heavy shuffles; the
    # locus/alleles 4-tuple is rebuilt at the end via _remap_vidx_pairs_to_loci.
    vp_heavy = vp_heavy.select("v1_idx", "v2_idx").cache()

    result = _compute_counts_for_subset(
        vp_heavy, encoded_gt_ht, n_samples, "heavy", n_partitions,
        heavy_variants, emit_phase=emit_phase, **pop_kwargs,
    )
    return _remap_vidx_pairs_to_loci(result, var_idx_ht)


def _empty_heavy_variants() -> hl.Table:
    """Empty ``v_idx``-keyed heavy-variants table (contribution / split_count).

    Passed to :func:`compute_counts_light` to make it treat EVERY pair as light
    (co-partitioned counting), i.e. count all pairs with no light/heavy split.
    """
    return hl.Table.parallelize(
        [],
        hl.tstruct(v_idx=hl.tint64, contribution=hl.tint64, split_count=hl.tint32),
        key=["v_idx"],
    )


def compute_counts_by_pop(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    encoded_gt_ht: hl.Table,
    pop_ht: hl.Table,
    pops: Optional[list] = None,
    *,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """Per-population genotype counts by restrict-per-group + flat count.

    For each genetic-ancestry group, restricts the encoded table to that group's
    samples (:func:`restrict_encoded_to_pops`, re-indexed) and runs a
    co-partitioned flat count (:func:`compute_counts_light` with an empty heavy
    set, so every pair is counted with no light/heavy split). Each pass thus
    shuffles only the small per-group sets and emits a small flat expression —
    avoiding both the full-cohort shuffle (disk blow-out) and the
    ``ClassTooLargeException`` that the inlined multi-pop
    :func:`_count_from_sets_by_pop` hits at 3+ groups.

    Assembles ``gt_counts_by_pop`` = ``dict<pop, struct{raw, adj}>``.
    ``GLOBAL_POP`` ("all") is the element-wise sum over the counted groups
    (matching v2's fold-sum; equals the full cohort when the groups partition
    it, i.e. every sample has a group label), and is also written as the flat
    ``gt_counts_raw`` / ``gt_counts_adj``.

    :param vp_ht: Variant pair list Table.
    :param var_idx_ht: ``(locus, alleles) → v_idx`` lookup.
    :param encoded_gt_ht: Encoded GT table (pre-projection).
    :param pop_ht: ``s → pop`` meta table (see ``resources.get_sample_pop_ht``).
    :param pops: Optional subset of groups; ``None`` = all groups present.
    :return: Counts Table keyed by the pair, with ``gt_counts_raw`` /
        ``gt_counts_adj`` (= "all") and ``gt_counts_by_pop``.
    """
    strat_pops, _sets, _sizes = _pop_stratification_for(encoded_gt_ht, pop_ht, pops)
    specific = [p for p in strat_pops if p != GLOBAL_POP]
    if not specific:
        raise ValueError("compute_counts_by_pop: no genetic-ancestry groups to count.")
    logger.info("compute_counts_by_pop: counting groups %s", specific)

    empty_heavy = _empty_heavy_variants()
    per_pop = {}
    for p in specific:
        restricted = restrict_encoded_to_pops(encoded_gt_ht, pop_ht, [p]).checkpoint(
            hl.utils.new_temp_file(f"encoded_pop_{p}", "ht")
        )
        per_pop[p] = compute_counts_light(
            vp_ht, var_idx_ht, restricted, empty_heavy,
            max_join_partitions=max_join_partitions,
        ).checkpoint(hl.utils.new_temp_file(f"counts_pop_{p}", "ht"))

    # Assemble: one representative pair-keyed base, join each group's counts,
    # sum for "all". Each join / the array-sum chain is a small expression.
    base = per_pop[specific[0]].select()
    pp = {p: per_pop[p][base.key] for p in specific}
    raw = {p: pp[p].gt_counts_raw for p in specific}
    adj = {p: pp[p].gt_counts_adj for p in specific}

    def _sum(arrs):
        total = arrs[0]
        for a in arrs[1:]:
            total = hl.zip(total, a).map(lambda x: x[0] + x[1])
        return total

    all_raw = _sum([raw[p] for p in specific])
    all_adj = _sum([adj[p] for p in specific])
    by_pop = hl.dict(
        [(GLOBAL_POP, hl.struct(raw=all_raw, adj=all_adj))]
        + [(p, hl.struct(raw=raw[p], adj=adj[p])) for p in specific]
    )
    return base.annotate(
        gt_counts_raw=all_raw, gt_counts_adj=all_adj, gt_counts_by_pop=by_pop
    )


def main(args):
    """Create variant pair matrix from gnomAD v4 VDS."""
    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix
    data_type = args.data_type
    min_an_pct = args.min_an_pct
    scope_flags = [bool(args.gene), bool(args.test_genes), bool(args.interval)]
    if sum(scope_flags) > 1:
        raise ValueError(
            "--gene, --test-genes, and --interval are mutually exclusive."
        )
    test = args.test or any(scope_flags)
    if args.gene:
        test_intervals = {args.gene: TEST_INTERVALS[args.gene]}
    elif args.test_genes:
        wanted = [g.strip() for g in args.test_genes.split(",") if g.strip()]
        missing = [g for g in wanted if g not in TEST_INTERVALS]
        if missing:
            raise ValueError(
                f"--test-genes contains unknown TEST_INTERVALS keys: {missing}. "
                f"Available: {sorted(TEST_INTERVALS)}"
            )
        test_intervals = {g: TEST_INTERVALS[g] for g in wanted}
    elif args.interval:
        test_intervals = {args.interval: args.interval}
    else:
        test_intervals = TEST_INTERVALS
    # Normalize --test-chrom (e.g. "5" -> "chr5"); None when not set.
    test_chrom = args.test_chrom
    if test_chrom and not test_chrom.startswith("chr"):
        test_chrom = f"chr{test_chrom}"

    # Scoped output postfix (mirrors create_vp_list.py) so a --gene /
    # --test-genes / --interval run reads + writes the same `{gene}_test` postfix
    # the upstream create_vp_list step used. Without it a test run falls back to
    # the `pcnt_test` default and can't find the gene's pair list / intermediates.
    if test and output_postfix is None:
        if args.gene:
            output_postfix = f"{args.gene}_test"
        elif args.test_genes:
            output_postfix = "_".join(sorted(test_intervals)) + "_test"
        elif args.interval:
            output_postfix = (
                args.interval.replace(":", "_").replace("-", "_") + "_test"
            )
        else:
            output_postfix = "all_test"

    # --vds-subset: densify + encode + count one by-sample subset (non_ukb /
    # ukb.<group>) of the cohort. The count-group OUTPUTS are qualified with the
    # subset (via get_variant_pair_resources(subset=...) + count_output_dir
    # below) so per-subset runs are isolated; the pair-list INPUT stays the
    # full-dataset artifact. --merge-subset-counts later sums the per-subset
    # count HTs back to the full cohort.
    subset = args.vds_subset
    if subset is not None:
        valid_subsets = get_count_subsets(data_type)
        if subset not in valid_subsets:
            raise ValueError(
                f"--vds-subset {subset!r} is not a valid subset for {data_type}; "
                f"choose one of {valid_subsets}."
            )
        if data_type != "exomes":
            raise ValueError("--vds-subset is only defined for exomes.")
        if args.stratify_by_pop or args.pops:
            raise ValueError(
                "--vds-subset is incompatible with --stratify-by-pop / --pops "
                "(the subset is itself the stratification; merge sums full-cohort "
                "counts only)."
            )
        if args.emit_no_pbt_counts:
            raise ValueError(
                "--vds-subset is incompatible with --emit-no-pbt-counts "
                "(run the no-PBT subtraction on the merged full-cohort counts)."
            )

    hl.init(
        log=os.path.join(tempfile.gettempdir(), "compute_vp_counts.log"),
        tmp_dir=tmp_dir,
        backend=args.backend,
    )

    logger.info(
        f"""
        Running script with the following parameters:

            Data type: {data_type}
            Backend: {args.backend}
            Test: {test}
            Gene: {args.gene or 'all test intervals'}
            Output postfix: {output_postfix}
            Overwrite: {overwrite}
            Tmp dir: {tmp_dir}
        """
    )

    # Get variant co-occurrence pipeline resources. `subset` qualifies only the
    # count-group outputs; the pair-list input stays on the plain postfix.
    resources = get_variant_pair_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir if test else None,
        output_postfix=output_postfix,
        overwrite=overwrite,
        subset=subset,
    )
    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )

    # --- Genotype count steps (4 phases, can run on different clusters) ---
    # Subset-qualified postfix for the count intermediates dir (matches the
    # subset-qualified output resources above) so per-subset encodes/counts
    # never collide.
    count_postfix = (
        (f"{output_postfix}.{subset}" if output_postfix is not None else subset)
        if subset is not None
        else output_postfix
    )
    count_output_dir = f"{tmp_dir}/genotype_count_intermediates{_get_output_postfix(count_postfix, test)}"
    # --stratify-by-pop / --pops derive per-pop counts at count time by joining
    # the meta pop label to the (existing) encoded `samples` global — no
    # re-encode. --pops restricts to a subset of groups (and implies
    # stratification); otherwise all groups present are used.
    requested_pops = None
    if args.pops:
        requested_pops = [p.strip() for p in args.pops.split(",") if p.strip()]
        valid = set(get_pops(data_type))
        unknown = [p for p in requested_pops if p not in valid]
        if unknown:
            raise ValueError(
                f"--pops has unknown group(s) {unknown}; valid groups for "
                f"{data_type}: {sorted(valid)}"
            )
        # GLOBAL_POP ("all") is always included; keep only specific groups here.
        requested_pops = [p for p in requested_pops if p != GLOBAL_POP]
    stratify_by_pop = args.stratify_by_pop or requested_pops is not None
    pop_ht = get_sample_pop_ht(data_type) if stratify_by_pop else None

    if args.encode_genotypes:
        logger.info("Densifying pair-list variants and encoding genotypes...")
        res = resources.create_variant_pair_genotype_counts_ht
        counts_release_only = args.counts_release_only

        if test_chrom:
            logger.info(
                "Single-chromosome test mode: restricting --encode-genotypes "
                "to %s and reading the vp_list_ht from the test-chrom path.",
                test_chrom,
            )
            chrom_interval = hl.parse_locus_interval(
                test_chrom, reference_genome="GRCh38"
            )
            filter_intervals = [chrom_interval]
            vp_path = (
                f"{DEFAULT_TMP_DIR}/exomes.variant_pairs.{test_chrom}_test.ht"
            )
            vp_ht = hl.read_table(vp_path)
        elif test:
            res.check_resource_existence()
            filter_intervals = list(test_intervals.values())
            vp_ht = res.vp_list_ht.ht()
        else:
            res.check_resource_existence()
            filter_intervals = None
            vp_ht = res.vp_list_ht.ht()

        # Densify only the (AN-floored) pair-list variants (via the shared
        # densify_encode_input_mt: unsplit read + phase-keeping split + af /
        # fixed_homalt_model joins), then encode. The dense MT is transient —
        # checkpointed to scratch so the two encode passes (var_idx + gt sets)
        # don't recompute the densify — and never persisted. --min-an-pct is
        # stamped onto the encoded table's globals so the count steps can refuse
        # to lower it. Per-pop stratification is a count-time concern (the pop
        # label is joined to the encoded `samples` global at count time), so the
        # encoding is pop-agnostic and reusable across stratified / unstratified
        # counts.
        ht = create_variant_pair_filter_ht(
            filter_pairs_by_an_pct(vp_ht, min_an_pct)
        )
        # --vds-subset: read the by-sample subset VDS (ALL cohort samples, only
        # interval-restricted) and hand it to the shared densify, which restricts
        # to the release / high-quality cohort AFTER densify (a column filter, so
        # a variant monomorphic within the subset keeps its hom-ref row). Full
        # run: densify reads + release-filters the VDS via get_vds_func itself.
        subset_vds = None
        subset_restrict_ht = None
        if subset is not None:
            subset_vds = _read_count_subset_vds(
                subset, filter_intervals=filter_intervals,
            )
            meta_ht = meta(data_type=data_type).ht()
            subset_restrict_ht = meta_ht.filter(
                meta_ht.release if counts_release_only else meta_ht.high_quality
            )
        mt = densify_encode_input_mt(
            get_vds_func, ht, data_type,
            release_only=counts_release_only,
            filter_intervals=filter_intervals,
            vds=subset_vds,
            restrict_samples_ht=subset_restrict_ht,
        )
        mt = mt.checkpoint(
            hl.utils.new_temp_file("encode_genotypes.dense", "mt"),
            overwrite=True,
        )
        encode_genotypes(
            mt,
            vp_ht,
            output_dir=count_output_dir,
            min_an_pct=min_an_pct,
        )
        logger.info("Encoded genotypes written to %s", count_output_dir)

    # --- Excluded-genes HT: build it FIRST if the user passed
    # --exclude-gene-ids, since downstream steps read it. ---
    excluded_genes_res = resources.build_excluded_genes_ht.excluded_genes_ht
    excluded_genes_ht = None
    if args.exclude_gene_ids:
        gene_ids = [
            g.strip() for g in args.exclude_gene_ids.split(",") if g.strip()
        ]
        logger.info(
            "Writing excluded-genes HT (%d gene_ids) → %s",
            len(gene_ids), excluded_genes_res.path,
        )
        hl.Table.parallelize(
            [hl.struct(gene_id=g) for g in gene_ids],
            hl.tstruct(gene_id=hl.tstr),
            key=["gene_id"],
        ).write(excluded_genes_res.path, overwrite=overwrite)
    # If the file exists from a prior --exclude-gene-ids invocation, use it.
    try:
        excluded_genes_ht = excluded_genes_res.ht()
    except Exception:
        excluded_genes_ht = None

    if args.build_variant_size_info:
        logger.info("Building variant size-info HT...")
        res = resources.create_variant_pair_genotype_counts_ht
        size_info_res = resources.build_variant_size_info_ht.variant_size_info_ht

        encoded_gt_ht = hl.read_table(
            f"{count_output_dir}/encoded_gt_sets_by_var_idx.ht"
        )
        encoded_floor = _read_min_an_pct(encoded_gt_ht)
        _assert_min_an_pct_not_lowered(
            min_an_pct, encoded_floor, "encoded genotype intermediates"
        )
        var_idx_ht = hl.read_table(f"{count_output_dir}/var_idx.ht")
        vp_ht = filter_pairs_by_an_pct(res.vp_list_ht.ht(), min_an_pct)
        # gene_id source for the size-info HT — normally the variant
        # filter HT; --gene-ids-from-pair-list opts into deriving it
        # from the pair list (one-off fallback when the filter HT for
        # this run hasn't been built).
        if args.gene_ids_from_pair_list:
            logger.info(
                "Synthesizing (locus, alleles) → gene_id mapping from "
                "the variant pair list (--gene-ids-from-pair-list)..."
            )
            # Re-bind to local variables after each key_by — referring to
            # vp_ht.foo after vp_ht.key_by(...) trips Hail's source
            # tracker (per CLAUDE.md).
            keyed_v1 = vp_ht.key_by(
                locus=vp_ht.locus1, alleles=vp_ht.alleles1,
            )
            v1_genes = keyed_v1.select(_gene_id=keyed_v1.gene_id)
            keyed_v2 = vp_ht.key_by(
                locus=vp_ht.locus2, alleles=vp_ht.alleles2,
            )
            v2_genes = keyed_v2.select(_gene_id=keyed_v2.gene_id)
            per_site = v1_genes.union(v2_genes)
            variant_filter_ht = per_site.group_by(
                per_site.locus, per_site.alleles,
            ).aggregate(
                gene_id=hl.array(
                    hl.agg.explode(
                        lambda g: hl.agg.collect_as_set(g),
                        per_site._gene_id,
                    )
                )
            )
        else:
            try:
                variant_filter_ht = (
                    resources.create_variant_filter_ht.variant_filter_ht.ht()
                )
            except Exception as e:
                logger.warning(
                    "Variant filter HT not available (%s); size-info HT "
                    "will omit gene_id. Gene exclusion will be a no-op. "
                    "(Pass --gene-ids-from-pair-list to derive gene_id "
                    "from the pair list itself.)", e,
                )
                variant_filter_ht = None
        size_info_ht = build_variant_size_info_ht(
            _project_count_fields(encoded_gt_ht), vp_ht, var_idx_ht,
            variant_filter_ht=variant_filter_ht,
        )
        size_info_ht.write(size_info_res.path, overwrite=overwrite)
        logger.info("Variant size-info HT written to %s", size_info_res.path)

    if args.build_size_info_report:
        size_info_path = (
            resources.build_variant_size_info_ht.variant_size_info_ht.path
        )
        cutoff = (
            args.heavy_contribution_cutoff
            if args.heavy_contribution_cutoff is not None
            else TARGET_HEAVY_PARTITION_BYTES
        )
        # Co-locate the report with the size-info HT's parent so it's
        # easy to find when inspecting outputs from the same run.
        report_parent = os.path.dirname(size_info_path.rstrip("/"))
        logger.info(
            "Building size-info report under %s/size_info_report/ "
            "(cutoff=%.2f GB)...",
            report_parent, cutoff / 1024**3,
        )
        report_md = build_report(
            size_info_path=size_info_path,
            output_dir=report_parent,
            heavy_contribution_cutoff=cutoff,
        )
        logger.info("Size-info report written to %s", report_md)

    if args.compute_counts_light or args.compute_counts_heavy:
        # Shared input loading — both count steps need the same HTs, and these
        # functions are pure transforms (all I/O lives here in main()).
        res = resources.create_variant_pair_genotype_counts_ht
        size_info_path = (
            resources.build_variant_size_info_ht.variant_size_info_ht.path
        )
        var_idx_ht = hl.read_table(f"{count_output_dir}/var_idx.ht")
        encoded_gt_ht = hl.read_table(
            f"{count_output_dir}/encoded_gt_sets_by_var_idx.ht"
        )
        vp_ht = filter_pairs_by_an_pct(res.vp_list_ht.ht(), min_an_pct)
        _assert_min_an_pct_not_lowered(
            min_an_pct, _read_min_an_pct(encoded_gt_ht),
            "encoded genotype intermediates",
        )
        if excluded_genes_ht is not None:
            vp_ht = _filter_pairs_by_excluded_genes(
                vp_ht, var_idx_ht, hl.read_table(size_info_path), excluded_genes_ht,
            )

        # PBT-member set for the no-PBT counts (release \ PBT via subtraction,
        # in the same sweep as the full-release counts). Full-cohort only.
        pbt_members_ht = None
        if args.emit_no_pbt_counts:
            if stratify_by_pop:
                raise ValueError(
                    "--emit-no-pbt-counts is not supported with --stratify-by-pop "
                    "(no-PBT counts are full-cohort only)."
                )
            ped_resource = (
                pedigree(finalized=True) if args.trio_set == "pedigree" else trios()
            )
            trio_samples = complete_trio_samples(ped_resource.pedigree())
            pbt_members_ht = samples_ht(trio_samples)
            logger.info(
                "no-PBT counts enabled: excluding %d PBT members (trio_set=%s) via "
                "count-time subtraction.",
                len(trio_samples), args.trio_set,
            )

        if stratify_by_pop:
            # Per-pop: loop-restrict-per-group + co-partitioned flat count +
            # assemble dict (compute_counts_by_pop). Each group is counted over
            # its own small restricted encoding, so there is no light/heavy
            # split (no size_info needed) and no full-cohort shuffle. The
            # complete per-pop result is written to counts_light.ht (an empty
            # counts_heavy.ht keeps --combine-counts a no-op union).
            logger.info(
                "Computing per-pop counts (compute_counts_by_pop) for groups: %s",
                requested_pops if requested_pops is not None else "all present",
            )
            ht = compute_counts_by_pop(
                vp_ht, var_idx_ht, encoded_gt_ht, pop_ht, requested_pops,
            )
            ht.write(f"{count_output_dir}/counts_light.ht", overwrite=overwrite)
            if args.compute_counts_heavy:
                _empty_counts_ht(stratify_by_pop=True).write(
                    f"{count_output_dir}/counts_heavy.ht", overwrite=overwrite,
                )
            logger.info("Per-pop counts written.")
        else:
            # Physical-phase refinement of the AaBb cell (n_phased_cis /
            # n_phased_trans) — full-cohort only; the per-pop path above does
            # not carry it yet.
            emit_phase = args.emit_phase_counts
            cutoff = (
                args.heavy_contribution_cutoff
                if args.heavy_contribution_cutoff is not None
                else TARGET_HEAVY_PARTITION_BYTES
            )
            heavy_variants = _size_info_to_heavy_variants(
                hl.read_table(size_info_path),
                heavy_contribution_cutoff=cutoff,
                excluded_genes_ht=excluded_genes_ht,
            )
            # No variants heavy at this cutoff → count every pair via the
            # indexed-lookup plan; heavy = empty.
            n_heavy = heavy_variants.count()
            if n_heavy == 0:
                logger.info(
                    "No variants heavy at cutoff → count_all_pairs_via_index "
                    "for all pairs, empty heavy table."
                )
                if args.compute_counts_light:
                    ht = count_all_pairs_via_index(
                        vp_ht, var_idx_ht, encoded_gt_ht, emit_phase=emit_phase,
                    )
                    ht.write(
                        f"{count_output_dir}/counts_light.ht", overwrite=overwrite
                    )
                    logger.info("Light counts written.")
                if args.compute_counts_heavy:
                    _empty_counts_ht(emit_phase=emit_phase).write(
                        f"{count_output_dir}/counts_heavy.ht", overwrite=overwrite,
                    )
                    logger.info("Empty heavy counts written.")
            else:
                if args.compute_counts_light:
                    logger.info("Computing counts for light split...")
                    ht = compute_counts_light(
                        vp_ht=vp_ht,
                        var_idx_ht=var_idx_ht,
                        encoded_gt_ht=encoded_gt_ht,
                        heavy_variants=heavy_variants,
                        size_info_ht=hl.read_table(size_info_path),
                        emit_phase=emit_phase,
                    )
                    ht.write(
                        f"{count_output_dir}/counts_light.ht", overwrite=overwrite
                    )
                    logger.info("Light counts written.")

                if args.compute_counts_heavy:
                    logger.info("Computing counts for heavy split...")
                    ht = compute_counts_heavy(
                        vp_ht=vp_ht,
                        var_idx_ht=var_idx_ht,
                        encoded_gt_ht=encoded_gt_ht,
                        heavy_variants=heavy_variants,
                        emit_phase=emit_phase,
                    )
                    ht.write(
                        f"{count_output_dir}/counts_heavy.ht", overwrite=overwrite
                    )
                    logger.info("Heavy counts written.")

            # no-PBT: count the PBT∩cohort sub-population on a re-indexed
            # restriction of the SAME encode (small — no per-pair literal), and
            # write it to counts_pbt.ht so --combine-counts can subtract it from
            # the full-cohort counts (release \ PBT). Runs alongside the light
            # count (once). See restrict_encoded_to_samples / _subtract_pbt_counts.
            if args.emit_no_pbt_counts and args.compute_counts_light:
                logger.info(
                    "no-PBT: restricting the encode to PBT∩cohort and counting..."
                )
                pbt_encoded, n_pbt = restrict_encoded_to_samples(
                    encoded_gt_ht, pbt_members_ht
                )
                pbt_encoded = pbt_encoded.checkpoint(
                    f"{count_output_dir}/pbt_encoded.ht", overwrite=overwrite
                )
                pbt_counts = count_all_pairs_via_index(
                    vp_ht, var_idx_ht, pbt_encoded,
                )
                pbt_counts.annotate_globals(no_pbt_trio_set=args.trio_set).write(
                    f"{count_output_dir}/counts_pbt.ht", overwrite=overwrite
                )
                logger.info(
                    "no-PBT: PBT∩cohort counts written (%d samples).", n_pbt
                )

    if args.combine_counts:
        logger.info("Combining light + heavy counts...")
        res = resources.create_variant_pair_genotype_counts_ht

        light_path = f"{count_output_dir}/counts_light.ht"
        heavy_path = f"{count_output_dir}/counts_heavy.ht"
        tables = []
        for path in [light_path, heavy_path]:
            try:
                tables.append(hl.read_table(path))
            except Exception:
                logger.info("Skipping %s (not found).", path)
        if tables:
            # Light and heavy keep different pipeline-internal extras
            # (light: v1_idx,v2_idx; heavy: v_idx,_split_idx). Normalize
            # both to the common schema before union so the row types match.
            common = [
                "locus1", "alleles1", "locus2", "alleles2",
                "gt_counts_raw", "gt_counts_adj",
            ]
            # Preserve optional payloads (per-pop dict, physical-phase counts)
            # when every table carries them — light + heavy from one run share
            # a schema, so intersect across tables to stay safe.
            common += [
                c for c in ("gt_counts_by_pop", "n_phased_cis", "n_phased_trans")
                if all(c in t.row for t in tables)
            ]
            tables = [
                t.key_by().select(*common).key_by(
                    "locus1", "alleles1", "locus2", "alleles2"
                )
                for t in tables
            ]
            ht = tables[0] if len(tables) == 1 else tables[0].union(tables[1])

            # no-PBT: subtract the PBT∩cohort counts (counts_pbt.ht, from the
            # count step) to add gt_counts_*_no_pbt = release \ PBT, and carry the
            # no_pbt_trio_set stamp so trio_phasing can verify its --trio-set.
            try:
                pbt_ht = hl.read_table(f"{count_output_dir}/counts_pbt.ht")
            except Exception:
                pbt_ht = None
            if pbt_ht is not None:
                ht = _subtract_pbt_counts(ht, pbt_ht)
                ht = ht.annotate_globals(
                    no_pbt_trio_set=hl.eval(pbt_ht.index_globals().no_pbt_trio_set)
                )

            ht = ht.naive_coalesce(1000).checkpoint(res.vp_gt_counts_ht.path, overwrite=overwrite)
            logger.info("The variant pair genotype counts Table has been written...")

    if args.merge_subset_counts:
        # Sum the per-subset (--vds-subset) genotype-count HTs back to the full
        # cohort. Each subset's counts live at the subset-qualified path; the
        # merged full-cohort result is written to the plain (un-subset-qualified)
        # genotype-counts path — i.e. what a single full-cohort run would write.
        if args.merge_subsets in (None, "all"):
            merge_subsets = get_count_subsets(data_type)
        else:
            merge_subsets = [
                s.strip() for s in args.merge_subsets.split(",") if s.strip()
            ]
            valid_subsets = set(get_count_subsets(data_type))
            unknown = [s for s in merge_subsets if s not in valid_subsets]
            if unknown:
                raise ValueError(
                    f"--merge-subsets has unknown subset(s) {unknown}; valid "
                    f"subsets for {data_type}: {sorted(valid_subsets)}."
                )
        logger.info("Merging %d subset count HTs: %s", len(merge_subsets), merge_subsets)
        subset_hts = []
        for s in merge_subsets:
            s_postfix = (
                f"{output_postfix}.{s}" if output_postfix is not None else s
            )
            s_path = get_variant_pair_genotype_counts_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir if test else None,
                output_postfix=s_postfix,
            ).path
            logger.info("Reading subset %s counts from %s", s, s_path)
            subset_hts.append(hl.read_table(s_path))
        merged = merge_subset_counts(subset_hts)
        merged = merged.annotate_globals(merged_from_subsets=merge_subsets)
        out_path = get_variant_pair_genotype_counts_ht(
            data_type=data_type,
            test=test,
            tmp_dir=tmp_dir if test else None,
            output_postfix=output_postfix,
        ).path
        merged.naive_coalesce(1000).write(out_path, overwrite=overwrite)
        logger.info("Merged full-cohort genotype counts written to %s", out_path)

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
        "--backend",
        default="spark",
        choices=("spark", "batch"),
        help=(
            "Hail Query backend: 'spark' (default, uses local/Dataproc Spark) or "
            "'batch' (Hail Query-on-Batch). Use 'batch' to run on Hail Batch instead "
            "of locally or on Dataproc."
        ),
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to test intervals (all genes in TEST_INTERVALS) for testing.",
    )
    parser.add_argument(
        "--gene",
        choices=list(TEST_INTERVALS),
        help=(
            "Run on a single gene; uses that gene's interval from TEST_INTERVALS and "
            "implies --test."
        ),
    )
    parser.add_argument(
        "--test-genes",
        help=(
            "Comma-separated subset of TEST_INTERVALS keys (e.g. "
            "'SGCA,CAPN3,HFE'). Restricts the test intervals to just those "
            "genes; implies --test. Mutually exclusive with --gene."
        ),
    )
    parser.add_argument(
        "--interval",
        help=(
            "Explicit locus interval (e.g. 'chr19:1-58617616') for chromosome-"
            "scale scans. Bypasses the TEST_INTERVALS lookup; implies --test. "
            "Mutually exclusive with --gene and --test-genes."
        ),
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
            f"is {DEFAULT_DATA_TYPE}."
        ),
    )
    parser.add_argument(
        "--min-an-pct",
        type=int,
        default=0,
        help=(
            "Exclusive AN_percent floor required on both sides of a pair, "
            "applied when a pair list is consumed (dense-MT / encode / "
            "count steps). Default 0 drops only uncallable (an_pct==0) "
            "endpoints; raise to e.g. 10 for a stricter floor. Pass a "
            "negative value to keep every pair. Must be passed consistently "
            "across the steps that read the pair list."
        ),
    )
    parser.add_argument(
        "--counts-release-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Use release-only samples for the dense MT densified in genotype "
            "counting (--encode-genotypes). When False, uses high-quality "
            "samples instead. Default is True."
        ),
    )
    parser.add_argument(
        "--test-chrom",
        nargs="?",
        const="chr19",
        default=None,
        help=(
            "Restrict --encode-genotypes to a single chromosome (useful for "
            "full-genome runtime / cost estimation). Pass without a value to "
            "default to chr19, or specify e.g. '--test-chrom 5'. Reads the "
            "chrom-postfixed variant-pair list produced by create_vp_list.py; "
            "production locations are untouched."
        ),
    )
    parser.add_argument(
        "--encode-genotypes",
        action="store_true",
        help=(
            "Step A: Densify the pair-list variants out of the VDS (transient, "
            "not persisted) and encode them into per-variant sample sets. Run "
            "once; reuse for light/heavy with any threshold."
        ),
    )
    parser.add_argument(
        "--build-variant-size-info",
        action="store_true",
        help=(
            "Step A2: Build the per-variant size-info HT (one row per "
            "variant with _contribution, _cum_before, _bytes, "
            "split_count, locus, alleles, gene_id). Written to the "
            "variant_size_info resource path. Downstream "
            "--compute-counts-light and --compute-counts-heavy read this "
            "HT, applying the runtime --heavy-contribution-cutoff and "
            "--exclude-gene-ids filters. Inspect (e.g. exported to TSV "
            "or via --build-size-info-report) to retune the cutoff "
            "without rebuilding."
        ),
    )
    parser.add_argument(
        "--build-size-info-report",
        action="store_true",
        help=(
            "Step A2b: Generate a per-variant + per-gene report (figures "
            "+ stats.json + report.md) from the variant size-info HT. "
            "Uses --heavy-contribution-cutoff to decide is_heavy at "
            "report time. Output goes to "
            "{tmp_dir}/size_info_report{postfix}/."
        ),
    )
    parser.add_argument(
        "--heavy-contribution-cutoff",
        type=int,
        default=None,
        help=(
            "Bytes. Heavy iff _contribution >= this value. The size-info "
            "HT carries no heavy/light decision itself; the cutoff lives "
            "here so it can be retuned after inspecting the HT without a "
            "rebuild. Defaults to TARGET_HEAVY_PARTITION_BYTES (500 MB) "
            "— the tightest sensible value (anything past it needs "
            "split_count>1 anyway)."
        ),
    )
    parser.add_argument(
        "--exclude-gene-ids",
        type=str,
        default=None,
        help=(
            "Comma-separated Ensembl gene_ids to exclude from BOTH light "
            "and heavy counts. Writes a small excluded_genes.ht resource "
            "and applies it as a filter in both count steps (pairs are "
            "dropped when either side's gene_id array intersects this "
            "set). Use to defer chr19 mega-genes (MUC16, RYR1, ...) to "
            "per-gene jobs."
        ),
    )
    parser.add_argument(
        "--gene-ids-from-pair-list",
        action="store_true",
        help=(
            "One-off fallback for --build-variant-size-info: derive the "
            "(locus, alleles) → gene_id mapping from the variant pair "
            "list itself (which carries gene_id per pair) instead of "
            "the variant_filter HT. Use when the variant_filter HT for "
            "this run hasn't been built. The normal path is to use "
            "variant_filter_ht."
        ),
    )
    parser.add_argument(
        "--compute-counts-light",
        action="store_true",
        help=(
            "Step B: Compute genotype counts for light pairs, using the "
            "size-info HT from --build-variant-size-info plus any "
            "--heavy-contribution-cutoff and --exclude-gene-ids "
            "overrides. The light path auto-selects a fast indexed "
            "lookup when no heavy variants remain, and a co-partition + "
            "semi_join plan otherwise."
        ),
    )
    parser.add_argument(
        "--compute-counts-heavy",
        action="store_true",
        help=(
            "Step C: Compute genotype counts for heavy pairs (at least "
            "one variant in the heavy set after applying "
            "--heavy-contribution-cutoff). Returns an empty Table when "
            "no variants are heavy."
        ),
    )
    parser.add_argument(
        "--combine-counts",
        action="store_true",
        help="Step D: Union light + heavy counts into the final output Table.",
    )
    parser.add_argument(
        "--stratify-by-pop",
        action="store_true",
        help=(
            "On the count step (--compute-counts-light/heavy), also emit "
            "per-genetic-ancestry-group "
            "genotype counts (gt_counts_by_pop dict, keyed by pop incl. 'all') "
            "alongside the flat full-cohort gt_counts_raw/adj. The pop label "
            "(meta.population_inference.pop) is joined to the EXISTING encoded "
            "`samples` global at count time — no re-encode needed. The pop "
            "breakdown intersects each variant's sets with per-pop sample-index "
            "sets, so it multiplies per-pair set work by the number of groups. "
            "Use --pops to restrict to a subset of groups."
        ),
    )
    parser.add_argument(
        "--pops",
        help=(
            "Comma-separated genetic-ancestry groups to stratify by (e.g. "
            "'nfe,afr,eas'). Implies --stratify-by-pop and restricts the "
            "per-pop breakdown to these groups; 'all' is always included. "
            "Group names must be valid for the data type (see GEN_ANC_GROUPS); "
            "requested groups with no samples in the cohort are skipped. "
            "Default (with --stratify-by-pop, no --pops): all groups present."
        ),
    )
    parser.add_argument(
        "--emit-phase-counts",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "On the full-cohort count step (--compute-counts-light/heavy), also "
            "emit n_phased_cis / n_phased_trans — the physically-phased (PGT + "
            "PID) refinement of the adj double-het (AaBb) cell, matching the "
            "gnomAD MNV pipeline's same-haplotype call. Requires the encoded "
            "intermediates to carry `phased_het` (produced by --encode-genotypes "
            "here; empty on phase-less inputs). Not yet emitted on the per-pop "
            "(--stratify-by-pop) path. Default True."
        ),
    )
    parser.add_argument(
        "--emit-no-pbt-counts",
        action="store_true",
        help=(
            "On the full-cohort count step (--compute-counts-light), also emit "
            "gt_counts_raw_no_pbt / gt_counts_adj_no_pbt at --combine-counts — the "
            "full-cohort 9-cell counts with the PBT (trio) members subtracted out "
            "(release \\ PBT). Computed by restricting the SAME encode to "
            "PBT∩cohort (restrict_encoded_to_samples — re-indexed once, no "
            "per-pair literal), counting it to counts_pbt.ht, and subtracting "
            "element-wise at combine (exact for all 9 cells; disjoint-cohort "
            "additivity), so the trio-comparison pipeline needs no separate "
            "gnomAD-no-PBT densify / count pass. Full-cohort only (incompatible "
            "with --stratify-by-pop)."
        ),
    )
    parser.add_argument(
        "--trio-set",
        default="pedigree",
        choices=("pedigree", "trios"),
        help=(
            "Which gnomad_qc v4 pedigree resource defines the PBT members for "
            "--emit-no-pbt-counts: 'pedigree' (all trios) or 'trios' (one random "
            "trio per family). Stamped onto the output as the `no_pbt_trio_set` "
            "global and enforced by trio_phasing.py (which refuses a comparison "
            "whose --trio-set differs). Default 'pedigree'. Ignored without "
            "--emit-no-pbt-counts."
        ),
    )

    parser.add_argument(
        "--vds-subset",
        default=None,
        help=(
            "Densify + encode + count a single by-sample subset of the cohort "
            "(one of: non_ukb, ukb.<group> — see resources.get_count_subsets) "
            "instead of the full VDS. Reads the subset VDS by path, applies the "
            "same release/high-quality + interval filtering as the full run, and "
            "writes subset-qualified count-group outputs "
            "(genotype_count_intermediates.{postfix}.{subset}/, "
            "variant_pairs.genotype_counts.{postfix}.{subset}.ht). The pair-list "
            "input stays the full-dataset artifact. Run once per subset (on "
            "smaller clusters), then --merge-subset-counts to sum back to the "
            "full cohort. Exomes only; incompatible with --stratify-by-pop / "
            "--pops / --emit-no-pbt-counts."
        ),
    )
    parser.add_argument(
        "--merge-subset-counts",
        action="store_true",
        help=(
            "Sum the per-subset (--vds-subset) genotype-count HTs element-wise "
            "into the full-cohort counts. Reads each subset's "
            "variant_pairs.genotype_counts.{postfix}.{subset}.ht and writes the "
            "full-cohort variant_pairs.genotype_counts.{postfix}.ht (what a "
            "single full-cohort run would produce). Exact for all 9 cells "
            "(disjoint-cohort additivity). Use --merge-subsets to pick which "
            "subsets; default is all."
        ),
    )
    parser.add_argument(
        "--merge-subsets",
        default=None,
        help=(
            "Comma-separated subsets to merge with --merge-subset-counts (e.g. "
            "'non_ukb,ukb.nfe,ukb.afr'), or 'all' (default) for every subset in "
            "resources.get_count_subsets."
        ),
    )

    args = parser.parse_args()

    # Copy the Hail log off the driver's local disk to GCS so a finished (or
    # crashed) run can be inspected after the cluster is torn down — mirrors
    # gnomad_qc's get_logging_path / copy_log idiom.
    _log_label = (
        args.output_postfix or args.gene or args.test_genes or args.interval or "all"
    )
    _log_label = _log_label.replace(":", "_").replace("-", "_").replace(",", "_")
    _gcs_log_path = os.path.join(
        args.tmp_dir, "logs", f"compute_vp_counts.{_log_label}.log"
    )
    try:
        main(args)
    finally:
        try:
            logger.info("Copying Hail log to %s", _gcs_log_path)
            hl.copy_log(_gcs_log_path)
        except Exception as e:  # noqa: BLE001
            logger.warning("Could not copy Hail log to %s: %s", _gcs_log_path, e)

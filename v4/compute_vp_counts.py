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
   as the gnomAD MNV pipeline.

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
from gnomad.utils.annotations import get_adj_expr
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds

from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_TMP_DIR,
    GLOBAL_POP,
    TEST_INTERVALS,
    _get_output_postfix,
    get_pops,
    get_sample_pop_ht,
    get_variant_filter_ht,
    get_variant_pair_resources,
)
from gnomad_chets.v4.size_info_report import build_report
from gnomad_chets.v4.utils import (
    calculate_partitions_by_size,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("compute_vp_counts")
logger.setLevel(logging.INFO)


def filter_pairs_by_an_pct(ht: hl.Table, min_an_pct: int) -> hl.Table:
    """Drop pairs whose AN_percent is at or below ``min_an_pct`` on either side.

    With ~no callable samples at a locus (``an_pct == 0``) the AABB cell
    (``min(n_callable_v1, n_callable_v2)``) collapses and the haplotype EM
    degenerates, so those pairs carry no co-occurrence signal; higher floors
    additionally trade power for AN quality. Applied at pair-list
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
    :return: Split variant-data MT keyed by ``(locus, alleles)``.
    """
    split_entries = [
        "L" + e if e in _LOCAL_ENTRY_REMAP else e
        for e in (entries_to_keep + ["LA"])
    ]
    variant_mt = variant_mt.select_entries(*split_entries)
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


_ENCODED_SET_FIELDS = (
    "all_samples", "raw_hr_adj_missing",
    "raw_het", "raw_hv", "adj_het", "adj_hv",
    "raw_callable", "adj_callable",
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
    shuffle disk). Complement flags are preserved — a set stored as ``N∖A``
    becomes ``keep∖A_keep`` within the kept cohort, so the ``_count_from_sets``
    decode identity still holds with ``n_samples = n_keep``.

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
    n_keep = len(keep)
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

    def _n_pos(field, is_comp_field):
        # Positive within-kept size from the re-indexed stored set.
        return hl.if_else(
            e[is_comp_field], n_keep - reidx[field].length(), reidx[field].length()
        )

    n_a = _n_pos("all_samples", "all_samples_is_complement")
    n_f = _n_pos("raw_hr_adj_missing", "raw_hr_adj_missing_is_complement")
    out = e.annotate(
        **reidx,
        n_with_data=hl.int32(n_a + n_f),
        n_raw_hr_adj_missing=hl.int32(n_f),
        n_raw_callable=hl.int32(_n_pos("raw_callable", "raw_callable_is_complement")),
        n_adj_callable=hl.int32(_n_pos("adj_callable", "adj_callable_is_complement")),
    )
    return out.annotate_globals(samples=hl.literal(keep_samples, samples_dtype))


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

    ``all_samples`` in **primary form** (use_complement=False) stores
    cats 1, 3-6 and is disjoint from ``raw_hr_adj_missing`` (cat 2). In
    **complement form** (use_complement=True), ``all_samples`` stores
    the proper complement of A_pos = cats 1, 3-6, which equals
    ``cat 2 ∪ cat 7`` — i.e. it overlaps with ``raw_hr_adj_missing``.
    The decoder identity ``|pos ∩ A_pos| = |pos| − |pos ∩ stored|`` only
    holds when ``stored`` is the proper complement; storing only cat 7
    silently over-counts ``|pos ∩ A_pos|`` by ``|pos ∩ cat 2|``, which
    drives downstream cells negative on pairs where ``cat 2`` is
    non-trivial.

    The count algebra in :func:`_count_from_sets` reconstructs
    ``D'_v = A_v ∪ F_v`` (= cats 1-6) by summing the four cross terms;
    each decode uses the proper-complement identity above.

    Two of the per-category sets can still be large at low-coverage
    variants (raw-hom-ref carriers below adj threshold), so two
    complement-form stored sets save space; each carries a flag + size
    so downstream can resolve without materializing the complement.

    Output schema:

        ----------------------------------------
        Global fields:
            'samples': array<struct { s: str }>
        ----------------------------------------
        Row fields:
            'v_idx': int64
            'all_samples': set<int32>      # primary: cats 1, 3-6
                                           # complement: cats 2 ∪ cat 7
                                           # (proper complement of A_pos in N).
                                           # Disjoint from raw_hr_adj_missing
                                           # only in primary form.
            'all_samples_is_complement': bool
            'n_with_data': int32           # |cats 1-6| = samples NOT in
                                           # adj-PASS-0/0 majority
                                           # (= |all_samples positive|
                                           # + |raw_hr_adj_missing|)
            'raw_hr_adj_missing': set<int32>  # smaller of (cat 2) or its
                                           # complement
            'raw_hr_adj_missing_is_complement': bool
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
    adj_gt_count_expr = hl.if_else(
        adj_pass_expr, gt_count_expr, 0, missing_false=True
    )

    # Physical (read-backed) phase, when the dense MT carries it (PGT + PID
    # kept via _split_variant_data_keeping_phase). For a phased het-ref call
    # we record the phase-set id (PID) and which haplotype carries the alt
    # (PGT[0]); at count time two het variants of the same sample are cis iff
    # they share a PID and agree on PGT[0]. Reference-block (hom-ref) samples
    # have no PGT/PID, so phase is missing for them. Absent on MTs without
    # phase (e.g. the exploded PBT trio MT) → an empty phased_het dict.
    has_phase = "PGT" in mt.entry and "PID" in mt.entry
    select_entry_exprs = dict(raw_gt=gt_count_expr, adj_gt=adj_gt_count_expr)
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
    # (cat 2) so the same sample never appears in both. Complement form
    # is still cat 7 (adj-PASS-0/0 majority). raw_hr_adj_missing (cat 2)
    # gets independent complement-form treatment because at low-coverage
    # variants it can be large.
    gt = hl.enumerate(ht._entries)
    # ``not_adj_hom_ref`` = cats 1-6 (used only to derive n_with_data;
    # not stored). The set actually stored as ``all_samples`` is the
    # cats-1,3-6 subset below, which is disjoint from
    # ``raw_hr_adj_missing`` (cat 2) so the same sample never appears
    # in both stored sets.
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
    adj_hom_ref = gt.filter(
        lambda x: hl.is_defined(x[1])
        & hl.is_missing(x[1].raw_gt)
        & hl.is_missing(x[1].adj_gt)
    )
    raw_hr_adj_missing = gt.filter(
        lambda x: hl.is_defined(x[1])
        & hl.is_missing(x[1].raw_gt)
        & (x[1].adj_gt == 0)
    )
    raw_hr_adj_missing_complement = gt.filter(
        lambda x: hl.is_missing(x[1])
        | hl.is_defined(x[1].raw_gt)
        | hl.is_missing(x[1].adj_gt)
        | (x[1].adj_gt != 0)
    )
    n_with = not_adj_hom_ref.length()
    n_raw_hr_adj_missing = raw_hr_adj_missing.length()
    # use_complement compares the positive form (cats 1, 3-6) to its
    # proper complement in N. Since A_pos = cats 1, 3-6, the complement
    # of A_pos in N is cats 2 ∪ cat 7 = raw_hr_adj_missing ∪ adj_hom_ref.
    # So the complement-form storage size is |cat 2| + |cat 7|.
    n_all_samples_positive = n_with - n_raw_hr_adj_missing  # |cats 1, 3-6|
    use_complement = n_all_samples_positive > (
        adj_hom_ref.length() + n_raw_hr_adj_missing
    )
    F_use_complement = n_raw_hr_adj_missing > (
        raw_hr_adj_missing_complement.length()
    )
    # Callability sets for exact co-callable hom-ref counts in the
    # per-sample path. raw-callable = cats 2-7 (any GT call = present
    # entry); adj-callable = cats 4,6,7 (adj-PASS). Each stored as the
    # smaller of itself / its complement. ``implicit_homref`` records which
    # bulk category (cat 7 = 0/0 adj-PASS, or cat 1 = no-entry) is the
    # omitted majority in the transpose's gt_info, so the per-sample step
    # can resolve an "absent" partner side.
    n_total = hl.len(ht._entries)
    present = gt.filter(lambda x: hl.is_defined(x[1]))
    n_present = present.length()
    n_cat1 = n_total - n_present
    cat1_samples = gt.filter(lambda x: hl.is_missing(x[1]))
    raw_callable_use_complement = n_present > n_cat1
    adj_callable = present.filter(
        lambda x: hl.is_missing(x[1].adj_gt) | (x[1].adj_gt != 0)
    )
    n_adj_callable = adj_callable.length()
    adj_uncallable = gt.filter(
        lambda x: hl.is_missing(x[1])
        | (hl.is_defined(x[1].adj_gt) & (x[1].adj_gt == 0))
    )
    adj_callable_use_complement = n_adj_callable > (n_total - n_adj_callable)
    implicit_homref = adj_hom_ref.length() >= n_cat1
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
        # Stored set:
        #   - primary (use_complement=False): cats 1, 3-6 (disjoint from
        #     raw_hr_adj_missing / cat 2).
        #   - complement (use_complement=True): the proper complement of
        #     A_pos in N = cats 2 ∪ cat 7 = adj_hom_ref ∪ raw_hr_adj_missing.
        #
        # Storing the proper complement is required for the decoder's
        # ``|pos ∩ A_pos| = |pos| − |pos ∩ stored|`` identity to be
        # correct. Storing only cat 7 (the old buggy form) silently
        # over-counts ``|pos ∩ A_pos|`` by ``|pos ∩ cat 2|``, which
        # produces negative-valued cells downstream when ``cat 2`` is
        # non-trivial.
        #
        # In complement form, ``all_samples`` and ``raw_hr_adj_missing``
        # are NO LONGER disjoint — cat 2 is in both. The dedup-storage
        # benefit (saved bytes) only applies to the primary-form branch.
        # In complement form the additional |cat 2| samples don't
        # meaningfully grow storage because cat 2 is small whenever the
        # complement branch was chosen.
        all_samples=hl.if_else(
            use_complement,
            hl.set(adj_hom_ref.map(lambda x: x[0])).union(
                hl.set(raw_hr_adj_missing.map(lambda x: x[0]))
            ),
            hl.set(not_adj_hom_ref_no_F.map(lambda x: x[0])),
        ),
        all_samples_is_complement=use_complement,
        n_with_data=hl.int32(n_with),
        raw_hr_adj_missing=hl.if_else(
            F_use_complement,
            hl.set(raw_hr_adj_missing_complement.map(lambda x: x[0])),
            hl.set(raw_hr_adj_missing.map(lambda x: x[0])),
        ),
        raw_hr_adj_missing_is_complement=F_use_complement,
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
        implicit_homref=implicit_homref,
        raw_callable=hl.if_else(
            raw_callable_use_complement,
            hl.set(cat1_samples.map(lambda x: x[0])),
            hl.set(present.map(lambda x: x[0])),
        ),
        raw_callable_is_complement=raw_callable_use_complement,
        n_raw_callable=hl.int32(n_present),
        adj_callable=hl.if_else(
            adj_callable_use_complement,
            hl.set(adj_uncallable.map(lambda x: x[0])),
            hl.set(adj_callable.map(lambda x: x[0])),
        ),
        adj_callable_is_complement=adj_callable_use_complement,
        n_adj_callable=hl.int32(n_adj_callable),
    )

    # Rekey by var_idx and drop the locus/alleles fields to shrink row size.
    ht = ht.annotate(v_idx=var_idx_ht[ht.locus, ht.alleles].var_idx)
    return ht.key_by("v_idx").drop("locus", "alleles")



def _count_from_sets(
    v1_het: hl.expr.SetExpression,
    v1_hv: hl.expr.SetExpression,
    v1_all: hl.expr.SetExpression,
    v1_n: hl.expr.Int32Expression,
    v1_is_complement: hl.expr.BooleanExpression,
    v1_F: hl.expr.SetExpression,
    v1_n_F: hl.expr.Int32Expression,
    v1_F_is_complement: hl.expr.BooleanExpression,
    v2_het: hl.expr.SetExpression,
    v2_hv: hl.expr.SetExpression,
    v2_all: hl.expr.SetExpression,
    v2_n: hl.expr.Int32Expression,
    v2_is_complement: hl.expr.BooleanExpression,
    v2_F: hl.expr.SetExpression,
    v2_n_F: hl.expr.Int32Expression,
    v2_F_is_complement: hl.expr.BooleanExpression,
    n_samples: hl.expr.Int32Expression,
    *,
    include_raw_hr_adj_missing: bool,
) -> hl.expr.ArrayExpression:
    """
    Compute 9-element genotype count array from per-variant sample sets.

    Per-variant storage:

      - ``v_all`` stores either ``A_v`` (= cats 1, 3-6) in primary form,
        or the proper complement ``N − A_v = cats 2 ∪ cat 7`` in
        complement form. The decoder identity
        ``|pos ∩ A_v| = |pos| − |pos ∩ stored|`` requires storing the
        proper complement; storing only cat 7 silently over-counts by
        ``|pos ∩ cat 2|`` (drives downstream cells negative). In primary
        form ``v_all`` is disjoint from ``F_v``; in complement form
        cat 2 is in both ``v_all`` and ``F_v``.
        ``v_n`` is always ``|cats 1-6| = |A_v| + |F_v|`` — i.e.
        n_with_data, the "samples with any data" count. The positive
        A_v size is then ``v_n - v_n_F``.
      - ``v_F`` stores either ``raw_hr_adj_missing`` (cat 2) or its
        complement. ``v_n_F`` is always |cat 2|.

    The "real" not-adj-hom-ref set is ``D'_v = A_v ∪ F_v = cats 1-6``;
    we reconstruct intersections over D' from intersections over A and
    F at read time.

    Two modes (Python-level branch via ``include_raw_hr_adj_missing``):

      - ``False`` (ADJ cells): hom-ref = adj-PASS-0/0 = N \\ D'. Pass
        ``v_het = adj_het``, ``v_hv = adj_hv``.
      - ``True`` (RAW cells): hom-ref = raw-0/0 = (adj-PASS-0/0) ∪ F.
        Pass ``v_het = raw_het`` (all het), ``v_hv = raw_hv`` (all hv).
        F adds the adj-fail-hom-ref samples to the hom-ref pool.

    Count array layout: ``[AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]``
    where A/a = v1 ref/alt, B/b = v2 ref/alt. All branching on Hail
    complement flags operates on ints (set lengths and intersection
    lengths), never on sets themselves.

    :param v1_het, v1_hv: Sample sets for v1 het / hom-var.
    :param v1_all, v1_n, v1_is_complement: v1 ``all_samples`` triple.
    :param v1_F, v1_n_F, v1_F_is_complement: v1 ``raw_hr_adj_missing`` triple.
    :param v2_...: same for v2.
    :param n_samples: Total number of samples in the cohort.
    :param include_raw_hr_adj_missing: ``True`` for raw cells, ``False``
        for adj cells.
    :return: 9-element genotype count array.
    """
    # |A_pos ∩ B_pos| with optional complement-form storage. See the
    # docstring for the case derivation.
    def _isect_pos(a, a_n_pos, a_is_comp, b, b_n_pos, b_is_comp):
        raw = a.intersection(b).length()
        return (
            hl.case()
            .when(~a_is_comp & ~b_is_comp, raw)
            .when(~a_is_comp & b_is_comp, a_n_pos - raw)
            .when(a_is_comp & ~b_is_comp, b_n_pos - raw)
            .default(a_n_pos + b_n_pos - n_samples + raw)
        )

    # |positive_set ∩ S| where positive_set is stored as-is and S is
    # complement-aware: returns |pos ∩ S_pos|.
    def _pos_in_comp_aware(pos_set, comp_aware, comp_aware_is_comp):
        raw = pos_set.intersection(comp_aware).length()
        return hl.if_else(
            comp_aware_is_comp, pos_set.length() - raw, raw,
        )

    # Positive-form sizes for the disjoint A_v and F_v sets:
    #   v_all in primary form contains cats 1, 3-6, with size
    #     v_n - v_n_F (= |cats 1-6| - |cat 2|).
    #   v_F  in primary form contains cat 2, with size v_n_F.
    v1_n_A = v1_n - v1_n_F
    v2_n_A = v2_n - v2_n_F

    # D'_v = A_v ∪ F_v (disjoint per variant). Decompose
    # |D'_v1 ∩ D'_v2| into the four cross terms; reuse the A-F and F-F
    # pieces in the raw-cells branch below.
    a1_isect_a2 = _isect_pos(
        v1_all, v1_n_A, v1_is_complement,
        v2_all, v2_n_A, v2_is_complement,
    )
    a1_isect_f2 = _isect_pos(
        v1_all, v1_n_A, v1_is_complement,
        v2_F, v2_n_F, v2_F_is_complement,
    )
    f1_isect_a2 = _isect_pos(
        v1_F, v1_n_F, v1_F_is_complement,
        v2_all, v2_n_A, v2_is_complement,
    )
    f1_isect_f2 = _isect_pos(
        v1_F, v1_n_F, v1_F_is_complement,
        v2_F, v2_n_F, v2_F_is_complement,
    )
    d1_isect_d2 = a1_isect_a2 + a1_isect_f2 + f1_isect_a2 + f1_isect_f2

    # |H_v1 ∩ H_v2| where H = N \ D'  (adj-PASS-0/0 at both).
    #   = N - |D'_v1 ∪ D'_v2| = N - |D'_v1| - |D'_v2| + |D'_v1 ∩ D'_v2|
    # v_n = |D'_v| still — that semantics is unchanged.
    h1_isect_h2 = n_samples - v1_n - v2_n + d1_isect_d2

    # Carrier-carrier cells (het/hv sets always positive form).
    het_het = v1_het.intersection(v2_het).length()
    het_hv = v1_het.intersection(v2_hv).length()
    hv_het = v1_hv.intersection(v2_het).length()
    hv_hv = v1_hv.intersection(v2_hv).length()

    # Edge cells (adj component): |carrier_v ∩ H_other|
    #   = |carrier_v| - |carrier_v ∩ D'_other|
    # where |carrier_v ∩ D'_other| = |carrier_v ∩ A_other|
    #                              + |carrier_v ∩ F_other|.
    v1_het_in_d2 = (
        _pos_in_comp_aware(v1_het, v2_all, v2_is_complement)
        + _pos_in_comp_aware(v1_het, v2_F, v2_F_is_complement)
    )
    v1_hv_in_d2 = (
        _pos_in_comp_aware(v1_hv, v2_all, v2_is_complement)
        + _pos_in_comp_aware(v1_hv, v2_F, v2_F_is_complement)
    )
    v2_het_in_d1 = (
        _pos_in_comp_aware(v2_het, v1_all, v1_is_complement)
        + _pos_in_comp_aware(v2_het, v1_F, v1_F_is_complement)
    )
    v2_hv_in_d1 = (
        _pos_in_comp_aware(v2_hv, v1_all, v1_is_complement)
        + _pos_in_comp_aware(v2_hv, v1_F, v1_F_is_complement)
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
        v1_het_in_f2 = _pos_in_comp_aware(v1_het, v2_F, v2_F_is_complement)
        v1_hv_in_f2 = _pos_in_comp_aware(v1_hv, v2_F, v2_F_is_complement)
        v2_het_in_f1 = _pos_in_comp_aware(v2_het, v1_F, v1_F_is_complement)
        v2_hv_in_f1 = _pos_in_comp_aware(v2_hv, v1_F, v1_F_is_complement)

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


def _pop_restrict_variant(v, het, hv, pop_set, pop_size):
    """Restrict one variant's :func:`_count_from_sets` inputs to ``pop_set``.

    ``het`` / ``hv`` are the raw *or* adj carrier sets (passed explicitly so the
    same helper serves both counts). Returns the 8-tuple of pop-restricted args
    in the order :func:`_count_from_sets` consumes per variant:
    ``(het, hv, all, n_with_data, all_is_complement, F, n_F, F_is_complement)``.

    All sets are intersected with the pop's sample-index set; the complement
    flags are unchanged (a set stored as ``N∖A`` restricted to ``pop`` becomes
    ``pop∖A_pop``, still the complement — now within ``pop``). The positive
    within-pop sizes use the proper-complement identity
    ``|A ∩ pop| = |pop| − |stored_complement ∩ pop|``, matching how
    :func:`_count_from_sets` interprets a complement-form set against
    ``n_samples = pop_size``.
    """
    all_p = v.all_samples.intersection(pop_set)
    f_p = v.raw_hr_adj_missing.intersection(pop_set)
    n_a_p = hl.if_else(
        v.all_samples_is_complement, pop_size - all_p.length(), all_p.length()
    )
    n_f_p = hl.if_else(
        v.raw_hr_adj_missing_is_complement,
        pop_size - f_p.length(),
        f_p.length(),
    )
    return (
        het.intersection(pop_set),
        hv.intersection(pop_set),
        all_p,
        hl.int32(n_a_p + n_f_p),
        v.all_samples_is_complement,
        f_p,
        hl.int32(n_f_p),
        v.raw_hr_adj_missing_is_complement,
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
            *_pop_restrict_variant(v1, v1.raw_het, v1.raw_hv, ps, sz),
            *_pop_restrict_variant(v2, v2.raw_het, v2.raw_hv, ps, sz),
            sz,
            include_raw_hr_adj_missing=True,
        )
        adj = _count_from_sets(
            *_pop_restrict_variant(v1, v1.adj_het, v1.adj_hv, ps, sz),
            *_pop_restrict_variant(v2, v2.adj_het, v2.adj_hv, ps, sz),
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

    :param vp_subset: Pair Table with ``v1_idx, v2_idx, locus1, alleles1,
        locus2, alleles2``.
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
            pairs=hl.agg.collect(
                hl.struct(
                    v2_idx=vp_subset.v2_idx,
                    _v2_split_idx=vp_subset._v2_split_idx,
                    locus1=vp_subset.locus1,
                    alleles1=vp_subset.alleles1,
                    locus2=vp_subset.locus2,
                    alleles2=vp_subset.alleles2,
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
        locus1=vp_exploded.pairs.locus1,
        alleles1=vp_exploded.pairs.alleles1,
        locus2=vp_exploded.pairs.locus2,
        alleles2=vp_exploded.pairs.alleles2,
    )
    vp_exploded = vp_exploded.key_by().drop("v_idx", "_split_idx").cache()
    vp_exploded = vp_exploded.key_by(
        v_idx=vp_exploded.v2_idx,
        _split_idx=vp_exploded._v2_split_idx,
    ).drop("v2_idx", "_v2_split_idx")

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
        v1.all_samples_is_complement,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v1.raw_hr_adj_missing_is_complement,
        v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
        v2.all_samples_is_complement,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        v2.raw_hr_adj_missing_is_complement,
        n_samples,
        include_raw_hr_adj_missing=True,
    )
    gt_counts_adj = _count_from_sets(
        v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
        v1.all_samples_is_complement,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v1.raw_hr_adj_missing_is_complement,
        v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
        v2.all_samples_is_complement,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        v2.raw_hr_adj_missing_is_complement,
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
        "locus1",
        "alleles1",
        "locus2",
        "alleles2",
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
heavy step."""

_RESUME_LIGHT_GT_PATH: Optional[str] = None
_RESUME_LIGHT_VP_PATH: Optional[str] = None
"""ONE-SHOT resume hooks for :func:`compute_counts_light`. Both intermediates
(``gt_light``, ``vp_light``) are written by the light step just before the
co-partition + zip-join. When these are set, the function skips the
v_idx-annotation, heavy-filter, ``semi_join``, and both writes; it reads
the two paths directly and jumps to co-partitioning. Reset to ``None`` after
the resume run completes — stale paths from a different postfix would inject
the wrong intermediates."""


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
    "all_samples", "n_with_data", "all_samples_is_complement",
    "raw_hr_adj_missing", "n_raw_hr_adj_missing",
    "raw_hr_adj_missing_is_complement",
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
    """Drop encoder fields the ``_count_from_sets`` paths don't read.

    The per-sample rework added ``implicit_homref`` + raw/adj callable sets
    to the encoded table; for low-coverage variants those sets can be large
    (~315k integers per row). The light/heavy paths don't use them but
    would otherwise carry them through every shuffle. Projecting them away
    here is purely a row-size optimisation.
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
    ``n_phased_trans`` so it unions with a phase-annotated light result.
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

    When ``pops`` is provided (with ``pop_index_sets`` / ``pop_sizes`` from
    :func:`_build_pop_stratification`), additionally emits a
    ``gt_counts_by_pop`` dict<pop, struct{raw, adj}> alongside the flat
    ``gt_counts_raw`` / ``gt_counts_adj`` (which remain the full-cohort counts).
    """
    v1 = encoded_gt_ht[vp.v1_idx]
    v2 = encoded_gt_ht[vp.v2_idx]
    gt_counts_raw = _count_from_sets(
        v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
        v1.all_samples_is_complement,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v1.raw_hr_adj_missing_is_complement,
        v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
        v2.all_samples_is_complement,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        v2.raw_hr_adj_missing_is_complement,
        n_samples,
        include_raw_hr_adj_missing=True,
    )
    gt_counts_adj = _count_from_sets(
        v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
        v1.all_samples_is_complement,
        v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
        v1.raw_hr_adj_missing_is_complement,
        v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
        v2.all_samples_is_complement,
        v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
        v2.raw_hr_adj_missing_is_complement,
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
    vp = vp.select(
        "locus1", "alleles1", "locus2", "alleles2", **count_fields
    ).cache()
    return vp.key_by("locus1", "alleles1", "locus2", "alleles2").cache()


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
    return _count_pairs_via_index(
        vp_ht, encoded_gt_ht, n_samples, emit_phase=emit_phase, **pop_kwargs
    )


def compute_counts_light(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    encoded_gt_ht: hl.Table,
    heavy_variants: hl.Table,
    max_join_partitions: int = 10000,
    *,
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
        vp_light = vp_light.key_by("v1_idx", "v2_idx").select(
            "locus1", "alleles1", "locus2", "alleles2"
        ).cache()

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

    # Co-partition on v1_idx for the v1 zip-join. Weight ``n_with_data`` by
    # per-variant v1 pair degree so a genomic cluster of high-degree variants
    # gets subdivided across partitions instead of collapsing into one hot
    # spot. Without this weight, a partition with 6k variants of moderate
    # ``n_with_data`` but very high pair counts holds ~30% of the total
    # work and OOMs the executor even though its GT-set storage looks
    # balanced by the plain ``n_with_data`` metric.
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

    return _count_pairs_via_index(
        vp_light, gt_light, n_samples, emit_phase=emit_phase, **pop_kwargs
    )


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
    vp_heavy = vp_heavy.select(
        "v1_idx", "v2_idx", "locus1", "alleles1", "locus2", "alleles2",
    ).cache()

    result = _compute_counts_for_subset(
        vp_heavy, encoded_gt_ht, n_samples, "heavy", n_partitions,
        heavy_variants, emit_phase=emit_phase, **pop_kwargs,
    )
    return result.key_by("locus1", "alleles1", "locus2", "alleles2")


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
    :param encoded_gt_ht: Encoded GT table (pre-projection; complement form must
        be the PROPER complement — repair a buggy encoding first).
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

    # Get variant co-occurrence pipeline resources.
    resources = get_variant_pair_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir if test else None,
        output_postfix=output_postfix,
        overwrite=overwrite,
    )
    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )

    # --- Genotype count steps (4 phases, can run on different clusters) ---
    count_output_dir = f"{tmp_dir}/genotype_count_intermediates{_get_output_postfix(output_postfix, test)}"
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

        # Densify only the (AN-floored) pair-list variants straight out of the
        # VDS, then encode. The dense MT is transient — checkpointed to scratch
        # so the two encode passes (var_idx + gt sets) don't recompute the
        # densify — and never persisted. --min-an-pct is stamped onto the
        # encoded table's globals so the count steps can refuse to lower it.
        # Per-pop stratification is a count-time concern (the pop label is
        # joined to the encoded `samples` global at count time), so the encoding
        # is pop-agnostic and reusable across stratified / unstratified counts.
        ht = create_variant_pair_filter_ht(
            filter_pairs_by_an_pct(vp_ht, min_an_pct)
        )
        # Read unsplit so the local phased-GT field (LPGT) — which the
        # gnomad_qc split loader drops — survives. We split it ourselves with
        # _split_variant_data_keeping_phase so physical (read-backed) phase
        # (PGT + PID) reaches the encoder. get_gnomad_v4_vds rejects
        # filter_variant_ht on unsplit reads, so the variant restriction is
        # applied inside the helper instead.
        vds = get_vds_func(
            release_only=counts_release_only,
            high_quality_only=not counts_release_only,
            split=False,
            filter_intervals=filter_intervals,
            split_reference_blocks=False,
        )
        variant_mt = _split_variant_data_keeping_phase(
            vds.variant_data, ht, ["GT", "GQ", "DP", "AD", "PGT", "PID"],
        )
        vds = hl.vds.VariantDataset(vds.reference_data, variant_mt)
        mt = hl.vds.to_dense_mt(vds)
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
                    ht.write(f"{count_output_dir}/counts_light.ht", overwrite=overwrite)
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
                        emit_phase=emit_phase,
                    )
                    ht.write(f"{count_output_dir}/counts_light.ht", overwrite=overwrite)
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
                    ht.write(f"{count_output_dir}/counts_heavy.ht", overwrite=overwrite)
                    logger.info("Heavy counts written.")

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
            ht = ht.naive_coalesce(1000).checkpoint(res.vp_gt_counts_ht.path, overwrite=overwrite)
            logger.info("The variant pair genotype counts Table has been written...")

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
            f"is {DEFAULT_DATA_TYPE}.",
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

    args = parser.parse_args()
    main(args)

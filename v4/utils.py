"""Utility functions for variant co-occurrence pipeline."""

import logging
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad_chets.v4.resources import (
    CLINVAR_CATEGORIES,
    CLINVAR_CATEGORY_BLB,
    CLINVAR_CATEGORY_PLP,
    CLINVAR_CATEGORY_VUS,
    TEST_INTERVALS,
)

logger = logging.getLogger(__name__)


########################################################################################
### Shared per-variant QC / AF annotation
########################################################################################


def annotate_filters_af(
    ht: hl.Table,
    filter_ht: hl.Table,
    freq_ht: hl.Table,
    *,
    include_an: bool = False,
) -> hl.Table:
    """Annotate ``ht`` with ``filters`` and ``af`` (optionally ``an``).

    Every ``create_*_variant_filter_ht`` builder needs the same per-variant
    QC-pass + AF lookup; this consolidates it so the orchestrator can
    annotate (and pre-filter to QC-PASS + AF > 0) once and pass the result
    to each builder. ``an`` is needed only by the partner-set creation
    path (for the freq-based ``n_pair`` estimate).

    :param ht: Table keyed by ``(locus, alleles)``.
    :param filter_ht: gnomAD final-filter Table.
    :param freq_ht: gnomAD frequency Table.
    :param include_an: Whether to also annotate ``an`` (from
        ``freq[0].AN``).
    :return: ``ht`` with ``filters``, ``af`` (and optionally ``an``) added.
    """
    annotations = {
        "filters": filter_ht[ht.locus, ht.alleles].filters,
        "af": freq_ht[ht.locus, ht.alleles].freq[0].AF,
    }
    if include_an:
        annotations["an"] = freq_ht[ht.locus, ht.alleles].freq[0].AN
    return ht.annotate(**annotations)


########################################################################################
### ClinVar category filter
########################################################################################
# TODO: upstream to gnomad_methods/gnomad/utils/filtering.py as a
# generalization of filter_to_clinvar_pathogenic.

_CLINVAR_NO_STAR_ASSERTIONS = hl.literal(
    {
        "no_assertion_provided",
        "no_assertion_criteria_provided",
        "no_interpretation_for_the_individual_variant",
    }
)


def clinvar_category_match_expr(
    clnsig: hl.expr.ArrayExpression,
    category: str,
    clnrevstat: Optional[hl.expr.ArrayExpression] = None,
    clnsigconf: Optional[hl.expr.StringExpression] = None,
    remove_no_assertion: bool = True,
    remove_conflicting: bool = True,
) -> hl.expr.BooleanExpression:
    """Row-level expression for ClinVar significance-category membership.

    Categories (matched against ``CLNSIG`` substrings, case-insensitive):

    * ``"plp"``: any entry contains ``"pathogenic"`` (matches
      ``Pathogenic`` / ``Likely_pathogenic``).
    * ``"blb"``: at least one entry contains ``"benign"`` AND no entry
      contains ``"pathogenic"`` (excludes mixed records).
    * ``"vus"``: any entry contains ``"uncertain_significance"``.

    When ``remove_no_assertion`` is ``True``, ``clnrevstat`` is required
    and zero-star records are excluded; when ``remove_conflicting`` is
    ``True``, ``clnsigconf`` is required and rows with a conflicting-
    interpretation entry are excluded.

    :param clnsig: ``info.CLNSIG`` array (case-insensitive substring match).
    :param category: One of :data:`CLINVAR_CATEGORIES`.
    :param clnrevstat: ``info.CLNREVSTAT`` array; required when
        ``remove_no_assertion`` is ``True``.
    :param clnsigconf: ``info.CLNSIGCONF`` string; required when
        ``remove_conflicting`` is ``True``.
    :return: BooleanExpression that is ``True`` for rows in ``category``.
    :raises ValueError: If ``category`` is not in :data:`CLINVAR_CATEGORIES`,
        or a required field is missing.
    """
    if remove_no_assertion and clnrevstat is None:
        raise ValueError("clnrevstat is required when remove_no_assertion is True.")
    if remove_conflicting and clnsigconf is None:
        raise ValueError("clnsigconf is required when remove_conflicting is True.")
    sigs_lower = clnsig.map(lambda x: x.lower())
    if category == CLINVAR_CATEGORY_PLP:
        match_expr = sigs_lower.any(lambda x: x.contains("pathogenic"))
    elif category == CLINVAR_CATEGORY_BLB:
        has_benign = sigs_lower.any(lambda x: x.contains("benign"))
        has_pathogenic = sigs_lower.any(lambda x: x.contains("pathogenic"))
        match_expr = has_benign & ~has_pathogenic
    elif category == CLINVAR_CATEGORY_VUS:
        match_expr = sigs_lower.any(
            lambda x: x.contains("uncertain_significance")
        )
    else:
        raise ValueError(
            f"Unknown ClinVar category: {category!r}. Valid: "
            f"{CLINVAR_CATEGORIES}."
        )
    if remove_no_assertion:
        match_expr = match_expr & (
            hl.set(clnrevstat)
            .intersection(_CLINVAR_NO_STAR_ASSERTIONS)
            .length()
            == 0
        )
    if remove_conflicting:
        match_expr = match_expr & hl.is_missing(clnsigconf)
    return match_expr


def filter_clinvar_by_category(
    ht: hl.Table,
    category: str,
    clnrevstat_field: str = "CLNREVSTAT",
    clnsig_field: str = "CLNSIG",
    clnsigconf_field: str = "CLNSIGCONF",
    remove_no_assertion: bool = True,
    remove_conflicting: bool = True,
) -> hl.Table:
    """Filter a ClinVar Table to one significance category.

    Thin wrapper over :func:`clinvar_category_match_expr` operating on a
    ClinVar Table with the standard ``info`` schema.

    :param ht: ClinVar Table (e.g. from
        ``gnomad.resources.grch38.reference_data.clinvar.ht()``).
    :param category: One of :data:`CLINVAR_CATEGORIES`.
    :return: Filtered Table.
    """
    match_expr = clinvar_category_match_expr(
        clnsig=ht.info[clnsig_field],
        category=category,
        clnrevstat=ht.info[clnrevstat_field] if remove_no_assertion else None,
        clnsigconf=ht.info[clnsigconf_field] if remove_conflicting else None,
        remove_no_assertion=remove_no_assertion,
        remove_conflicting=remove_conflicting,
    )
    return ht.filter(match_expr)


########################################################################################
### AN_percent annotation (per-locus, X/Y-aware)
########################################################################################

ADJ_FREQ_META = {"group": "adj"}
"""Strata-meta filter for the ``adj`` (all-samples) entry in the AN HT."""

AN_CUTOFFS: List[int] = [50, 75, 80, 85, 90, 95]
"""Standard AN_percent cutoffs (percent of max AN) used for variant/pair QC."""


def _find_strata_index(strata_meta: List[dict], criteria: dict) -> int:
    """Find the index of the first strata-meta entry matching ``criteria``."""
    for i, m in enumerate(strata_meta):
        if all(m.get(k) == v for k, v in criteria.items()):
            return i
    raise ValueError(
        f"No strata-meta entry matches {criteria}; "
        f"available entries: {strata_meta}"
    )


def get_an_percent_expr(
    an_ht: hl.Table,
    locus_expr: hl.expr.LocusExpression,
) -> hl.expr.Int32Expression:
    """
    Build an expression that returns AN_percent at ``locus_expr``.

    AN_percent is ``int((AN / total_AN) * 100)`` where ``total_AN`` is the
    maximum possible AN for the locus given the sample-set composition:

    * Autosomes: ``2 * N`` (where ``N = strata_sample_count[adj]``).
    * X non-PAR: ``2 * XX + XY``.
    * Y non-PAR: ``XY`` (only XY samples have Y).

    The strata indexing mirrors
    :func:`gnomad_constraint.utils.constraint.get_exome_coverage_expr`.

    :param an_ht: gnomAD all-sites AN Table (keyed by ``locus``) with row
        field ``AN`` and global fields ``strata_sample_count`` /
        ``strata_meta``. Typically returned by
        :func:`gnomad.resources.grch38.gnomad.all_sites_an`.
    :param locus_expr: Locus expression to look up against ``an_ht``.
    :return: Int32 expression with AN_percent (typically 0-100).
    """
    an_globals = hl.eval(an_ht.index_globals())
    strata_sample_count = list(an_globals.strata_sample_count)
    strata_meta = [dict(m) for m in an_globals.strata_meta]

    adj_idx = _find_strata_index(strata_meta, ADJ_FREQ_META)
    xx_count = strata_sample_count[
        _find_strata_index(strata_meta, {**ADJ_FREQ_META, "sex": "XX"})
    ]
    xy_count = strata_sample_count[
        _find_strata_index(strata_meta, {**ADJ_FREQ_META, "sex": "XY"})
    ]
    adj_count = strata_sample_count[adj_idx]

    an_count = (
        hl.case()
        .when(locus_expr.in_x_nonpar(), (xx_count * 2) + xy_count)
        .when(locus_expr.in_y_nonpar(), xy_count)
        .default(adj_count * 2)
    )

    # AN is an array indexed by strata; pull the adj-group value.
    return hl.int((an_ht[locus_expr].AN[adj_idx] / an_count) * 100)


def calculate_partitions_by_size(
    ht: hl.Table,
    n_partitions: int,
    size_field: str,
    weight_ht: Optional[hl.Table] = None,
    weight_field: Optional[str] = None,
) -> List:
    """
    Calculate partition intervals that balance by estimated row size.

    Unlike ``_calculate_new_partitions`` which splits by row count, this uses
    the length of an array field as a size proxy to create intervals where each
    partition holds roughly equal total data.

    Works with any key type (e.g., ``(locus, alleles)``, ``int64``, etc.).

    When ``weight_field`` is provided, the size for each key is multiplied by
    the length of the weight field. This is useful when a downstream join
    replicates data — e.g., a variant appearing in many variant pairs will have
    its genotype array replicated for each pair, so the effective size is
    ``len(gt_info) * len(variant_pairs)``. The weight field can live on ``ht``
    itself or on a separate ``weight_ht`` (joined by key).

    :param ht: Sorted Table with a size-proxy array field.
    :param n_partitions: Desired number of partitions.
    :param size_field: Name of the array field on ``ht`` to use as size proxy.
    :param weight_ht: Optional second Table with the same key as ``ht`` holding
        the weight field. If ``None`` and ``weight_field`` is provided, the
        weight field is read from ``ht``.
    :param weight_field: Name of the array field to use as a multiplier on
        row size. Looked up on ``weight_ht`` if provided, otherwise on ``ht``.
    :return: List of partition intervals suitable for
        ``hl.read_table(_intervals=...)``.
    """
    def _as_int64(expr):
        """Coerce a scalar or array expression to an int64 size proxy."""
        if isinstance(expr, hl.expr.CollectionExpression):
            return hl.or_else(hl.int64(hl.len(expr)), hl.int64(0))
        return hl.or_else(hl.int64(expr), hl.int64(0))

    row_size = _as_int64(ht[size_field])
    if weight_field is not None:
        weight_expr = (
            weight_ht[ht.key][weight_field]
            if weight_ht is not None
            else ht[weight_field]
        )
        row_size = row_size * _as_int64(weight_expr)

    # Compute row size and exclusive cumulative sum in a single pass, then
    # checkpoint. Total is derived from the last row (scan.sum is exclusive,
    # so total = last _cumsum + last _row_size). Boundary detection is done
    # with a filter on the checkpointed table, which only reads — no
    # recomputation of the scan.
    ht = ht.select(
        _row_size=row_size,
        _cumsum=hl.scan.sum(row_size),
    )
    ht = ht.checkpoint(
        hl.utils.new_temp_file("calc_partitions_by_size", "ht")
    )

    key_fields = list(ht.key)
    key_dtype = ht.key.dtype

    # Total = last row's cumsum + last row's size.
    last_row = ht.tail(1).collect()[0]
    total_size = last_row._cumsum + last_row._row_size
    target_per_partition = total_size / n_partitions

    logger.info(
        "Calculating size-balanced partitions: total_size=%d, n_partitions=%d, "
        "target_per_partition=%d",
        total_size,
        n_partitions,
        target_per_partition,
    )

    # Partition index for each row = floor(_cumsum / target). A boundary is
    # where the partition index exceeds the max seen so far (via scan.max on
    # the checkpointed _cumsum field — no nested scans since _cumsum is a
    # stored field). Annotate + filter in a single pipeline; the scan runs
    # over the checkpoint's stored data.
    pidx = hl.int32(ht._cumsum / target_per_partition)
    ht = ht.annotate(_pidx=pidx, _prev_max_pidx=hl.scan.max(pidx))
    boundary_ht = ht.filter(
        hl.is_missing(ht._prev_max_pidx) | (ht._pidx > ht._prev_max_pidx)
    )
    boundary_keys = boundary_ht.key_by().select(*key_fields).collect()
    last_key = {f: getattr(last_row, f) for f in key_fields}

    def _to_struct(row):
        return hl.Struct(**{f: row[f] for f in key_fields})

    # Build intervals: each boundary starts a partition; successive boundaries
    # are the end of the previous partition. Final interval extends (inclusive)
    # to the last key of the table.
    intervals = []
    structs = [_to_struct(k) for k in boundary_keys]
    for i, start in enumerate(structs):
        if i + 1 < len(structs):
            end, includes_end = structs[i + 1], False
        else:
            end, includes_end = _to_struct(last_key), True
        intervals.append(
            hl.Interval(
                start=start,
                end=end,
                includes_end=includes_end,
                point_type=key_dtype,
            )
        )

    logger.info("Created %d size-balanced partition intervals.", len(intervals))
    return intervals


def compute_v2_independent_set(
    variants: List[tuple],
    edges: List[tuple],
) -> set:
    """
    Compute a lightweight independent set (V2) via greedy algorithm.

    Variants assigned to V2 will always be on the v2 (shuffle) side of pairs.
    The algorithm greedily selects the cheapest variants for V2 (low
    ``n_with_data`` relative to their degree), ensuring no two V2 variants
    are paired together.

    :param variants: List of ``(v_idx, n_with_data)`` tuples.
    :param edges: List of ``(v1_idx, v2_idx)`` tuples (the pair graph edges).
    :return: Set of ``v_idx`` values assigned to V2.
    """
    from collections import defaultdict

    # Build adjacency list.
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)

    # Score = n_with_data / max(degree, 1).
    # Lower score → better V2 candidate (cheap to shuffle, covers many edges).
    scored = []
    for v_idx, n in variants:
        deg = len(adj.get(v_idx, set()))
        score = n / max(deg, 1)
        scored.append((score, n, v_idx))
    scored.sort()

    v2_set = set()
    forbidden = set()
    for _, _, v_idx in scored:
        if v_idx not in forbidden:
            v2_set.add(v_idx)
            forbidden.update(adj[v_idx])

    logger.info(
        "V2 independent set: %d of %d variants (%.1f%%), covering %d of %d edges",
        len(v2_set),
        len(variants),
        100 * len(v2_set) / max(len(variants), 1),
        sum(1 for a, b in edges if a in v2_set or b in v2_set),
        len(edges),
    )
    return v2_set


def filter_for_testing(
    data: Union[hl.Table, hl.MatrixTable],
    test_intervals: Dict[str, str] = TEST_INTERVALS,
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter a Table or MatrixTable to a test interval.

    :param data: Hail Table or MatrixTable to filter.
    :param test_intervals: Dictionary of gene names and intervals to filter to.
        Default is TEST_INTERVALS.
    :return: Filtered Table or MatrixTable.
    """
    intervals = [hl.parse_locus_interval(interval, reference_genome="GRCh38") for interval in test_intervals.values()]
    if isinstance(data, hl.Table) or isinstance(data, hl.MatrixTable):
        return hl.filter_intervals(data, intervals)
    else:
        raise ValueError(
            f"Unsupported type: {type(data)}. Must be hl.Table or hl.MatrixTable."
        )

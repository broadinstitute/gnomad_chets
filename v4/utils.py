"""Utility functions for variant co-occurrence pipeline."""

import logging
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad_chets.v4.resources import TEST_INTERVALS

logger = logging.getLogger(__name__)


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

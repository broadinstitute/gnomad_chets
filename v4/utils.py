"""Utility functions for variant co-occurrence pipeline."""

import logging
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad_chets.v4.resources import TEST_INTERVALS

logger = logging.getLogger(__name__)


def calculate_partitions_by_size(
    ht: hl.Table,
    n_partitions: int,
    size_field: str = "gt_info",
    weight_ht: Optional[hl.Table] = None,
    weight_field: Optional[str] = None,
) -> List:
    """
    Calculate partition intervals that balance by estimated row size.

    Unlike ``_calculate_new_partitions`` which splits by row count, this uses
    the length of an array field as a size proxy to create intervals where each
    partition holds roughly equal total data.

    When ``weight_ht`` and ``weight_field`` are provided, the size for each key
    is multiplied by the length of the weight field from the second table. This
    is useful when a downstream join replicates data — e.g., a variant appearing
    in many variant pairs will have its genotype array replicated for each pair,
    so the effective size is ``len(gt_info) * len(variant_pairs)``.

    The approach:

        1. Annotate each row with its estimated size (with optional weight).
        2. Aggregate to get total estimated size across all rows.
        3. Annotate each row with a cumulative size via ``hl.scan.sum``.
        4. Assign each row a partition index = floor(cumsum / target).
        5. Filter to only the first row of each new partition bucket (where
           partition index > previous partition index, detected via a second
           scan tracking the running max of partition indices seen so far).
        6. Collect only those boundary keys (~n_partitions rows).
        7. Build ``hl.Interval`` list from boundary keys.

    :param ht: Table keyed by (locus, alleles) with a size-proxy array field.
    :param n_partitions: Desired number of partitions.
    :param size_field: Name of the array field on ``ht`` to use as size proxy.
    :param weight_ht: Optional second Table keyed by (locus, alleles) with a
        multiplier array field for row size.
    :param weight_field: Name of the array field on ``weight_ht`` to use as a
        multiplier. Required if ``weight_ht`` is provided.
    :return: List of partition intervals suitable for
        ``hl.read_table(_intervals=...)``.
    """
    row_size = hl.int64(hl.len(ht[size_field]))
    if weight_ht is not None:
        weight = weight_ht[ht.key][weight_field]
        row_size = row_size * hl.or_else(hl.len(weight), 0)

    ht = ht.annotate(_row_size=row_size)
    total_size = ht.aggregate(hl.agg.sum(ht._row_size))
    target_per_partition = total_size / n_partitions

    logger.info(
        "Calculating size-balanced partitions: total_size=%d, n_partitions=%d, "
        "target_per_partition=%d",
        total_size,
        n_partitions,
        target_per_partition,
    )

    # Cumulative size via scan (exclusive: value is sum of all preceding rows).
    ht = ht.select(
        _partition_idx=hl.int32(
            hl.scan.sum(ht._row_size) / target_per_partition
        ),
    )
    # The max partition index seen *before* this row (via scan.max) tells us
    # the previous row's partition index. When _partition_idx > _prev_max,
    # this row is the first in a new bucket — a boundary.
    ht = ht.annotate(_prev_max=hl.scan.max(ht._partition_idx))
    boundary_ht = ht.filter(ht._partition_idx > ht._prev_max)
    boundary_keys = boundary_ht.key_by().select("locus", "alleles").collect()

    # Build intervals from boundary keys.
    ref = ht.locus.dtype.reference_genome
    contig_info = ref.lengths
    first_contig = ref.contigs[0]
    last_contig = ref.contigs[-1]

    intervals = []
    start = hl.Struct(
        locus=hl.Locus(first_contig, 1, reference_genome=ref),
        alleles=["", ""],
    )
    for key in boundary_keys:
        end = hl.Struct(locus=key.locus, alleles=key.alleles)
        intervals.append(hl.Interval(start=start, end=end, includes_end=False))
        start = end
    # Final interval to the end.
    intervals.append(
        hl.Interval(
            start=start,
            end=hl.Struct(
                locus=hl.Locus(
                    last_contig, contig_info[last_contig], reference_genome=ref
                ),
                alleles=["zzz", "zzz"],
            ),
            includes_end=True,
        )
    )

    logger.info("Created %d size-balanced partition intervals.", len(intervals))
    return intervals


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

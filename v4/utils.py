"""Utility functions for variant co-occurrence pipeline."""

from typing import Union

import hail as hl

from gnomad_chets.v4.resources import TEST_INTERVAL


def filter_for_testing(
    data: Union[hl.Table, hl.MatrixTable], test_interval: str = TEST_INTERVAL
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter a Table or MatrixTable to a test interval.

    :param data: Hail Table or MatrixTable to filter.
    :param test_interval: Interval string to filter to (e.g., "chr21:46324141-46445769").
        Default is TEST_INTERVAL.
    :return: Filtered Table or MatrixTable.
    """
    interval = hl.parse_locus_interval(test_interval, reference_genome="GRCh38")
    if isinstance(data, hl.Table) or isinstance(data, hl.MatrixTable):
        return hl.filter_intervals(data, [interval])
    else:
        raise ValueError(
            f"Unsupported type: {type(data)}. Must be hl.Table or hl.MatrixTable."
        )

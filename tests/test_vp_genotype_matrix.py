"""Semantics tests for vp_genotype_matrix's genotype counting.

Self-contained: builds tiny MatrixTables in-process, no GCS or gnomAD access, so it
runs anywhere Hail runs. From the repo root:

    python -m pytest tests/test_vp_genotype_matrix.py
"""
import os
import sys

import hail as hl
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp_genotype_matrix import (  # noqa: E402
    GENOTYPE_CLASSES,
    create_variant_pair_genotype_counts_ht,
)

# The encoder's 7 per-variant categories, as (GT, adj):
#   0: no call        -> uncallable; must NOT land in the hom-ref/hom-ref cell
#   1: 0/0, adj FAIL  -> hom-ref for raw, excluded from adj
#   2: 0/1, adj FAIL  -> het for raw, excluded from adj
#   3: 0/1, adj PASS  -> het for both
#   4: 1/1, adj FAIL  -> hom-var for raw, excluded from adj
#   5: 1/1, adj PASS  -> hom-var for both
#   6: 0/0, adj PASS  -> hom-ref for both
CATEGORIES = [
    (None, False),
    ("0/0", False),
    ("0/1", False),
    ("0/1", True),
    ("1/1", False),
    ("1/1", True),
    ("0/0", True),
]
N_CATEGORIES = len(CATEGORIES)


def _call(gt):
    if gt is None:
        return hl.missing(hl.tcall)
    a, b = gt.split("/")
    return hl.call(int(a), int(b))


def _locus(pos):
    return hl.locus("1", pos, "GRCh37")


@pytest.fixture(scope="module")
def category_grid():
    """One sample per (category at v1, category at v2) combination.

    With 49 samples covering all 7x7 combinations exactly once, every genotype class
    is reachable and the expected counts are analytic: each of the 9 raw cells draws
    from a 2x2 block of the grid (4 samples), and each adj cell from exactly 1.
    """
    positions = {"v1": 1000, "v2": 2000}
    rows = []
    for sample_idx in range(N_CATEGORIES**2):
        cat_v1, cat_v2 = divmod(sample_idx, N_CATEGORIES)
        for name, cat in (("v1", cat_v1), ("v2", cat_v2)):
            gt, adj = CATEGORIES[cat]
            rows.append(
                {
                    "locus": _locus(positions[name]),
                    "alleles": ["A", "C"],
                    "s": f"s{sample_idx:02d}",
                    "GT": _call(gt),
                    "adj": adj,
                }
            )

    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh37"),
            alleles=hl.tarray(hl.tstr),
            s=hl.tstr,
            GT=hl.tcall,
            adj=hl.tbool,
        ),
    )
    mt = ht.to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])

    vp_ht = hl.Table.parallelize(
        [
            {
                "locus1": _locus(positions["v1"]),
                "alleles1": ["A", "C"],
                "locus2": _locus(positions["v2"]),
                "alleles2": ["A", "C"],
            }
        ],
        hl.tstruct(
            locus1=hl.tlocus("GRCh37"),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh37"),
            alleles2=hl.tarray(hl.tstr),
        ),
    ).key_by("locus1", "alleles1", "locus2", "alleles2")

    counts = create_variant_pair_genotype_counts_ht(mt, vp_ht, mt.adj, run_tag=None)
    return counts.collect()[0]


def test_raw_counts_over_category_grid(category_grid):
    """Raw: hom-ref = {no-adj 0/0, adj 0/0}, het = both 0/1, hom-var = both 1/1.

    That's a 2-category block per genotype per variant, so 4 samples per cell.
    """
    assert list(category_grid.gt_counts_raw) == [4] * len(GENOTYPE_CLASSES)


def test_adj_counts_over_category_grid(category_grid):
    """Adj: only the adj-PASS category counts for each genotype, so 1 per cell."""
    assert list(category_grid.gt_counts_adj) == [1] * len(GENOTYPE_CLASSES)


def test_uncallable_samples_excluded(category_grid):
    """No-call samples must not inflate any cell.

    49 samples, but the 13 with a no-call at either variant (row 0 and column 0 of
    the grid) carry no genotype at all, so they belong to no cell.
    """
    assert sum(category_grid.gt_counts_raw) == N_CATEGORIES**2 - (
        2 * N_CATEGORIES - 1
    )


def test_counts_are_zero_not_missing_when_nothing_survives():
    """A pair with no adj-passing sample yields zeros, not NA.

    gnomAD v2 falls back to a zero array here; returning NA instead would make such
    pairs read as absent rather than uncounted.
    """
    rows = []
    for sample_idx in range(4):
        for pos in (1000, 2000):
            rows.append(
                {
                    "locus": _locus(pos),
                    "alleles": ["A", "C"],
                    "s": f"s{sample_idx}",
                    "GT": _call("0/1"),
                    "adj": False,
                }
            )
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh37"),
            alleles=hl.tarray(hl.tstr),
            s=hl.tstr,
            GT=hl.tcall,
            adj=hl.tbool,
        ),
    )
    mt = ht.to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
    vp_ht = hl.Table.parallelize(
        [
            {
                "locus1": _locus(1000),
                "alleles1": ["A", "C"],
                "locus2": _locus(2000),
                "alleles2": ["A", "C"],
            }
        ],
        hl.tstruct(
            locus1=hl.tlocus("GRCh37"),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh37"),
            alleles2=hl.tarray(hl.tstr),
        ),
    ).key_by("locus1", "alleles1", "locus2", "alleles2")

    row = create_variant_pair_genotype_counts_ht(
        mt, vp_ht, mt.adj, run_tag=None
    ).collect()[0]

    assert list(row.gt_counts_adj) == [0] * len(GENOTYPE_CLASSES)
    # All four samples are het/het in raw, and nothing is hom-ref.
    assert row.gt_counts_raw[GENOTYPE_CLASSES.index("AaBb")] == 4
    assert row.gt_counts_raw[GENOTYPE_CLASSES.index("AABB")] == 0

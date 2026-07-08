"""
Unit tests for v4/in_trans_oe.py.

Run as:

    cd /Users/jgoodric/PycharmProjects && python -m gnomad_chets.v4.test_in_trans_oe

Tests are runnable scripts (matching the existing
``v4/test_create_vp_matrix.py`` style) — they spin up a local Hail backend,
build small in-memory fixture Tables, and assert on per-pair math
(``annotate_pair_oe_terms``).
"""

import hail as hl

from gnomad_chets.v4.in_trans_oe import (
    PARTNER_SET_CLINVAR_PLP,
    aggregate_oe_per_candidate,
    annotate_pair_oe_terms,
)


def _build_pair_ht(rows):
    return hl.Table.parallelize(
        rows,
        schema=hl.tstruct(
            locus1=hl.tlocus("GRCh38"),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh38"),
            alleles2=hl.tarray(hl.tstr),
            gt_counts_raw=hl.tarray(hl.tint32),
            gt_counts_adj=hl.tarray(hl.tint32),
        ),
        key=["locus1", "alleles1", "locus2", "alleles2"],
    )


def _build_freq_ht(rows):
    return hl.Table.parallelize(
        rows,
        schema=hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            freq=hl.tarray(hl.tstruct(AF=hl.tfloat64)),
        ),
        key=["locus", "alleles"],
    )


def test_annotate_pair_oe_terms_basic():
    """Per-pair math wiring: AF lookup, double_carriers, e_pair, o_pair."""
    pair_ht = _build_pair_ht([
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"),
            alleles1=["A", "T"],
            locus2=hl.locus("chr21", 2000, "GRCh38"),
            alleles2=["G", "C"],
            gt_counts_raw=[80, 5, 0, 5, 8, 1, 0, 1, 0],
            gt_counts_adj=[80, 5, 0, 5, 8, 1, 0, 1, 0],
        ),
    ])
    freq_ht = _build_freq_ht([
        hl.struct(
            locus=hl.locus("chr21", 1000, "GRCh38"),
            alleles=["A", "T"],
            freq=[hl.struct(AF=0.1)],
        ),
        hl.struct(
            locus=hl.locus("chr21", 2000, "GRCh38"),
            alleles=["G", "C"],
            freq=[hl.struct(AF=0.05)],
        ),
    ])

    n_samples = 100
    annotated = annotate_pair_oe_terms(pair_ht, freq_ht, n_samples=n_samples)
    rows = annotated.collect()

    assert len(rows) == 1
    row = rows[0]

    assert abs(row.af1 - 0.1) < 1e-9, f"af1: {row.af1}"
    assert abs(row.af2 - 0.05) < 1e-9, f"af2: {row.af2}"

    assert row.double_carriers == 8, f"double_carriers: {row.double_carriers}"

    expected_e = 2 * n_samples * 0.1 * 0.05
    assert abs(row.e_pair - expected_e) < 1e-9, f"e_pair: {row.e_pair}"

    assert 0.0 <= row.p_chet <= 1.0, f"p_chet out of [0,1]: {row.p_chet}"

    assert abs(row.o_pair - row.p_chet * row.double_carriers) < 1e-9, (
        f"o_pair {row.o_pair} != p_chet {row.p_chet} * dc {row.double_carriers}"
    )

    print(
        f"PASS test_annotate_pair_oe_terms_basic: "
        f"af1={row.af1:.4f}, af2={row.af2:.4f}, "
        f"double_carriers={row.double_carriers}, "
        f"e_pair={row.e_pair:.4f}, p_chet={row.p_chet:.4f}, "
        f"o_pair={row.o_pair:.4f}"
    )


def test_annotate_pair_oe_terms_use_raw_vs_adj():
    """``use_adj`` toggles which 9-array drives double_carriers and p_chet."""
    pair_ht = _build_pair_ht([
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"),
            alleles1=["A", "T"],
            locus2=hl.locus("chr21", 2000, "GRCh38"),
            alleles2=["G", "C"],
            gt_counts_raw=[70, 4, 0, 4, 12, 1, 0, 1, 0],
            gt_counts_adj=[80, 5, 0, 5, 7, 1, 0, 1, 0],
        ),
    ])
    freq_ht = _build_freq_ht([
        hl.struct(
            locus=hl.locus("chr21", 1000, "GRCh38"),
            alleles=["A", "T"],
            freq=[hl.struct(AF=0.1)],
        ),
        hl.struct(
            locus=hl.locus("chr21", 2000, "GRCh38"),
            alleles=["G", "C"],
            freq=[hl.struct(AF=0.05)],
        ),
    ])

    adj_row = annotate_pair_oe_terms(
        pair_ht, freq_ht, n_samples=100, use_adj=True
    ).collect()[0]
    raw_row = annotate_pair_oe_terms(
        pair_ht, freq_ht, n_samples=100, use_adj=False
    ).collect()[0]

    assert adj_row.double_carriers == 7, f"adj double_carriers: {adj_row.double_carriers}"
    assert raw_row.double_carriers == 12, f"raw double_carriers: {raw_row.double_carriers}"

    # e_pair depends only on n_samples and AFs, not on counts.
    assert abs(adj_row.e_pair - raw_row.e_pair) < 1e-9

    print(
        f"PASS test_annotate_pair_oe_terms_use_raw_vs_adj: "
        f"adj.double_carriers={adj_row.double_carriers}, "
        f"raw.double_carriers={raw_row.double_carriers}"
    )


def test_annotate_pair_oe_terms_n_scaling():
    """e_pair scales linearly with n_samples; o_pair does not."""
    pair_ht = _build_pair_ht([
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"),
            alleles1=["A", "T"],
            locus2=hl.locus("chr21", 2000, "GRCh38"),
            alleles2=["G", "C"],
            gt_counts_raw=[80, 5, 0, 5, 8, 1, 0, 1, 0],
            gt_counts_adj=[80, 5, 0, 5, 8, 1, 0, 1, 0],
        ),
    ])
    freq_ht = _build_freq_ht([
        hl.struct(
            locus=hl.locus("chr21", 1000, "GRCh38"),
            alleles=["A", "T"],
            freq=[hl.struct(AF=0.2)],
        ),
        hl.struct(
            locus=hl.locus("chr21", 2000, "GRCh38"),
            alleles=["G", "C"],
            freq=[hl.struct(AF=0.1)],
        ),
    ])

    r100 = annotate_pair_oe_terms(pair_ht, freq_ht, n_samples=100).collect()[0]
    r1000 = annotate_pair_oe_terms(pair_ht, freq_ht, n_samples=1000).collect()[0]

    assert abs(r100.e_pair - 2 * 100 * 0.2 * 0.1) < 1e-9
    assert abs(r1000.e_pair - 2 * 1000 * 0.2 * 0.1) < 1e-9
    assert abs(r1000.e_pair - 10 * r100.e_pair) < 1e-9
    assert abs(r100.o_pair - r1000.o_pair) < 1e-9, "o_pair must not depend on n_samples"

    print(
        f"PASS test_annotate_pair_oe_terms_n_scaling: "
        f"e@100={r100.e_pair:.4f}, e@1000={r1000.e_pair:.4f} "
        f"(ratio {r1000.e_pair / r100.e_pair:.2f}); o invariant={r100.o_pair:.4f}"
    )


def _build_annotated_pair_ht(rows):
    return hl.Table.parallelize(
        rows,
        schema=hl.tstruct(
            locus1=hl.tlocus("GRCh38"),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh38"),
            alleles2=hl.tarray(hl.tstr),
            af1=hl.tfloat64,
            af2=hl.tfloat64,
            double_carriers=hl.tint32,
            p_chet=hl.tfloat64,
            e_pair=hl.tfloat64,
            o_pair=hl.tfloat64,
        ),
        key=["locus1", "alleles1", "locus2", "alleles2"],
    )


def _build_variant_ht(rows):
    return hl.Table.parallelize(
        rows,
        schema=hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            gene_id=hl.tarray(hl.tstr),
        ),
        key=["locus", "alleles"],
    )


def test_aggregate_oe_per_candidate_basic():
    """Two candidates × one partner in the same gene; one extra pair where the
    "partner" is not in partner_ht should be excluded."""
    # Variants A=chr21:1000:A>T, B=chr21:2000:G>C, P=chr21:3000:T>A — all in gene G1.
    # Pairs (v1<=v2): (A,B), (A,P), (B,P).
    # P is in partner_ht; A and B are not.
    # Expected output:
    #   (A, G1): partners=[P]
    #   (B, G1): partners=[P]
    annotated_pair_ht = _build_annotated_pair_ht([
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"), alleles1=["A", "T"],
            locus2=hl.locus("chr21", 2000, "GRCh38"), alleles2=["G", "C"],
            af1=0.01, af2=0.005, double_carriers=0, p_chet=0.5, e_pair=0.1, o_pair=0.0,
        ),
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"), alleles1=["A", "T"],
            locus2=hl.locus("chr21", 3000, "GRCh38"), alleles2=["T", "A"],
            af1=0.01, af2=0.001, double_carriers=0, p_chet=0.5, e_pair=0.02, o_pair=0.0,
        ),
        hl.struct(
            locus1=hl.locus("chr21", 2000, "GRCh38"), alleles1=["G", "C"],
            locus2=hl.locus("chr21", 3000, "GRCh38"), alleles2=["T", "A"],
            af1=0.005, af2=0.001, double_carriers=0, p_chet=0.5, e_pair=0.01, o_pair=0.0,
        ),
    ])
    partner_ht = _build_variant_ht([
        hl.struct(locus=hl.locus("chr21", 3000, "GRCh38"), alleles=["T", "A"], gene_id=["G1"]),
    ])
    candidate_ht = _build_variant_ht([
        hl.struct(locus=hl.locus("chr21", 1000, "GRCh38"), alleles=["A", "T"], gene_id=["G1"]),
        hl.struct(locus=hl.locus("chr21", 2000, "GRCh38"), alleles=["G", "C"], gene_id=["G1"]),
        hl.struct(locus=hl.locus("chr21", 3000, "GRCh38"), alleles=["T", "A"], gene_id=["G1"]),
    ])

    out = aggregate_oe_per_candidate(
        annotated_pair_ht, partner_ht, candidate_ht, partner_set=PARTNER_SET_CLINVAR_PLP
    )
    rows = sorted(out.collect(), key=lambda r: r.locus.position)

    assert len(rows) == 2, f"expected 2 rows, got {len(rows)}"

    a, b = rows
    assert a.locus.position == 1000 and a.alleles == ["A", "T"]
    assert a.gene_id == "G1"
    assert a.partner_set == "CLINVAR_PLP"
    assert a.n_partners == 1
    assert abs(a.candidate_af - 0.01) < 1e-9
    assert abs(a.total_expected_in_trans - 0.02) < 1e-9
    assert a.total_observed_in_trans == 0.0
    assert len(a.partners) == 1
    assert a.partners[0].locus.position == 3000

    assert b.locus.position == 2000 and b.alleles == ["G", "C"]
    assert b.n_partners == 1
    assert b.partners[0].locus.position == 3000

    print(f"PASS test_aggregate_oe_per_candidate_basic: 2 rows")


def test_aggregate_oe_per_candidate_gene_intersection_and_warnings():
    """Multi-gene candidate × multi-gene partner produces one row per shared gene;
    pairs with no shared gene are excluded; low-AF candidate triggers a warning."""
    # A in {G1, G2}; P in {G1, G2}; Q in {G3} (no shared gene with anyone).
    # Pairs: (A, P), (A, Q).
    # Expected: (A, G1) with partner P; (A, G2) with partner P; nothing for Q.
    # A's AF=0.001 → low_candidate_af warning. e_pair=0.05 → low_total_expected too.
    annotated_pair_ht = _build_annotated_pair_ht([
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"), alleles1=["A", "T"],
            locus2=hl.locus("chr21", 2000, "GRCh38"), alleles2=["G", "C"],
            af1=0.001, af2=0.01, double_carriers=0, p_chet=0.5, e_pair=0.05, o_pair=0.0,
        ),
        hl.struct(
            locus1=hl.locus("chr21", 1000, "GRCh38"), alleles1=["A", "T"],
            locus2=hl.locus("chr21", 4000, "GRCh38"), alleles2=["C", "G"],
            af1=0.001, af2=0.02, double_carriers=0, p_chet=0.5, e_pair=0.1, o_pair=0.0,
        ),
    ])
    partner_ht = _build_variant_ht([
        hl.struct(locus=hl.locus("chr21", 2000, "GRCh38"), alleles=["G", "C"], gene_id=["G1", "G2"]),
        hl.struct(locus=hl.locus("chr21", 4000, "GRCh38"), alleles=["C", "G"], gene_id=["G3"]),
    ])
    candidate_ht = _build_variant_ht([
        hl.struct(locus=hl.locus("chr21", 1000, "GRCh38"), alleles=["A", "T"], gene_id=["G1", "G2"]),
        hl.struct(locus=hl.locus("chr21", 2000, "GRCh38"), alleles=["G", "C"], gene_id=["G1", "G2"]),
        hl.struct(locus=hl.locus("chr21", 4000, "GRCh38"), alleles=["C", "G"], gene_id=["G3"]),
    ])

    out = aggregate_oe_per_candidate(
        annotated_pair_ht, partner_ht, candidate_ht, partner_set=PARTNER_SET_CLINVAR_PLP
    )
    rows = sorted(out.collect(), key=lambda r: (r.locus.position, r.gene_id))

    # Two rows from candidate A: one in G1, one in G2. Q-partner is filtered (no shared gene).
    assert len(rows) == 2, f"expected 2 rows, got {len(rows)}: {[(r.locus.position, r.gene_id) for r in rows]}"

    g1_row, g2_row = rows
    assert g1_row.gene_id == "G1"
    assert g2_row.gene_id == "G2"

    for row in (g1_row, g2_row):
        assert row.locus.position == 1000
        assert row.n_partners == 1
        assert row.partners[0].locus.position == 2000
        # Both warnings fire: AF=0.001 < 0.005, and total_expected=0.05 < 1
        warnings = list(row.warnings)
        assert "low_candidate_af" in warnings
        assert "low_total_expected" in warnings

    print(
        f"PASS test_aggregate_oe_per_candidate_gene_intersection_and_warnings: "
        f"2 rows, both with low_candidate_af + low_total_expected warnings"
    )


if __name__ == "__main__":
    hl.init(quiet=True, log="/tmp/test_in_trans_oe.log")
    test_annotate_pair_oe_terms_basic()
    test_annotate_pair_oe_terms_use_raw_vs_adj()
    test_annotate_pair_oe_terms_n_scaling()
    test_aggregate_oe_per_candidate_basic()
    test_aggregate_oe_per_candidate_gene_intersection_and_warnings()
    print("All tests passed.")

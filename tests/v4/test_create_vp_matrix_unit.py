"""Pytest unit tests for the smaller helpers in
``gnomad_chets.v4.compute_vp_counts``.

The existing ``v4/test_create_vp_matrix.py`` is an integration script that
runs the full pipeline against a fixture dense MT; this file covers the
pure-transform helpers and recently-added bug fixes with small in-memory
inputs that pytest can run in a few seconds.

Run with::

    pytest v4/test_create_vp_matrix_unit.py -v
"""
import hail as hl
import pytest

from gnomad_chets.v4.compute_vp_counts import (
    DEFAULT_SHUFFLE_BUDGET_BYTES,
    MIN_HEAVY_PARTITIONS,
    TARGET_HEAVY_PARTITION_BYTES,
    _build_pop_stratification,
    _COUNT_FROM_SETS_FIELDS,
    _count_from_sets,
    _count_from_sets_by_pop,
    _create_var_idx_ht,
    _drop_pairs_missing_v_idx,
    _empty_counts_ht,
    _heavy_filter_by_contribution,
    _project_count_fields,
    _read_min_an_pct,
    filter_pairs_by_an_pct,
)
from gnomad_chets.v4.resources import GLOBAL_POP


# ---------------------------------------------------------------------------
# Session-scoped Hail init
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def _hail_session():
    hl.init(idempotent=True, quiet=True)
    yield


# ---------------------------------------------------------------------------
# Helpers for building tiny test fixtures
# ---------------------------------------------------------------------------

def _vp_table(rows):
    """Build a variant-pair Table with the standard 4-field key.

    Each row dict needs keys ``c1, p1, a1, c2, p2, a2`` for v1 and v2.
    """
    structs = [
        hl.Struct(
            locus1=hl.locus(r["c1"], r["p1"], "GRCh38"),
            alleles1=r["a1"],
            locus2=hl.locus(r["c2"], r["p2"], "GRCh38"),
            alleles2=r["a2"],
        )
        for r in rows
    ]
    return hl.Table.parallelize(
        structs,
        hl.tstruct(
            locus1=hl.tlocus("GRCh38"),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh38"),
            alleles2=hl.tarray(hl.tstr),
        ),
        key=["locus1", "alleles1", "locus2", "alleles2"],
    )


def _var_idx_table(variants):
    """Build a var_idx Table keyed by (locus, alleles)."""
    structs = [
        hl.Struct(
            locus=hl.locus(c, p, "GRCh38"),
            alleles=a,
            var_idx=hl.int64(i),
        )
        for i, (c, p, a) in enumerate(variants)
    ]
    return hl.Table.parallelize(
        structs,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            var_idx=hl.tint64,
        ),
        key=["locus", "alleles"],
    )


def _annotate_vidx(vp, var_idx_ht):
    return vp.annotate(
        v1_idx=var_idx_ht[vp.locus1, vp.alleles1].var_idx,
        v2_idx=var_idx_ht[vp.locus2, vp.alleles2].var_idx,
    )


# ===========================================================================
# _empty_counts_ht
# ===========================================================================

class TestEmptyCountsHt:

    def test_returns_empty(self):
        ht = _empty_counts_ht("GRCh38")
        assert ht.count() == 0

    def test_keyed_by_4_pair_fields(self):
        ht = _empty_counts_ht("GRCh38")
        assert list(ht.key) == ["locus1", "alleles1", "locus2", "alleles2"]

    def test_value_fields_are_gt_counts(self):
        ht = _empty_counts_ht("GRCh38")
        non_key = [f for f in ht.row if f not in ht.key]
        assert "gt_counts_raw" in non_key
        assert "gt_counts_adj" in non_key

    def test_no_by_pop_field_by_default(self):
        ht = _empty_counts_ht("GRCh38")
        assert "gt_counts_by_pop" not in ht.row

    def test_by_pop_schema_when_stratified(self):
        ht = _empty_counts_ht("GRCh38", stratify_by_pop=True)
        assert "gt_counts_by_pop" in ht.row
        assert ht.gt_counts_by_pop.dtype == hl.tdict(
            hl.tstr,
            hl.tstruct(raw=hl.tarray(hl.tint32), adj=hl.tarray(hl.tint32)),
        )

    def test_union_with_typed_result_works(self):
        # Real call site: empty result unioned with computed counts.
        empty = _empty_counts_ht("GRCh38")
        other = hl.Table.parallelize(
            [
                hl.Struct(
                    locus1=hl.locus("chr1", 100, "GRCh38"),
                    alleles1=["A", "T"],
                    locus2=hl.locus("chr1", 200, "GRCh38"),
                    alleles2=["G", "C"],
                    gt_counts_raw=[hl.int64(1)] * 9,
                    gt_counts_adj=[hl.int64(0)] * 9,
                )
            ],
            empty.row.dtype,
            key=list(empty.key),
        )
        out = empty.union(other)
        assert out.count() == 1


# ===========================================================================
# _project_count_fields
# ===========================================================================

class TestProjectCountFields:

    def _make_encoded_ht(self):
        rows = [
            hl.Struct(
                v_idx=hl.int64(0),
                raw_het=hl.empty_set(hl.tint32),
                raw_hv=hl.empty_set(hl.tint32),
                adj_het=hl.empty_set(hl.tint32),
                adj_hv=hl.empty_set(hl.tint32),
                all_samples=hl.empty_set(hl.tint32),
                raw_hr_adj_missing=hl.empty_set(hl.tint32),
                n_with_data=hl.int32(0),
                n_raw_hr_adj_missing=hl.int32(0),
                all_samples_is_complement=False,
                raw_hr_adj_missing_is_complement=False,
                # Encoder-only extras that should be dropped:
                implicit_homref=hl.empty_set(hl.tint32),
                _extra_diagnostics=hl.int64(123),
            )
        ]
        return hl.Table.parallelize(
            rows,
            hl.tstruct(
                v_idx=hl.tint64,
                raw_het=hl.tset(hl.tint32),
                raw_hv=hl.tset(hl.tint32),
                adj_het=hl.tset(hl.tint32),
                adj_hv=hl.tset(hl.tint32),
                all_samples=hl.tset(hl.tint32),
                raw_hr_adj_missing=hl.tset(hl.tint32),
                n_with_data=hl.tint32,
                n_raw_hr_adj_missing=hl.tint32,
                all_samples_is_complement=hl.tbool,
                raw_hr_adj_missing_is_complement=hl.tbool,
                implicit_homref=hl.tset(hl.tint32),
                _extra_diagnostics=hl.tint64,
            ),
            key=["v_idx"],
        )

    def test_keeps_only_count_from_sets_fields(self):
        encoded = self._make_encoded_ht()
        projected = _project_count_fields(encoded)
        kept = set(projected.row) - set(projected.key)
        assert kept == set(_COUNT_FROM_SETS_FIELDS)

    def test_extras_dropped(self):
        encoded = self._make_encoded_ht()
        projected = _project_count_fields(encoded)
        assert "implicit_homref" not in projected.row
        assert "_extra_diagnostics" not in projected.row

    def test_rows_preserved(self):
        encoded = self._make_encoded_ht()
        assert _project_count_fields(encoded).count() == encoded.count()


# ===========================================================================
# _drop_pairs_missing_v_idx
# ===========================================================================

class TestDropPairsMissingVIdx:

    def _build_pair_with_nulls(self):
        """Pair table where some pairs reference missing v_idxs."""
        # var_idx has only v_idx 0 (chr1:100 A>T). Pair table references
        # both that variant and a pair with chr1:999 (missing).
        var_idx = _var_idx_table([
            ("chr1", 100, ["A", "T"]),
            ("chr1", 200, ["A", "G"]),
        ])
        vp = _vp_table([
            # Both sides resolve.
            {"c1": "chr1", "p1": 100, "a1": ["A", "T"],
             "c2": "chr1", "p2": 200, "a2": ["A", "G"]},
            # v1 missing.
            {"c1": "chr1", "p1": 999, "a1": ["A", "T"],
             "c2": "chr1", "p2": 200, "a2": ["A", "G"]},
            # v2 missing.
            {"c1": "chr1", "p1": 100, "a1": ["A", "T"],
             "c2": "chr1", "p2": 998, "a2": ["A", "G"]},
            # Both missing.
            {"c1": "chr1", "p1": 999, "a1": ["A", "T"],
             "c2": "chr1", "p2": 998, "a2": ["A", "G"]},
        ])
        return _annotate_vidx(vp, var_idx)

    def test_no_drop_when_all_resolved(self):
        var_idx = _var_idx_table([
            ("chr1", 100, ["A", "T"]),
            ("chr1", 200, ["A", "G"]),
        ])
        vp = _vp_table([
            {"c1": "chr1", "p1": 100, "a1": ["A", "T"],
             "c2": "chr1", "p2": 200, "a2": ["A", "G"]},
        ])
        vp = _annotate_vidx(vp, var_idx)
        out = _drop_pairs_missing_v_idx(vp, "test")
        assert out.count() == 1

    def test_drops_v1_null(self):
        vp = self._build_pair_with_nulls()
        out = _drop_pairs_missing_v_idx(vp, "test")
        # Of 4 input pairs, only 1 (both resolved) survives.
        assert out.count() == 1

    def test_survivor_has_both_v_idxs(self):
        vp = self._build_pair_with_nulls()
        out = _drop_pairs_missing_v_idx(vp, "test")
        row = out.collect()[0]
        assert row.v1_idx is not None
        assert row.v2_idx is not None

    def test_logs_warning(self, caplog):
        import logging
        vp = self._build_pair_with_nulls()
        with caplog.at_level(logging.WARNING, logger="compute_vp_counts"):
            _drop_pairs_missing_v_idx(vp, "test_caller")
        assert any(
            "test_caller" in rec.message and "3 pairs" in rec.message
            for rec in caplog.records
        )


# ===========================================================================
# _create_var_idx_ht
# ===========================================================================

class TestCreateVarIdxHt:

    def test_indexes_match_row_order(self):
        # Build a 3-row, 1-col MatrixTable and check var_idx is 0, 1, 2.
        mt = hl.balding_nichols_model(n_populations=1, n_samples=1, n_variants=3)
        var_idx_ht = _create_var_idx_ht(mt)
        # In row key order, var_idx must be a permutation of [0, 1, 2].
        idxs = sorted(var_idx_ht.var_idx.collect())
        assert idxs == [0, 1, 2]

    def test_keyed_by_locus_alleles(self):
        mt = hl.balding_nichols_model(n_populations=1, n_samples=1, n_variants=2)
        var_idx_ht = _create_var_idx_ht(mt)
        assert list(var_idx_ht.key) == ["locus", "alleles"]


# ===========================================================================
# _heavy_filter_by_contribution
# ===========================================================================

class TestHeavyFilterByContribution:

    def _build_inputs(self, payloads_per_var, n_pairs_per_var):
        """Build a (encoded, vp, var_idx) triple with controlled per-variant
        payload size and degree.

        ``payloads_per_var``: list of per-variant payload bytes (= 4 * sum
        of set lengths). We inflate the ``all_samples`` set to reach that
        target.

        ``n_pairs_per_var``: list of per-variant degrees. Each variant is
        paired with var 0 ``n_pairs_per_var[v]`` times (variant 0 acts as a
        hub) for simple controlled degree.
        """
        n_vars = len(payloads_per_var)
        # var_idx: variants at chr1:100, chr1:200, chr1:300, ...
        variants = [("chr1", 100 + 100 * i, ["A", "T"]) for i in range(n_vars)]
        var_idx_ht = _var_idx_table(variants)

        # Encoded: all_samples = range(0, payload/4) so len * 4 == payload.
        encoded_rows = []
        for v in range(n_vars):
            set_len = max(1, payloads_per_var[v] // 4)
            encoded_rows.append(
                hl.Struct(
                    v_idx=hl.int64(v),
                    raw_het=hl.empty_set(hl.tint32),
                    raw_hv=hl.empty_set(hl.tint32),
                    adj_het=hl.empty_set(hl.tint32),
                    adj_hv=hl.empty_set(hl.tint32),
                    raw_hr_adj_missing=hl.empty_set(hl.tint32),
                    all_samples=hl.set(hl.range(0, set_len)),
                    n_with_data=hl.int32(set_len),
                    n_raw_hr_adj_missing=hl.int32(0),
                    all_samples_is_complement=False,
                    raw_hr_adj_missing_is_complement=False,
                )
            )
        encoded_ht = hl.Table.parallelize(
            encoded_rows,
            hl.tstruct(
                v_idx=hl.tint64,
                raw_het=hl.tset(hl.tint32),
                raw_hv=hl.tset(hl.tint32),
                adj_het=hl.tset(hl.tint32),
                adj_hv=hl.tset(hl.tint32),
                raw_hr_adj_missing=hl.tset(hl.tint32),
                all_samples=hl.tset(hl.tint32),
                n_with_data=hl.tint32,
                n_raw_hr_adj_missing=hl.tint32,
                all_samples_is_complement=hl.tbool,
                raw_hr_adj_missing_is_complement=hl.tbool,
            ),
            key=["v_idx"],
        )

        # Pair table: for each v > 0, create n_pairs_per_var[v] copies of
        # (v0, vN) — multiple rows to simulate v0's degree.
        pair_rows = []
        for v in range(1, n_vars):
            for _ in range(n_pairs_per_var[v]):
                pair_rows.append({
                    "c1": "chr1", "p1": 100, "a1": ["A", "T"],
                    "c2": "chr1", "p2": 100 + 100 * v, "a2": ["A", "T"],
                })
        vp_ht = _vp_table(pair_rows)
        return encoded_ht, vp_ht, var_idx_ht

    def test_empty_when_total_below_budget(self):
        # Tiny inputs: 100 bytes payload * 2 pairs = 200 bytes total << 10 GB.
        encoded, vp, var_idx = self._build_inputs([100, 100, 100], [0, 1, 1])
        result = _heavy_filter_by_contribution(
            encoded, vp, var_idx, DEFAULT_SHUFFLE_BUDGET_BYTES,
        )
        assert result.count() == 0

    def test_returns_split_count_field(self):
        # Push the budget low so the filter pulls something.
        encoded, vp, var_idx = self._build_inputs([400, 400, 400], [0, 1, 1])
        result = _heavy_filter_by_contribution(
            encoded, vp, var_idx, shuffle_budget_bytes=1,
        )
        # Schema must contain split_count.
        assert "split_count" in result.row
        assert "contribution" in result.row
        assert "v_idx" in result.row

    def test_split_count_at_least_one(self):
        encoded, vp, var_idx = self._build_inputs([400, 400, 400], [0, 1, 1])
        result = _heavy_filter_by_contribution(
            encoded, vp, var_idx, shuffle_budget_bytes=1,
        )
        sc = result.split_count.collect()
        assert all(s >= 1 for s in sc)

    def test_split_count_scales_with_contribution(self):
        # Variant 1 is paired with var 0 many times → high degree → high
        # contribution. Its split_count should exceed 1 when the
        # contribution exceeds TARGET_HEAVY_PARTITION_BYTES.
        big_payload = TARGET_HEAVY_PARTITION_BYTES // 2  # 250 MB per variant
        # var 1 gets degree 4 → contribution ≈ 4 * 250 MB = 1 GB → 2 splits.
        encoded, vp, var_idx = self._build_inputs(
            [big_payload, big_payload, big_payload],
            [0, 4, 0],
        )
        result = _heavy_filter_by_contribution(
            encoded, vp, var_idx, shuffle_budget_bytes=1,
        )
        rows = {r.v_idx: r for r in result.collect()}
        assert 1 in rows, "Expected var 1 (high degree) to be pulled"
        # 4 × 250 MB = 1 GB, target 500 MB → ceil(1024/500) = 3 splits.
        assert rows[1].split_count >= 2


# ===========================================================================
# filter_pairs_by_an_pct
# ===========================================================================

class TestFilterPairsByAnPct:

    def _build_an_annotated_pairs(self, rows):
        """Build a pair table with v1_an_pct and v2_an_pct fields."""
        structs = [
            hl.Struct(
                locus1=hl.locus("chr1", r["p1"], "GRCh38"),
                alleles1=["A", "T"],
                locus2=hl.locus("chr1", r["p2"], "GRCh38"),
                alleles2=["A", "G"],
                v1_an_pct=hl.int32(r["v1_pct"]),
                v2_an_pct=hl.int32(r["v2_pct"]),
            )
            for r in rows
        ]
        ht = hl.Table.parallelize(
            structs,
            hl.tstruct(
                locus1=hl.tlocus("GRCh38"),
                alleles1=hl.tarray(hl.tstr),
                locus2=hl.tlocus("GRCh38"),
                alleles2=hl.tarray(hl.tstr),
                v1_an_pct=hl.tint32,
                v2_an_pct=hl.tint32,
            ),
            key=["locus1", "alleles1", "locus2", "alleles2"],
        )
        return ht

    def test_keeps_when_both_above_floor(self):
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 95, "v2_pct": 90},
            {"p1": 101, "p2": 201, "v1_pct": 80, "v2_pct": 85},
        ])
        assert filter_pairs_by_an_pct(ht, 80).count() == 2

    def test_drops_when_v1_below_floor(self):
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 80, "v2_pct": 90},
            {"p1": 101, "p2": 201, "v1_pct": 70, "v2_pct": 95},  # v1 fails
        ])
        out = filter_pairs_by_an_pct(ht, 80)
        positions = out.locus1.position.collect()
        assert 101 not in positions

    def test_drops_when_v2_below_floor(self):
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 95, "v2_pct": 75},  # v2 fails
            {"p1": 101, "p2": 201, "v1_pct": 95, "v2_pct": 80},
        ])
        out = filter_pairs_by_an_pct(ht, 80)
        positions = out.locus1.position.collect()
        assert 100 not in positions
        assert 101 in positions

    def test_no_op_when_min_is_negative(self):
        # Negative floor = "no filter".
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 0, "v2_pct": 0},
            {"p1": 101, "p2": 201, "v1_pct": 0, "v2_pct": 0},
        ])
        assert filter_pairs_by_an_pct(ht, -1).count() == 2


# ===========================================================================
# _read_min_an_pct
# ===========================================================================

class TestReadMinAnPct:

    def test_reads_global(self):
        ht = hl.utils.range_table(1).annotate_globals(min_an_pct=42)
        assert _read_min_an_pct(ht) == 42

    def test_returns_negative_one_when_missing(self):
        ht = hl.utils.range_table(1)
        assert _read_min_an_pct(ht) == -1


# ===========================================================================
# _count_from_sets — real-data regression (chr19, partition 63)
# ===========================================================================

# Pre-extracted fixture (see scratchpad/extract_tier3_fixtures.py): a small
# set of validated (encoded_v1, encoded_v2, expected_gt_counts) tuples
# pulled from the chr19 test pipeline. The all_samples field has been
# rewritten to the FIXED proper-complement form so the test exercises the
# post-fix decoder; pairs whose endpoints couldn't be fixed without
# enumerating the full sample space (complement-A + complement-F) were
# excluded at extraction time.
#
# The fixture is read once per session via _real_data_fixture; per-row tests
# parameterise over the resulting list of Struct rows so a failing pair is
# named clearly in pytest output.

_REAL_DATA_FIXTURE_PATH = (
    "gs://gnomad-tmp-30day/test_fixtures/count_from_sets_chr19_real.ht"
)


@pytest.fixture(scope="session")
def _real_data_fixture():
    """Load the chr19 real-data fixture HT once per test session.

    Returns a tuple ``(rows, n_samples)`` where ``rows`` is the collected
    list of fixture rows and ``n_samples`` is the cohort size stored as
    the HT's ``n_samples`` global.
    """
    ht = hl.read_table(_REAL_DATA_FIXTURE_PATH)
    n_samples = hl.eval(ht.index_globals().n_samples)
    rows = ht.collect()
    return rows, int(n_samples)


def _row_id(row):
    """Stable pytest id for a fixture row: 'chr19:pos:ref:alt|chr19:pos:ref:alt'."""
    def _v(locus, alleles):
        return f"{locus.contig}:{locus.position}:{'/'.join(alleles)}"
    return f"{_v(row.locus1, row.alleles1)}|{_v(row.locus2, row.alleles2)}"


def _real_data_rows():
    """Module-level loader so pytest can parameterise over fixture rows.

    pytest.mark.parametrize evaluates at collection time, before the
    session-scoped fixture runs, so the load lives here instead.
    """
    try:
        ht = hl.read_table(_REAL_DATA_FIXTURE_PATH)
        return ht.collect()
    except Exception:
        # If the fixture isn't reachable (e.g. running offline without GCS
        # creds) skip the whole class rather than hard-fail collection.
        return []


_REAL_DATA_ROWS = _real_data_rows()


class TestCountFromSetsRealData:
    """Regression tests for the ``_count_from_sets`` complement-form bug.

    Each fixture row is a real ``(v1_encoded, v2_encoded, expected_counts)``
    triple pulled from the chr19 test pipeline. The bug class — buggy
    ``all_samples`` in complement form producing negative / over-large
    cells — is caught both by the exact-match per-row tests and by the
    ``sum(adj) <= n_samples`` invariant test.
    """

    @staticmethod
    def _call(v1, v2, n_samples, include_raw_hr_adj_missing):
        return hl.eval(
            _count_from_sets(
                v1.raw_het if include_raw_hr_adj_missing else v1.adj_het,
                v1.raw_hv if include_raw_hr_adj_missing else v1.adj_hv,
                v1.all_samples,
                hl.int32(v1.n_with_data),
                hl.bool(v1.all_samples_is_complement),
                v1.raw_hr_adj_missing,
                hl.int32(v1.n_raw_hr_adj_missing),
                hl.bool(v1.raw_hr_adj_missing_is_complement),
                v2.raw_het if include_raw_hr_adj_missing else v2.adj_het,
                v2.raw_hv if include_raw_hr_adj_missing else v2.adj_hv,
                v2.all_samples,
                hl.int32(v2.n_with_data),
                hl.bool(v2.all_samples_is_complement),
                v2.raw_hr_adj_missing,
                hl.int32(v2.n_raw_hr_adj_missing),
                hl.bool(v2.raw_hr_adj_missing_is_complement),
                hl.int32(n_samples),
                include_raw_hr_adj_missing=include_raw_hr_adj_missing,
            )
        )

    @pytest.mark.skipif(
        not _REAL_DATA_ROWS,
        reason="chr19 real-data fixture not reachable (GCS or fixture missing)",
    )
    @pytest.mark.parametrize(
        "row",
        _REAL_DATA_ROWS,
        ids=[_row_id(r) for r in _REAL_DATA_ROWS] if _REAL_DATA_ROWS else None,
    )
    def test_raw_counts_match(self, row, _real_data_fixture):
        _rows, n_samples = _real_data_fixture
        result = self._call(
            row.v1_encoded, row.v2_encoded, n_samples,
            include_raw_hr_adj_missing=True,
        )
        expected = list(row.expected_gt_counts_raw)
        assert list(result) == expected, (
            f"raw counts mismatch for pair {_row_id(row)}: "
            f"got {list(result)}, expected {expected}"
        )

    @pytest.mark.skipif(
        not _REAL_DATA_ROWS,
        reason="chr19 real-data fixture not reachable (GCS or fixture missing)",
    )
    @pytest.mark.parametrize(
        "row",
        _REAL_DATA_ROWS,
        ids=[_row_id(r) for r in _REAL_DATA_ROWS] if _REAL_DATA_ROWS else None,
    )
    def test_adj_counts_match(self, row, _real_data_fixture):
        _rows, n_samples = _real_data_fixture
        result = self._call(
            row.v1_encoded, row.v2_encoded, n_samples,
            include_raw_hr_adj_missing=False,
        )
        expected = list(row.expected_gt_counts_adj)
        assert list(result) == expected, (
            f"adj counts mismatch for pair {_row_id(row)}: "
            f"got {list(result)}, expected {expected}"
        )

    @pytest.mark.skipif(
        not _REAL_DATA_ROWS,
        reason="chr19 real-data fixture not reachable (GCS or fixture missing)",
    )
    @pytest.mark.parametrize(
        "row",
        _REAL_DATA_ROWS,
        ids=[_row_id(r) for r in _REAL_DATA_ROWS] if _REAL_DATA_ROWS else None,
    )
    def test_adj_sum_below_n_samples(self, row, _real_data_fixture):
        """Bug-class invariant: sum(adj) <= n_samples.

        The pre-fix decoder over-counted ``|pos ∩ A_pos|`` for complement-
        stored variants, which inflated edge cells past the cohort size.
        This invariant catches the same bug class without depending on the
        exact ground-truth values.
        """
        _rows, n_samples = _real_data_fixture
        result = self._call(
            row.v1_encoded, row.v2_encoded, n_samples,
            include_raw_hr_adj_missing=False,
        )
        total = sum(result)
        assert total <= n_samples, (
            f"adj sum {total} exceeds n_samples {n_samples} for pair "
            f"{_row_id(row)} (bug-class signal)"
        )
        # Cells must be non-negative; pre-fix the complement-bug drove
        # specific cells negative.
        assert all(c >= 0 for c in result), (
            f"adj cells contain negative entries for pair {_row_id(row)}: "
            f"{list(result)}"
        )


# ===========================================================================
# _build_pop_stratification — driver-side per-pop index sets
# ===========================================================================

class TestBuildPopStratification:
    """The samples-global → (pops, index sets, sizes) helper for per-pop counts."""

    def test_groups_indices_by_pop(self):
        samples = [
            hl.Struct(s="a", pop="nfe"),
            hl.Struct(s="b", pop="afr"),
            hl.Struct(s="c", pop="nfe"),
            hl.Struct(s="d", pop="afr"),
        ]
        pops, index_sets, sizes = _build_pop_stratification(samples)
        # "all" first, then the specific groups sorted.
        assert pops == [GLOBAL_POP, "afr", "nfe"]
        # "all" is handled via the flat counts, so it is NOT in the dicts.
        assert set(index_sets) == {"afr", "nfe"}
        assert sizes == {"afr": 2, "nfe": 2}
        assert hl.eval(index_sets["nfe"]) == {0, 2}
        assert hl.eval(index_sets["afr"]) == {1, 3}

    def test_missing_pop_excluded_from_specific(self):
        samples = [
            hl.Struct(s="a", pop="nfe"),
            hl.Struct(s="b", pop=None),
            hl.Struct(s="c", pop="nfe"),
        ]
        pops, index_sets, sizes = _build_pop_stratification(samples)
        assert pops == [GLOBAL_POP, "nfe"]
        assert sizes == {"nfe": 2}
        # index 1 (no pop) is in no specific-pop set.
        assert hl.eval(index_sets["nfe"]) == {0, 2}

    def test_raises_without_pop_source(self):
        # No pop_map and no `pop` field on the samples → cannot stratify.
        with pytest.raises(ValueError, match="needs a pop_map"):
            _build_pop_stratification([hl.Struct(s="a")])

    def test_pop_map_path(self):
        # The count-time path: samples carry only `s`, pop comes from the
        # meta-derived {s: pop} map — no re-encode / no `pop` on the encoding.
        samples = [hl.Struct(s="a"), hl.Struct(s="b"), hl.Struct(s="c")]
        pop_map = {"a": "nfe", "b": "afr", "c": "nfe"}
        pops, index_sets, sizes = _build_pop_stratification(samples, pop_map)
        assert pops == [GLOBAL_POP, "afr", "nfe"]
        assert sizes == {"afr": 1, "nfe": 2}
        assert hl.eval(index_sets["nfe"]) == {0, 2}
        assert hl.eval(index_sets["afr"]) == {1}

    def test_pop_map_missing_sample_excluded(self):
        samples = [hl.Struct(s="a"), hl.Struct(s="b")]
        pop_map = {"a": "nfe"}  # 'b' absent from the map
        pops, index_sets, sizes = _build_pop_stratification(samples, pop_map)
        assert pops == [GLOBAL_POP, "nfe"]
        assert sizes == {"nfe": 1}
        assert hl.eval(index_sets["nfe"]) == {0}

    def test_requested_pops_subset(self):
        samples = [
            hl.Struct(s="a", pop="nfe"),
            hl.Struct(s="b", pop="afr"),
            hl.Struct(s="c", pop="eas"),
        ]
        # Restrict to a subset; requested order is preserved, others dropped.
        pops, index_sets, sizes = _build_pop_stratification(
            samples, requested_pops=["eas", "nfe"]
        )
        assert pops == [GLOBAL_POP, "eas", "nfe"]
        assert set(index_sets) == {"eas", "nfe"}
        assert "afr" not in sizes

    def test_requested_pop_absent_is_skipped(self):
        samples = [hl.Struct(s="a", pop="nfe")]
        # 'sas' requested but has no samples → silently skipped.
        pops, index_sets, sizes = _build_pop_stratification(
            samples, requested_pops=["nfe", "sas"]
        )
        assert pops == [GLOBAL_POP, "nfe"]
        assert sizes == {"nfe": 1}


# ===========================================================================
# _count_from_sets_by_pop — per-pop counts additivity
# ===========================================================================

def _enc_variant(het, hv, all_samples, raw_hr_adj_missing, *, adj_het, adj_hv):
    """Build a primary-form encoded-variant struct for _count_from_sets_by_pop.

    Primary form (both is_complement flags False) so membership is explicit.
    n_with_data = |all_samples| + |raw_hr_adj_missing| (= |cats 1-6|); the
    stored ``all_samples`` holds cats 1,3-6 (disjoint from cat 2).
    """
    def _s(xs):
        return hl.literal(set(xs), hl.tset(hl.tint32))

    return hl.struct(
        raw_het=_s(het),
        raw_hv=_s(hv),
        adj_het=_s(adj_het),
        adj_hv=_s(adj_hv),
        all_samples=_s(all_samples),
        n_with_data=hl.int32(len(all_samples) + len(raw_hr_adj_missing)),
        all_samples_is_complement=hl.bool(False),
        raw_hr_adj_missing=_s(raw_hr_adj_missing),
        n_raw_hr_adj_missing=hl.int32(len(raw_hr_adj_missing)),
        raw_hr_adj_missing_is_complement=hl.bool(False),
    )


class TestCountFromSetsByPop:
    """Per-pop counts reuse _count_from_sets over pop-restricted sample sets.

    Every 9-cell entry is a count of samples, so partitioning the cohort into
    disjoint groups that cover all N indices must make the per-pop cells sum
    (element-wise) to the flat full-cohort cells. This is method-independent
    and catches restriction / size / complement-handling bugs.
    """

    N = 10
    # Two variants over a 10-sample cohort (indices 0-9), primary form.
    V1 = _enc_variant(
        het=[1, 2], hv=[3], all_samples=[1, 2, 3, 4],
        raw_hr_adj_missing=[5], adj_het=[1], adj_hv=[3],
    )
    V2 = _enc_variant(
        het=[2, 6], hv=[7], all_samples=[2, 6, 7, 8],
        raw_hr_adj_missing=[9], adj_het=[6], adj_hv=[7],
    )

    def _flat(self, include_raw):
        return _count_from_sets(
            self.V1.raw_het if include_raw else self.V1.adj_het,
            self.V1.raw_hv if include_raw else self.V1.adj_hv,
            self.V1.all_samples, self.V1.n_with_data,
            self.V1.all_samples_is_complement,
            self.V1.raw_hr_adj_missing, self.V1.n_raw_hr_adj_missing,
            self.V1.raw_hr_adj_missing_is_complement,
            self.V2.raw_het if include_raw else self.V2.adj_het,
            self.V2.raw_hv if include_raw else self.V2.adj_hv,
            self.V2.all_samples, self.V2.n_with_data,
            self.V2.all_samples_is_complement,
            self.V2.raw_hr_adj_missing, self.V2.n_raw_hr_adj_missing,
            self.V2.raw_hr_adj_missing_is_complement,
            hl.int32(self.N),
            include_raw_hr_adj_missing=include_raw,
        )

    def _by_pop(self, index_sets, sizes, pops):
        return hl.eval(
            _count_from_sets_by_pop(
                self.V1, self.V2, pops, index_sets, sizes,
                self._flat(True), self._flat(False),
            )
        )

    def test_all_entry_equals_flat(self):
        result = self._by_pop(
            {"even": hl.literal({0, 2, 4, 6, 8}, hl.tset(hl.tint32)),
             "odd": hl.literal({1, 3, 5, 7, 9}, hl.tset(hl.tint32))},
            {"even": 5, "odd": 5},
            [GLOBAL_POP, "even", "odd"],
        )
        flat_raw = list(hl.eval(self._flat(True)))
        flat_adj = list(hl.eval(self._flat(False)))
        assert list(result[GLOBAL_POP].raw) == flat_raw
        assert list(result[GLOBAL_POP].adj) == flat_adj

    def test_partition_sums_to_global(self):
        # even ∪ odd = all 10 indices, disjoint → per-pop cells sum to flat.
        result = self._by_pop(
            {"even": hl.literal({0, 2, 4, 6, 8}, hl.tset(hl.tint32)),
             "odd": hl.literal({1, 3, 5, 7, 9}, hl.tset(hl.tint32))},
            {"even": 5, "odd": 5},
            [GLOBAL_POP, "even", "odd"],
        )
        for field in ("raw", "adj"):
            e = list(result["even"][field])
            o = list(result["odd"][field])
            a = list(result[GLOBAL_POP][field])
            assert [x + y for x, y in zip(e, o)] == a, (
                f"{field}: even {e} + odd {o} != all {a}"
            )

    def test_single_pop_covering_all_equals_global(self):
        result = self._by_pop(
            {"whole": hl.literal(set(range(self.N)), hl.tset(hl.tint32))},
            {"whole": self.N},
            [GLOBAL_POP, "whole"],
        )
        assert list(result["whole"].raw) == list(result[GLOBAL_POP].raw)
        assert list(result["whole"].adj) == list(result[GLOBAL_POP].adj)

    def test_per_pop_adj_sum_below_pop_size(self):
        sizes = {"even": 5, "odd": 5}
        result = self._by_pop(
            {"even": hl.literal({0, 2, 4, 6, 8}, hl.tset(hl.tint32)),
             "odd": hl.literal({1, 3, 5, 7, 9}, hl.tset(hl.tint32))},
            sizes,
            [GLOBAL_POP, "even", "odd"],
        )
        for pop, size in sizes.items():
            cells = list(result[pop].adj)
            assert all(c >= 0 for c in cells), f"{pop} has negative cells: {cells}"
            assert sum(cells) <= size, (
                f"{pop} adj sum {sum(cells)} exceeds pop size {size}"
            )

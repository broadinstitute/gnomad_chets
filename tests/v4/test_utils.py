"""Pytest unit tests for ``gnomad_chets.v4.utils``.

Run with::

    pytest v4/test_utils.py -v

Tests intentionally use tiny in-memory Hail Tables built from
``hl.Table.parallelize`` so the suite runs in a few seconds against a
local Spark backend (no Dataproc required).
"""
import hail as hl
import pytest

from gnomad_chets.v4.resources import (
    CLINVAR_CATEGORY_BLB,
    CLINVAR_CATEGORY_PLP,
    CLINVAR_CATEGORY_VUS,
)
from gnomad_chets.v4.utils import (
    _find_strata_index,
    annotate_filters_af,
    calculate_partitions_by_size,
    clinvar_category_match_expr,
    clinvar_review_flags_expr,
    compute_v2_independent_set,
    filter_clinvar_by_category,
    filter_for_testing,
    get_an_percent_expr,
)


# ---------------------------------------------------------------------------
# Session-scoped Hail init
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def _hail_session():
    """Initialize Hail once per pytest session."""
    hl.init(idempotent=True, quiet=True)
    yield
    # No teardown — Hail handles its own.


# ===========================================================================
# _find_strata_index
# ===========================================================================

class TestFindStrataIndex:
    META = [
        {"group": "raw"},
        {"group": "adj"},
        {"group": "adj", "sex": "XX"},
        {"group": "adj", "sex": "XY"},
        {"group": "adj", "pop": "afr"},
    ]

    def test_finds_first_full_match(self):
        assert _find_strata_index(self.META, {"group": "adj"}) == 1

    def test_partial_match_returns_first(self):
        # {"group": "adj", "sex": "XX"} matches index 2 (full match), not 1.
        assert _find_strata_index(
            self.META, {"group": "adj", "sex": "XX"},
        ) == 2

    def test_match_with_extra_criteria_fields(self):
        assert _find_strata_index(self.META, {"pop": "afr"}) == 4

    def test_no_match_raises(self):
        with pytest.raises(ValueError, match="No strata-meta entry matches"):
            _find_strata_index(self.META, {"pop": "nfe"})

    def test_empty_meta_raises(self):
        with pytest.raises(ValueError):
            _find_strata_index([], {"group": "adj"})

    def test_empty_criteria_returns_first(self):
        # Every entry matches the empty criteria; the first is index 0.
        assert _find_strata_index(self.META, {}) == 0


# ===========================================================================
# calculate_partitions_by_size
# ===========================================================================

def _range_ht_with_size(n_rows, size_per_row):
    """Build a small int64-keyed Table whose ``_size`` is constant."""
    ht = hl.utils.range_table(n_rows).key_by("idx")
    ht = ht.annotate(_size=hl.int64(size_per_row))
    return ht


class TestCalculatePartitionsBySize:

    def test_uniform_rows_produce_n_partitions(self):
        ht = _range_ht_with_size(20, 1)
        intervals = calculate_partitions_by_size(ht, 4, "_size")
        # Uniform row sizes → boundaries placed evenly; first row always
        # produces a boundary, so we get exactly `n_partitions` intervals.
        assert len(intervals) == 4

    def test_total_size_reachable_in_intervals_covers_all_rows(self):
        ht = _range_ht_with_size(10, 100)
        intervals = calculate_partitions_by_size(ht, 3, "_size")
        # Apply the intervals back to the table; every row must survive.
        path = hl.utils.new_temp_file("test_calc_parts", "ht")
        ht.write(path, overwrite=True)
        kept = hl.read_table(path, _intervals=intervals).count()
        assert kept == 10

    def test_array_size_field_uses_length(self):
        # When size_field is an array, the function uses hl.len.
        ht = hl.utils.range_table(6).key_by("idx")
        ht = ht.annotate(_arr=hl.range(0, ht.idx + 1))  # lengths 1,2,3,4,5,6
        intervals = calculate_partitions_by_size(ht, 3, "_arr")
        assert 1 <= len(intervals) <= 3

    def test_single_huge_row_only_one_boundary_at_that_key(self):
        # A row whose size exceeds total/n_partitions cannot be split — it
        # gets one boundary, and the algorithm cannot achieve `n_partitions`.
        ht = hl.utils.range_table(5).key_by("idx")
        ht = ht.annotate(
            _size=hl.if_else(ht.idx == 2, hl.int64(1000), hl.int64(1))
        )
        intervals = calculate_partitions_by_size(ht, 10, "_size")
        # Documented limitation: cannot exceed the number of distinct keys.
        assert len(intervals) <= 5

    def test_weight_field_multiplies_size(self):
        # Row size is multiplied by weight_field length, increasing boundaries.
        ht = hl.utils.range_table(10).key_by("idx")
        ht = ht.annotate(
            _size=hl.int64(1),
            _weight=hl.range(0, hl.if_else(ht.idx == 5, 100, 1)),
        )
        intervals_no_w = calculate_partitions_by_size(ht, 4, "_size")
        intervals_w = calculate_partitions_by_size(
            ht, 4, "_size", weight_field="_weight"
        )
        # With weight, row 5 dominates → boundary placement shifts.
        # Both runs must still cover every row.
        path = hl.utils.new_temp_file("test_calc_parts_weight", "ht")
        ht.write(path, overwrite=True)
        assert hl.read_table(path, _intervals=intervals_no_w).count() == 10
        assert hl.read_table(path, _intervals=intervals_w).count() == 10

    def test_weight_from_separate_table(self):
        ht = hl.utils.range_table(8).key_by("idx")
        ht = ht.annotate(_size=hl.int64(2))
        weight_ht = hl.utils.range_table(8).key_by("idx")
        weight_ht = weight_ht.annotate(
            _w=hl.range(0, hl.if_else(weight_ht.idx == 3, 50, 1))
        )
        intervals = calculate_partitions_by_size(
            ht, 4, "_size", weight_ht=weight_ht, weight_field="_w",
        )
        path = hl.utils.new_temp_file("test_calc_parts_weight_ht", "ht")
        ht.write(path, overwrite=True)
        assert hl.read_table(path, _intervals=intervals).count() == 8


# ===========================================================================
# compute_v2_independent_set
# ===========================================================================

class TestComputeV2IndependentSet:

    def test_empty(self):
        assert compute_v2_independent_set([], []) == set()

    def test_no_edges_picks_everything(self):
        # With no edges, every variant can go to V2 (no neighbors → no
        # conflicts). The greedy algorithm picks them all.
        result = compute_v2_independent_set(
            [(0, 100), (1, 100), (2, 100)], [],
        )
        assert result == {0, 1, 2}

    def test_one_edge_picks_one_endpoint(self):
        # Edge between 0 and 1 forbids putting both in V2; greedy picks one.
        result = compute_v2_independent_set(
            [(0, 100), (1, 100)], [(0, 1)],
        )
        assert len(result) == 1
        assert result.issubset({0, 1})

    def test_triangle_picks_one(self):
        # K3: edges (0-1), (1-2), (0-2). Any one vertex is independent.
        result = compute_v2_independent_set(
            [(0, 100), (1, 100), (2, 100)],
            [(0, 1), (1, 2), (0, 2)],
        )
        assert len(result) == 1

    def test_greedy_prefers_cheap_variant(self):
        # Variant 0 has lowest n_with_data → smallest score → picked first.
        result = compute_v2_independent_set(
            [(0, 1), (1, 100), (2, 100)], [(0, 1), (0, 2)],
        )
        # 0 is the cheapest and a hub; picking it forbids 1, 2.
        assert result == {0}

    def test_degree_breaks_tie(self):
        # Same n_with_data but variant 0 has higher degree → lower score
        # (n / max(degree, 1)) → picked first.
        result = compute_v2_independent_set(
            [(0, 100), (1, 100), (2, 100)],
            [(0, 1), (0, 2)],  # 0 has degree 2; 1, 2 have degree 1.
        )
        # 0 dominates → picked first → forbids 1, 2.
        assert result == {0}


# ===========================================================================
# clinvar_category_match_expr / filter_clinvar_by_category
# ===========================================================================

def _make_clinvar_ht(rows):
    """Build a tiny ClinVar-shaped Table from a list of dicts."""
    typed_rows = []
    for r in rows:
        typed_rows.append(
            hl.Struct(
                locus=hl.locus(r["chr"], r["pos"], "GRCh38"),
                alleles=r["alleles"],
                info=hl.Struct(
                    CLNSIG=r.get("clnsig", []),
                    CLNREVSTAT=r.get("clnrevstat", []),
                    CLNSIGCONF=r.get("clnsigconf", hl.missing(hl.tstr)),
                ),
            )
        )
    ht = hl.Table.parallelize(
        typed_rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            info=hl.tstruct(
                CLNSIG=hl.tarray(hl.tstr),
                CLNREVSTAT=hl.tarray(hl.tstr),
                CLNSIGCONF=hl.tstr,
            ),
        ),
        key=["locus", "alleles"],
    )
    return ht


_REVSTAT_GOOD = ["reviewed_by_expert_panel"]
_REVSTAT_BAD = ["no_assertion_provided"]


class TestClinvarCategoryMatchExpr:

    def test_unknown_category_raises(self):
        ht = _make_clinvar_ht([{
            "chr": "chr1", "pos": 100, "alleles": ["A", "T"],
            "clnsig": ["Pathogenic"], "clnrevstat": _REVSTAT_GOOD,
        }])
        with pytest.raises(ValueError, match="Unknown ClinVar category"):
            clinvar_category_match_expr(
                clnsig=ht.info.CLNSIG, category="bogus",
                clnrevstat=ht.info.CLNREVSTAT, clnsigconf=ht.info.CLNSIGCONF,
            )

    def test_missing_clnrevstat_raises_when_required(self):
        ht = _make_clinvar_ht([{
            "chr": "chr1", "pos": 100, "alleles": ["A", "T"],
            "clnsig": ["Pathogenic"], "clnrevstat": _REVSTAT_GOOD,
        }])
        with pytest.raises(ValueError, match="clnrevstat is required"):
            clinvar_category_match_expr(
                clnsig=ht.info.CLNSIG, category=CLINVAR_CATEGORY_PLP,
                clnrevstat=None, clnsigconf=ht.info.CLNSIGCONF,
                remove_no_assertion=True,
            )

    def test_missing_clnsigconf_raises_when_required(self):
        ht = _make_clinvar_ht([{
            "chr": "chr1", "pos": 100, "alleles": ["A", "T"],
            "clnsig": ["Pathogenic"], "clnrevstat": _REVSTAT_GOOD,
        }])
        with pytest.raises(ValueError, match="clnsigconf is required"):
            clinvar_category_match_expr(
                clnsig=ht.info.CLNSIG, category=CLINVAR_CATEGORY_PLP,
                clnrevstat=ht.info.CLNREVSTAT, clnsigconf=None,
                remove_no_assertion=False, remove_conflicting=True,
            )


class TestFilterClinvarByCategory:

    @pytest.fixture(scope="class")
    def fixture_ht(self):
        return _make_clinvar_ht([
            # 0: PLP, good review.
            {"chr": "chr1", "pos": 100, "alleles": ["A", "T"],
             "clnsig": ["Pathogenic"], "clnrevstat": _REVSTAT_GOOD},
            # 1: Likely_pathogenic, good review.
            {"chr": "chr1", "pos": 101, "alleles": ["A", "C"],
             "clnsig": ["Likely_pathogenic"], "clnrevstat": _REVSTAT_GOOD},
            # 2: BLB.
            {"chr": "chr1", "pos": 102, "alleles": ["A", "G"],
             "clnsig": ["Benign"], "clnrevstat": _REVSTAT_GOOD},
            # 3: VUS.
            {"chr": "chr1", "pos": 103, "alleles": ["A", "T"],
             "clnsig": ["Uncertain_significance"], "clnrevstat": _REVSTAT_GOOD},
            # 4: Mixed B/P (excluded from BLB).
            {"chr": "chr1", "pos": 104, "alleles": ["A", "T"],
             "clnsig": ["Benign", "Pathogenic"], "clnrevstat": _REVSTAT_GOOD},
            # 5: PLP but no-star review (excluded by default).
            {"chr": "chr1", "pos": 105, "alleles": ["A", "T"],
             "clnsig": ["Pathogenic"], "clnrevstat": _REVSTAT_BAD},
            # 6: PLP with conflicting interpretation (excluded by default).
            {"chr": "chr1", "pos": 106, "alleles": ["A", "T"],
             "clnsig": ["Pathogenic"], "clnrevstat": _REVSTAT_GOOD,
             "clnsigconf": "Pathogenic(2),Likely_benign(1)"},
        ])

    def test_plp_excludes_no_star_and_conflicting(self, fixture_ht):
        result = filter_clinvar_by_category(fixture_ht, CLINVAR_CATEGORY_PLP)
        positions = set(result.locus.position.collect())
        assert positions == {100, 101, 104}

    def test_blb_excludes_mixed(self, fixture_ht):
        result = filter_clinvar_by_category(fixture_ht, CLINVAR_CATEGORY_BLB)
        positions = set(result.locus.position.collect())
        # 102 only; 104 has both Benign + Pathogenic.
        assert positions == {102}

    def test_vus(self, fixture_ht):
        result = filter_clinvar_by_category(fixture_ht, CLINVAR_CATEGORY_VUS)
        positions = set(result.locus.position.collect())
        assert positions == {103}

    def test_keep_no_assertion(self, fixture_ht):
        result = filter_clinvar_by_category(
            fixture_ht, CLINVAR_CATEGORY_PLP, remove_no_assertion=False,
        )
        positions = set(result.locus.position.collect())
        # Adds back the no-star row (105).
        assert positions == {100, 101, 104, 105}

    def test_keep_conflicting(self, fixture_ht):
        result = filter_clinvar_by_category(
            fixture_ht, CLINVAR_CATEGORY_PLP, remove_conflicting=False,
        )
        positions = set(result.locus.position.collect())
        # Adds back the conflicting row (106).
        assert positions == {100, 101, 104, 106}


# ===========================================================================
# annotate_filters_af
# ===========================================================================

def _make_locus_alleles_ht(rows):
    """Build a tiny key-only Table keyed by (locus, alleles).

    Each input row dict needs ``chr``, ``pos``, ``alleles``.
    """
    typed = [
        hl.Struct(
            locus=hl.locus(r["chr"], r["pos"], "GRCh38"),
            alleles=r["alleles"],
        )
        for r in rows
    ]
    return hl.Table.parallelize(
        typed,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
        ),
        key=["locus", "alleles"],
    )


class TestClinvarReviewFlags:

    def _flags(self, clnrevstat, clnsigconf):
        return set(hl.eval(clinvar_review_flags_expr(
            hl.literal(clnrevstat, hl.tarray(hl.tstr)),
            hl.literal(clnsigconf, hl.tstr) if clnsigconf is not None
            else hl.missing(hl.tstr),
        )))

    def test_clean_record_no_flags(self):
        assert self._flags(_REVSTAT_GOOD, None) == set()

    def test_no_assertion_flag(self):
        assert self._flags(_REVSTAT_BAD, None) == {"no_assertion"}

    def test_conflicting_flag(self):
        assert self._flags(_REVSTAT_GOOD, "Pathogenic(2),Likely_benign(1)") == {
            "conflicting"
        }

    def test_both_flags(self):
        assert self._flags(_REVSTAT_BAD, "Pathogenic(1),Benign(1)") == {
            "no_assertion",
            "conflicting",
        }


class TestAnnotateFiltersAf:

    def test_annotates_filters_and_af(self):
        base = _make_locus_alleles_ht([
            {"chr": "chr1", "pos": 100, "alleles": ["A", "T"]},
            {"chr": "chr1", "pos": 101, "alleles": ["A", "C"]},
        ])
        filter_ht = hl.Table.parallelize(
            [
                hl.Struct(
                    locus=hl.locus("chr1", 100, "GRCh38"),
                    alleles=["A", "T"],
                    filters=hl.empty_set(hl.tstr),
                ),
                hl.Struct(
                    locus=hl.locus("chr1", 101, "GRCh38"),
                    alleles=["A", "C"],
                    filters={"InbreedingCoeff"},
                ),
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                filters=hl.tset(hl.tstr),
            ),
            key=["locus", "alleles"],
        )
        freq_ht = hl.Table.parallelize(
            [
                hl.Struct(
                    locus=hl.locus("chr1", 100, "GRCh38"),
                    alleles=["A", "T"],
                    freq=[hl.Struct(AF=0.001, AN=1000)],
                ),
                hl.Struct(
                    locus=hl.locus("chr1", 101, "GRCh38"),
                    alleles=["A", "C"],
                    freq=[hl.Struct(AF=0.01, AN=2000)],
                ),
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                freq=hl.tarray(hl.tstruct(AF=hl.tfloat64, AN=hl.tint32)),
            ),
            key=["locus", "alleles"],
        )
        out = annotate_filters_af(base, filter_ht, freq_ht).collect()
        out_by_pos = {r.locus.position: r for r in out}
        assert out_by_pos[100].af == pytest.approx(0.001)
        assert out_by_pos[100].filters == set()
        assert out_by_pos[101].af == pytest.approx(0.01)
        assert out_by_pos[101].filters == {"InbreedingCoeff"}
        # `an` not included by default.
        assert not hasattr(out_by_pos[100], "an")

    def test_include_an(self):
        base = _make_locus_alleles_ht([
            {"chr": "chr1", "pos": 100, "alleles": ["A", "T"]},
        ])
        filter_ht = hl.Table.parallelize(
            [hl.Struct(
                locus=hl.locus("chr1", 100, "GRCh38"),
                alleles=["A", "T"],
                filters=hl.empty_set(hl.tstr),
            )],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                filters=hl.tset(hl.tstr),
            ),
            key=["locus", "alleles"],
        )
        freq_ht = hl.Table.parallelize(
            [hl.Struct(
                locus=hl.locus("chr1", 100, "GRCh38"),
                alleles=["A", "T"],
                freq=[hl.Struct(AF=0.005, AN=1500)],
            )],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                freq=hl.tarray(hl.tstruct(AF=hl.tfloat64, AN=hl.tint32)),
            ),
            key=["locus", "alleles"],
        )
        out = annotate_filters_af(
            base, filter_ht, freq_ht, include_an=True,
        ).collect()[0]
        assert out.an == 1500


# ===========================================================================
# get_an_percent_expr
# ===========================================================================

class TestGetAnPercentExpr:

    def _make_an_ht(self, an_value_adj, strata_sample_count, strata_meta):
        """Tiny AN HT with one autosomal locus."""
        rows = [
            hl.Struct(
                locus=hl.locus("chr1", 100, "GRCh38"),
                AN=hl.literal(an_value_adj),
            )
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                AN=hl.tarray(hl.tint32),
            ),
            key=["locus"],
        )
        return ht.annotate_globals(
            strata_sample_count=strata_sample_count,
            strata_meta=strata_meta,
        )

    def _eval_pct(self, an_ht, chrom, pos):
        # get_an_percent_expr indexes an_ht by locus; index with a locus
        # *column* (a scalar hl.locus expression can't index a Table).
        q = hl.utils.range_table(1)
        q = q.annotate(locus=hl.locus(chrom, pos, "GRCh38"))
        return q.annotate(pct=get_an_percent_expr(an_ht, q.locus)).pct.collect()[0]

    def test_autosome_full_an(self):
        # adj_count = 100 samples; max_AN = 200; AN=200 → 100%.
        meta = [
            {"group": "raw"},
            {"group": "adj"},
            {"group": "adj", "sex": "XX"},
            {"group": "adj", "sex": "XY"},
        ]
        ht = self._make_an_ht(
            an_value_adj=[400, 200, 90, 110],
            strata_sample_count=[100, 100, 45, 55],
            strata_meta=meta,
        )
        pct = self._eval_pct(ht, "chr1", 100)
        # AN[adj=1] = 200, total = adj_count * 2 = 200 → 100%
        assert pct == 100

    def test_autosome_half_an(self):
        meta = [
            {"group": "raw"},
            {"group": "adj"},
            {"group": "adj", "sex": "XX"},
            {"group": "adj", "sex": "XY"},
        ]
        ht = self._make_an_ht(
            an_value_adj=[200, 100, 50, 50],
            strata_sample_count=[100, 100, 50, 50],
            strata_meta=meta,
        )
        pct = self._eval_pct(ht, "chr1", 100)
        # AN[adj] = 100, total = 200 → 50%
        assert pct == 50

    def test_missing_adj_strata_raises(self):
        # No adj group at all -> the adj strata lookup raises at build time,
        # before any indexing. (A meta with only sex-split adj entries would
        # partial-match {group: adj}, so it must omit adj entirely.)
        meta = [{"group": "raw"}]
        ht = self._make_an_ht(
            an_value_adj=[200],
            strata_sample_count=[100],
            strata_meta=meta,
        )
        with pytest.raises(ValueError, match="No strata-meta entry matches"):
            get_an_percent_expr(ht, hl.locus("chr1", 100, "GRCh38"))


# ===========================================================================
# filter_for_testing
# ===========================================================================

class TestFilterForTesting:

    def test_table_filtered_to_interval(self):
        # Build a small Table spanning multiple loci; filter to one CAPN3
        # interval defined in TEST_INTERVALS.
        from gnomad_chets.v4.resources import TEST_INTERVALS
        # Use any gene listed in TEST_INTERVALS; pick CAPN3 from
        # resources to make the test self-contained.
        capn3_interval = TEST_INTERVALS.get("CAPN3")
        if capn3_interval is None:
            pytest.skip("TEST_INTERVALS doesn't define CAPN3 — skipping")
        capn3_loc = hl.parse_locus_interval(capn3_interval, reference_genome="GRCh38")
        start_pos = hl.eval(capn3_loc.start.position)
        chrom = hl.eval(capn3_loc.start.contig)

        rows = [
            # Inside CAPN3.
            hl.Struct(
                locus=hl.locus(chrom, start_pos + 100, "GRCh38"),
                alleles=["A", "T"],
            ),
            # Outside CAPN3 (different chrom).
            hl.Struct(
                locus=hl.locus("chr1", 100, "GRCh38"),
                alleles=["A", "T"],
            ),
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
            ),
            key=["locus", "alleles"],
        )
        # Limit to just CAPN3 to make the test deterministic regardless of
        # other genes in TEST_INTERVALS.
        result = filter_for_testing(ht, test_intervals={"CAPN3": capn3_interval})
        positions = result.locus.position.collect()
        assert len(positions) == 1
        assert positions[0] == start_pos + 100

    def test_unsupported_type_raises(self):
        with pytest.raises(ValueError, match="Unsupported type"):
            filter_for_testing("not a table")  # type: ignore[arg-type]

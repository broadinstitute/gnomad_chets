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
    _build_pop_stratification,
    _COUNT_FROM_SETS_FIELDS,
    _count_from_sets,
    _count_from_sets_by_pop,
    _count_phase_from_sets,
    _create_var_idx_ht,
    _encode_genotype_sets_by_var_idx,
    _drop_pairs_missing_v_idx,
    _empty_counts_ht,
    _no_pbt_count_fields,
    _pbt_index_set_for,
    _pop_restrict_variant,
    _project_count_fields,
    _read_min_an_pct,
    filter_pairs_by_an_pct,
    restrict_encoded_to_pops,
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
                phased_het=hl.empty_dict(
                    hl.tint32, hl.tstruct(pid=hl.tstr, gt0=hl.tint32)
                ),
                # Encoder-only extra that should be dropped:
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
                phased_het=hl.tdict(
                    hl.tint32, hl.tstruct(pid=hl.tstr, gt0=hl.tint32)
                ),
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
        assert "_extra_diagnostics" not in projected.row

    def test_rows_preserved(self):
        encoded = self._make_encoded_ht()
        assert _project_count_fields(encoded).count() == encoded.count()


# ===========================================================================
# _count_phase_from_sets — physical-phase refinement of the AaBb cell
# ===========================================================================

class TestCountPhaseFromSets:
    """cis / trans classification of double-het samples from phase sidecars."""

    @staticmethod
    def _set(xs):
        return hl.set([hl.int32(x) for x in xs])

    @staticmethod
    def _pdict(d):
        # d: {sample_idx: (pid, gt0)}
        if not d:
            return hl.empty_dict(
                hl.tint32, hl.tstruct(pid=hl.tstr, gt0=hl.tint32)
            )
        return hl.dict(
            [
                (hl.int32(k), hl.struct(pid=pid, gt0=hl.int32(gt0)))
                for k, (pid, gt0) in d.items()
            ]
        )

    def _count(self, v1_het, v1_phase, v2_het, v2_phase):
        return hl.eval(
            _count_phase_from_sets(
                self._set(v1_het), self._pdict(v1_phase),
                self._set(v2_het), self._pdict(v2_phase),
            )
        )

    def test_cis_same_pid_same_gt0(self):
        r = self._count([5], {5: ("p1", 0)}, [5], {5: ("p1", 0)})
        assert (r.n_phased_cis, r.n_phased_trans) == (1, 0)

    def test_trans_same_pid_diff_gt0(self):
        r = self._count([5], {5: ("p1", 0)}, [5], {5: ("p1", 1)})
        assert (r.n_phased_cis, r.n_phased_trans) == (0, 1)

    def test_different_pid_is_unphased(self):
        # Same sample phased on both sides but in different phase sets → neither.
        r = self._count([5], {5: ("p1", 0)}, [5], {5: ("p2", 0)})
        assert (r.n_phased_cis, r.n_phased_trans) == (0, 0)

    def test_absent_from_dict_is_unphased(self):
        # Double-het but no phase recorded → excluded from cis/trans.
        r = self._count([5], {}, [5], {})
        assert (r.n_phased_cis, r.n_phased_trans) == (0, 0)

    def test_not_double_het_excluded(self):
        # Het at v1 only (disjoint het sets) → no double-het to classify.
        r = self._count([5], {5: ("p1", 0)}, [6], {6: ("p1", 0)})
        assert (r.n_phased_cis, r.n_phased_trans) == (0, 0)

    def test_mixed_cis_trans_unphased(self):
        v1p = {1: ("p", 0), 2: ("p", 0), 3: ("p", 0)}
        v2p = {1: ("p", 0), 2: ("p", 1)}  # 1 cis, 2 trans, 3 unphased (absent)
        r = self._count([1, 2, 3], v1p, [1, 2, 3], v2p)
        assert (r.n_phased_cis, r.n_phased_trans) == (1, 1)

    def test_one_side_phased_other_not(self):
        # Phased at v1, not at v2 → not counted (needs both sides).
        r = self._count([5], {5: ("p1", 0)}, [5], {})
        assert (r.n_phased_cis, r.n_phased_trans) == (0, 0)


# ===========================================================================
# restrict_encoded_to_pops — per-pop restriction (re-index incl. phased_het)
# ===========================================================================

class TestRestrictEncodedToPops:
    """Per-pop restriction must re-index every set AND keep phased_het present
    and re-indexed, so the downstream _project_count_fields (which selects
    phased_het) doesn't crash. Regression for the --stratify-by-pop path.
    """

    _PHASE_T = hl.tstruct(pid=hl.tstr, gt0=hl.tint32)

    def _encoded_ht(self):
        # 4-sample cohort s0..s3 (indices 0-3), one variant. Sample 0 is
        # no-entry (in all_samples), 1/2 het, 3 hom-var; 1 & 2 phased.
        row = hl.Struct(
            v_idx=hl.int64(0),
            raw_het=hl.set([hl.int32(1), hl.int32(2)]),
            raw_hv=hl.set([hl.int32(3)]),
            adj_het=hl.set([hl.int32(1), hl.int32(2)]),
            adj_hv=hl.set([hl.int32(3)]),
            all_samples=hl.set([hl.int32(i) for i in (0, 1, 2, 3)]),
            n_with_data=hl.int32(4),
            raw_hr_adj_missing=hl.empty_set(hl.tint32),
            n_raw_hr_adj_missing=hl.int32(0),
            phased_het=hl.dict([
                (hl.int32(1), hl.struct(pid="p", gt0=hl.int32(0))),
                (hl.int32(2), hl.struct(pid="p", gt0=hl.int32(1))),
            ]),
        )
        ht = hl.Table.parallelize(
            [row],
            hl.tstruct(
                v_idx=hl.tint64,
                raw_het=hl.tset(hl.tint32), raw_hv=hl.tset(hl.tint32),
                adj_het=hl.tset(hl.tint32), adj_hv=hl.tset(hl.tint32),
                all_samples=hl.tset(hl.tint32), n_with_data=hl.tint32,
                raw_hr_adj_missing=hl.tset(hl.tint32),
                n_raw_hr_adj_missing=hl.tint32,
                phased_het=hl.tdict(hl.tint32, self._PHASE_T),
            ),
            key=["v_idx"],
        )
        return ht.annotate_globals(
            samples=hl.literal(
                [hl.Struct(s=f"s{i}") for i in range(4)],
                hl.tarray(hl.tstruct(s=hl.tstr)),
            )
        )

    def _pop_ht(self):
        # s1, s3 → nfe (kept); s0, s2 → afr (dropped when restricting to nfe).
        return hl.Table.parallelize(
            [hl.Struct(s=f"s{i}", pop=p)
             for i, p in enumerate(["afr", "nfe", "afr", "nfe"])],
            hl.tstruct(s=hl.tstr, pop=hl.tstr),
            key=["s"],
        )

    def test_phased_het_reindexed_and_kept(self):
        out = restrict_encoded_to_pops(self._encoded_ht(), self._pop_ht(), ["nfe"])
        r = out.collect()[0]
        # keep=[1,3] → remap {1:0, 3:1}. Set/dict keys re-indexed; kept-out
        # keys dropped. phased_het: key 1 kept→0, key 2 dropped.
        assert dict(r.phased_het) == {0: hl.Struct(pid="p", gt0=0)}
        assert set(r.all_samples) == {0, 1}        # {1,3} → {0,1}
        assert set(r.adj_het) == {0}               # {1} kept (2 dropped) → {0}
        assert set(r.adj_hv) == {1}                # {3} → {1}
        assert r.n_with_data == 2
        assert hl.eval(out.index_globals().samples) == [
            hl.Struct(s="s1"), hl.Struct(s="s3"),
        ]

    def test_project_count_fields_survives_restriction(self):
        # The regression: restricted encoding must still carry every
        # _COUNT_FROM_SETS_FIELDS (incl. phased_het) so projection doesn't fail.
        out = restrict_encoded_to_pops(self._encoded_ht(), self._pop_ht(), ["nfe"])
        projected = _project_count_fields(out)
        kept = set(projected.row) - set(projected.key)
        assert kept == set(_COUNT_FROM_SETS_FIELDS)


# ===========================================================================
# high-AB het -> hom-alt correction in _encode_genotype_sets_by_var_idx
# ===========================================================================

class TestHighAbHetCorrection:
    """v4 GATK<4.1.4.1 high-AB het correction: an adj het-ref (AB>0.9), not
    het-non-ref, unfixed-model, at AF>0.01, is reclassified hom-var in BOTH raw
    and adj (matching the release's ab_adjusted_freq, which corrects freq[1]=raw
    too). The raw reclassification is gated on adj-pass. Exemptions
    (het_non_ref / fixed_homalt_model / AF) must leave it a het.
    """

    @staticmethod
    def _encode_one_variant(af, samples):
        # samples: list of dict(gt, ad, het_nr, fixed[, gq, dp]); one variant.
        # gt=None -> missing GT; a bare "0"/"1" -> haploid call. GQ/DP default
        # 50/20 (adj-passing) unless overridden per sample.
        n = len(samples)
        mt = hl.utils.range_matrix_table(1, n)
        mt = mt.annotate_rows(
            locus=hl.locus("chr1", 100, "GRCh38"),
            alleles=["A", "T"],
            af=hl.float64(af),
        ).key_rows_by("locus", "alleles").drop("row_idx")
        gt = hl.literal([s["gt"] for s in samples], hl.tarray(hl.tstr))
        ad = hl.literal([s["ad"] for s in samples])
        hnr = hl.literal([s["het_nr"] for s in samples])
        fixed = hl.literal([s["fixed"] for s in samples])
        gq = hl.literal([s.get("gq", 50) for s in samples])
        dp = hl.literal([s.get("dp", 20) for s in samples])
        mt = mt.annotate_cols(
            s=hl.str(mt.col_idx), fixed_homalt_model=fixed[mt.col_idx],
        )
        mt = mt.annotate_entries(
            GT=hl.parse_call(gt[mt.col_idx]),
            GQ=hl.int32(gq[mt.col_idx]),
            DP=hl.int32(dp[mt.col_idx]),
            AD=ad[mt.col_idx].map(hl.int32),
            _het_non_ref=hnr[mt.col_idx],
        ).key_cols_by("s")
        enc = _encode_genotype_sets_by_var_idx(mt, _create_var_idx_ht(mt))
        return enc.collect()[0]

    def test_correction_and_exemptions(self):
        # 0: high-AB het, eligible -> hom-var in raw AND adj
        # 1: high-AB het but het_non_ref -> het (exempt)
        # 2: high-AB het but fixed_homalt_model -> het (exempt)
        # 3: normal het (AB 0.5) -> het
        # 4: hom-var -> hom-var
        r = self._encode_one_variant(0.02, [
            {"gt": "0/1", "ad": [1, 19], "het_nr": False, "fixed": False},
            {"gt": "0/1", "ad": [1, 19], "het_nr": True, "fixed": False},
            {"gt": "0/1", "ad": [1, 19], "het_nr": False, "fixed": True},
            {"gt": "0/1", "ad": [10, 10], "het_nr": False, "fixed": False},
            {"gt": "1/1", "ad": [0, 20], "het_nr": False, "fixed": False},
        ])
        # Sample 0 is now corrected out of the het set and into the hom set in
        # BOTH raw and adj (adj-passing, so the raw gate lets it through).
        assert set(r.raw_het) == {1, 2, 3}
        assert set(r.raw_hv) == {0, 4}
        assert set(r.adj_het) == {1, 2, 3}
        assert set(r.adj_hv) == {0, 4}

    def test_raw_correction_gated_on_adj_pass(self):
        # High-AB het that FAILS adj (low GQ): not in gnomAD's (adj-gated)
        # correction set, so it stays a het in RAW and is dropped from the adj
        # carrier sets entirely (adj_gt=0).
        r = self._encode_one_variant(0.02, [
            {"gt": "0/1", "ad": [1, 19], "het_nr": False, "fixed": False, "gq": 5},
        ])
        assert set(r.raw_het) == {0}      # stays raw het (adj-fail => not corrected)
        assert set(r.raw_hv) == set()
        assert set(r.adj_het) == set()    # adj-fail => not an adj carrier
        assert set(r.adj_hv) == set()

    def test_af_gate_below_threshold_not_corrected(self):
        # Same high-AB het but the variant's AF <= threshold -> not corrected
        # (het in both raw and adj).
        r = self._encode_one_variant(0.005, [
            {"gt": "0/1", "ad": [1, 19], "het_nr": False, "fixed": False},
        ])
        assert set(r.adj_het) == {0}
        assert set(r.adj_hv) == set()
        assert set(r.raw_het) == {0}
        assert set(r.raw_hv) == set()


class TestSexPloidyClassification:
    """The encoder's het/hom classifier handles the outputs of
    ``adjusted_sex_ploidy_expr`` (applied upstream in densify_encode_input_mt)
    with no special-casing: a hemizygous alt (haploid ``Call(1)``) -> hom cell,
    a haploid ref -> hom-ref majority (in no stored set), and a dropped het
    (missing GT on a defined entry) -> uncallable (in all_samples/n_with_data
    but in no genotype cell, so excluded from AABB).
    """

    @staticmethod
    def _encode(samples):
        return TestHighAbHetCorrection._encode_one_variant(0.02, samples)

    def test_haploid_alt_is_hom(self):
        # 0: hemizygous alt (haploid) -> hom; 1: diploid hom-var (control).
        r = self._encode([
            {"gt": "1", "ad": [0, 20], "het_nr": False, "fixed": False},
            {"gt": "1/1", "ad": [0, 20], "het_nr": False, "fixed": False},
        ])
        assert set(r.raw_hv) == {0, 1}
        assert set(r.adj_hv) == {0, 1}
        assert set(r.raw_het) == set()

    def test_haploid_ref_and_missing_gt(self):
        # 0: haploid ref -> hom-ref majority (not a carrier, NOT in all_samples).
        # 1: missing GT, defined entry (dropped chrX XY het) -> uncallable:
        #    in all_samples + n_with_data but in no genotype cell.
        # 2: diploid het (control) -> raw/adj het, in all_samples.
        r = self._encode([
            {"gt": "0", "ad": [20, 0], "het_nr": False, "fixed": False},
            {"gt": None, "ad": [0, 0], "het_nr": False, "fixed": False},
            {"gt": "0/1", "ad": [10, 10], "het_nr": False, "fixed": False},
        ])
        assert set(r.raw_het) == {2}
        # haploid ref: hom-ref majority, in no stored set.
        assert 0 not in set(r.all_samples)
        assert 0 not in set(r.raw_hv) and 0 not in set(r.adj_hv)
        # missing GT: uncallable -> in all_samples but in no genotype cell.
        assert 1 in set(r.all_samples)
        assert 1 not in set(r.raw_het) and 1 not in set(r.raw_hv)
        assert 1 not in set(r.adj_het) and 1 not in set(r.adj_hv)


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
# filter_pairs_by_an_pct
# ===========================================================================

class TestFilterPairsByAnPct:

    def _build_an_annotated_pairs(self, rows):
        """Build a pair table with an_pct1 / an_pct2 fields (as create_vp_list
        produces them; filter_pairs_by_an_pct reads those names)."""
        structs = [
            hl.Struct(
                locus1=hl.locus("chr1", r["p1"], "GRCh38"),
                alleles1=["A", "T"],
                locus2=hl.locus("chr1", r["p2"], "GRCh38"),
                alleles2=["A", "G"],
                an_pct1=hl.int32(r["v1_pct"]),
                an_pct2=hl.int32(r["v2_pct"]),
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
                an_pct1=hl.tint32,
                an_pct2=hl.tint32,
            ),
            key=["locus1", "alleles1", "locus2", "alleles2"],
        )
        return ht

    def test_keeps_when_both_above_floor(self):
        # Floor is exclusive (an_pct > min); both sides strictly above 80.
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 95, "v2_pct": 90},
            {"p1": 101, "p2": 201, "v1_pct": 85, "v2_pct": 85},
        ])
        assert filter_pairs_by_an_pct(ht, 80).count() == 2

    def test_drops_when_v1_below_floor(self):
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 90, "v2_pct": 90},
            {"p1": 101, "p2": 201, "v1_pct": 70, "v2_pct": 95},  # v1 fails
        ])
        out = filter_pairs_by_an_pct(ht, 80)
        positions = out.locus1.position.collect()
        assert 100 in positions
        assert 101 not in positions

    def test_drops_when_v2_below_floor(self):
        ht = self._build_an_annotated_pairs([
            {"p1": 100, "p2": 200, "v1_pct": 95, "v2_pct": 75},  # v2 fails
            {"p1": 101, "p2": 201, "v1_pct": 95, "v2_pct": 90},
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
# set of validated (v1_encoded, v2_encoded, expected_gt_counts) rows pulled
# from the chr19 test pipeline. The sets are stored complement-aware (with the
# *_is_complement flags); _count_from_sets now takes positive form only, so
# _real_data_results normalizes each complement-stored set to positive form
# (N \ stored) — done as a DATA-plane transform over the fixture Table's
# columns (never hl.literal-ing the big sets, which would blow up codegen for
# low-AN rows whose positive all_samples ≈ n_samples), running the kernel over
# the whole table in one job and collecting only the tiny 9-int result arrays.
# Expected counts are storage-form-independent, so exact matches still
# regression-test the chr19 counts.

_REAL_DATA_FIXTURE_PATH = (
    "gs://gnomad-tmp-30day/test_fixtures/count_from_sets_chr19_real.ht"
)


def _row_id(row):
    """Stable id for a fixture row: 'chr19:pos:ref:alt|chr19:pos:ref:alt'."""
    def _v(locus, alleles):
        return f"{locus.contig}:{locus.position}:{'/'.join(alleles)}"
    return f"{_v(row.locus1, row.alleles1)}|{_v(row.locus2, row.alleles2)}"


def _real_data_reachable():
    """Cheap reachability probe (keys only) so the class skips when the GCS
    fixture isn't available, without collecting the big sample sets."""
    try:
        hl.read_table(_REAL_DATA_FIXTURE_PATH).select()._force_count()
        return True
    except Exception:
        return False


_REAL_DATA_AVAILABLE = _real_data_reachable()


@pytest.fixture(scope="session")
def _real_data_results():
    """Run _count_from_sets on every fixture row in ONE data-plane Hail job.

    Returns ``(rows, n_samples)`` where each collected row carries the computed
    ``_got_raw`` / ``_got_adj`` 9-cell arrays plus the fixture's expected
    counts + pair key. The big sample sets stay in the data plane (streamed),
    so complement-stored low-AN rows don't trigger literal-codegen blowup.
    """
    ht = hl.read_table(_REAL_DATA_FIXTURE_PATH)
    n_samples = int(hl.eval(ht.index_globals().n_samples))
    ns = hl.int32(n_samples)
    full = hl.set(hl.range(0, ns))

    def _pos(stored, is_comp):
        return hl.if_else(is_comp, full.difference(stored), stored)

    def _row_counts(v1, v2, include_raw):
        return _count_from_sets(
            v1.raw_het if include_raw else v1.adj_het,
            v1.raw_hv if include_raw else v1.adj_hv,
            _pos(v1.all_samples, v1.all_samples_is_complement),
            hl.int32(v1.n_with_data),
            _pos(v1.raw_hr_adj_missing, v1.raw_hr_adj_missing_is_complement),
            hl.int32(v1.n_raw_hr_adj_missing),
            v2.raw_het if include_raw else v2.adj_het,
            v2.raw_hv if include_raw else v2.adj_hv,
            _pos(v2.all_samples, v2.all_samples_is_complement),
            hl.int32(v2.n_with_data),
            _pos(v2.raw_hr_adj_missing, v2.raw_hr_adj_missing_is_complement),
            hl.int32(v2.n_raw_hr_adj_missing),
            ns,
            include_raw_hr_adj_missing=include_raw,
        )

    ht = ht.annotate(
        _got_raw=_row_counts(ht.v1_encoded, ht.v2_encoded, True),
        _got_adj=_row_counts(ht.v1_encoded, ht.v2_encoded, False),
    )
    rows = ht.key_by().select(
        "locus1", "alleles1", "locus2", "alleles2",
        "_got_raw", "_got_adj",
        "expected_gt_counts_raw", "expected_gt_counts_adj",
    ).collect()
    return rows, n_samples


@pytest.mark.skipif(
    not _REAL_DATA_AVAILABLE,
    reason="chr19 real-data fixture not reachable (GCS or fixture missing)",
)
class TestCountFromSetsRealData:
    """Exact-count regression for _count_from_sets on real chr19 data.

    Exact per-row matches plus the ``sum(adj) <= n_samples`` / non-negativity
    invariants guard against off-by-one / over-counting — the class of failure
    the historical complement-storage bug produced.
    """

    def test_raw_counts_match(self, _real_data_results):
        rows, _ = _real_data_results
        for r in rows:
            assert list(r._got_raw) == list(r.expected_gt_counts_raw), (
                f"raw counts mismatch for pair {_row_id(r)}: got "
                f"{list(r._got_raw)}, expected {list(r.expected_gt_counts_raw)}"
            )

    def test_adj_counts_match(self, _real_data_results):
        rows, _ = _real_data_results
        for r in rows:
            assert list(r._got_adj) == list(r.expected_gt_counts_adj), (
                f"adj counts mismatch for pair {_row_id(r)}: got "
                f"{list(r._got_adj)}, expected {list(r.expected_gt_counts_adj)}"
            )

    def test_adj_sum_below_n_samples(self, _real_data_results):
        rows, n_samples = _real_data_results
        for r in rows:
            got = list(r._got_adj)
            assert sum(got) <= n_samples, (
                f"adj sum {sum(got)} exceeds n_samples {n_samples} for pair "
                f"{_row_id(r)} (bug-class signal)"
            )
            assert all(c >= 0 for c in got), (
                f"adj cells contain negatives for pair {_row_id(r)}: {got}"
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
    """Build an encoded-variant struct for _count_from_sets_by_pop.

    Sets are stored positive form (membership explicit). n_with_data =
    |all_samples| + |raw_hr_adj_missing| (= |cats 1-6|); ``all_samples``
    holds cats 1,3-6 (disjoint from cat 2 = raw_hr_adj_missing).
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
        raw_hr_adj_missing=_s(raw_hr_adj_missing),
        n_raw_hr_adj_missing=hl.int32(len(raw_hr_adj_missing)),
    )


class TestCountFromSetsByPop:
    """Per-pop counts reuse _count_from_sets over pop-restricted sample sets.

    Every 9-cell entry is a count of samples, so partitioning the cohort into
    disjoint groups that cover all N indices must make the per-pop cells sum
    (element-wise) to the flat full-cohort cells. This is method-independent
    and catches restriction / size bugs.
    """

    N = 10
    # Two variants over a 10-sample cohort (indices 0-9), positive form.
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
            self.V1.raw_hr_adj_missing, self.V1.n_raw_hr_adj_missing,
            self.V2.raw_het if include_raw else self.V2.adj_het,
            self.V2.raw_hv if include_raw else self.V2.adj_hv,
            self.V2.all_samples, self.V2.n_with_data,
            self.V2.raw_hr_adj_missing, self.V2.n_raw_hr_adj_missing,
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


# ===========================================================================
# _no_pbt_count_fields — release \ PBT via count-time subtraction
# ===========================================================================

class TestNoPbtCountFields:
    """No-PBT counts subtract the PBT∩cohort sub-population from the full cohort.

    The correctness property (disputed then confirmed): because the counted
    cohort partitions disjointly into ``PBT∩cohort`` and ``cohort \\ PBT``, every
    one of the 9 cells — AABB included — is an additive per-sample count, so
    ``full − pbt`` computed element-wise EQUALS counting ``cohort \\ PBT``
    directly. These tests pin that equivalence (esp. for AABB, which is derived
    as ``n − |D'_1 ∪ D'_2|`` rather than stored) so a future change that broke
    the additivity would fail here.
    """

    N = 10
    V1 = _enc_variant(
        het=[1, 2], hv=[3], all_samples=[1, 2, 3, 4],
        raw_hr_adj_missing=[5], adj_het=[1], adj_hv=[3],
    )
    V2 = _enc_variant(
        het=[2, 6], hv=[7], all_samples=[2, 6, 7, 8],
        raw_hr_adj_missing=[9], adj_het=[6], adj_hv=[7],
    )
    PBT = {2, 3, 9}
    KEEP = set(range(N)) - PBT  # cohort \ PBT

    def _flat(self, include_raw):
        return _count_from_sets(
            self.V1.raw_het if include_raw else self.V1.adj_het,
            self.V1.raw_hv if include_raw else self.V1.adj_hv,
            self.V1.all_samples, self.V1.n_with_data,
            self.V1.raw_hr_adj_missing, self.V1.n_raw_hr_adj_missing,
            self.V2.raw_het if include_raw else self.V2.adj_het,
            self.V2.raw_hv if include_raw else self.V2.adj_hv,
            self.V2.all_samples, self.V2.n_with_data,
            self.V2.raw_hr_adj_missing, self.V2.n_raw_hr_adj_missing,
            hl.int32(self.N),
            include_raw_hr_adj_missing=include_raw,
        )

    def _direct_restricted(self, keep, include_raw):
        """Count a sub-cohort directly by restricting to its index set."""
        ks = hl.literal(keep, hl.tset(hl.tint32))
        return list(hl.eval(_count_from_sets(
            *_pop_restrict_variant(
                self.V1,
                self.V1.raw_het if include_raw else self.V1.adj_het,
                self.V1.raw_hv if include_raw else self.V1.adj_hv,
                ks,
            ),
            *_pop_restrict_variant(
                self.V2,
                self.V2.raw_het if include_raw else self.V2.adj_het,
                self.V2.raw_hv if include_raw else self.V2.adj_hv,
                ks,
            ),
            hl.int32(len(keep)),
            include_raw_hr_adj_missing=include_raw,
        )))

    def _no_pbt(self):
        pbt_set = hl.literal(self.PBT, hl.tset(hl.tint32))
        fields = _no_pbt_count_fields(
            self.V1, self.V2, pbt_set, len(self.PBT),
            self._flat(True), self._flat(False),
        )
        return hl.eval(hl.struct(**fields))

    def test_subtract_equals_direct_complement_count(self):
        # The crux: full − (PBT∩cohort) == count(cohort \ PBT) directly, all
        # 9 cells including AABB.
        got = self._no_pbt()
        assert list(got.gt_counts_raw_no_pbt) == self._direct_restricted(
            self.KEEP, include_raw=True
        )
        assert list(got.gt_counts_adj_no_pbt) == self._direct_restricted(
            self.KEEP, include_raw=False
        )

    def test_pbt_plus_no_pbt_equals_full(self):
        # Disjoint partition ⇒ per-cohort cells sum to the full-cohort cells.
        got = self._no_pbt()
        pbt_raw = self._direct_restricted(self.PBT, include_raw=True)
        pbt_adj = self._direct_restricted(self.PBT, include_raw=False)
        flat_raw = list(hl.eval(self._flat(True)))
        flat_adj = list(hl.eval(self._flat(False)))
        assert [a + b for a, b in zip(pbt_raw, got.gt_counts_raw_no_pbt)] == flat_raw
        assert [a + b for a, b in zip(pbt_adj, got.gt_counts_adj_no_pbt)] == flat_adj

    def test_no_pbt_cells_nonneg_and_bounded(self):
        got = self._no_pbt()
        for arr in (got.gt_counts_raw_no_pbt, got.gt_counts_adj_no_pbt):
            cells = list(arr)
            assert all(c >= 0 for c in cells), f"negative no-PBT cells: {cells}"
            assert sum(cells) <= len(self.KEEP), (
                f"no-PBT sum {sum(cells)} exceeds cohort\\PBT size {len(self.KEEP)}"
            )

    def test_empty_pbt_is_noop(self):
        # No PBT members ⇒ no_pbt == full cohort.
        fields = _no_pbt_count_fields(
            self.V1, self.V2,
            hl.empty_set(hl.tint32), 0,
            self._flat(True), self._flat(False),
        )
        got = hl.eval(hl.struct(**fields))
        assert list(got.gt_counts_raw_no_pbt) == list(hl.eval(self._flat(True)))
        assert list(got.gt_counts_adj_no_pbt) == list(hl.eval(self._flat(False)))


# ===========================================================================
# _pbt_index_set_for — driver-side PBT∩cohort index set from the samples global
# ===========================================================================

class TestPbtIndexSetFor:
    """Match the encoded ``samples`` global against a PBT-member HT by ``s``."""

    def _encoded(self, sample_ids):
        return hl.utils.range_table(1).annotate_globals(
            samples=hl.array([hl.struct(s=hl.str(s)) for s in sample_ids])
        )

    def test_indices_and_size(self):
        enc = self._encoded([f"S{i}" for i in range(5)])  # S0..S4 at idx 0..4
        members = hl.Table.parallelize(
            [{"s": "S1"}, {"s": "S3"}, {"s": "SX"}],  # SX not in cohort
            schema=hl.tstruct(s=hl.tstr), key="s",
        )
        idx_set, size = _pbt_index_set_for(enc, members)
        assert size == 2
        assert hl.eval(idx_set) == {1, 3}

    def test_no_members_present(self):
        enc = self._encoded(["A", "B", "C"])
        members = hl.Table.parallelize(
            [{"s": "Z"}], schema=hl.tstruct(s=hl.tstr), key="s",
        )
        idx_set, size = _pbt_index_set_for(enc, members)
        assert size == 0
        assert hl.eval(idx_set) == set()

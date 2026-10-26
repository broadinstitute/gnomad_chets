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
    _count_cohort_vds_kwargs,
    _count_from_sets,
    _count_from_sets_by_pop,
    _count_phase_from_sets,
    _create_var_idx_ht,
    _encode_genotype_sets_by_var_idx,
    _drop_pairs_missing_v_idx,
    _empty_counts_ht,
    _pop_restrict_variant,
    _project_count_fields,
    _read_min_an_pct,
    _subtract_pbt_counts,
    count_all_pairs_via_index,
    filter_pairs_by_an_pct,
    merge_subset_counts,
    restrict_dense_mt_for_encode,
    restrict_encoded_to_pops,
    restrict_encoded_to_samples,
    split_encoded_by_pop,
)
from gnomad_chets.v4.resources import (
    GLOBAL_POP,
    get_count_strata,
    get_count_subsets,
    get_variant_pair_resources,
)


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

    def test_restrict_to_samples_reindexes(self):
        # The no-PBT restriction: keep an explicit sample set (here s1, s3) and
        # re-index — same re-index core as the pop path, so the keep-set literal
        # stays out of the per-pair count.
        keep_ht = hl.Table.parallelize(
            [{"s": "s1"}, {"s": "s3"}, {"s": "sX"}],  # sX not in cohort
            schema=hl.tstruct(s=hl.tstr), key="s",
        )
        out, n_keep = restrict_encoded_to_samples(self._encoded_ht(), keep_ht)
        assert n_keep == 2
        r = out.collect()[0]
        assert set(r.all_samples) == {0, 1}   # {1,3} → {0,1}
        assert set(r.adj_het) == {0}          # {1} kept → {0}
        assert set(r.adj_hv) == {1}           # {3} → {1}
        assert r.n_with_data == 2
        assert dict(r.phased_het) == {0: hl.Struct(pid="p", gt0=0)}
        assert hl.eval(out.index_globals().samples) == [
            hl.Struct(s="s1"), hl.Struct(s="s3"),
        ]


# ===========================================================================
# split_encoded_by_pop — one encode, one count cohort per group
# ===========================================================================

class TestSplitEncodedByPop:
    """--split-encode-by-pop: restricting one encode to each group and counting
    each restriction must reproduce the cohort's counts exactly once summed,
    all 9 cells (AABB from each group's own sample count) and the phase counts,
    and the per-group breakdown must be each group's own counts.
    """

    _PHASE_T = hl.tstruct(pid=hl.tstr, gt0=hl.tint32)
    POPS = ["afr", "nfe", "afr", "nfe", "sas", "nfe", "afr", "sas"]
    VARIANTS = [
        ("chr1", 100, ["A", "C"]),
        ("chr1", 200, ["G", "T"]),
        ("chr1", 300, ["C", "A"]),
    ]

    def _encoded_ht(self, n_samples=8):
        # s0 has no entry at v0 (in all_samples, no call) and is raw hom-ref
        # failing adj at v1; s4 is a raw-only het at v0; s5 is raw hom-ref
        # failing adj at v0. Phase: at (v0, v1), s1 is cis, s2 and s6 trans.
        rows = [
            dict(
                v_idx=0, raw_het={1, 2, 4, 6}, raw_hv={3}, adj_het={1, 2, 6},
                adj_hv={3}, all_samples={0, 1, 2, 3, 4, 6},
                raw_hr_adj_missing={5},
                phased_het={1: ("p1", 0), 2: ("p2", 1), 6: ("p6", 0)},
            ),
            dict(
                v_idx=1, raw_het={1, 2, 5, 6}, raw_hv={7}, adj_het={1, 2, 6},
                adj_hv={7}, all_samples={1, 2, 5, 6, 7},
                raw_hr_adj_missing={0},
                phased_het={1: ("p1", 0), 2: ("p2", 0), 6: ("p6", 1)},
            ),
            dict(
                v_idx=2, raw_het={3}, raw_hv=set(), adj_het={3}, adj_hv=set(),
                all_samples={3, 4}, raw_hr_adj_missing=set(), phased_het={},
            ),
        ]
        structs = [
            hl.Struct(
                v_idx=r["v_idx"],
                raw_het=r["raw_het"], raw_hv=r["raw_hv"],
                adj_het=r["adj_het"], adj_hv=r["adj_hv"],
                all_samples=r["all_samples"],
                n_with_data=len(r["all_samples"]) + len(r["raw_hr_adj_missing"]),
                raw_hr_adj_missing=r["raw_hr_adj_missing"],
                n_raw_hr_adj_missing=len(r["raw_hr_adj_missing"]),
                phased_het={
                    k: hl.Struct(pid=pid, gt0=gt0)
                    for k, (pid, gt0) in r["phased_het"].items()
                },
            )
            for r in rows
        ]
        ht = hl.Table.parallelize(
            structs,
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
                [hl.Struct(s=f"s{i}") for i in range(n_samples)],
                hl.tarray(hl.tstruct(s=hl.tstr)),
            )
        )

    def _pop_ht(self, pops):
        return hl.Table.parallelize(
            [hl.Struct(s=f"s{i}", pop=p) for i, p in enumerate(pops)],
            hl.tstruct(s=hl.tstr, pop=hl.tstr),
            key=["s"],
        )

    def _vp_ht(self):
        v = self.VARIANTS
        return _vp_table([
            dict(c1=v[i][0], p1=v[i][1], a1=v[i][2],
                 c2=v[j][0], p2=v[j][1], a2=v[j][2])
            for i, j in ((0, 1), (0, 2), (1, 2))
        ])

    @staticmethod
    def _by_pos(ht):
        return {
            (r.locus1.position, r.locus2.position): r for r in ht.collect()
        }

    def _count(self, encoded_ht):
        return count_all_pairs_via_index(
            self._vp_ht(), _var_idx_table(self.VARIANTS), encoded_ht,
            emit_phase=True,
        )

    def test_partitions_cohort_and_reindexes(self):
        parts = split_encoded_by_pop(self._encoded_ht(), self._pop_ht(self.POPS))
        assert list(parts) == ["afr", "nfe", "sas"]
        members = {
            p: [smp.s for smp in hl.eval(t.index_globals().samples)]
            for p, t in parts.items()
        }
        assert members == {
            "afr": ["s0", "s2", "s6"],
            "nfe": ["s1", "s3", "s5"],
            "sas": ["s4", "s7"],
        }
        # afr keeps s0, s2, s6 → dense 0, 1, 2.
        v0 = {r.v_idx: r for r in parts["afr"].collect()}[0]
        assert set(v0.all_samples) == {0, 1, 2}
        assert set(v0.adj_het) == {1, 2}
        assert dict(v0.phased_het) == {
            1: hl.Struct(pid="p2", gt0=1), 2: hl.Struct(pid="p6", gt0=0),
        }

    def test_group_counts_sum_to_cohort_counts(self):
        encoded = self._encoded_ht()
        parts = split_encoded_by_pop(encoded, self._pop_ht(self.POPS))
        per_pop = {p: self._count(t) for p, t in parts.items()}
        merged = self._by_pos(
            merge_subset_counts(list(per_pop.values()), pops=list(per_pop))
        )
        full = self._by_pos(self._count(encoded))
        assert set(merged) == set(full) == {(100, 200), (100, 300), (200, 300)}
        for key, f in full.items():
            m = merged[key]
            assert list(m.gt_counts_raw) == list(f.gt_counts_raw)
            assert list(m.gt_counts_adj) == list(f.gt_counts_adj)
            assert (m.n_phased_cis, m.n_phased_trans) == (
                f.n_phased_cis, f.n_phased_trans,
            )
            assert m.gt_counts_by_pop[GLOBAL_POP] == hl.Struct(
                raw=m.gt_counts_raw, adj=m.gt_counts_adj,
                n_phased_cis=m.n_phased_cis, n_phased_trans=m.n_phased_trans,
            )
        # The fixture exercises the phase counts and every group's AABB.
        assert (full[(100, 200)].n_phased_cis,
                full[(100, 200)].n_phased_trans) == (1, 2)
        for p, t in per_pop.items():
            n_pop = self.POPS.count(p)
            for key, r in self._by_pos(t).items():
                by_pop = merged[key].gt_counts_by_pop[p]
                assert list(by_pop.raw) == list(r.gt_counts_raw)
                assert list(by_pop.adj) == list(r.gt_counts_adj)
                assert (by_pop.n_phased_cis, by_pop.n_phased_trans) == (
                    r.n_phased_cis, r.n_phased_trans,
                )
                assert sum(r.gt_counts_adj) <= n_pop

    def test_no_pbt_by_group_equals_direct_count(self):
        # --count-pbt-members on each group + --merge-no-pbt: every group's
        # release-minus-trio counts, and the total's, must equal counting the
        # cohort without its trio members directly, phase included. s1 (cis)
        # and s6 (trans) are trio members; sas holds none, so its PBT counts
        # are empty and its no-PBT counts equal its counts.
        encoded = self._encoded_ht()
        pops = self._pop_ht(self.POPS)
        pbt = {"s1", "s6"}

        def _samples(names):
            return hl.Table.parallelize(
                [hl.Struct(s=x) for x in sorted(names)], hl.tstruct(s=hl.tstr),
                key="s",
            )

        parts = split_encoded_by_pop(encoded, pops)
        counts, pbt_counts = [], []
        for t in parts.values():
            counts.append(self._count(t))
            members = {x.s for x in hl.eval(t.index_globals().samples)} & pbt
            pbt_counts.append(
                self._count(restrict_encoded_to_samples(t, _samples(members))[0])
                if members
                else _empty_counts_ht(emit_phase=True)
            )
        labels = list(parts)
        got = self._by_pos(_subtract_pbt_counts(
            merge_subset_counts(counts, pops=labels),
            merge_subset_counts(pbt_counts, pops=labels),
        ))

        others = {f"s{i}" for i in range(len(self.POPS))} - pbt
        direct_encoded, _ = restrict_encoded_to_samples(encoded, _samples(others))
        direct_all = self._by_pos(self._count(direct_encoded))
        direct = {
            p: self._by_pos(self._count(t))
            for p, t in split_encoded_by_pop(direct_encoded, pops).items()
        }
        assert set(got) == set(direct_all)
        for key, r in got.items():
            d = direct_all[key]
            assert list(r.gt_counts_raw_no_pbt) == list(d.gt_counts_raw)
            assert list(r.gt_counts_adj_no_pbt) == list(d.gt_counts_adj)
            assert (r.n_phased_cis_no_pbt, r.n_phased_trans_no_pbt) == (
                d.n_phased_cis, d.n_phased_trans,
            )
            assert r.gt_counts_by_pop_no_pbt[GLOBAL_POP] == hl.Struct(
                raw=r.gt_counts_raw_no_pbt, adj=r.gt_counts_adj_no_pbt,
                n_phased_cis=r.n_phased_cis_no_pbt,
                n_phased_trans=r.n_phased_trans_no_pbt,
            )
            for p, by_key in direct.items():
                g, dp = r.gt_counts_by_pop_no_pbt[p], by_key[key]
                assert list(g.raw) == list(dp.gt_counts_raw)
                assert list(g.adj) == list(dp.gt_counts_adj)
                assert (g.n_phased_cis, g.n_phased_trans) == (
                    dp.n_phased_cis, dp.n_phased_trans,
                )
            assert r.gt_counts_by_pop_no_pbt["sas"] == r.gt_counts_by_pop["sas"]
        # The trio members carried the (100, 200) pair's cis and one trans call.
        assert (got[(100, 200)].n_phased_cis_no_pbt,
                got[(100, 200)].n_phased_trans_no_pbt) == (0, 1)

    def test_unlabeled_sample_rejected(self):
        with pytest.raises(ValueError, match="no genetic-ancestry group"):
            split_encoded_by_pop(
                self._encoded_ht(), self._pop_ht(self.POPS[:-1])
            )

    def test_strata_partition_by_group(self):
        # Every group appears once in UKB and once outside it, and the UKB
        # strata are exactly the UKB subset VDSes.
        strata = get_count_strata("exomes")
        ukb = {s for s in strata if s.startswith("ukb.")}
        assert ukb == set(get_count_subsets("exomes")) - {"non_ukb"}
        groups = sorted(set(strata.values()))
        assert sorted(strata.values()) == sorted(groups * 2)
        assert all(s.endswith(f".{p}") for s, p in strata.items())


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

    def test_genomes_style_all_false_is_corrected(self):
        # Genomes carry no fixed_homalt_model; the encode input sets it False for
        # every sample, reproducing hom_alt_depletion_fix (no per-sample gate).
        # With every sample False, an eligible high-AB het must be corrected.
        r = self._encode_one_variant(0.02, [
            {"gt": "0/1", "ad": [1, 19], "het_nr": False, "fixed": False},
            {"gt": "0/1", "ad": [10, 10], "het_nr": False, "fixed": False},
        ])
        assert set(r.adj_hv) == {0}
        assert set(r.adj_het) == {1}

    def test_absent_fixed_homalt_model_skips_correction(self):
        # The trap the genomes fix avoids: the encoder feature-detects
        # fixed_homalt_model and SKIPS the whole correction when the column is
        # absent -- so the eligible high-AB het below stays a het. Leaving the
        # field off genomes, rather than setting it False, would run cleanly and
        # silently leave genomes counts hom-alt depleted.
        mt = hl.utils.range_matrix_table(1, 1)
        mt = mt.annotate_rows(
            locus=hl.locus("chr1", 100, "GRCh38"), alleles=["A", "T"],
            af=hl.float64(0.02),
        ).key_rows_by("locus", "alleles").drop("row_idx")
        mt = mt.annotate_cols(s=hl.str(mt.col_idx))
        mt = mt.annotate_entries(
            GT=hl.parse_call("0/1"), GQ=hl.int32(50), DP=hl.int32(20),
            AD=hl.literal([1, 19]), _het_non_ref=False,
        ).key_cols_by("s")
        assert "fixed_homalt_model" not in mt.col
        r = _encode_genotype_sets_by_var_idx(mt, _create_var_idx_ht(mt)).collect()[0]
        assert set(r.adj_het) == {0}
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


class TestCountCohort:
    """Selecting a count cohort for either data type, encoding from an existing
    dense MT, and qualifying outputs per cohort (--count-cohort)."""

    def test_exomes_loader_kwargs(self):
        assert _count_cohort_vds_kwargs("exomes", release_only=False) == {
            "release_only": False, "high_quality_only": True,
        }
        assert _count_cohort_vds_kwargs("exomes", release_only=True) == {
            "release_only": True, "high_quality_only": False,
        }

    def test_genomes_release_kwargs(self):
        # The genomes loader rejects high_quality_only, and its
        # remove_hard_filtered_samples default reads a v3.1 HT that no longer
        # exists -- so neither may reach it. (The high-quality branch reads the
        # meta HT and is exercised on the cluster, not here.)
        kw = _count_cohort_vds_kwargs("genomes", release_only=True)
        assert "high_quality_only" not in kw
        assert kw["remove_hard_filtered_samples"] is False
        assert kw["release_only"] is True

    def test_dense_mt_missing_fields_rejected(self):
        # A GT-only MT (what a phaser VCF export holds) cannot feed the encoder:
        # adj needs GQ/DP/AD and the high-AB correction needs _het_non_ref.
        # Rejected before any metadata is read.
        mt = hl.utils.range_matrix_table(1, 1)
        mt = mt.annotate_rows(
            locus=hl.locus("chr1", 100, "GRCh38"), alleles=["A", "T"],
        ).key_rows_by("locus", "alleles").drop("row_idx")
        mt = mt.annotate_cols(s=hl.str(mt.col_idx)).key_cols_by("s")
        mt = mt.annotate_entries(GT=hl.parse_call("0/1"))
        with pytest.raises(ValueError, match="missing entry fields"):
            restrict_dense_mt_for_encode(
                mt, mt.rows().select(), "genomes",
                release_only=False, filter_intervals=None,
            )

    def test_count_cohort_qualifies_outputs_not_pair_list(self):
        # Two cohorts counted against one pair list must not collide: the count
        # outputs carry the cohort, the shared pair-list input does not.
        paths = {}
        for cohort in ("s1", "s2"):
            r = get_variant_pair_resources(
                data_type="genomes", test=True, tmp_dir="gs://tmp",
                output_postfix="chr21_vp_upper", subset=cohort,
            )
            step = r.create_variant_pair_genotype_counts_ht
            paths[cohort] = (step.vp_list_ht.path, step.vp_gt_counts_ht.path)
        assert paths["s1"][0] == paths["s2"][0]
        assert ".s1" not in paths["s1"][0]
        assert paths["s1"][1] != paths["s2"][1]
        assert "chr21_vp_upper.s1" in paths["s1"][1]
        assert "chr21_vp_upper.s2" in paths["s2"][1]


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
# No-PBT: release \ PBT = full − (PBT∩cohort restriction), subtracted at combine
# ===========================================================================

class TestPbtSubtractionExactness:
    """The property the no-PBT feature relies on: because the counted cohort
    partitions disjointly into ``PBT∩cohort`` and ``cohort \\ PBT``, every one of
    the 9 cells — AABB included — is an additive per-sample count, so ``full −
    pbt`` element-wise EQUALS counting ``cohort \\ PBT`` directly. In the pipeline
    the PBT counts come from :func:`restrict_encoded_to_samples` + a normal
    count; here we exercise the same math via :func:`_count_from_sets` +
    :func:`_pop_restrict_variant` (the restriction primitive), so a change that
    broke the additivity (esp. for AABB) fails here.
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
        """full − (PBT∩cohort), the pipeline's `_subtract_pbt_counts` at the
        array level."""
        out = {}
        for key, inc in (("raw", True), ("adj", False)):
            flat = list(hl.eval(self._flat(inc)))
            pbt = self._direct_restricted(self.PBT, inc)
            out[key] = [a - b for a, b in zip(flat, pbt)]
        return out

    def test_subtract_equals_direct_complement_count(self):
        # The crux: full − (PBT∩cohort) == count(cohort \ PBT) directly, all
        # 9 cells including AABB.
        got = self._no_pbt()
        assert got["raw"] == self._direct_restricted(self.KEEP, include_raw=True)
        assert got["adj"] == self._direct_restricted(self.KEEP, include_raw=False)

    def test_pbt_plus_no_pbt_equals_full(self):
        # Disjoint partition ⇒ per-cohort cells sum to the full-cohort cells.
        got = self._no_pbt()
        pbt_raw = self._direct_restricted(self.PBT, include_raw=True)
        pbt_adj = self._direct_restricted(self.PBT, include_raw=False)
        assert [a + b for a, b in zip(pbt_raw, got["raw"])] == list(
            hl.eval(self._flat(True))
        )
        assert [a + b for a, b in zip(pbt_adj, got["adj"])] == list(
            hl.eval(self._flat(False))
        )

    def test_no_pbt_cells_nonneg_and_bounded(self):
        got = self._no_pbt()
        for cells in (got["raw"], got["adj"]):
            assert all(c >= 0 for c in cells), f"negative no-PBT cells: {cells}"
            assert sum(cells) <= len(self.KEEP), (
                f"no-PBT sum {sum(cells)} exceeds cohort\\PBT size {len(self.KEEP)}"
            )

    def test_empty_pbt_is_noop(self):
        # No PBT members ⇒ PBT counts all zero ⇒ no_pbt == full cohort.
        assert self._direct_restricted(set(), include_raw=True) == [0] * 9
        got = self._no_pbt()
        # (recompute against full to make the no-op explicit)
        assert got["raw"] != [0] * 9  # sanity: the fixture is non-trivial


class TestSubtractPbtCounts:
    """`_subtract_pbt_counts` annotates gt_counts_*_no_pbt = full − pbt,
    element-wise per pair, joining full and PBT count HTs on the pair key."""

    @staticmethod
    def _counts_ht(rows):
        return hl.Table.parallelize(
            rows,
            hl.tstruct(
                locus1=hl.tlocus("GRCh38"), alleles1=hl.tarray(hl.tstr),
                locus2=hl.tlocus("GRCh38"), alleles2=hl.tarray(hl.tstr),
                gt_counts_raw=hl.tarray(hl.tint32),
                gt_counts_adj=hl.tarray(hl.tint32),
            ),
            key=["locus1", "alleles1", "locus2", "alleles2"],
        )

    def test_elementwise_subtract(self):
        key = dict(
            locus1=hl.locus("chr1", 100, "GRCh38"), alleles1=["A", "T"],
            locus2=hl.locus("chr1", 200, "GRCh38"), alleles2=["C", "G"],
        )
        full = self._counts_ht([{
            **key,
            "gt_counts_raw": [100, 5, 4, 3, 2, 1, 0, 0, 0],
            "gt_counts_adj": [90, 4, 3, 2, 1, 0, 0, 0, 0],
        }])
        pbt = self._counts_ht([{
            **key,
            "gt_counts_raw": [10, 1, 0, 1, 0, 0, 0, 0, 0],
            "gt_counts_adj": [9, 1, 0, 0, 0, 0, 0, 0, 0],
        }])
        out = _subtract_pbt_counts(full, pbt).collect()[0]
        assert list(out.gt_counts_raw_no_pbt) == [90, 4, 4, 2, 2, 1, 0, 0, 0]
        assert list(out.gt_counts_adj_no_pbt) == [81, 3, 3, 2, 1, 0, 0, 0, 0]


class TestMergeSubsetCounts:
    """`merge_subset_counts` sums per-subset (--vds-subset) genotype-count HTs
    element-wise into the full-cohort counts, over the shared pair key."""

    @staticmethod
    def _counts_ht(rows, *, phase=True):
        fields = dict(
            locus1=hl.tlocus("GRCh38"), alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh38"), alleles2=hl.tarray(hl.tstr),
            gt_counts_raw=hl.tarray(hl.tint32),
            gt_counts_adj=hl.tarray(hl.tint32),
        )
        if phase:
            fields["n_phased_cis"] = hl.tint32
            fields["n_phased_trans"] = hl.tint32
        return hl.Table.parallelize(
            rows, hl.tstruct(**fields),
            key=["locus1", "alleles1", "locus2", "alleles2"],
        )

    @staticmethod
    def _row(pos2, raw, adj, cis, trans):
        return dict(
            locus1=hl.locus("chr1", 100, "GRCh38"), alleles1=["A", "C"],
            locus2=hl.locus("chr1", pos2, "GRCh38"), alleles2=["G", "T"],
            gt_counts_raw=raw, gt_counts_adj=adj,
            n_phased_cis=cis, n_phased_trans=trans,
        )

    def test_elementwise_sum_and_dtype(self):
        # Pair at pos2=200 present in both subsets; pos2=300 only in the second
        # (full-outer union keeps it, summing just the present subset).
        a = self._counts_ht([
            self._row(200, [10, 1, 0, 1, 2, 0, 0, 0, 1],
                      [9, 1, 0, 1, 2, 0, 0, 0, 1], 3, 1),
        ])
        b = self._counts_ht([
            self._row(200, [20, 0, 0, 2, 1, 0, 0, 0, 0],
                      [20, 0, 0, 2, 1, 0, 0, 0, 0], 1, 0),
            self._row(300, [5, 0, 0, 0, 0, 0, 0, 0, 1],
                      [5, 0, 0, 0, 0, 0, 0, 0, 1], 0, 0),
        ])
        merged = merge_subset_counts([a, b])
        # Sums stay int32 (matches the per-subset schema).
        assert merged.gt_counts_raw.dtype == hl.tarray(hl.tint32)
        by_pos = {r.locus2.position: r for r in merged.collect()}
        both = by_pos[200]
        assert list(both.gt_counts_raw) == [30, 1, 0, 3, 3, 0, 0, 0, 1]
        assert list(both.gt_counts_adj) == [29, 1, 0, 3, 3, 0, 0, 0, 1]
        assert both.n_phased_cis == 4 and both.n_phased_trans == 1
        only_b = by_pos[300]
        assert list(only_b.gt_counts_raw) == [5, 0, 0, 0, 0, 0, 0, 0, 1]
        assert only_b.n_phased_cis == 0 and only_b.n_phased_trans == 0

    def test_no_phase_columns_dropped(self):
        # When any input lacks phase columns, the merge omits them.
        a = self._counts_ht([
            self._row(200, [1, 0, 0, 0, 0, 0, 0, 0, 0],
                      [1, 0, 0, 0, 0, 0, 0, 0, 0], 0, 0),
        ], phase=True)
        b = self._counts_ht([
            {k: v for k, v in self._row(
                200, [2, 0, 0, 0, 0, 0, 0, 0, 0],
                [2, 0, 0, 0, 0, 0, 0, 0, 0], 0, 0
            ).items() if k not in ("n_phased_cis", "n_phased_trans")},
        ], phase=False)
        merged = merge_subset_counts([a, b])
        assert "n_phased_cis" not in merged.row
        assert list(merged.collect()[0].gt_counts_raw) == [3, 0, 0, 0, 0, 0, 0, 0, 0]

    def test_by_pop_breakdown(self):
        # Two afr strata and one nfe stratum: afr sums its two, "all" is the
        # flat total, and every value stays int32.
        a = self._counts_ht([
            self._row(200, [10, 1, 0, 1, 2, 0, 0, 0, 1],
                      [9, 1, 0, 1, 2, 0, 0, 0, 1], 1, 1),
        ])
        b = self._counts_ht([
            self._row(200, [5, 0, 0, 0, 1, 0, 0, 0, 0],
                      [5, 0, 0, 0, 1, 0, 0, 0, 0], 1, 0),
        ])
        c = self._counts_ht([
            self._row(200, [20, 0, 0, 2, 1, 0, 0, 0, 0],
                      [20, 0, 0, 2, 1, 0, 0, 0, 0], 0, 1),
        ])
        merged = merge_subset_counts([a, b, c], pops=["afr", "nfe", "afr"])
        assert merged.gt_counts_by_pop.dtype == hl.tdict(
            hl.tstr,
            hl.tstruct(
                raw=hl.tarray(hl.tint32), adj=hl.tarray(hl.tint32),
                n_phased_cis=hl.tint32, n_phased_trans=hl.tint32,
            ),
        )
        r = merged.collect()[0]
        assert set(r.gt_counts_by_pop) == {GLOBAL_POP, "afr", "nfe"}
        afr = r.gt_counts_by_pop["afr"]
        assert list(afr.raw) == [30, 1, 0, 3, 3, 0, 0, 0, 1]
        assert (afr.n_phased_cis, afr.n_phased_trans) == (1, 2)
        assert list(r.gt_counts_by_pop["nfe"].adj) == [5, 0, 0, 0, 1, 0, 0, 0, 0]
        total = r.gt_counts_by_pop[GLOBAL_POP]
        assert list(total.raw) == list(r.gt_counts_raw) == [
            35, 1, 0, 3, 4, 0, 0, 0, 1,
        ]
        assert (total.n_phased_cis, total.n_phased_trans) == (2, 2)

    def test_by_pop_label_count_must_match(self):
        a = self._counts_ht([
            self._row(200, [1, 0, 0, 0, 0, 0, 0, 0, 0],
                      [1, 0, 0, 0, 0, 0, 0, 0, 0], 0, 0),
        ])
        with pytest.raises(ValueError, match="group labels"):
            merge_subset_counts([a, a], pops=["afr"])

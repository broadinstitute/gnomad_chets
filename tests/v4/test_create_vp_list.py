"""Unit tests for ``gnomad_chets.v4.create_vp_list.assemble_sites_ht``.

Run with::

    pytest v4/test_create_vp_list.py -v

Covers the ``--preprocess-sites-ht`` assemble step: the multi-source join,
the adj/raw callstats, the optional-source schema, and — most importantly —
the QC gate that drops a variant only when *both* filter sources agree it
hard-fails (dropping AS_VQSR and both-agree InbreedingCoeff, retaining AC0
and InbreedingCoeff disagreements).

The per-variant helpers this step also calls (``get_an_percent_expr``,
``clinvar_category_match_expr``, ``filter_for_testing``, ``_find_strata_index``)
are covered in ``test_utils.py``.

Tables are built with ``hl.Table.parallelize`` and run against the local
Spark backend; the Hail session comes from ``v4/conftest.py``.
"""

import hail as hl
import pytest

from gnomad_chets.v4.create_vp_list import assemble_sites_ht

REF = "GRCh38"
CHR = "chr1"


# ---------------------------------------------------------------------------
# Synthetic-input builders
# ---------------------------------------------------------------------------
def _filters_ht(rows):
    """(pos, [filter tokens]) rows -> HT keyed (locus, alleles) with ``filters``."""
    typed = [
        hl.Struct(
            locus=hl.locus(CHR, pos, REF),
            alleles=["A", "T"],
            filters=hl.set(hl.literal(list(toks), hl.tarray(hl.tstr))),
        )
        for pos, toks in rows
    ]
    return hl.Table.parallelize(
        typed,
        hl.tstruct(
            locus=hl.tlocus(REF),
            alleles=hl.tarray(hl.tstr),
            filters=hl.tset(hl.tstr),
        ),
        key=["locus", "alleles"],
    )


def _freq_ht(specs):
    """(pos, (AC,AF,AN)_adj, (AC,AF,AN)_raw) rows -> freq HT (freq[0]=adj, [1]=raw)."""
    rows = [
        hl.Struct(
            locus=hl.locus(CHR, pos, REF),
            alleles=["A", "T"],
            freq=[
                hl.Struct(AC=adj[0], AF=adj[1], AN=adj[2]),
                hl.Struct(AC=raw[0], AF=raw[1], AN=raw[2]),
            ],
        )
        for pos, adj, raw in specs
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus(REF),
            alleles=hl.tarray(hl.tstr),
            freq=hl.tarray(
                hl.tstruct(AC=hl.tint32, AF=hl.tfloat64, AN=hl.tint32)
            ),
        ),
        key=["locus", "alleles"],
    )


def _freq_ht_dummy(positions):
    """Freq HT with placeholder adj/raw callstats for each position."""
    return _freq_ht([(p, (1, 0.1, 10), (2, 0.2, 10)) for p in positions])


def _vep_ht(positions):
    """VEP HT carrying a minimal ``vep`` struct."""
    rows = [
        hl.Struct(
            locus=hl.locus(CHR, p, REF),
            alleles=["A", "T"],
            vep=hl.Struct(most_severe_consequence="missense_variant"),
        )
        for p in positions
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus(REF),
            alleles=hl.tarray(hl.tstr),
            vep=hl.tstruct(most_severe_consequence=hl.tstr),
        ),
        key=["locus", "alleles"],
    )


def _spliceai_ht(specs):
    """(pos, ds) rows -> SpliceAI HT with ``spliceai_ds_max``."""
    rows = [
        hl.Struct(locus=hl.locus(CHR, p, REF), alleles=["A", "T"], spliceai_ds_max=ds)
        for p, ds in specs
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus(REF),
            alleles=hl.tarray(hl.tstr),
            spliceai_ds_max=hl.tfloat32,
        ),
        key=["locus", "alleles"],
    )


def _pangolin_ht(specs):
    """(pos, ds) rows -> Pangolin HT with ``pangolin_largest_ds``."""
    rows = [
        hl.Struct(
            locus=hl.locus(CHR, p, REF), alleles=["A", "T"], pangolin_largest_ds=ds
        )
        for p, ds in specs
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus(REF),
            alleles=hl.tarray(hl.tstr),
            pangolin_largest_ds=hl.tfloat64,
        ),
        key=["locus", "alleles"],
    )


def _clinvar_ht(specs):
    """(pos, clnsig, geneinfo) rows -> ClinVar-shaped HT with ``info``."""
    rows = [
        hl.Struct(
            locus=hl.locus(CHR, p, REF),
            alleles=["A", "T"],
            info=hl.Struct(
                CLNSIG=clnsig,
                CLNREVSTAT=["reviewed_by_expert_panel"],
                CLNSIGCONF=hl.missing(hl.tstr),
                GENEINFO=geneinfo,
            ),
        )
        for p, clnsig, geneinfo in specs
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus(REF),
            alleles=hl.tarray(hl.tstr),
            info=hl.tstruct(
                CLNSIG=hl.tarray(hl.tstr),
                CLNREVSTAT=hl.tarray(hl.tstr),
                CLNSIGCONF=hl.tstr,
                GENEINFO=hl.tstr,
            ),
        ),
        key=["locus", "alleles"],
    )


def _an_ht(positions):
    """All-sites AN HT (locus-keyed) with adj=50% globals for ``get_an_percent_expr``."""
    meta = [
        {"group": "raw"},
        {"group": "adj"},
        {"group": "adj", "sex": "XX"},
        {"group": "adj", "sex": "XY"},
    ]
    rows = [
        hl.Struct(locus=hl.locus(CHR, p, REF), AN=[200, 100, 50, 50]) for p in positions
    ]
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(locus=hl.tlocus(REF), AN=hl.tarray(hl.tint32)),
        key=["locus"],
    )
    return ht.annotate_globals(
        strata_sample_count=[100, 100, 50, 50], strata_meta=meta
    )


def _assemble(variants, with_release=True, **kwargs):
    """Assemble from a list of ``dict(pos=, of=[...], rel=[...])`` and return the HT.

    ``of`` is the only_filters token list; include ``rel`` to place the variant
    in ``release_filter_ht`` with that token list (omit to leave it absent).
    """
    positions = [v["pos"] for v in variants]
    filter_ht = _filters_ht([(v["pos"], v["of"]) for v in variants])
    release_ht = (
        _filters_ht([(v["pos"], v["rel"]) for v in variants if "rel" in v])
        if with_release
        else None
    )
    return assemble_sites_ht(
        filter_ht=filter_ht,
        freq_ht=_freq_ht_dummy(positions),
        vep_ht=_vep_ht(positions),
        release_filter_ht=release_ht,
        **kwargs,
    )


def _surviving(variants, **kwargs):
    """Set of positions that survive the assemble gate."""
    return set(_assemble(variants, **kwargs).locus.position.collect())


# ===========================================================================
# Gate: which variants survive
# ===========================================================================
class TestAssembleSitesHtGate:
    def test_pass_and_ac0_kept_hard_filter_dropped(self):
        survivors = _surviving(
            [
                {"pos": 100, "of": [], "rel": []},  # PASS
                {"pos": 101, "of": ["AC0"], "rel": ["AC0"]},  # AC0-only
                {"pos": 102, "of": ["AS_VQSR"], "rel": ["AS_VQSR"]},  # hard fail
                {"pos": 103, "of": ["AC0", "AS_VQSR"], "rel": ["AC0", "AS_VQSR"]},
            ]
        )
        assert survivors == {100, 101}

    def test_both_agree_inbreeding_coeff_dropped(self):
        survivors = _surviving(
            [{"pos": 200, "of": ["InbreedingCoeff"], "rel": ["InbreedingCoeff"]}]
        )
        assert survivors == set()

    def test_both_agree_ic_with_ac0_dropped(self):
        survivors = _surviving(
            [
                {
                    "pos": 201,
                    "of": ["AC0", "InbreedingCoeff"],
                    "rel": ["AC0", "InbreedingCoeff"],
                }
            ]
        )
        assert survivors == set()

    @pytest.mark.parametrize(
        "of, rel",
        [
            (["InbreedingCoeff"], []),  # only_filters flags IC, release doesn't
            ([], ["InbreedingCoeff"]),  # release flags IC, only_filters doesn't
            (["AC0", "InbreedingCoeff"], ["AC0"]),  # disagreement alongside AC0
            (["AC0"], ["AC0", "InbreedingCoeff"]),
        ],
    )
    def test_inbreeding_coeff_disagreement_kept(self, of, rel):
        survivors = _surviving([{"pos": 300, "of": of, "rel": rel}])
        assert survivors == {300}

    def test_as_vqsr_dropped_regardless_of_ic_disagreement(self):
        # AS_VQSR present in both; IC differs. AS_VQSR still drops it.
        survivors = _surviving(
            [{"pos": 301, "of": ["AS_VQSR"], "rel": ["AS_VQSR", "InbreedingCoeff"]}]
        )
        assert survivors == set()

    def test_release_missing_falls_back_to_only_filters(self):
        # No ``rel`` key => variant absent from release_filter_ht.
        survivors = _surviving(
            [
                {"pos": 400, "of": []},  # PASS, unreleased -> kept
                {"pos": 401, "of": ["AC0"]},  # AC0-only, unreleased -> kept
                {"pos": 402, "of": ["InbreedingCoeff"]},  # IC, no release -> dropped
                {"pos": 403, "of": ["AS_VQSR"]},  # hard fail -> dropped
            ]
        )
        assert survivors == {400, 401}

    def test_no_release_ht_gates_on_only_filters(self):
        # release_filter_ht=None => InbreedingCoeff is a drop (single source).
        survivors = _surviving(
            [
                {"pos": 500, "of": []},
                {"pos": 501, "of": ["AC0"]},
                {"pos": 502, "of": ["InbreedingCoeff"]},
                {"pos": 503, "of": ["AS_VQSR"]},
            ],
            with_release=False,
        )
        assert survivors == {500, 501}


# ===========================================================================
# Annotations / schema
# ===========================================================================
class TestAssembleSitesHtAnnotations:
    def test_adj_and_raw_callstats(self):
        filter_ht = _filters_ht([(100, [])])
        freq_ht = _freq_ht([(100, (3, 0.03, 1000), (7, 0.07, 1200))])
        ht = assemble_sites_ht(
            filter_ht=filter_ht, freq_ht=freq_ht, vep_ht=_vep_ht([100])
        )
        row = ht.collect()[0]
        assert (row.ac, row.af, row.an) == (3, 0.03, 1000)
        assert (row.ac_raw, row.af_raw, row.an_raw) == (7, 0.07, 1200)

    def test_missing_from_freq_yields_null_callstats(self):
        filter_ht = _filters_ht([(100, [])])
        empty_freq = _freq_ht([])  # variant absent from freq
        ht = assemble_sites_ht(
            filter_ht=filter_ht, freq_ht=empty_freq, vep_ht=_vep_ht([100])
        )
        row = ht.collect()[0]
        assert row.ac is None and row.af is None and row.an is None
        assert row.ac_raw is None and row.an_raw is None

    def test_vep_carried_through(self):
        ht = _assemble([{"pos": 100, "of": [], "rel": []}])
        assert ht.collect()[0].vep.most_severe_consequence == "missense_variant"

    def test_filters_release_present_when_provided(self):
        ht = _assemble([{"pos": 100, "of": ["AC0"], "rel": ["AC0", "InbreedingCoeff"]}])
        assert "filters_release" in ht.row
        row = ht.collect()[0]
        assert row.filters == {"AC0"}
        assert row.filters_release == {"AC0", "InbreedingCoeff"}

    def test_filters_release_absent_when_no_release_ht(self):
        ht = _assemble([{"pos": 100, "of": []}], with_release=False)
        assert "filters_release" not in ht.row

    def test_filters_release_null_for_unreleased_variant(self):
        # Present in only_filters (PASS) but absent from release_filter_ht.
        ht = _assemble([{"pos": 100, "of": []}])  # no rel key
        assert "filters_release" in ht.row
        assert ht.collect()[0].filters_release is None

    def test_optional_sources_omitted_when_none(self):
        ht = _assemble([{"pos": 100, "of": [], "rel": []}])
        for f in ("an_pct", "spliceai_ds_max", "pangolin_largest_ds", "clinvar"):
            assert f not in ht.row

    def test_all_optional_sources_present(self):
        filter_ht = _filters_ht([(100, [])])
        ht = assemble_sites_ht(
            filter_ht=filter_ht,
            freq_ht=_freq_ht_dummy([100]),
            vep_ht=_vep_ht([100]),
            release_filter_ht=_filters_ht([(100, [])]),
            an_ht=_an_ht([100]),
            spliceai_ht=_spliceai_ht([(100, 0.9)]),
            pangolin_ht=_pangolin_ht([(100, 0.8)]),
            clinvar_ht=_clinvar_ht([(100, ["Pathogenic"], "SGCA:6442")]),
        )
        row = ht.collect()[0]
        assert row.an_pct == 50  # AN[adj]=100 / (adj_count*2=200) -> 50%
        assert abs(row.spliceai_ds_max - 0.9) < 1e-6
        assert abs(row.pangolin_largest_ds - 0.8) < 1e-9
        assert row.clinvar.is_plp is True
        assert row.clinvar.is_blb is False
        assert row.clinvar.GENEINFO == "SGCA:6442"

    def test_clinvar_missing_variant_is_null(self):
        filter_ht = _filters_ht([(100, []), (101, [])])
        ht = assemble_sites_ht(
            filter_ht=filter_ht,
            freq_ht=_freq_ht_dummy([100, 101]),
            vep_ht=_vep_ht([100, 101]),
            clinvar_ht=_clinvar_ht([(100, ["Pathogenic"], "G:1")]),  # 101 absent
        )
        by_pos = {r.locus.position: r for r in ht.collect()}
        assert by_pos[100].clinvar is not None
        assert by_pos[101].clinvar is None

    def test_row_key_is_locus_alleles(self):
        ht = _assemble([{"pos": 100, "of": [], "rel": []}])
        assert list(ht.key) == ["locus", "alleles"]

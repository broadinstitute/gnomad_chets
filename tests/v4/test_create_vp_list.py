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

from gnomad_chets.v4.create_vp_list import (
    _build_intronic_padding_interval_ht,
    _get_clinvar_gene_id_expr,
    _get_hc_lof_gene_id_expr,
    _get_intronic_padding_gene_id_expr,
    _get_vep_gene_id_expr,
    _validate_clinvar_categories,
    _validate_least_consequence,
    _validate_pathogenic_splice,
    assemble_sites_ht,
    create_variant_filter_ht,
    create_variant_pair_ht,
)
from gnomad_chets.v4.resources import (
    CLINVAR_CATEGORY_FIELD_FMT,
    CLINVAR_CATEGORY_PLP,
    CLINVAR_CATEGORIES,
    SITES_FIELD_CLINVAR,
    SITES_FIELD_PANGOLIN,
    SITES_FIELD_SPLICEAI,
    SOURCE_CLINVAR_PLP,
    SOURCE_HC_LOF,
    SOURCE_IN_TRANS_OE_CANDIDATE,
    SOURCE_SPLICE_PATH,
)

REF = "GRCh38"
CHR = "chr1"

# Transcript-consequence element type used to build synthetic VEP structs.
_TC_TYPE = hl.tstruct(
    transcript_id=hl.tstr,
    biotype=hl.tstr,
    consequence_terms=hl.tarray(hl.tstr),
    gene_id=hl.tstr,
    gene_symbol=hl.tstr,
    lof=hl.tstr,
)


def _tc(
    gene_id,
    gene_symbol="SYM",
    terms=("missense_variant",),
    lof="",
    transcript_id="ENST0",
    biotype="protein_coding",
):
    """One protein-coding Ensembl transcript_consequences struct."""
    return hl.Struct(
        transcript_id=transcript_id,
        biotype=biotype,
        consequence_terms=list(terms),
        gene_id=gene_id,
        gene_symbol=gene_symbol,
        lof=lof,
    )


def _csqs(tcs):
    """Literal transcript_consequences array expression from _tc structs."""
    return hl.literal(tcs, hl.tarray(_TC_TYPE))


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
        assert row.clinvar.clinvar_review == set()  # reviewed, non-conflicting

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


# ---------------------------------------------------------------------------
# create_variant_filter_ht inputs
# ---------------------------------------------------------------------------
def _clinvar_struct(plp=False, blb=False, vus=False, geneinfo="", review=()):
    """A sites-HT ``clinvar`` struct (is_<cat> bools + clinvar_review + GENEINFO)."""
    flags = dict(zip(CLINVAR_CATEGORIES, (plp, blb, vus)))
    return hl.struct(
        **{CLINVAR_CATEGORY_FIELD_FMT.format(category=c): flags[c] for c in CLINVAR_CATEGORIES},
        clinvar_review=hl.set(hl.literal(list(review), hl.tarray(hl.tstr))),
        GENEINFO=geneinfo,
    )


def _cvf_ht(rows, with_clinvar=False, with_splice=False):
    """Minimal sites HT for create_variant_filter_ht.

    Each row dict: pos, af, tcs (list of _tc structs), filters (list), and —
    when the flag is set — clinvar (dict of plp/blb/vus/geneinfo) and
    spliceai/pangolin floats.
    """
    fields = {
        "locus": hl.tlocus(REF),
        "alleles": hl.tarray(hl.tstr),
        "af": hl.tfloat64,
        "filters": hl.tset(hl.tstr),
        "vep": hl.tstruct(transcript_consequences=hl.tarray(_TC_TYPE)),
    }
    cv_type = hl.tstruct(
        **{CLINVAR_CATEGORY_FIELD_FMT.format(category=c): hl.tbool for c in CLINVAR_CATEGORIES},
        clinvar_review=hl.tset(hl.tstr),
        GENEINFO=hl.tstr,
    )
    if with_clinvar:
        fields[SITES_FIELD_CLINVAR] = cv_type
    if with_splice:
        fields[SITES_FIELD_SPLICEAI] = hl.tfloat32
        fields[SITES_FIELD_PANGOLIN] = hl.tfloat64
    typed = []
    for r in rows:
        _tcs = r.get("tcs", [])
        row = {
            "locus": hl.locus(CHR, r["pos"], REF),
            "alleles": ["A", "T"],
            "af": r["af"],
            "filters": set(r.get("filters", [])),
            # tcs=None -> a NULL vep struct (variant absent from the VEP HT).
            "vep": (
                hl.missing(fields["vep"])
                if _tcs is None
                else hl.Struct(transcript_consequences=_tcs)
            ),
        }
        if with_clinvar:
            cv = r.get("clinvar")
            row[SITES_FIELD_CLINVAR] = (
                hl.missing(cv_type) if cv is None else _clinvar_struct(**cv)
            )
        if with_splice:
            row[SITES_FIELD_SPLICEAI] = r.get("spliceai", 0.0)
            row[SITES_FIELD_PANGOLIN] = r.get("pangolin", 0.0)
        typed.append(hl.Struct(**row))
    return hl.Table.parallelize(typed, hl.tstruct(**fields), key=["locus", "alleles"])


# ===========================================================================
# Filter-step validators
# ===========================================================================
class TestFilterStepValidators:
    def test_clinvar_missing_field_raises(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}])
        with pytest.raises(ValueError, match="missing the 'clinvar'"):
            _validate_clinvar_categories([CLINVAR_CATEGORY_PLP], ht)

    def test_clinvar_unknown_category_raises(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}], with_clinvar=True)
        with pytest.raises(ValueError, match="Unknown ClinVar categories"):
            _validate_clinvar_categories(["not_a_category"], ht)

    def test_clinvar_none_or_empty_ok(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}])
        _validate_clinvar_categories(None, ht)
        _validate_clinvar_categories([], ht)

    def test_clinvar_valid_ok(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}], with_clinvar=True)
        _validate_clinvar_categories([CLINVAR_CATEGORY_PLP], ht)

    def test_pathogenic_splice_missing_raises(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}])
        with pytest.raises(ValueError, match="rebuild it with"):
            _validate_pathogenic_splice(True, ht)

    def test_pathogenic_splice_present_ok(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}], with_splice=True)
        _validate_pathogenic_splice(True, ht)

    def test_pathogenic_splice_false_ok(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01}])
        _validate_pathogenic_splice(False, ht)

    def test_least_consequence_unknown_raises(self):
        with pytest.raises(ValueError, match="not in CSQ_ORDER"):
            _validate_least_consequence("not_a_consequence")

    def test_least_consequence_valid_ok(self):
        _validate_least_consequence("missense_variant")


# ===========================================================================
# _get_vep_gene_id_expr
# ===========================================================================
class TestVepGeneIdExpr:
    def test_keeps_at_or_above_least_consequence_and_dedups(self):
        csqs = _csqs([
            _tc(gene_id="G1", terms=["missense_variant"]),  # kept
            _tc(gene_id="G2", terms=["3_prime_UTR_variant"]),  # below -> dropped
            _tc(gene_id="G1", terms=["stop_gained"]),  # more severe, dup gene
        ])
        genes = hl.eval(_get_vep_gene_id_expr(csqs, "missense_variant"))
        assert set(genes) == {"G1"}

    def test_lowering_threshold_admits_utr(self):
        csqs = _csqs([_tc(gene_id="G2", terms=["3_prime_UTR_variant"])])
        genes = hl.eval(_get_vep_gene_id_expr(csqs, "3_prime_UTR_variant"))
        assert set(genes) == {"G2"}

    def test_no_qualifying_transcript_empty(self):
        csqs = _csqs([_tc(gene_id="G2", terms=["3_prime_UTR_variant"])])
        genes = hl.eval(_get_vep_gene_id_expr(csqs, "missense_variant"))
        assert list(genes) == []


# ===========================================================================
# _get_hc_lof_gene_id_expr
# ===========================================================================
class TestHcLofGeneIdExpr:
    def test_keeps_hc_only_and_dedups(self):
        csqs = _csqs([
            _tc(gene_id="G1", lof="HC"),
            _tc(gene_id="G2", lof="LC"),
            _tc(gene_id="G1", lof="HC"),
        ])
        genes = hl.eval(_get_hc_lof_gene_id_expr(csqs))
        assert set(genes) == {"G1"}

    def test_no_hc_empty(self):
        csqs = _csqs([_tc(gene_id="G1", lof="LC"), _tc(gene_id="G2", lof="")])
        assert list(hl.eval(_get_hc_lof_gene_id_expr(csqs))) == []


# ===========================================================================
# _get_clinvar_gene_id_expr
# ===========================================================================
class TestClinvarGeneIdExpr:
    def test_in_category_matches_geneinfo_symbols(self):
        cv = _clinvar_struct(plp=True, geneinfo="SYMA:111|SYMB:222")
        csqs = _csqs([
            _tc(gene_id="G_A", gene_symbol="SYMA"),
            _tc(gene_id="G_X", gene_symbol="SYMX"),
        ])
        genes, fallback = _get_clinvar_gene_id_expr(cv, csqs, CLINVAR_CATEGORY_PLP)
        assert set(hl.eval(genes)) == {"G_A"}
        assert hl.eval(fallback) is False

    def test_fallback_to_all_vep_when_no_symbol_match(self):
        cv = _clinvar_struct(plp=True, geneinfo="SYMZ:999")  # no csq match
        csqs = _csqs([
            _tc(gene_id="G_A", gene_symbol="SYMA"),
            _tc(gene_id="G_X", gene_symbol="SYMX"),
        ])
        genes, fallback = _get_clinvar_gene_id_expr(cv, csqs, CLINVAR_CATEGORY_PLP)
        assert set(hl.eval(genes)) == {"G_A", "G_X"}
        assert hl.eval(fallback) is True

    def test_not_in_category_is_empty(self):
        cv = _clinvar_struct(plp=False, blb=True, geneinfo="SYMA:111")
        csqs = _csqs([_tc(gene_id="G_A", gene_symbol="SYMA")])
        genes, fallback = _get_clinvar_gene_id_expr(cv, csqs, CLINVAR_CATEGORY_PLP)
        assert list(hl.eval(genes)) == []
        assert hl.eval(fallback) is False


# ===========================================================================
# create_variant_filter_ht
# ===========================================================================
class TestCreateVariantFilterHt:
    def test_base_vep_csq_af_filter_and_drop_empty(self):
        ht = _cvf_ht([
            {"pos": 100, "af": 0.01, "tcs": [_tc("G1", terms=["missense_variant"])]},
            {"pos": 101, "af": 0.9, "tcs": [_tc("G2", terms=["missense_variant"])]},
            {"pos": 102, "af": 0.01, "tcs": [_tc("G3", terms=["3_prime_UTR_variant"])]},
        ])
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05
        )
        by_pos = {r.locus.position: r for r in out.collect()}
        # 100: af ok + qualifying csq -> vep_csq; 101: af too high; 102: no
        # qualifying consequence. Only 100 survives (source non-empty).
        assert set(by_pos) == {100}
        assert by_pos[100].source == {"vep_csq"}
        assert by_pos[100].gene_id == ["G1"]
        assert by_pos[100].af == 0.01  # af carried through for downstream re-thresholding

    def test_null_af_is_kept(self):
        # af=None (variant absent from release freq HT) must NOT be dropped:
        # it qualifies via its VEP consequence and the NULL-tolerant AF gate.
        ht = _cvf_ht([{"pos": 100, "af": None, "tcs": [_tc("G1")]}])
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05
        )
        rows = out.collect()
        assert len(rows) == 1
        assert rows[0].source == {"vep_csq"}
        assert rows[0].af is None

    def test_max_freq_boundary_inclusive(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.05, "tcs": [_tc("G1")]}])
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05
        )
        assert out.count() == 1  # af == max_freq is kept (<=)

    def test_in_trans_oe_candidate_af_threshold(self):
        # af=0.4: above max_freq (0.05) so no vep_csq, but <= in_trans_oe_max_af.
        ht = _cvf_ht([{"pos": 100, "af": 0.4, "tcs": [_tc("G1")]}])
        out = create_variant_filter_ht(
            ht,
            least_consequence="missense_variant",
            max_freq=0.05,
            include_in_trans_oe_candidates=True,
            in_trans_oe_max_af=0.5,
        )
        row = out.collect()[0]
        assert SOURCE_IN_TRANS_OE_CANDIDATE in row.source
        assert "vep_csq" not in row.source

    def test_hc_lof_source(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01, "tcs": [_tc("G1", lof="HC")]}])
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_hc_lof=True,
        )
        assert SOURCE_HC_LOF in out.collect()[0].source

    def test_clinvar_plp_source(self):
        ht = _cvf_ht(
            [{
                "pos": 100, "af": 0.01,
                "tcs": [_tc("G1", gene_symbol="SYMA", terms=["3_prime_UTR_variant"])],
                "clinvar": {"plp": True, "geneinfo": "SYMA:1"},
            }],
            with_clinvar=True,
        )
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_clinvar_categories=[CLINVAR_CATEGORY_PLP],
        )
        row = out.collect()[0]
        # ClinVar has no consequence filter, so a UTR-only variant still tags.
        assert SOURCE_CLINVAR_PLP in row.source
        assert row.gene_id == ["G1"]
        assert row.clinvar_review == set()  # clean record -> no flags

    def test_clinvar_uses_in_trans_oe_af_cap(self):
        # ClinVar tags cap at in_trans_oe_max_af (0.2), not max_freq (0.05):
        # a P/LP at af=0.1 is kept; at af=0.3 it's dropped.
        def _run(af):
            ht = _cvf_ht(
                [{"pos": 100, "af": af,
                  "tcs": [_tc("G1", gene_symbol="SYMA", terms=["3_prime_UTR_variant"])],
                  "clinvar": {"plp": True, "geneinfo": "SYMA:1"}}],
                with_clinvar=True,
            )
            return create_variant_filter_ht(
                ht, least_consequence="missense_variant", max_freq=0.05,
                include_clinvar_categories=[CLINVAR_CATEGORY_PLP],
                in_trans_oe_max_af=0.2,
            ).collect()
        kept = _run(0.1)
        assert len(kept) == 1 and SOURCE_CLINVAR_PLP in kept[0].source
        assert _run(0.3) == []  # above the 0.2 cap, no other source -> dropped

    def test_clinvar_review_flag_carried(self):
        # A conflicting P/LP is still included (clinvar_plp) and its
        # clinvar_review flag is carried through so it's droppable downstream.
        ht = _cvf_ht(
            [{"pos": 100, "af": 0.01,
              "tcs": [_tc("G1", gene_symbol="SYMA")],
              "clinvar": {"plp": True, "geneinfo": "SYMA:1", "review": ["conflicting"]}}],
            with_clinvar=True,
        )
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_clinvar_categories=[CLINVAR_CATEGORY_PLP],
        ).collect()
        assert len(out) == 1
        assert SOURCE_CLINVAR_PLP in out[0].source
        assert out[0].clinvar_review == {"conflicting"}

    def test_records_variant_filter_params_global(self):
        ht = _cvf_ht([{"pos": 100, "af": 0.01, "tcs": [_tc("G1")]}])
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_hc_lof=True,
        )
        g = hl.eval(out.index_globals().variant_filter_params)
        assert g.least_consequence == "missense_variant"
        assert g.max_freq == 0.05
        assert g.include_hc_lof is True
        assert g.include_pathogenic_splice is False


# ===========================================================================
# GENCODE intronic-padding helpers (now take a gencode HT -> unit-testable)
# ===========================================================================
def _gencode_ht(exons):
    """Synthetic GENCODE HT. Each exon: (gene_id, strand, start, end[, feature, transcript_type])."""
    rows = []
    for e in exons:
        rows.append(
            hl.Struct(
                interval=hl.interval(
                    hl.locus(CHR, e[2], REF),
                    hl.locus(CHR, e[3], REF),
                    includes_end=True,
                ),
                feature=e[4] if len(e) > 4 else "exon",
                transcript_type=e[5] if len(e) > 5 else "protein_coding",
                gene_id=e[0],
                strand=e[1],
            )
        )
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            interval=hl.tinterval(hl.tlocus(REF)),
            feature=hl.tstr,
            transcript_type=hl.tstr,
            gene_id=hl.tstr,
            strand=hl.tstr,
        ),
        key=["interval"],
    )


def _loci_ht(positions):
    return hl.Table.parallelize(
        [hl.Struct(locus=hl.locus(CHR, p, REF)) for p in positions],
        hl.tstruct(locus=hl.tlocus(REF)),
        key=["locus"],
    )


class TestNullHandling:
    """Missing-value inputs must never NA-poison the `source` set and drop a
    variant that qualifies for a source (the class of bug behind the NULL-`af`
    drop)."""

    def test_null_vep_not_poisoned_via_splice(self):
        # NULL vep (variant absent from the VEP HT) but flagged by SpliceAI.
        # The vep-derived tags must not NA-poison the source and drop it.
        ht = _cvf_ht(
            [{"pos": 100, "af": 0.01, "tcs": None, "spliceai": 0.9, "pangolin": 0.0}],
            with_splice=True,
        )
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_pathogenic_splice=True,
        ).collect()
        assert len(out) == 1  # must be kept, not dropped
        assert SOURCE_SPLICE_PATH in out[0].source

    def test_null_clinvar_struct_not_poisoned(self):
        # NULL clinvar struct + clinvar categories enabled: no clinvar tag, but
        # must still be kept via vep_csq (clinvar tag False, not NA).
        ht = _cvf_ht(
            [{"pos": 100, "af": 0.01, "tcs": [_tc("G1")], "clinvar": None}],
            with_clinvar=True,
        )
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_clinvar_categories=[CLINVAR_CATEGORY_PLP],
        ).collect()
        assert len(out) == 1
        assert out[0].source == {"vep_csq"}

    def test_null_spliceai_pangolin_not_poisoned(self):
        # NULL splice scores + pathogenic-splice enabled: splice_path False (not
        # NA), variant kept via vep_csq.
        ht = _cvf_ht(
            [{"pos": 100, "af": 0.01, "tcs": [_tc("G1")],
              "spliceai": None, "pangolin": None}],
            with_splice=True,
        )
        out = create_variant_filter_ht(
            ht, least_consequence="missense_variant", max_freq=0.05,
            include_pathogenic_splice=True,
        ).collect()
        assert len(out) == 1
        assert out[0].source == {"vep_csq"}


class TestIntronicPadding:
    def test_build_flanking_intervals_plus_strand(self):
        # + strand gene, exons [100,200] and [300,400]. acceptor=3, donor=8.
        gc = _gencode_ht([("G1", "+", 100, 200), ("G1", "+", 300, 400)])
        iv = _build_intronic_padding_interval_ht(gc, 3, 8).collect()
        pairs = sorted((r.interval.start.position, r.interval.end.position) for r in iv)
        # exon1 donor (right) flank [200,208]; exon2 acceptor (left) flank [297,300).
        # Genomic-edge sides (exon1 left, exon2 right) clip to the gene body.
        assert pairs == [(200, 208), (297, 300)]
        assert {r.gene_id for r in iv} == {"G1"}

    def test_build_minus_strand_swaps_acceptor_donor(self):
        gc = _gencode_ht([("G1", "-", 100, 200), ("G1", "-", 300, 400)])
        iv = _build_intronic_padding_interval_ht(gc, 3, 8).collect()
        pairs = sorted((r.interval.start.position, r.interval.end.position) for r in iv)
        # On '-' strand donor/acceptor swap sides: exon1 right flank = acceptor(3)
        # -> [200,203]; exon2 left flank = donor(8) -> [292,300).
        assert pairs == [(200, 203), (292, 300)]

    def test_build_excludes_non_exon_and_non_protein_coding(self):
        gc = _gencode_ht([
            ("G1", "+", 100, 200),
            ("G1", "+", 300, 400),
            ("G2", "+", 500, 600, "CDS", "protein_coding"),  # not an exon
            ("G3", "+", 700, 800, "exon", "lncRNA"),  # not protein-coding
        ])
        iv = _build_intronic_padding_interval_ht(gc, 3, 8).collect()
        assert {r.gene_id for r in iv} == {"G1"}

    def test_get_intronic_overlap(self):
        gc = _gencode_ht([("G1", "+", 100, 200), ("G1", "+", 300, 400)])
        q = _loci_ht([205, 250, 298])
        q = q.annotate(genes=_get_intronic_padding_gene_id_expr(q, gc, 3, 8))
        by = {r.locus.position: set(r.genes) for r in q.collect()}
        assert by[205] == {"G1"}  # in donor flank [200,208]
        assert by[250] == set()  # deep intron, no interval
        assert by[298] == {"G1"}  # in acceptor flank [297,300)

    def test_get_intronic_donor_expansion(self):
        # donor=12 captures +9..+12 that donor=8 misses (position 210).
        gc = _gencode_ht([("G1", "+", 100, 200), ("G1", "+", 300, 400)])
        q = _loci_ht([210])
        d8 = q.annotate(g=_get_intronic_padding_gene_id_expr(q, gc, 3, 8)).g.collect()[0]
        d12 = q.annotate(g=_get_intronic_padding_gene_id_expr(q, gc, 3, 12)).g.collect()[0]
        assert list(d8) == []
        assert set(d12) == {"G1"}


def _pair_mt(carriers, samples):
    """Tiny GT MatrixTable: 3 variants at chr1:1-3, cols=``samples``.

    ``carriers`` = set of (row_idx, col_idx) getting a het GT (rest hom-ref).
    """
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=len(samples))
    mt = mt.annotate_rows(locus=hl.locus(CHR, mt.row_idx + 1, REF), alleles=["A", "C"])
    mt = mt.annotate_cols(s=hl.array(list(samples))[mt.col_idx])
    car = hl.literal(set(carriers))
    mt = mt.annotate_entries(
        GT=hl.if_else(
            car.contains((mt.row_idx, mt.col_idx)), hl.call(0, 1), hl.call(0, 0)
        )
    )
    return mt.key_rows_by("locus", "alleles").key_cols_by("s").select_entries("GT")


def _pair_filter_ht(mt):
    return mt.rows().annotate(
        gene_id=hl.set(["G1"]), an_pct=100.0, source=hl.set(["vep_csq"])
    )


class TestCreateVariantPairHtFlags:
    # s1 = release+trio (carries v1,v2); s2 = release only (v2,v3);
    # s3 = neither (v1,v3). s3 is absent from the subset HT -> both flags False.
    _CARRIERS = {(0, 0), (1, 0), (1, 1), (2, 1), (0, 2), (2, 2)}
    _SAMPLES = ["s1", "s2", "s3"]

    def _subset_ht(self):
        return hl.Table.parallelize(
            [
                {"s": "s1", "is_release": True, "is_trio": True},
                {"s": "s2", "is_release": True, "is_trio": False},
            ],
            hl.tstruct(s=hl.tstr, is_release=hl.tbool, is_trio=hl.tbool),
            key="s",
        )

    def test_in_release_in_trios_flags(self):
        mt = _pair_mt(self._CARRIERS, self._SAMPLES)
        ht = create_variant_pair_ht(
            mt, _pair_filter_ht(mt), sample_subset_ht=self._subset_ht()
        )
        got = {
            (r.locus1.position, r.locus2.position): (r.in_release, r.in_trios)
            for r in ht.collect()
        }
        assert got == {
            (1, 2): (True, True),
            (2, 3): (True, False),
            (1, 3): (False, False),
        }

    def test_flags_absent_without_subset_ht(self):
        mt = _pair_mt(self._CARRIERS, self._SAMPLES)
        ht = create_variant_pair_ht(mt, _pair_filter_ht(mt))
        assert "in_release" not in set(ht.row)
        assert "in_trios" not in set(ht.row)
        assert ht.count() == 3

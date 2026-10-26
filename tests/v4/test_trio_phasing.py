"""Pytest units for ``gnomad_chets.v4.trio_phasing`` on tiny in-memory inputs."""
import hail as hl
import pytest

from gnomad_chets.v4.trio_phasing import call_trio_chet, call_trio_chet_within_genes


@pytest.fixture(scope="session", autouse=True)
def _hail_session():
    hl.init(idempotent=True, quiet=True)
    yield


def _pbt_mt():
    # Two probands, three variants. p0: 0|1 at v0 and v1 (cis). p1: 1|0 at v0,
    # 0|1 at v1 (trans), and an unphased het at v2, which no proband carries
    # as a phased het.
    calls = [
        ["0|1", "1|0"],
        ["0|1", "0|1"],
        ["0|0", "0/1"],
    ]
    mt = hl.utils.range_matrix_table(3, 2)
    mt = mt.key_rows_by(
        locus=hl.locus("chr1", 100 + 100 * mt.row_idx, "GRCh38"),
        alleles=hl.literal(["A", "T"]),
    )
    mt = mt.key_cols_by(s=hl.str("p") + hl.str(mt.col_idx))
    mt = mt.annotate_entries(
        PBT_GT=hl.parse_call(hl.literal(calls)[mt.row_idx][mt.col_idx]),
        trio_adj=True,
    )
    return mt.select_entries("PBT_GT", "trio_adj")


def _pairs():
    def v(i):
        return hl.locus("chr1", 100 + 100 * i, "GRCh38"), ["A", "T"]

    rows = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        (l1, a1), (l2, a2) = v(i), v(j)
        rows.append(hl.Struct(locus1=l1, alleles1=a1, locus2=l2, alleles2=a2))
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus1=hl.tlocus("GRCh38"), alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh38"), alleles2=hl.tarray(hl.tstr),
        ),
        key=["locus1", "alleles1", "locus2", "alleles2"],
    )


def test_call_trio_chet_counts_shared_phased_probands():
    # Only (v0, v1) has a proband phased het at both: p0 cis, p1 trans. The
    # pairs touching v2 have no phased-het carrier there and are dropped.
    out = call_trio_chet(_pbt_mt(), _pairs()).collect()
    assert len(out) == 1
    r = out[0]
    assert (r.locus1.position, r.locus2.position) == (100, 200)
    assert (r.raw.n_same_hap, r.raw.n_chet) == (1, 1)
    assert (r.adj.n_same_hap, r.adj.n_chet) == (1, 1)


def _pops():
    return hl.Table.parallelize(
        [hl.Struct(s="p0", pop="afr"), hl.Struct(s="p1", pop="nfe")],
        hl.tstruct(s=hl.tstr, pop=hl.tstr),
        key="s",
    )


def _genes():
    rows = [
        hl.Struct(locus=hl.Locus("chr1", 100 + 100 * i, "GRCh38"), alleles=["A", "T"], gene_id={"G1"})
        for i in range(3)
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr), gene_id=hl.tset(hl.tstr)),
        key=["locus", "alleles"],
    )


def test_call_trio_chet_by_pop_splits_probands_by_group():
    # p0 (afr) is cis and p1 (nfe) trans at (v0, v1); "all" keeps both, and the
    # within-gene self-join gives the same per-group counts.
    for out in (
        call_trio_chet(_pbt_mt(), _pairs(), pop_ht=_pops()).collect(),
        call_trio_chet_within_genes(_pbt_mt(), _pairs(), _genes(), n_partitions=2, pop_ht=_pops()).collect(),
    ):
        assert len(out) == 1
        by = out[0].adj_by_pop
        assert (by["afr"].n_same_hap, by["afr"].n_chet) == (1, 0)
        assert (by["nfe"].n_same_hap, by["nfe"].n_chet) == (0, 1)
        assert (by["all"].n_same_hap, by["all"].n_chet) == (1, 1)

"""Unit tests for ``gnomad_chets.v4.phase_from_vcf``.

Run with::

    pytest tests/v4/test_phase_from_vcf.py -v

The script's whole premise is that statistical phase can be fed to the existing
read-backed-phase machinery by setting ``PGT = GT`` and a constant ``PID``. That
premise is load-bearing, so these tests drive the real
``encode_genotypes`` / ``count_all_pairs_via_index`` path rather than
re-implementing the counting, and check the resulting ``n_phased_cis`` /
``n_phased_trans`` against hand-computed values.

Tables are built with ``hl.Table.parallelize`` / ``range_matrix_table`` and run
against the local Spark backend; the Hail session comes from ``tests/conftest.py``.
"""

import hail as hl
import pytest

from gnomad_chets.v4.phase_from_vcf import (
    annotate_vcf_chet,
    call_vcf_chet,
    prepare_phased_mt,
)

REFERENCE_GENOME = "GRCh37"
CONTIG = "21"


def _phased_mt(genotypes) -> hl.MatrixTable:
    """Build a phased MatrixTable from per-variant lists of ``(hap0, hap1)``.

    :param genotypes: ``genotypes[variant][sample]`` giving the two alleles.
    :return: MatrixTable keyed by ``(locus, alleles)`` with a phased ``GT``.
    """
    n_variants = len(genotypes)
    n_samples = len(genotypes[0])
    mt = hl.utils.range_matrix_table(n_variants, n_samples)
    lit = hl.literal(genotypes)
    mt = mt.annotate_entries(
        GT=hl.call(
            lit[mt.row_idx][mt.col_idx][0],
            lit[mt.row_idx][mt.col_idx][1],
            phased=True,
        )
    )
    mt = mt.annotate_cols(s=hl.str(mt.col_idx)).key_cols_by("s")
    mt = mt.annotate_rows(
        locus=hl.locus(
            CONTIG, 1000 + mt.row_idx * 100, reference_genome=REFERENCE_GENOME
        ),
        alleles=["A", "C"],
    ).key_rows_by("locus", "alleles")
    return mt.drop("row_idx", "col_idx")


def _pair_ht(pos1: int, pos2: int) -> hl.Table:
    """Build a one-row variant pair Table for the given positions."""
    return hl.Table.parallelize(
        [
            hl.struct(
                locus1=hl.locus(CONTIG, pos1, REFERENCE_GENOME),
                alleles1=["A", "C"],
                locus2=hl.locus(CONTIG, pos2, REFERENCE_GENOME),
                alleles2=["A", "C"],
            )
        ],
        schema=hl.tstruct(
            locus1=hl.tlocus(REFERENCE_GENOME),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus(REFERENCE_GENOME),
            alleles2=hl.tarray(hl.tstr),
        ),
    ).key_by("locus1", "alleles1", "locus2", "alleles2")


class TestPreparePhasedMt:
    """``prepare_phased_mt`` must present statistical phase as physical phase."""

    def test_adds_pgt_pid_and_adj(self):
        """PGT mirrors GT, PID is the constant phase set, adj is set."""
        mt = _phased_mt([[(0, 1), (1, 0)]])
        out = prepare_phased_mt(mt, "chr21")
        assert "PGT" in out.entry
        assert "PID" in out.entry
        assert "adj" in out.entry

        rows = out.entries().collect()
        assert all(r.PGT == r.GT for r in rows)
        assert all(r.PID == "chr21" for r in rows)
        assert all(r.adj for r in rows)

    def test_pid_missing_for_unphased_calls(self):
        """An unphased call gets no PID, so it cannot be paired across variants."""
        mt = _phased_mt([[(0, 1)]])
        mt = mt.annotate_entries(GT=hl.call(0, 1, phased=False))
        out = prepare_phased_mt(mt, "chr21")
        assert out.entries().collect()[0].PID is None


class TestAnnotateVcfChet:
    """The tri-state call must match ``trio_phasing._chet_call`` semantics."""

    @pytest.mark.parametrize(
        "n_cis,n_trans,expected",
        [
            (0, 3, True),  # trans only
            (3, 0, False),  # cis only
            (2, 2, None),  # contradictory across samples
            (0, 0, None),  # no double-het carrier
        ],
    )
    def test_call(self, n_cis, n_trans, expected):
        """True is trans, False is cis, missing is no or mixed evidence."""
        ht = hl.Table.parallelize(
            [hl.struct(n_phased_cis=n_cis, n_phased_trans=n_trans)],
            schema=hl.tstruct(n_phased_cis=hl.tint32, n_phased_trans=hl.tint32),
        )
        row = annotate_vcf_chet(ht).collect()[0]
        assert row.vcf_chet == expected
        assert row.n_vcf_phased == n_cis + n_trans


class TestCallVcfChetEndToEnd:
    """Drive the real encoder and counter, not a re-implementation."""

    def test_cis_and_trans_counts(self, tmp_path):
        """Double-het samples split into cis and trans by haplotype agreement.

        s0 carries both alts on hap1 and s2 carries both on hap0, so both are
        cis. s1 and s3 carry one alt on each haplotype, so both are trans. s4 is
        not a double het and must not be counted either way.
        """
        genotypes = [
            # v1 at 1000: s0    s1      s2      s3      s4
            [(0, 1), (0, 1), (1, 0), (1, 0), (0, 0)],
            # v2 at 1100
            [(0, 1), (1, 0), (1, 0), (0, 1), (0, 1)],
        ]
        mt = _phased_mt(genotypes)
        vp_ht = _pair_ht(1000, 1100)

        ht = call_vcf_chet(mt, vp_ht, str(tmp_path), "chr21")
        row = ht.collect()[0]

        assert row.n_phased_cis == 2
        assert row.n_phased_trans == 2
        assert row.n_vcf_phased == 4
        # Mixed evidence, so no single call for the pair.
        assert row.vcf_chet is None

    def test_reused_encode_gives_the_same_counts(self, tmp_path):
        """``encoded_dir`` skips the encode and counts from an earlier run's."""
        genotypes = [
            [(0, 1), (0, 1), (1, 0), (1, 0), (0, 0)],
            [(0, 1), (1, 0), (1, 0), (0, 1), (0, 1)],
        ]
        mt = _phased_mt(genotypes)
        vp_ht = _pair_ht(1000, 1100)
        first = call_vcf_chet(mt, vp_ht, str(tmp_path / "first"), "chr21").collect()
        (encode_dir,) = (tmp_path / "first").glob("phase_from_vcf_encoded_*")

        again = call_vcf_chet(
            mt, vp_ht, str(tmp_path / "again"), "chr21", encoded_dir=str(encode_dir)
        ).collect()

        assert again == first
        assert not (tmp_path / "again").exists()

    def test_all_cis_gives_a_cis_call(self, tmp_path):
        """A pair where every double het is cis resolves to a cis call."""
        genotypes = [
            [(0, 1), (1, 0), (0, 0)],
            [(0, 1), (1, 0), (0, 1)],
        ]
        mt = _phased_mt(genotypes)
        ht = call_vcf_chet(mt, _pair_ht(1000, 1100), str(tmp_path), "chr21")
        row = ht.collect()[0]

        assert row.n_phased_cis == 2
        assert row.n_phased_trans == 0
        assert row.vcf_chet is False

    def test_phase_counts_do_not_exceed_the_double_het_cell(self, tmp_path):
        """``n_phased_cis + n_phased_trans`` cannot exceed the ``AaBb`` cell.

        With every genotype phased the two should be equal; the invariant is
        checked as an inequality so it still holds if a caller supplies a
        partially phased VCF.
        """
        genotypes = [
            [(0, 1), (1, 0), (1, 1), (0, 0)],
            [(0, 1), (0, 1), (0, 1), (0, 1)],
        ]
        mt = _phased_mt(genotypes)
        ht = call_vcf_chet(mt, _pair_ht(1000, 1100), str(tmp_path), "chr21")
        row = ht.collect()[0]

        double_het = row.gt_counts_adj[4]
        assert row.n_phased_cis + row.n_phased_trans <= double_het
        assert row.n_phased_cis + row.n_phased_trans == double_het

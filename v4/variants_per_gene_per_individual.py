"""
Count variants per gene per individual in gnomAD individual-level data,
broken down by variant class (lof / missense / synonymous / noncoding /
other), using VEP's `most_severe_consequence` field. Works with gnomAD v2
(exomes only) or v4 (exomes) layouts.

Gene assignment uses ONLY the canonical transcript (per VEP's `canonical`
flag): each variant contributes to at most one gene, via that gene's
canonical transcript. Variants with no canonical-transcript hit (e.g.
purely intergenic) are dropped. This means each variant is counted once,
not fanned out across every gene any transcript happens to touch -- so
there's no row explode anywhere in this pipeline.

Output: a Hail Table keyed by (gene_symbol, variant_class, s) with
`n_variants` (distinct variant count).

v4 loading:
  - Genotypes come from `get_gnomad_v4_vds()` (broadinstitute/gnomad_qc),
    not the raw VDS path directly, so they're pre-filtered: drops the
    excessively multiallelic chr19:5787204 site, removes duplicate/
    withdrawn UKB samples, and removes hard-filtered samples. That site
    drop is kept -- it blows up split_multi, and the co-occurrence
    pipeline drops it too, so suppressing it would both cost compute and
    diverge from the table we compare against.
  - We do NOT densify. `to_dense_mt()` fills in explicit hom-ref calls for
    every sample at every site by merging in the VDS's reference-block
    data -- expensive, and unnecessary here, since we only ever count
    non-ref genotypes. `vds.variant_data` (after split=True) is already
    the sparse table of actual variant calls, which is exactly what's
    needed; densifying first and then filtering to non-ref later would
    do a huge amount of wasted work materializing hom-ref calls just to
    immediately throw them away.
  - VEP is NOT part of the raw v4 genotype VDS (confirmed against
    gnomad_qc.v4.resources.annotations -- get_vep() is a wholly separate
    VersionedTableResource, exactly analogous to v2). So v4 needs the
    same style of VEP join as v2 -- but that join is deferred to
    join_vep_late(), called from main() only after entry/row filtering,
    not done at load time. load_matrix_table() returns genotypes only.
  Requires a gnomad_qc recent enough that get_gnomad_v4_vds() accepts
  `filter_intervals` -- older releases lack it and the call raises TypeError.
  This requires the gnomad_qc package to be importable on the cluster:
    pip install git+https://github.com/broadinstitute/gnomad_qc.git
  (gnomad_qc itself depends on gnomad_methods, which pip pulls in as
  `gnomad`). It also requires gnomAD-team GCS bucket access -- the
  underlying data is not public.

v2 loading:
  - Restricted to exomes (gnomAD v2 chet/individual-level work is
    exomes-only in practice).
  - v2's hardcalls MT (gnomad.exomes.mt) also has NO `vep` field of its
    own -- genotypes only. VEP lives on the public release sites Table
    (get_gnomad_public_data("exomes")), same as v4's separate VEP Table --
    but again, that join is deferred to join_vep_late(), not done here.
  - There is no v2 equivalent of v4's hardcoded chr19:5787204 site drop in
    gnomad_qc (checked v2/resources/variant_qc.py, found none). v2's
    actual mechanism for excluding problematic sites is the release
    Table's `filters` field (RF hard/soft filtering: PASS vs
    AC0/RF/InbreedingCoeff/etc.). main() applies that filter once for BOTH
    builds from the release sites Table (--skip-filter-pass to disable).

Usage:
    hailctl dataproc submit <cluster> variants_per_gene_per_individual.py \
        --output-prefix gs://path/to/output_prefix \
        --gnomad-version v4 \
        [--mt-path gs://path/to/gnomad_data]   # overrides the version default \
        [--interval-path gs://path/to/genes.bed] \
        [--skip-filter-pass]                    # v2 only
"""

import argparse
from typing import List, Optional, Union

import hail as hl
from gnomad.assessment.summary_stats import freq_bin_expr
from gnomad.resources.grch37.gnomad import public_release as public_release_grch37
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.vep import (
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_LOW_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    CSQ_NON_CODING,
    filter_vep_transcript_csqs_expr,
    get_most_severe_consequence_expr,
)

# A couple of well-known genes for quick --gene validation runs, both
# builds. Add more as needed, or just pass --gene-interval directly.
COHORT_COUNT_INTERVAL = {"v2": "1:1-2", "v4": "chr1:1-2"}
"""Tiny interval used only to count cohort columns.

``--summary-stats-only`` needs the cohort size but not the genotypes. Column
count is independent of which rows are read, so loading this one-locus slice
is exact and cheap, where counting distinct samples in the per-individual table
would miss everyone carrying nothing in scope.
"""

KNOWN_GENE_INTERVALS = {
    "v2": {  # GRCh37
        "BRCA1": "17:41196312-41277500",
        "PCSK9": "1:55505221-55530525",
    },
    "v4": {  # GRCh38
        "BRCA1": "chr17:43044295-43125364",
        "PCSK9": "chr1:55039475-55064852",
        # The co-occurrence pipeline's 5-gene test bundle, so runs here can be
        # compared directly against the v4 numbers already computed for them at
        # gs://gnomad/v4.1/variant_cooccurrence/test-5-gene/. Intervals match
        # v4/resources.py TEST_INTERVALS exactly.
        "AHNAK2": "chr14:104937244-104978374",
        "ANO5": "chr11:21782659-22283567",
        "CAPN3": "chr15:42359498-42412949",
        "DYSF": "chr2:71453561-71686763",
        "SGCA": "chr17:50164214-50175928",
    },
}

# --- consequence_terms -> variant class (classified off the CANONICAL
# transcript's own terms, not the variant-wide most_severe_consequence
# field -- see canonical_transcript_annotations_expr below for why) ------

# Variant classes are derived from gnomad_methods' severity buckets rather than
# hand-listed, so every departure from that grouping is visible as an explicit
# delta below. Those buckets are current as of VEP v105
# (gnomad.utils.vep.CURRENT_VEP_VERSION) with some terms kept for backwards
# compatibility, and are "loosely based on VEP's categories but ... adjusted to
# better serve gnomAD's use cases" -- so a delta here may agree with VEP v105
# even where it disagrees with gnomad_methods. Drop the deltas to get the
# standard gnomAD grouping.
_LOF_EXTRA = {"start_lost"}  # Considered high impact in VEP v105, previously medium.
"""Terms added to ``CSQ_CODING_HIGH_IMPACT`` to form the lof class."""

_MISSENSE_EXCLUDE = {
    "start_lost",  # Considered high impact in VEP v105, previously medium.
}
"""``CSQ_CODING_MEDIUM_IMPACT`` terms not counted as missense.

Only ``start_lost``, which moves to lof via :data:`_LOF_EXTRA`. Every other
MEDIUM-impact term is kept, so no consequence falls through to "other".
"""

_LOW_IMPACT_AS_NONCODING = {
    "splice_region_variant",  # Considered low impact in VEP v105, previously medium.
    "coding_sequence_variant",  # Considered modifier/non-coding in VEP v105, but keeping as low.
    "incomplete_terminal_codon_variant",
}
"""``CSQ_CODING_LOW_IMPACT`` terms counted as noncoding instead of synonymous."""

_LOF_TERMS = hl.set(set(CSQ_CODING_HIGH_IMPACT) | _LOF_EXTRA)
"""Consequence terms counted as loss-of-function."""

_MISSENSE_TERMS = hl.set(set(CSQ_CODING_MEDIUM_IMPACT) - _MISSENSE_EXCLUDE)
"""Consequence terms counted as missense."""

_SYNONYMOUS_TERMS = hl.set(set(CSQ_CODING_LOW_IMPACT) - _LOW_IMPACT_AS_NONCODING)
"""Consequence terms counted as synonymous."""

_NONCODING_TERMS = hl.set(set(CSQ_NON_CODING) | _LOW_IMPACT_AS_NONCODING)
"""Consequence terms counted as noncoding."""


def classify_consequence_expr(
    consequence: hl.expr.StringExpression,
) -> hl.expr.StringExpression:
    """
    Map a single VEP consequence term to a variant class.

    Takes the term already chosen by
    :func:`gnomad.utils.vep.get_most_severe_consequence_expr`, so no priority
    arbitration happens here -- ``CSQ_ORDER`` has done the ranking. The four
    term sets are disjoint, so the lookup is unambiguous.

    A missing ``consequence`` classifies as "other" rather than propagating,
    because ``hl.set.contains`` on a missing value is False. That is the wanted
    behaviour here: :func:`gnomad.utils.vep.get_most_severe_consequence_expr`
    returns missing when NO term is in ``CSQ_ORDER``, i.e. for a novel VEP term,
    and "other" is the right bucket for one. A variant with no canonical
    transcript at all is already excluded upstream by the caller, so that case
    never reaches here.

    Every term in ``CSQ_ORDER`` is classified, which matters here: this takes
    the MOST SEVERE term, so an unclassified top term would send the whole
    transcript to "other" even when a lesser classified term is present. Only
    two deltas from gnomad_methods' buckets remain, both deliberate --
    :data:`_LOF_EXTRA` and :data:`_LOW_IMPACT_AS_NONCODING`.

    :param consequence: One VEP consequence term, typically from
        :func:`gnomad.utils.vep.get_most_severe_consequence_expr`.
    :return: One of "lof", "missense", "synonymous", "noncoding", "other".
    """
    return (
        hl.case()
        .when(_LOF_TERMS.contains(consequence), "lof")
        .when(_MISSENSE_TERMS.contains(consequence), "missense")
        .when(_SYNONYMOUS_TERMS.contains(consequence), "synonymous")
        .when(_NONCODING_TERMS.contains(consequence), "noncoding")
        .default("other")
    )


def canonical_transcript_annotations_expr(
    vep_struct: hl.expr.StructExpression,
) -> hl.expr.StructExpression:
    """gene_symbol AND variant_class, both derived from the SAME canonical
    transcript -- so variant_class always reflects what's actually driving
    the gene assignment. Using vep.most_severe_consequence instead would be
    wrong here: it's a variant-wide field reflecting the single worst
    consequence across ALL transcripts, which could come from a different,
    non-canonical transcript than the one whose gene we're attributing the
    variant to (e.g. 'stop_gained' on some other transcript but merely
    'intron_variant' on the gene's own canonical transcript).

    gene_symbol is missing if there's no canonical-transcript hit (e.g.
    purely intergenic) -- deliberately a single scalar, not a set over all
    transcripts, so each variant maps to at most one gene and there's no
    row explode anywhere in this pipeline.

    Transcripts are restricted to canonical AND protein-coding AND Ensembl
    (``transcript_id`` starting "ENST"). Each rules out a real failure mode:
    without ``protein_coding`` a variant can be attributed to a lncRNA or
    pseudogene, and without ``ensembl_only`` it can land on a RefSeq transcript
    whose ``gene_id`` is an Entrez id rather than an ENSG -- the cause of genes
    silently splitting across two id spaces. A variant with no transcript
    meeting all three gets a missing gene and is dropped downstream.

    NOTE most variants still have several qualifying transcripts, and
    ``[0]`` takes whichever VEP listed first, with no tie-break on
    ``mane_select``. That only affects which of several equally-canonical
    protein-coding transcripts supplies the gene symbol.

    :param vep_struct: The full VEP struct for a variant.
    :return: Struct with ``gene_symbol`` and ``variant_class``, both missing
        when the variant has no canonical protein-coding Ensembl transcript.
    """
    canonical_tcs = filter_vep_transcript_csqs_expr(
        vep_struct.transcript_consequences,
        canonical=True,
        ensembl_only=True,
        protein_coding=True,
    )
    # Taking [0] is a real choice, not a formality. Once the filter is canonical
    # AND protein-coding AND Ensembl, having several qualifying transcripts
    # almost always means the variant sits in the canonical transcript of two
    # OVERLAPPING protein-coding genes -- so [0] is picking which gene the
    # variant is attributed to, not merely which transcript of one gene. On
    # chr20: 265,409 of 1,640,361 variants (16.2%) have more than one qualifying
    # transcript, and for 264,232 of those (99.6%) the transcripts belong to
    # different genes.
    #
    # mane_select does not fix this: all 265,409 have a MANE Select transcript
    # and it differs from [0] in only 10,842 cases, because each overlapping
    # gene has its own MANE transcript. Resolving it properly needs a rule for
    # choosing BETWEEN genes (most severe consequence, or emitting both and
    # accepting the row explode this design avoids).
    primary_tc = hl.or_missing(hl.len(canonical_tcs) > 0, canonical_tcs[0])
    return hl.struct(
        gene_symbol=primary_tc.gene_symbol,
        variant_class=hl.or_missing(
            hl.is_defined(primary_tc),
            classify_consequence_expr(
                get_most_severe_consequence_expr(primary_tc.consequence_terms)
            ),
        ),
    )


# --- AF binning -----------------------------------------------------------

AF_CUTOFFS = [1e-4, 1e-3, 1e-2, 0.05, 0.1]
"""AF bin edges handed to :func:`gnomad.assessment.summary_stats.freq_bin_expr`."""

AF_UPPER = 0.5
"""Top AF cutoff. Variants above it get their own bin.

Binning is on the alt-allele AF exactly as released, so this bin is real and
holds variants whose alt allele is the major one.
"""

AF_BIN_ORDER = [
    "<0.01%",
    "0.01% - 0.1%",
    "0.1% - 1.0%",
    "1.0% - 5.0%",
    "5.0% - 10.0%",
    "10.0% - 50.0%",
    ">50.0%",
]
"""The bins :data:`AF_CUTOFFS` and :data:`AF_UPPER` produce, in increasing-AF order.

``freq_bin_expr`` labels do NOT sort lexically -- "10.0% - 50.0%" sorts before
"5.0% - 10.0%", and "<0.01%" sorts last -- so anything needing them in order
(TSV columns, plot axes, report tables) must use this list rather than sorting
the strings.
"""


def af_bin_expr(freq_expr: hl.expr.StructExpression) -> hl.expr.StringExpression:
    """
    Bucket a frequency struct into the reporting bins.

    A thin wrapper over :func:`gnomad.assessment.summary_stats.freq_bin_expr`
    pinned to this script's cutoffs. ``ac_cutoffs`` is passed an empty list, not
    None: upstream types it ``Optional`` but calls ``sorted()`` on it, so None
    raises. Empty disables the AC0/singleton/doubleton bins, which would
    otherwise subdivide the rarest AF bin.

    :param freq_expr: Frequency struct with ``AC`` and ``AF`` fields, e.g.
        ``freq[0]`` from a gnomAD release sites Table.
    :return: One of :data:`AF_BIN_ORDER`, or "Missing" when ``AC`` is missing.
    """
    return freq_bin_expr(
        freq_expr, ac_cutoffs=[], af_cutoffs=AF_CUTOFFS, upper_af=AF_UPPER
    )


# --- Loading -------------------------------------------------------------

V2_EXOMES_HARDCALLS_MT_PATH = (
    "gs://gnomad_v2/hardcalls/hail-0.2/mt/exomes/gnomad.exomes.mt"
)
V4_RAW_EXOMES_VDS_PATH = "gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds"
V4_VDS_INTERVALS_PATH = "gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.intervals.he"
"""Exact partition intervals of the raw v4 exomes VDS variant data.

Reading the sites Table with these boundaries makes its join against the
genotype MT a zip of aligned partitions rather than a shuffle. They are the
bounds Hail wrote to disk for the untouched VDS, so every caller aligns to the
same ones; see analysis/write_v4_vds_intervals.py, which produced this file.
"""


def load_matrix_table(
    mt_path: Optional[str],
    gnomad_version: str,
    release_only: bool,
    high_quality_only: bool,
    chrom: Optional[str] = None,
    skip_v4_qc_wrapper: bool = False,
    filter_intervals: Optional[List[Union[str, hl.Interval]]] = None,
) -> hl.MatrixTable:
    """Returns genotypes only -- NO vep annotation yet, and NO PASS filter for
    either build: main() applies that once from the release sites Table after
    this returns -- and neither is the VEP annotation, which main() picks up
    from that same sites Table in the same join. VEP annotation is the most
    expensive join in this pipeline (full transcript_consequences arrays
    per row), so it should touch the smallest possible set of rows --
    joining it here, before non-ref filtering, would mean annotating rows
    that might get dropped moments later for having zero carriers.

    filter_intervals (from --gene or --interval-path) MUST be applied INSIDE this
    function, not by the caller after it returns. For v4 specifically,
    get_gnomad_v4_vds() runs an eager vds.variant_data.count_cols() (a
    real Spark job, not lazy) partway through, and a full split_multi at
    the end -- both over the WHOLE genome-wide callset if no interval has
    been applied yet. Restricting afterward (e.g. in main()) is too late:
    by then get_gnomad_v4_vds() has already scanned everything. Passing
    the interval into get_gnomad_v4_vds()'s own filter_intervals param
    applies it near the top of that function, before both of those.
    :param mt_path: Explicit MatrixTable/VDS path. None uses the version default
        (v4: ``get_gnomad_v4_vds()``; v2: the exomes hardcalls MT).
    :param gnomad_version: "v2" or "v4".
    :param release_only: v4 only. Restrict to release samples.
    :param high_quality_only: v4 only. Restrict to high-quality samples.
    :param chrom: Restrict to one chromosome; either naming convention is accepted.
    :param skip_v4_qc_wrapper: v4 only, ignored when ``mt_path`` is set. Read the
        raw VDS directly, bypassing ALL of ``get_gnomad_v4_vds()``'s QC steps.
    :param filter_intervals: Intervals applied INSIDE ``get_gnomad_v4_vds()``,
        before its eager ``count_cols()`` and ``split_multi``. Must be passed here
        rather than applied to the return value.
    :return: MatrixTable of genotypes only, with no VEP annotation.
    """
    if gnomad_version == "v4":
        norm_chrom = f"chr{chrom.replace('chr', '')}" if chrom else None

        if mt_path is not None or skip_v4_qc_wrapper:
            # Bypass get_gnomad_v4_vds() -- either the user gave an
            # explicit path, or explicitly asked to skip it
            # (--skip-v4-qc-wrapper). NOTE this also skips get_gnomad_v4_vds's
            # OTHER QC steps, not just the chr19:5787204 multiallelic drop:
            # duplicate/withdrawn UKB sample removal and hard-filtered
            # sample removal are bundled into that same function and have
            # no separate on/off switch. Fine for a quick single-gene
            # smoke test; reconsider for a real production run.
            read_path = mt_path or V4_RAW_EXOMES_VDS_PATH
            if skip_v4_qc_wrapper and mt_path is None:
                print(
                    f"--skip-v4-qc-wrapper: reading raw VDS directly from {read_path}, "
                    "bypassing get_gnomad_v4_vds() entirely -- no chr19:5787204 drop, no "
                    "duplicate/withdrawn UKB removal, no hard-filtered sample removal, "
                    "and NO release_only/high_quality_only restriction: every sample in "
                    "the VDS is counted."
                )
            vds = hl.vds.read_vds(read_path)
            if norm_chrom:
                vds = hl.vds.filter_chromosomes(vds, keep=[norm_chrom])
            if filter_intervals:
                # Before split_multi, same reasoning as get_gnomad_v4_vds's
                # own filter_intervals-before-split ordering below.
                vds = hl.vds.filter_intervals(
                    vds,
                    parse_intervals(
                        filter_intervals,
                        vds.reference_data.locus.dtype.reference_genome,
                    ),
                )
            # A freshly-read VDS is unsplit (LGT/LA, not global GT), so
            # split it here.
            vds = hl.vds.split_multi(vds)
        else:
            from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

            print(
                "Loading gnomAD v4 exomes via get_gnomad_v4_vds() "
                f"(release_only={release_only}, high_quality_only={high_quality_only}, "
                f"filter_intervals={filter_intervals}); "
                "this removes duplicate/withdrawn UKB samples and hard-filtered "
                "samples, and drops the excessively multiallelic chr19:5787204 site."
            )
            vds = get_gnomad_v4_vds(
                split=True,
                remove_hard_filtered_samples=True,
                high_quality_only=high_quality_only,
                release_only=release_only,
                chrom=norm_chrom,
                # Applied inside get_gnomad_v4_vds BEFORE its eager
                # count_cols() and split_multi -- this is the whole
                # point of threading it through rather than filtering
                # the return value in main().
                filter_intervals=filter_intervals,
            )

        # No densify: variant_data (post split=True) is already the
        # sparse table of actual variant calls -- exactly what's needed
        # to count non-ref genotypes, without materializing hom-ref calls
        # for every sample at every site first.
        return vds.variant_data

    else:  # v2, exomes only
        mt_path = mt_path or V2_EXOMES_HARDCALLS_MT_PATH
        if mt_path != V2_EXOMES_HARDCALLS_MT_PATH:
            print(f"Using custom v2 path: {mt_path}")
        mt = hl.read_matrix_table(mt_path)
        if chrom:
            mt = restrict_to_chrom(mt, chrom, gnomad_version)
        if filter_intervals:
            # hl.read_matrix_table is lazy (unlike get_gnomad_v4_vds,
            # nothing eager happens above), but restricting here still
            # keeps the PASS-filter join below from dealing with more
            # rows than necessary.
            mt = hl.filter_intervals(
                mt,
                parse_intervals(filter_intervals, mt.locus.dtype.reference_genome),
            )
        return mt


def parse_intervals(
    intervals: List[Union[str, hl.Interval]], reference_genome: str
) -> List[hl.Interval]:
    """
    Parse any string entries in an interval list, passing the rest through.

    :param intervals: Intervals as strings, `hl.Interval`, or a mix.
    :param reference_genome: Reference genome used to parse the strings.
    :return: List of `hl.Interval`.
    """
    return [
        hl.parse_locus_interval(i, reference_genome=reference_genome)
        if isinstance(i, str)
        else i
        for i in intervals
    ]


def read_interval_list(interval_path: str, reference_genome: str) -> List[hl.Interval]:
    """
    Read a BED/interval file into a list of intervals.

    Returned as a list rather than applied directly so the caller can push the
    intervals INTO the data read -- see :func:`load_matrix_table`, where
    filtering after the fact costs a genome-wide split.

    :param interval_path: BED/interval file readable by
        ``hl.import_locus_intervals``.
    :param reference_genome: Reference genome for the loci, e.g. "GRCh38".
    :return: List of `hl.Interval` over the file's regions.
    """
    ht = hl.import_locus_intervals(interval_path, reference_genome=reference_genome)
    return ht.interval.collect()


def restrict_to_chrom(
    mt: Union[hl.MatrixTable, hl.Table], chrom: str, gnomad_version: str
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Restrict to a single chromosome, normalising the contig name.

    v2 is GRCh37 (contigs '1', '19', 'X'); v4 is GRCh38 ('chr1', 'chr19',
    'chrX'). Either input form is accepted for either build.

    :param mt: MatrixTable or Table keyed by locus.
    :param chrom: Chromosome, with or without the "chr" prefix.
    :param gnomad_version: "v2" or "v4"; picks the contig naming convention.
    :return: Input restricted to that chromosome.
    """
    reference_genome = mt.locus.dtype.reference_genome
    chrom = chrom.replace("chr", "")
    if gnomad_version == "v4":
        chrom = f"chr{chrom}"
    interval = hl.parse_locus_interval(chrom, reference_genome=reference_genome)
    return hl.filter_intervals(mt, [interval])


def resolve_gene_interval(
    gene: str, gene_interval: Optional[str], gnomad_version: str
) -> str:
    """
    Resolve a gene to the locus interval used to pre-filter the data.

    :param gene: Gene symbol.
    :param gene_interval: Explicit interval, which wins when supplied.
    :param gnomad_version: "v2" (GRCh37) or "v4" (GRCh38); selects which
        ``KNOWN_GENE_INTERVALS`` table to fall back to.
    :return: Locus interval string.
    :raises ValueError: If no interval is supplied and the gene is not in
        ``KNOWN_GENE_INTERVALS`` for this build.
    """
    if gene_interval:
        return gene_interval
    interval = KNOWN_GENE_INTERVALS.get(gnomad_version, {}).get(gene)
    if interval is None:
        raise ValueError(
            f"No known interval for {gene} ({gnomad_version}); pass --gene-interval explicitly, "
            f"e.g. --gene-interval chr17:43044295-43125364 (v4/GRCh38) or 17:41196312-41277500 (v2/GRCh37)."
        )
    return interval


def get_release_sites_ht(
    gnomad_version: str,
    push_down_interval: Optional[str] = None,
    chrom: Optional[str] = None,
) -> hl.Table:
    """
    Read the public release sites Table for a build, narrowed to the run scope.

    The sites Table carries ``vep`` alongside ``freq``, so one read supplies the
    PASS filters, the frequencies the AF bins key off, AND the VEP annotation --
    no separate VEP join, and no need for gnomad_qc's ``get_vep()`` on v4.

    On v4 it is read co-partitioned with the VDS the genotypes come from, so the
    downstream join is a zip of aligned partitions rather than a shuffle. v2 has
    no equivalent: the persisted intervals are the v4 VDS's, and GRCh38.

    :param gnomad_version: "v2" (GRCh37) or "v4" (GRCh38).
    :param push_down_interval: Restrict to this locus interval, if given.
    :param chrom: Restrict to this chromosome when no interval is given.
    :return: Table keyed by (locus, alleles) with ``freq``, ``filters``, ``vep``.
    """
    if gnomad_version == "v4":
        ht = public_release("exomes").ht(
            read_args={
                "_intervals": hl.eval(
                    hl.experimental.read_expression(V4_VDS_INTERVALS_PATH)
                )
            }
        )
        reference_genome = "GRCh38"
    else:
        ht = public_release_grch37("exomes").ht()
        reference_genome = "GRCh37"

    ht = ht.select("freq", "filters", "vep")
    if push_down_interval:
        return hl.filter_intervals(
            ht,
            [
                hl.parse_locus_interval(
                    push_down_interval, reference_genome=reference_genome
                )
            ],
        )
    if chrom:
        # Otherwise a whole-chromosome run joins the genome-wide sites Table.
        return restrict_to_chrom(ht, chrom, gnomad_version)
    return ht


def annotate_sites_and_filter_pass(
    mt: hl.MatrixTable,
    sites: hl.Table,
    skip_pass_filter: bool = False,
    verbose_counts: bool = False,
) -> hl.MatrixTable:
    """
    Annotate a genotype MT from the release sites Table and filter to PASS.

    An empty ``filters`` set is what PASS means in the release sites Table, so
    the filter drops everything the release flagged (AC0, RF/VQSR,
    InbreedingCoeff...). The ``is_defined(af)`` half additionally drops variants
    absent from the sites Table entirely, and applies even under
    ``skip_pass_filter`` because the AF bins need an AF.

    :param mt: Genotype MatrixTable keyed by (locus, alleles).
    :param sites: Release sites Table from :func:`get_release_sites_ht`.
    :param skip_pass_filter: Keep variants the release flagged.
    :param verbose_counts: Print kept/total site counts, at the cost of an extra
        full execution of the join.
    :return: ``mt`` with ``af``, ``filters``, ``af_bin`` and ``vep`` row
        annotations, filtered as described.
    """
    s = sites[mt.row_key]
    mt = mt.annotate_rows(
        af=s.freq[0].AF,
        filters=s.filters,
        af_bin=af_bin_expr(s.freq[0]),
        vep=s.vep,
    )
    is_pass = hl.is_defined(mt.af) & (hl.len(mt.filters) == 0)
    if verbose_counts:
        counts = mt.aggregate_rows(
            hl.struct(n_total=hl.agg.count(), n_pass=hl.agg.count_where(is_pass))
        )
        print(f"PASS filter: {counts.n_pass:,}/{counts.n_total:,} sites kept.")
    return mt.filter_rows(hl.is_defined(mt.af) if skip_pass_filter else is_pass)


def write_summary(ht: hl.Table, path: str) -> None:
    """
    Checkpoint a summary Table and export it as a TSV alongside.

    :param ht: Summary Table to write.
    :param path: Destination ``.ht`` path; the TSV replaces the extension.
    :return: None.
    """
    ht = ht.checkpoint(path, overwrite=True)
    ht.export(f"{path[:-3]}.tsv.bgz")
    print(f"Wrote {ht.count()} summary rows to {path}")


def carried_variants_ht(mt: hl.MatrixTable) -> hl.Table:
    """
    Flatten an annotated MatrixTable to one row per (carried variant, individual).

    Expects hom-ref entries to have been filtered out already, so every
    surviving entry is a raw carrier. That is what lets the caller count raw
    carriers with a plain ``count()`` rather than a ``count_where()``; only the
    adj flag has to be carried per row.

    adj is computed from GT/GQ/DP/AD. When ``AD`` is absent -- as on the v2
    hardcalls MT, which has no allele-depth field -- it falls back to True, so
    on v2 the adj counts equal the raw counts rather than being wrong.

    :param mt: MatrixTable with hom-ref entries already removed, row fields
        ``gene_symbol``, ``variant_class``, ``af`` and ``af_bin``, and entry
        field ``GT`` (plus ``GQ``/``DP``/``AD`` where available).
    :return: Table with one row per (carried variant, individual), carrying the
        sample id ``s``, the row annotations above, and ``carrier_adj``.
    """
    adj = (
        get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD)
        if "AD" in set(mt.entry)
        else hl.bool(True)
    )
    mt = mt.select_entries(carrier_adj=adj)
    # key_cols_by() BEFORE entries(): with columns keyed, entries() sorts to
    # produce row-major order, which is a full shuffle of the carrier table.
    # Unkeying first makes it a partition-local expansion. The result is unkeyed
    # either way, so this only removes the sort.
    et = mt.key_cols_by().entries().key_by()
    return et.select("s", "gene_symbol", "variant_class", "af_bin", "af", "carrier_adj")


def per_individual_counts_ht(ht: hl.Table) -> hl.Table:
    """
    Count carried variants per (gene, consequence class, AF bin, individual).

    Raw counts are a plain ``count()`` because every row of the input is already
    a raw carrier -- see :func:`carried_variants_ht`. adj counts additionally
    require the genotype to pass adj.

    This is the table that does not scale: it has one row per (individual, gene,
    class, AF bin), so a whole-chromosome run produces hundreds of millions.
    The pair grid built from the same input stays small regardless.

    :param ht: Long-form carried-variant Table from :func:`carried_variants_ht`.
    :return: Table keyed by (gene_symbol, variant_class, af_bin, s) with
        ``n_variants_raw`` and ``n_variants_adj``.
    """
    return ht.group_by(
        gene_symbol=ht.gene_symbol,
        variant_class=ht.variant_class,
        af_bin=ht.af_bin,
        s=ht.s,
    ).aggregate(
        n_variants_raw=hl.agg.count(),
        n_variants_adj=hl.agg.count_where(ht.carrier_adj),
    )


def pair_grid_ht(ht: hl.Table, individual_counts: bool = True) -> hl.Table:
    """
    Cross each individual's carried variants in a gene into the pair grid.

    For each individual, every combination of two DISTINCT variants they carry
    in a gene: k variants give k(k-1)/2 pairs, {A, B} and {B, A} are the same
    pair counted once, and a variant is never paired with itself. Within a pair,
    var1 is whichever member has the lower AF, so the grid is fully crossed on
    both axes and a single-class view is a roll-up of it.

    Counts are per individual, not per variant pair: one variant pair carried by
    500 people contributes 500. "Carried" means at least one alt allele, so a
    hom-alt variant is one variant, not two.

    ``individual_counts`` picks between two shapes of the same numbers:

    - True groups through (grid cell, individual) first, so ``n_individuals`` is
      exact. That intermediate is one row per (individual, gene, cell) and does
      not scale past a handful of genes.
    - False groups straight to the grid cell. The key space is then at most
      classes^2 x bins^2 per gene, so the shuffle is bounded by the number of
      cells rather than the number of pairs, which is what makes
      whole-chromosome scope feasible. Hail 0.2.134 has no cheap distinct-count
      aggregator, so ``n_individuals_raw``/``n_individuals_adj`` are left
      missing rather than estimated, keeping the schema stable either way.

    :param ht: Long-form carried-variant Table from :func:`carried_variants_ht`.
    :param individual_counts: Whether to compute exact ``n_individuals``.
    :return: Table keyed by (gene_symbol, class1, class2, af_bin1, af_bin2) with
        ``n_pairs_raw``/``n_pairs_adj`` and ``n_individuals_raw``/``_adj``.
    """
    carried = ht.group_by(gene_symbol=ht.gene_symbol, s=ht.s).aggregate(
        vs=hl.agg.collect(
            hl.struct(cls=ht.variant_class, bin=ht.af_bin, af=ht.af, adj=ht.carrier_adj)
        )
    )
    carried = carried.annotate(_n=hl.len(carried.vs))
    carried = carried.annotate(
        _pairs=hl.range(carried._n).flatmap(
            lambda i: hl.range(i + 1, carried._n).map(
                lambda j: hl.struct(x=carried.vs[i], y=carried.vs[j])
            )
        )
    )
    carried = carried.explode("_pairs")
    carried = carried.annotate(
        _lo=hl.if_else(
            carried._pairs.x.af <= carried._pairs.y.af,
            carried._pairs.x,
            carried._pairs.y,
        ),
        _hi=hl.if_else(
            carried._pairs.x.af <= carried._pairs.y.af,
            carried._pairs.y,
            carried._pairs.x,
        ),
    )
    both_adj = carried._lo.adj & carried._hi.adj
    grid_key = dict(
        gene_symbol=carried.gene_symbol,
        class1=carried._lo.cls,
        class2=carried._hi.cls,
        af_bin1=carried._lo.bin,
        af_bin2=carried._hi.bin,
    )
    if not individual_counts:
        return carried.group_by(**grid_key).aggregate(
            n_pairs_raw=hl.agg.count(),
            n_pairs_adj=hl.agg.count_where(both_adj),
            n_individuals_raw=hl.missing(hl.tint64),
            n_individuals_adj=hl.missing(hl.tint64),
        )

    cell = carried.group_by(s=carried.s, **grid_key).aggregate(
        n_pairs_raw=hl.agg.count(),
        n_pairs_adj=hl.agg.count_where(both_adj),
    )
    return cell.group_by(
        gene_symbol=cell.gene_symbol,
        class1=cell.class1,
        class2=cell.class2,
        af_bin1=cell.af_bin1,
        af_bin2=cell.af_bin2,
    ).aggregate(
        n_pairs_raw=hl.agg.sum(cell.n_pairs_raw),
        n_pairs_adj=hl.agg.sum(cell.n_pairs_adj),
        n_individuals_raw=hl.agg.count(),
        n_individuals_adj=hl.agg.count_where(cell.n_pairs_adj > 0),
    )


def summary_stats_ht(ht: hl.Table, n_samples: int) -> hl.Table:
    """
    Roll the per-individual counts up to per (gene, consequence class).

    Collapses the AF-bin axis first, so ``n_individuals`` counts a person once
    per gene and class no matter how many bins they carry variants in, and
    ``n_individuals_ge2`` is the compound-carrier count the in-trans work cares
    about.

    :param ht: Per-individual Table from :func:`per_individual_counts_ht`.
    :param n_samples: Cohort size, used as the denominator for the mean. Counts
        every sample, including those carrying nothing in the gene.
    :return: Table keyed by (gene_symbol, variant_class).
    """
    per = ht.group_by(
        gene_symbol=ht.gene_symbol, variant_class=ht.variant_class, s=ht.s
    ).aggregate(nr=hl.agg.sum(ht.n_variants_raw), na=hl.agg.sum(ht.n_variants_adj))
    return per.group_by(
        gene_symbol=per.gene_symbol, variant_class=per.variant_class
    ).aggregate(
        n_individuals_raw=hl.agg.count_where(per.nr > 0),
        n_individuals_adj=hl.agg.count_where(per.na > 0),
        n_individuals_ge2_adj=hl.agg.count_where(per.na >= 2),
        total_variants_raw=hl.agg.sum(per.nr),
        total_variants_adj=hl.agg.sum(per.na),
        max_variants_adj=hl.agg.max(per.na),
        mean_variants_per_individual_adj=hl.agg.sum(per.na) / n_samples,
    )


def pair_grid_summary_ht(ht: hl.Table) -> hl.Table:
    """
    Roll the pair grid up to per (gene, class1, class2), collapsing the AF axes.

    The mirror of :func:`summary_stats_ht` for the other output. Pair counts sum
    cleanly across AF bins because each (individual, pair) lands in exactly one
    cell.

    ``n_individuals`` does NOT sum: one person can carry pairs in several AF
    bins of the same class pair and would be counted once per bin. The largest
    single cell is reported instead, which is a lower bound on the distinct
    individuals, and is missing entirely when the grid was built with
    ``--no-pair-grid-individual-counts``.

    :param ht: Pair grid from :func:`pair_grid_ht`.
    :return: Table keyed by (gene_symbol, class1, class2).
    """
    return ht.group_by(
        gene_symbol=ht.gene_symbol, class1=ht.class1, class2=ht.class2
    ).aggregate(
        n_af_cells=hl.agg.count(),
        n_pairs_raw=hl.agg.sum(ht.n_pairs_raw),
        n_pairs_adj=hl.agg.sum(ht.n_pairs_adj),
        max_cell_individuals_adj=hl.agg.max(ht.n_individuals_adj),
    )


def build_variant_annotation_ht(
    push_down_interval: Optional[str] = None,
    chrom: Optional[str] = None,
) -> hl.Table:
    """
    Build a per-variant gene / class / frequency annotation Table.

    Built from the public release sites HT, which carries BOTH freq and vep, so
    the separate ``get_vep()`` join isn't needed here.

    :param push_down_interval: Restrict to this interval before annotating.
    :param chrom: Restrict to this chromosome before annotating.
    :return: Table keyed by (locus, alleles) with ``gene_symbol``,
        ``variant_class``, ``af`` and ``af_bin``, restricted to PASS variants
        that have a canonical-transcript gene assignment.
    """
    ht = public_release("exomes").ht()
    if push_down_interval:
        ht = hl.filter_intervals(
            ht,
            [hl.parse_locus_interval(push_down_interval, reference_genome="GRCh38")],
        )
    elif chrom:
        ht = restrict_to_chrom(ht, chrom, "v4")
    ht = ht.filter(hl.len(ht.filters) == 0)
    _ann = canonical_transcript_annotations_expr(ht.vep)
    ht = ht.select(
        gene_symbol=_ann.gene_symbol,
        variant_class=_ann.variant_class,
        af=ht.freq[0].AF,
        af_bin=af_bin_expr(ht.freq[0]),
    )
    return ht.filter(hl.is_defined(ht.gene_symbol) & hl.is_defined(ht.af))


def cooccurrence_pair_grid(
    counts_path: str, ann: hl.Table, gene: Optional[str] = None
) -> hl.Table:
    """Aggregate a co-occurrence genotype-counts HT onto this script's grid.

    gt_counts is [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb], with
    A/a = variant 1 and B/b = variant 2 (capital = ref). Individuals non-ref at
    BOTH variants are cells 4, 5, 7 and 8 -- the same `double_carriers`
    definition in_trans_oe uses -- and each is exactly one (individual, pair)
    co-occurrence, which is what this script's n_pairs counts. So the two
    tables are directly comparable cell by cell.
    :param counts_path: Path to a co-occurrence genotype-counts HT.
    :param ann: Per-variant annotation Table from
        :func:`build_variant_annotation_ht`.
    :param gene: Restrict to this gene symbol, if given.
    :return: Table keyed by (gene_symbol, class1, class2, af_bin1, af_bin2) with
        the pipeline's pair and phase counts, prefixed ``coocc_``.
    """
    ht = hl.read_table(counts_path)
    ht = ht.annotate(_a=ann[ht.locus1, ht.alleles1], _b=ann[ht.locus2, ht.alleles2])
    ht = ht.filter(
        hl.is_defined(ht._a)
        & hl.is_defined(ht._b)
        & (ht._a.gene_symbol == ht._b.gene_symbol)
    )
    if gene:
        ht = ht.filter(ht._a.gene_symbol == gene)

    def _dbl(c):
        return c[4] + c[5] + c[7] + c[8]

    ht = ht.annotate(
        _lo=hl.if_else(ht._a.af <= ht._b.af, ht._a, ht._b),
        _hi=hl.if_else(ht._a.af <= ht._b.af, ht._b, ht._a),
        _dr=_dbl(ht.gt_counts_raw),
        _da=_dbl(ht.gt_counts_adj),
    )
    aggs = dict(
        coocc_n_variant_pairs=hl.agg.count(),
        coocc_n_pairs_raw=hl.agg.sum(ht._dr),
        coocc_n_pairs_adj=hl.agg.sum(ht._da),
    )
    if "n_phased_cis" in set(ht.row):
        aggs["coocc_n_phased_cis"] = hl.agg.sum(ht.n_phased_cis)
        aggs["coocc_n_phased_trans"] = hl.agg.sum(ht.n_phased_trans)
    return ht.group_by(
        gene_symbol=ht._lo.gene_symbol,
        class1=ht._lo.variant_class,
        class2=ht._hi.variant_class,
        af_bin1=ht._lo.af_bin,
        af_bin2=ht._hi.af_bin,
    ).aggregate(**aggs)


def run_cooccurrence_comparison(
    output_prefix: str,
    counts_path: str,
    gene: Optional[str],
    push_down_interval: Optional[str],
) -> hl.Table:
    """
    Compare this script's pair grid against the co-occurrence pipeline's counts.

    Both are aggregated onto the same grid and outer-joined. Cells present on
    only one side are the point: they show where the two variant sets diverge.

    :param output_prefix: The ``--output-prefix`` of a previous run; its
        ``.pair_grid.ht`` sibling is read as this script's side of the comparison.
    :param counts_path: Path to a co-occurrence genotype-counts HT.
    :param gene: Restrict to this gene symbol, if given.
    :param push_down_interval: Restrict the annotation Table to this interval.
    :return: The outer-joined grid, also written alongside ``output_prefix``.
    """
    mine = hl.read_table(f"{output_prefix}.pair_grid.ht")
    coocc = cooccurrence_pair_grid(
        counts_path, build_variant_annotation_ht(push_down_interval), gene
    )
    j = mine.join(coocc, how="outer")
    j = j.annotate(
        ratio_pairs_raw=hl.or_missing(
            hl.is_defined(j.coocc_n_pairs_raw) & (j.coocc_n_pairs_raw > 0),
            hl.float64(j.n_pairs_raw) / j.coocc_n_pairs_raw,
        ),
        ratio_pairs_adj=hl.or_missing(
            hl.is_defined(j.coocc_n_pairs_adj) & (j.coocc_n_pairs_adj > 0),
            hl.float64(j.n_pairs_adj) / j.coocc_n_pairs_adj,
        ),
    )
    out = f"{output_prefix}.pair_grid_vs_cooccurrence"
    j = j.checkpoint(f"{out}.ht", overwrite=True)
    j.export(f"{out}.tsv.bgz")
    print(f"Wrote {j.count()} joined grid rows to {out}.ht")
    tot = j.aggregate(
        hl.struct(
            mine=hl.agg.sum(hl.or_else(j.n_pairs_adj, 0)),
            coocc=hl.agg.sum(hl.or_else(j.coocc_n_pairs_adj, 0)),
            only_mine=hl.agg.count_where(hl.is_missing(j.coocc_n_pairs_adj)),
            only_coocc=hl.agg.count_where(hl.is_missing(j.n_pairs_adj)),
            both=hl.agg.count_where(
                hl.is_defined(j.n_pairs_adj) & hl.is_defined(j.coocc_n_pairs_adj)
            ),
        )
    )
    print(
        f"  adj (individual,pair) co-occurrences -- this script: {tot.mine:,} | "
        f"co-occurrence pipeline: {tot.coocc:,}"
    )
    print(
        f"  grid cells: {tot.both} in both | {tot.only_mine} only here | "
        f"{tot.only_coocc} only in the pipeline"
    )
    return j


def main(args: argparse.Namespace) -> None:
    """
    Count variants per gene per individual, and cross them into a pair grid.

    Writes up to two tables: the per-individual counts at ``<prefix>.ht``, and
    the (class1, class2, af_bin1, af_bin2) pair grid at
    ``<prefix>.pair_grid.ht``. ``--compare-cooccurrence-ht`` adds a third,
    ``<prefix>.pair_grid_vs_cooccurrence.ht``, written after the grid so a
    single run produces both; ``--compare-only`` skips the pipeline and compares
    a grid an earlier run wrote, reading no genotypes.

    :param args: Parsed command-line arguments; see the parser at the bottom of
        this module for the full set and their meanings.
    :return: None. Results are written to GCS.
    """
    gene = args.gene
    gene_interval = args.gene_interval
    gnomad_version = args.gnomad_version
    chrom = args.chrom
    # Tolerate a prefix given with the extension, since that is what the old
    # --out-path took.
    output_prefix = args.output_prefix.rstrip("/")
    if output_prefix.endswith(".ht"):
        output_prefix = output_prefix[:-3]
    per_individual_path = f"{output_prefix}.ht"
    gcp_project = args.gcp_project
    tmp_dir = args.tmp_dir

    # gs://gnomad and gs://gnomad_v2 (raw genotypes, VEP annotations) are
    # requester-pays buckets -- reads fail with a 400 "Bucket is a
    # requester pays bucket but no user project provided" unless Hail is
    # told which GCP project to bill. Scoped to just the buckets this
    # script actually reads from, not blanket-enabled for all of GCS, so
    # it doesn't silently start billing reads elsewhere.
    init_kwargs = {"tmp_dir": tmp_dir} if tmp_dir else {}
    if gcp_project:
        hl.init(
            gcs_requester_pays_configuration=(
                gcp_project,
                ["gnomad", "gnomad_v2", "gnomad-tmp"],
            ),
            **init_kwargs,
        )
    else:
        hl.init(**init_kwargs)

    # Single-gene mode (e.g. for quickly validating the pipeline): resolve
    # the locus *before* loading anything. It has to be threaded into
    # load_matrix_table as push_down_interval rather than applied to mt
    # afterward -- for v4, get_gnomad_v4_vds() does an eager
    # vds.variant_data.count_cols() and a full split_multi() over whatever
    # it's handed, so restricting post-hoc still pays for a genome-wide
    # split/count first. Passed as filter_intervals, it's applied inside
    # get_gnomad_v4_vds() before either of those run.
    if args.compare_only and not args.compare_cooccurrence_ht:
        raise ValueError("--compare-only requires --compare-cooccurrence-ht.")

    push_down_interval = None
    if gene:
        push_down_interval = resolve_gene_interval(gene, gene_interval, gnomad_version)
        print(f"--gene {gene}: restricting to {push_down_interval} before loading.")

    if args.summary_stats_only:
        # Roll up a per-individual table an earlier run already wrote. Reads no
        # genotypes, so it is cheap to re-run against a finished output.
        result = hl.read_table(per_individual_path)
        if args.n_samples:
            n_samples = args.n_samples
        else:
            # Column count does not depend on which rows are read, so loading a
            # single-locus slice gives the exact cohort size -- including people
            # carrying nothing in scope, whom counting distinct `s` in the table
            # would miss -- without reading the genotypes this mode exists to
            # avoid.
            n_samples = load_matrix_table(
                args.mt_path,
                gnomad_version,
                args.release_only,
                args.high_quality_only,
                chrom,
                args.skip_v4_qc_wrapper,
                filter_intervals=[COHORT_COUNT_INTERVAL[gnomad_version]],
            ).count_cols()
        print(f"Cohort size: {n_samples:,} samples")
        write_summary(
            summary_stats_ht(result, n_samples), f"{output_prefix}.summary_stats.ht"
        )
        write_summary(
            pair_grid_summary_ht(hl.read_table(f"{output_prefix}.pair_grid.ht")),
            f"{output_prefix}.pair_grid_summary.ht",
        )
        return

    if args.compare_only:
        # Skip straight to the comparison against a pair grid a previous run
        # already wrote. Everything above is cheap setup; everything below reads
        # genotypes, so this returns at that boundary and never does.
        run_cooccurrence_comparison(
            output_prefix, args.compare_cooccurrence_ht, gene, push_down_interval
        )
        return

    # Genotypes only -- no VEP yet. v2 comes back already PASS-filtered
    # (a `filters`-field lookup, not VEP -- see load_matrix_table).
    # Push every interval restriction INTO the read. Filtering after the load
    # is what the loader's docstring warns against: get_gnomad_v4_vds() does an
    # eager count_cols() and a full split_multi first, over the whole genome if
    # nothing has narrowed it yet. --gene was already pushed down; --interval-path
    # was not, so a BED-scoped run used to split the genome before restricting.
    load_intervals = [push_down_interval] if push_down_interval else []
    if args.interval_path:
        load_intervals += read_interval_list(
            args.interval_path, "GRCh38" if gnomad_version == "v4" else "GRCh37"
        )

    n_samples = 0

    mt = load_matrix_table(
        args.mt_path,
        gnomad_version,
        args.release_only,
        args.high_quality_only,
        chrom,
        args.skip_v4_qc_wrapper,
        filter_intervals=load_intervals or None,
    )

    # Keep the fields adj needs; the original select_entries("GT") dropped them.
    _entry = set(mt.entry)
    mt = mt.select_entries(*[f for f in ("GT", "GQ", "DP", "AD") if f in _entry])

    if args.summary_stats:
        # Cohort size for the per-individual mean; taken before any entry
        # filtering so it counts everyone, not just carriers.
        n_samples = mt.count_cols()
        print(f"Cohort size: {n_samples:,} samples")

    # Drop hom-ref entries and then rows with no carriers left, BEFORE the VEP
    # join and the adj computation, so only rows guaranteed to reach the output
    # pay for a transcript_consequences array. It matters most on v2, whose
    # hardcalls MT is DENSE -- an entry per (variant, sample), hom-ref included
    # -- whereas v4's variant_data is already sparse, so there it is close to a
    # no-op.
    mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = mt.filter_rows(hl.agg.count() > 0)

    # PASS + AF + VEP from the release sites Table, for BOTH builds and in one
    # join. Previously only v2 was PASS-filtered, at load time, so v4 counted
    # non-PASS variants.
    sites = get_release_sites_ht(gnomad_version, push_down_interval, chrom)
    mt = annotate_sites_and_filter_pass(
        mt,
        sites,
        skip_pass_filter=args.skip_filter_pass,
        verbose_counts=args.verbose_counts,
    )

    # gene_symbol and variant_class both come from the SAME canonical
    # transcript (see canonical_transcript_annotations_expr) -- a single
    # scalar gene per variant, and a variant_class that's guaranteed
    # consistent with the transcript actually driving that gene call.
    mt = mt.annotate_rows(**canonical_transcript_annotations_expr(mt.vep))

    # A variant with no canonical protein-coding Ensembl transcript has no gene
    # to attribute it to -- drop it. No explode: gene_symbol is a single scalar,
    # so each variant contributes to at most one output row per sample.
    mt = mt.filter_rows(hl.is_defined(mt.gene_symbol))

    if gene:
        # Interval overlap alone can pull in a neighboring gene's
        # canonical-transcript variants too -- restrict precisely to the
        # requested gene now that gene_symbol is actually assigned.
        mt = mt.filter_rows(mt.gene_symbol == gene)

    mt = mt.select_rows("gene_symbol", "variant_class", "af", "af_bin")

    et = carried_variants_ht(mt)
    et = et.checkpoint(hl.utils.new_temp_file("carried", "ht"))

    # --- Output 1: per (gene, class, af_bin, individual) ------------------
    result = None
    if args.no_individual_variant_counts:
        print(
            "--no-individual-variant-counts: skipping the per-individual "
            "variant-count table."
        )
    else:
        # checkpoint is write-then-read-back, so the count and preview below
        # come off the written Table: a written Hail Table stores each
        # partition's row count in its metadata, making count() a metadata
        # lookup rather than a rescan of the whole upstream pipeline. (The
        # original called show() before write(), executing that pipeline twice
        # for the same rows.)
        result = per_individual_counts_ht(et).checkpoint(
            per_individual_path, overwrite=True
        )
        n_result_rows = result.count()

        if n_result_rows == 0:
            # An empty output is indistinguishable from a silent bug unless
            # it's called out explicitly -- e.g. a gene with no qualifying
            # canonical-transcript variants in this cohort/version, or an
            # overly narrow interval/filter combination.
            scope = (
                f"gene {gene}"
                if gene
                else f"chrom {chrom}"
                if chrom
                else "the requested scope"
            )
            print(
                f"NOTE: 0 rows written to {per_individual_path} -- no (gene, variant_class, "
                f"sample) combinations found for {scope} ({gnomad_version}). This is a real "
                "result (no qualifying variants survived filtering), not a write failure "
                "-- verify the gene/interval and QC flags (--release-only, "
                "--high-quality-only, --skip-filter-pass) are what you intended."
            )
        else:
            print(
                f"Wrote {n_result_rows} per-gene-per-class-per-individual variant count "
                f"rows to {per_individual_path}"
            )
            if gene:
                print(f"\n(variant_class, sample) -> counts for {gene}:")
                result.show(25)

    # --- Output 2: the class x class / AF x AF pair grid ------------------
    grid = pair_grid_ht(et, individual_counts=not args.no_pair_grid_individual_counts)

    grid_path = f"{output_prefix}.pair_grid.ht"
    # The grid's group_by shuffle fails on chromosome-scale input with the
    # default shuffler; the new one handles it. Scoped to this write and unset
    # again, since it is a compilation-time flag and everything else is fine
    # without it.
    hl._set_flags(use_new_shuffle="1")
    grid = grid.checkpoint(grid_path, overwrite=True)
    hl._set_flags(use_new_shuffle=None)
    grid.export(f"{output_prefix}.pair_grid.tsv.bgz")
    print(f"Wrote {grid.count()} pair-grid rows to {grid_path}")

    if args.summary_stats:
        hl._set_flags(use_new_shuffle="1")
        # Rolled up here rather than by a follow-up script, so the numbers that
        # get reported are reproducible from the same command. Both outputs get
        # one; the per-individual roll-up is skipped when its table was not
        # written.
        if result is not None:
            write_summary(
                summary_stats_ht(result, n_samples), f"{output_prefix}.summary_stats.ht"
            )
        write_summary(
            pair_grid_summary_ht(grid), f"{output_prefix}.pair_grid_summary.ht"
        )
        hl._set_flags(use_new_shuffle=None)

    if args.compare_cooccurrence_ht:
        # The grid was just written above, so this compares against fresh output
        # rather than needing a previous run.
        run_cooccurrence_comparison(
            output_prefix, args.compare_cooccurrence_ht, gene, push_down_interval
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    io_group = parser.add_argument_group(
        "Input and output",
        "Where the data comes from and where results are written.",
    )
    io_group.add_argument(
        "--gnomad-version",
        choices=["v2", "v4"],
        default="v4",
        help="Controls loading path/method and canonical-field typing (int in v2, bool in v4)",
    )
    io_group.add_argument(
        "--mt-path",
        default=None,
        help="Explicit MatrixTable/VDS path. If omitted: v4 loads via "
        "get_gnomad_v4_vds() (recommended, applies QC filtering); "
        "v2 falls back to the exomes hardcalls MT path.",
    )
    io_group.add_argument(
        "--output-prefix",
        required=True,
        help="Path prefix every output is derived from, e.g. "
        "gs://bucket/chr20_variants_per_gene. Writes <prefix>.ht (per-individual "
        "counts), <prefix>.pair_grid.ht, and, with --summary-stats, "
        "<prefix>.summary_stats.ht and <prefix>.pair_grid_summary.ht. A trailing "
        '".ht" is accepted and stripped. The read-only modes (--compare-only, '
        "--summary-stats-only) take the same prefix and read what an earlier run "
        "wrote there.",
    )
    io_group.add_argument(
        "--tmp-dir",
        default=None,
        help="GCS scratch dir for Hail (e.g. gs://your-tmp/). Strongly recommended on "
        "Dataproc: without it checkpoints land on the tiny HDFS /tmp and fail with "
        "'minReplication'/'Premature end of file'.",
    )
    io_group.add_argument(
        "--gcp-project",
        default=None,
        help="GCP project ID to bill for requester-pays reads of gs://gnomad (the v4 "
        "VDS and its persisted partition intervals) and gs://gnomad_v2. Required for "
        "v4 genotype access; the release sites Tables are public and need no billing "
        "project.",
    )

    scope = parser.add_argument_group(
        "Scope",
        "Which part of the genome to run on. Every one of these is pushed into "
        "the data read rather than applied afterwards.",
    )
    scope.add_argument(
        "--gene",
        default=None,
        help="Restrict to a single gene (by canonical-transcript gene_symbol), e.g. BRCA1 -- "
        "for quickly validating the pipeline on a cheap, eyeballable slice instead of a "
        "full chromosome/genome run. Restricts to the gene's locus interval immediately "
        "(before VEP, before entry filtering), then additionally filters to "
        "gene_symbol == --gene after canonical-transcript assignment (interval overlap "
        "alone can pull in a neighboring gene too). Requires --gene-interval unless the "
        "gene is in KNOWN_GENE_INTERVALS, which covers BRCA1, PCSK9 and the five "
        "co-occurrence test genes (AHNAK2, ANO5, CAPN3, DYSF, SGCA). Prints a "
        "preview of the result before writing.",
    )
    scope.add_argument(
        "--gene-interval",
        default=None,
        help="Locus interval for --gene, e.g. chr17:43044295-43125364 (v4/GRCh38) or "
        "17:41196312-41277500 (v2/GRCh37). Required for --gene unless the gene is in "
        "KNOWN_GENE_INTERVALS.",
    )
    scope.add_argument(
        "--interval-path",
        default=None,
        help="Optional BED/interval file to restrict to specific genes/regions before processing",
    )
    scope.add_argument(
        "--chrom",
        default=None,
        help="Restrict to a single chromosome, e.g. '19' or 'chr19' (either "
        "form works for both versions -- normalized internally: v2 is "
        "GRCh37/'19', v4 is GRCh38/'chr19'). Applied before splitting/"
        "reading full data, for efficiency.",
    )

    filters = parser.add_argument_group(
        "Sample and variant filtering",
        "Which samples and sites are counted.",
    )
    filters.add_argument(
        "--release-only",
        action="store_true",
        help="(v4 only, via get_gnomad_v4_vds) Restrict to release samples only",
    )
    filters.add_argument(
        "--high-quality-only",
        action="store_true",
        help="(v4 only, via get_gnomad_v4_vds) Restrict to high-quality samples only",
    )
    filters.add_argument(
        "--skip-filter-pass",
        action="store_true",
        help="Skip the PASS-site filter, keeping variants the release flagged (AC0, "
        "RF/VQSR, InbreedingCoeff...). Applies to both builds. Variants absent from "
        "the release sites Table are still dropped, since their AF is needed.",
    )
    filters.add_argument(
        "--skip-v4-qc-wrapper",
        action="store_true",
        help="(v4 only, ignored if --mt-path is set) Bypass get_gnomad_v4_vds() and read "
        "the raw VDS directly. This skips every filter that function applies: the "
        "chr19:5787204 multiallelic-site drop, duplicate/withdrawn UKB sample removal, "
        "and hard-filtered sample removal. It ALSO silently disables --release-only and "
        "--high-quality-only, which are only passed to get_gnomad_v4_vds() -- so the raw "
        "path counts every sample in the VDS, not the 730,947 release samples. Fine for "
        "a quick --gene smoke test; do not use it for numbers you intend to report.",
    )

    outputs = parser.add_argument_group(
        "What to compute",
        "The per-individual table and the pair grid are always written unless "
        "disabled here; the summary roll-up is opt-in.",
    )
    outputs.add_argument(
        "--no-individual-variant-counts",
        action="store_true",
        help="Skip the per-individual variant-count TABLE entirely, writing only the "
        "pair grid. That table is one row per (individual, gene, class, AF bin) -- "
        "333M rows for chr20 -- and is what blocks large scopes; the grid is a few "
        "hundred rows per gene either way. Distinct from "
        "--no-pair-grid-individual-counts, which drops two COLUMNS of the grid.",
    )
    outputs.add_argument(
        "--no-pair-grid-individual-counts",
        action="store_true",
        help="Leave the pair grid's n_individuals_raw/n_individuals_adj COLUMNS "
        "missing. Computing them exactly needs a group_by keyed on the individual, "
        "one row per (individual, gene, cell), which does not scale past a few "
        "genes; without them the grid groups straight to the cell. The grid is "
        "still written either way. Distinct from --no-individual-variant-counts, "
        "which drops a whole output TABLE.",
    )
    outputs.add_argument(
        "--summary-stats",
        action="store_true",
        help="Also write roll-ups of both outputs: <prefix>.summary_stats.ht "
        "(per-individual counts by gene and consequence class, with the number "
        "of individuals carrying at least one and at least two, totals and the "
        "mean per individual) and <prefix>.pair_grid_summary.ht (the grid "
        "collapsed over the AF axes, by gene and class pair). Each also gets a "
        ".tsv.bgz. The per-individual roll-up is skipped with "
        "--no-individual-variant-counts, which skips the table it rolls up.",
    )
    outputs.add_argument(
        "--n-samples",
        type=int,
        default=None,
        help="Override the cohort size used as the denominator for the "
        "per-individual mean. Rarely needed: a full run takes it from the loaded "
        "MatrixTable, and --summary-stats-only counts columns from a one-locus "
        "slice, which is exact and cheap.",
    )

    modes = parser.add_argument_group(
        "Alternate modes",
        "Each of these skips the genotype pipeline and works from output an earlier "
        "run already wrote.",
    )
    modes.add_argument(
        "--compare-cooccurrence-ht",
        default=None,
        help="Path to a co-occurrence pipeline genotype-counts HT (e.g. "
        "exomes.variant_pairs.genotype_counts.<postfix>.ht). Aggregates that "
        "table onto this script's (class1,class2,af_bin1,af_bin2) grid and "
        "outer-joins it against this run's pair grid, writing "
        "<prefix>.pair_grid_vs_cooccurrence.ht. Runs after the grid is "
        "written, so a normal run produces both; add --compare-only to skip the "
        "pipeline and compare a grid an earlier run already wrote.",
    )
    modes.add_argument(
        "--compare-only",
        action="store_true",
        help="Skip the pipeline and run only the comparison, against the pair grid "
        "an earlier --output-prefix run wrote. Requires --compare-cooccurrence-ht. "
        "Reads no genotypes, so it is cheap to re-run.",
    )
    modes.add_argument(
        "--summary-stats-only",
        action="store_true",
        help="Skip the pipeline and roll up the per-individual table and pair grid "
        "an earlier --output-prefix run wrote, producing the same two summaries "
        "as --summary-stats. Reads no genotypes; the cohort size comes from a "
        "one-locus column count unless --n-samples overrides it.",
    )

    debug = parser.add_argument_group(
        "Diagnostics",
        "Extra output, at the cost of extra work.",
    )
    debug.add_argument(
        "--verbose-counts",
        action="store_true",
        help="Print kept/total site counts for the PASS filter. Off by default "
        "because it forces an extra full execution of the join (Hail is lazy) -- "
        "only enable for debugging/small runs.",
    )

    main(parser.parse_args())

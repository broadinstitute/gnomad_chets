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
    withdrawn UKB samples, and removes hard-filtered samples.
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
    AC0/RF/InbreedingCoeff/etc.) -- this IS applied at load time (it's a
    site-validity check, not VEP, and only selects one small field), via
    filter_to_pass_v2() (--skip-filter-pass to disable).

Usage:
    hailctl dataproc submit <cluster> variants_per_gene_per_individual.py \
        --out-path gs://path/to/output_table.ht \
        --gnomad-version v4 \
        [--mt-path gs://path/to/gnomad_data]   # overrides the version default \
        [--interval-path gs://path/to/genes.bed] \
        [--skip-filter-pass]                    # v2 only
"""

import argparse
import contextlib
import hail as hl

CHR19_MULTIALLELIC_DROP_INTERVAL = "chr19:5787204-5787205"


@contextlib.contextmanager
def suppress_chr19_multiallelic_drop():
    """Surgical alternative to bypassing get_gnomad_v4_vds() entirely
    (--skip-v4-qc-wrapper): monkeypatch hl.vds.filter_intervals so ONLY
    the exact (interval, keep=False) call get_gnomad_v4_vds() hardcodes
    to drop the chr19:5787204 problematic multiallelic site becomes a
    no-op, while every other filter_intervals call (that function's own
    chrom filtering, anything else) passes through unchanged -- so UKB
    dedup and hard-filtered sample removal still happen normally.

    This is inherently fragile: it pattern-matches the exact interval
    gnomad_qc hardcodes today. If gnomad_qc ever changes that internal
    call, the patch just stops matching and silently does nothing (fails
    safe -- the drop still happens -- rather than silently breaking
    something else), but it's worth re-checking after any gnomad_qc
    upgrade. If you're not actually processing chr19, this entire patch
    is moot anyway: get_gnomad_v4_vds() applies its chrom filter BEFORE
    this drop, so on any other chromosome there's no chr19 data left for
    the drop to act on regardless.
    """
    orig_filter_intervals = hl.vds.filter_intervals
    problem_interval = hl.parse_locus_interval(
        CHR19_MULTIALLELIC_DROP_INTERVAL, reference_genome="GRCh38"
    )

    def patched(vds, intervals, keep=True, **kwargs):
        is_the_hardcoded_drop_call = (
            keep is False
            and len(intervals) == 1
            and hl.eval(intervals[0] == problem_interval)
        )
        if is_the_hardcoded_drop_call:
            print(
                f"suppress_chr19_multiallelic_drop: intercepted and skipped "
                f"get_gnomad_v4_vds()'s internal drop of {CHR19_MULTIALLELIC_DROP_INTERVAL}. "
                "This site can blow up memory/compute in split_multi if you're actually "
                "processing chr19 -- only suppress this if you specifically need that site."
            )
            return vds
        return orig_filter_intervals(vds, intervals, keep=keep, **kwargs)

    hl.vds.filter_intervals = patched
    try:
        yield
    finally:
        hl.vds.filter_intervals = orig_filter_intervals

# A couple of well-known genes for quick --gene validation runs, both
# builds. Add more as needed, or just pass --gene-interval directly.
KNOWN_GENE_INTERVALS = {
    "v2": {  # GRCh37
        "BRCA1": "17:41196312-41277500",
        "PCSK9": "1:55505221-55530525",
    },
    "v4": {  # GRCh38
        "BRCA1": "chr17:43044295-43125364",
        "PCSK9": "chr1:55039475-55064852",
    },
}

# --- consequence_terms -> variant class (classified off the CANONICAL
# transcript's own terms, not the variant-wide most_severe_consequence
# field -- see canonical_transcript_annotations_expr below for why) ------

_LOF_TERMS = hl.set([
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
])
_MISSENSE_TERMS = hl.set([
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant",
])
_SYNONYMOUS_TERMS = hl.set([
    "synonymous_variant", "stop_retained_variant", "start_retained_variant",
])
_NONCODING_TERMS = hl.set([
    "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
    "upstream_gene_variant", "downstream_gene_variant",
    "non_coding_transcript_exon_variant", "non_coding_transcript_variant",
    "intergenic_variant", "regulatory_region_variant",
    "TF_binding_site_variant", "mature_miRNA_variant",
    "splice_region_variant", "coding_sequence_variant",
    "incomplete_terminal_codon_variant",
])


def classify_terms_expr(consequence_terms):
    """consequence_terms is an array (a single transcript can have more
    than one term, e.g. ['missense_variant', 'splice_region_variant']) --
    pick the worst matching category, in lof > missense > synonymous >
    noncoding > other priority order."""
    terms = hl.set(consequence_terms)
    return (
        hl.case()
        .when(terms.intersection(_LOF_TERMS).size() > 0, "lof")
        .when(terms.intersection(_MISSENSE_TERMS).size() > 0, "missense")
        .when(terms.intersection(_SYNONYMOUS_TERMS).size() > 0, "synonymous")
        .when(terms.intersection(_NONCODING_TERMS).size() > 0, "noncoding")
        .default("other")
    )


# --- Canonical-transcript handling: int (1/0) in BOTH v2 and v4. Confirmed
# against a live v4 run (TypeError: 'filter': expected bool, found int32 --
# the previous version/bool assumption was wrong) and against this repo's
# own v4 code in gnomad_chets/v4/rf_ptrans_features.py, which already uses
# `tc.canonical == 1` for v4. No version branching needed. ------------------

def canonical_filter_expr(tc):
    return tc.canonical == 1


def canonical_transcript_annotations_expr(vep_struct, gnomad_version: str):
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
    """
    canonical_tcs = vep_struct.transcript_consequences.filter(canonical_filter_expr)
    primary_tc = hl.or_missing(hl.len(canonical_tcs) > 0, canonical_tcs[0])
    return hl.struct(
        gene_symbol=primary_tc.gene_symbol,
        variant_class=hl.or_missing(
            hl.is_defined(primary_tc), classify_terms_expr(primary_tc.consequence_terms)
        ),
    )


# --- Loading -------------------------------------------------------------

V2_EXOMES_HARDCALLS_MT_PATH = "gs://gnomad_v2/hardcalls/hail-0.2/mt/exomes/gnomad.exomes.mt"
V4_RAW_EXOMES_VDS_PATH = "gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds"


def load_matrix_table(
    mt_path: str,
    gnomad_version: str,
    release_only: bool,
    high_quality_only: bool,
    skip_filter_pass: bool,
    chrom: str = None,
    verbose_counts: bool = False,
    skip_v4_qc_wrapper: bool = False,
    push_down_interval: str = None,
) -> hl.MatrixTable:
    """Returns genotypes only -- NO vep annotation yet. v2 gets its
    site-validity PASS filter applied here (that's a `filters`-field
    lookup, not VEP), but the VEP join itself is deferred to
    join_vep_late(), called from main() only after entry/row filtering
    has already dropped everything it can. VEP annotation is the most
    expensive join in this pipeline (full transcript_consequences arrays
    per row), so it should touch the smallest possible set of rows --
    joining it here, before non-ref filtering, would mean annotating rows
    that might get dropped moments later for having zero carriers.

    push_down_interval (e.g. from --gene) MUST be applied INSIDE this
    function, not by the caller after it returns. For v4 specifically,
    get_gnomad_v4_vds() runs an eager vds.variant_data.count_cols() (a
    real Spark job, not lazy) partway through, and a full split_multi at
    the end -- both over the WHOLE genome-wide callset if no interval has
    been applied yet. Restricting afterward (e.g. in main()) is too late:
    by then get_gnomad_v4_vds() has already scanned everything. Passing
    the interval into get_gnomad_v4_vds()'s own filter_intervals param
    applies it near the top of that function, before both of those.
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
                    "bypassing get_gnomad_v4_vds() entirely -- no chr19:5787204 drop, "
                    "no duplicate/withdrawn UKB removal, no hard-filtered sample removal."
                )
            vds = hl.vds.read_vds(read_path)
            if norm_chrom:
                vds = hl.vds.filter_chromosomes(vds, keep=[norm_chrom])
            if push_down_interval:
                # Before split_multi, same reasoning as get_gnomad_v4_vds's
                # own filter_intervals-before-split ordering below.
                reference_genome = vds.reference_data.locus.dtype.reference_genome
                interval = hl.parse_locus_interval(push_down_interval, reference_genome=reference_genome)
                vds = hl.vds.filter_intervals(vds, [interval])
            # A freshly-read VDS is unsplit (LGT/LA, not global GT), so
            # split it here.
            vds = hl.vds.split_multi(vds)
        else:
            from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

            print(
                "Loading gnomAD v4 exomes via get_gnomad_v4_vds() "
                f"(release_only={release_only}, high_quality_only={high_quality_only}, "
                f"push_down_interval={push_down_interval}); "
                "this removes duplicate/withdrawn UKB samples and hard-filtered samples "
                f"(its internal {CHR19_MULTIALLELIC_DROP_INTERVAL} multiallelic-site drop "
                "is suppressed -- see suppress_chr19_multiallelic_drop)."
            )
            with suppress_chr19_multiallelic_drop():
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
                    filter_intervals=[push_down_interval] if push_down_interval else None,
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
        if push_down_interval:
            # hl.read_matrix_table is lazy (unlike get_gnomad_v4_vds,
            # nothing eager happens above), but restricting here still
            # keeps the PASS-filter join below from dealing with more
            # rows than necessary.
            reference_genome = mt.locus.dtype.reference_genome
            interval = hl.parse_locus_interval(push_down_interval, reference_genome=reference_genome)
            mt = hl.filter_intervals(mt, [interval])
        if not skip_filter_pass:
            mt = filter_to_pass_v2(mt, chrom, gnomad_version, verbose_counts, push_down_interval)
        return mt


def filter_to_pass_v2(
    mt: hl.MatrixTable, chrom: str, gnomad_version: str, verbose_counts: bool, push_down_interval: str = None,
) -> hl.MatrixTable:
    """v2's hardcalls MT has no site-quality annotation of its own; the
    RF-based PASS/AC0/RF/InbreedingCoeff `filters` field lives on the
    public release sites Table. This is a site-validity check, not VEP,
    so it's fine (and cheap, since only one small field is selected) to
    apply early -- unlike VEP, it doesn't carry per-row
    transcript_consequences arrays."""
    from gnomad_qc.v2.resources.basics import get_gnomad_public_data

    release_ht = get_gnomad_public_data("exomes", split=True).select("filters")
    if chrom:
        release_ht = restrict_to_chrom(release_ht, chrom, gnomad_version)
    if push_down_interval:
        reference_genome = release_ht.locus.dtype.reference_genome
        interval = hl.parse_locus_interval(push_down_interval, reference_genome=reference_genome)
        release_ht = hl.filter_intervals(release_ht, [interval])

    release_filters = release_ht[mt.row_key].filters

    if verbose_counts:
        counts = mt.aggregate_rows(
            hl.struct(
                n_total=hl.agg.count(),
                n_pass=hl.agg.count_where(hl.is_defined(release_filters) & (hl.len(release_filters) == 0)),
            )
        )
        print(f"v2 PASS-filter: kept {counts.n_pass}/{counts.n_total} sites.")

    return mt.filter_rows(hl.is_defined(release_filters) & (hl.len(release_filters) == 0))


def join_vep_late(
    mt: hl.MatrixTable, gnomad_version: str, chrom: str, push_down_interval: str = None
) -> hl.MatrixTable:
    """The one and only VEP join, called as late as possible (from
    main(), after entry-level non-ref filtering and after dropping rows
    with zero remaining carriers) -- so it only ever touches rows that
    are guaranteed to make it into the output. VEP lives in a completely
    separate Table from genotypes for BOTH versions (v2: the public
    release sites Table; v4: gnomad_qc.v4.resources.annotations.get_vep(),
    confirmed to be a standalone VersionedTableResource, not part of the
    raw genotype VDS).

    Restricting vep_ht is just as important here as restricting mt was in
    load_matrix_table: `vep_ht[mt.row_key]` is a join against whichever
    Table vep_ht is, and an unrestricted vep_ht is the FULL genome-wide
    VEP Table (every exome variant, full transcript_consequences arrays)
    regardless of how small mt already is. chrom alone doesn't cover
    --gene mode (that sets push_down_interval, not chrom), so both are
    applied here.
    """
    if gnomad_version == "v4":
        from gnomad_qc.v4.resources.annotations import get_vep

        vep_ht = get_vep(data_type="exomes").ht().select("vep")
    else:
        from gnomad_qc.v2.resources.basics import get_gnomad_public_data

        vep_ht = get_gnomad_public_data("exomes", split=True).select("vep")

    if chrom:
        vep_ht = restrict_to_chrom(vep_ht, chrom, gnomad_version)

    if push_down_interval:
        reference_genome = vep_ht.locus.dtype.reference_genome
        interval = hl.parse_locus_interval(push_down_interval, reference_genome=reference_genome)
        vep_ht = hl.filter_intervals(vep_ht, [interval])

    return mt.annotate_rows(vep=vep_ht[mt.row_key].vep)


def restrict_to_intervals(mt, interval_path: str):
    intervals = hl.import_locus_intervals(
        interval_path, reference_genome=mt.locus.dtype.reference_genome
    )
    return hl.filter_intervals(mt, intervals.interval.collect())


def restrict_to_chrom(mt, chrom: str, gnomad_version: str):
    """Works on a Table or MatrixTable keyed by locus. v2 is GRCh37
    (contigs named '1', '19', 'X', ...); v4 is GRCh38 (contigs named
    'chr1', 'chr19', 'chrX', ...). Normalize whatever the user passes
    (e.g. '19' or 'chr19') to the right convention."""
    reference_genome = mt.locus.dtype.reference_genome
    chrom = chrom.replace("chr", "")
    if gnomad_version == "v4":
        chrom = f"chr{chrom}"
    interval = hl.parse_locus_interval(chrom, reference_genome=reference_genome)
    return hl.filter_intervals(mt, [interval])


def resolve_gene_interval(gene: str, gene_interval: str, gnomad_version: str) -> str:
    if gene_interval:
        return gene_interval
    interval = KNOWN_GENE_INTERVALS.get(gnomad_version, {}).get(gene)
    if interval is None:
        raise ValueError(
            f"No known interval for {gene} ({gnomad_version}); pass --gene-interval explicitly, "
            f"e.g. --gene-interval chr17:43044295-43125364 (v4/GRCh38) or 17:41196312-41277500 (v2/GRCh37)."
        )
    return interval


def main(
    mt_path, out_path, gnomad_version, interval_path, chrom,
    release_only, high_quality_only, skip_filter_pass, verbose_counts,
    gene, gene_interval, gcp_project, skip_v4_qc_wrapper,
):
    # gs://gnomad and gs://gnomad_v2 (raw genotypes, VEP annotations) are
    # requester-pays buckets -- reads fail with a 400 "Bucket is a
    # requester pays bucket but no user project provided" unless Hail is
    # told which GCP project to bill. Scoped to just the buckets this
    # script actually reads from, not blanket-enabled for all of GCS, so
    # it doesn't silently start billing reads elsewhere.
    if gcp_project:
        hl.init(gcs_requester_pays_configuration=(gcp_project, ["gnomad", "gnomad_v2", "gnomad-tmp"]))
    else:
        hl.init()

    # Single-gene mode (e.g. for quickly validating the pipeline): resolve
    # the locus *before* loading anything. It has to be threaded into
    # load_matrix_table as push_down_interval rather than applied to mt
    # afterward -- for v4, get_gnomad_v4_vds() does an eager
    # vds.variant_data.count_cols() and a full split_multi() over whatever
    # it's handed, so restricting post-hoc still pays for a genome-wide
    # split/count first. Passed as filter_intervals, it's applied inside
    # get_gnomad_v4_vds() before either of those run.
    push_down_interval = None
    if gene:
        push_down_interval = resolve_gene_interval(gene, gene_interval, gnomad_version)
        print(f"--gene {gene}: restricting to {push_down_interval} before loading.")

    # Genotypes only -- no VEP yet. v2 comes back already PASS-filtered
    # (a `filters`-field lookup, not VEP -- see load_matrix_table).
    mt = load_matrix_table(
        mt_path, gnomad_version, release_only, high_quality_only, skip_filter_pass, chrom, verbose_counts,
        skip_v4_qc_wrapper, push_down_interval=push_down_interval,
    )

    mt = mt.select_entries("GT")

    if interval_path:
        mt = restrict_to_intervals(mt, interval_path)

    # Filter to non-ref entries and drop now-empty rows BEFORE touching
    # VEP at all -- this is the whole point of deferring the join: only
    # variants that are guaranteed to appear in the output ever get a
    # transcript_consequences array pulled in.
    mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = mt.filter_rows(hl.agg.count() > 0)

    # The one and only VEP join, as late as possible.
    mt = join_vep_late(mt, gnomad_version, chrom, push_down_interval)

    # gene_symbol and variant_class both come from the SAME canonical
    # transcript (see canonical_transcript_annotations_expr) -- a single
    # scalar gene per variant, and a variant_class that's guaranteed
    # consistent with the transcript actually driving that gene call.
    mt = mt.annotate_rows(**canonical_transcript_annotations_expr(mt.vep, gnomad_version))

    # A variant with no canonical-transcript hit (e.g. purely intergenic)
    # has no gene to attribute it to -- drop it. No explode: gene_symbol
    # is a single scalar, so each variant already contributes to at most
    # one output row per sample.
    mt = mt.filter_rows(hl.is_defined(mt.gene_symbol))

    if gene:
        # Interval overlap alone can pull in a neighboring gene's
        # canonical-transcript variants too -- restrict precisely to the
        # requested gene now that gene_symbol is actually assigned.
        mt = mt.filter_rows(mt.gene_symbol == gene)

    mt = mt.select_rows("variant_class", "gene_symbol")

    # Aggregate directly on the matrix table's column (sample) axis via
    # hl.agg.group_by, instead of materializing mt.entries() first. This
    # avoids the expensive global (row_key, col_key) sort entries()
    # triggers (the sort Hail explicitly warns about) -- annotate_cols
    # scans down each sample's column once and groups by (gene_symbol,
    # variant_class) as it goes, no shuffle needed.
    #
    # CORRECTNESS-CRITICAL: hl.agg.group_by's key (gene_symbol,
    # variant_class) is a ROW-level expression, constant across every
    # sample -- so without restricting the aggregation's scope, it groups
    # over EVERY row with that gene/class in the WHOLE dataset, not just
    # rows where THIS sample has a defined/non-ref GT. filter_entries()
    # earlier only masks the entry value to missing; it does NOT remove
    # that row from other columns' aggregation scope. Wrapping the whole
    # group_by in hl.agg.filter(is_defined(GT), ...) is what actually (a)
    # scopes hl.agg.count() to just this sample's real carrier rows
    # (otherwise it counts the cohort-wide total for that gene/class,
    # identically for every sample) and (b) keeps the resulting dict
    # sparse -- only (gene, class) pairs this sample actually carries a
    # variant in appear at all, rather than every combination in the
    # entire dataset with most counts sitting at 0.
    # Dict value is a bare count (hl.agg.count()), not a struct -- only
    # n_variants is needed, so there's no reason to pay for constructing
    # and later indexing into a one-field struct per group.
    mt = mt.annotate_cols(
        _gene_class_counts=hl.agg.filter(
            hl.is_defined(mt.GT),
            hl.agg.group_by(
                hl.tuple([mt.gene_symbol, mt.variant_class]),
                hl.agg.count(),
            ),
        )
    )

    # Un-nest the per-sample dict into the long-format (gene_symbol,
    # variant_class, s) -> n_variants Table. This explode happens on
    # cols() (one row per sample) AFTER aggregation has already collapsed
    # variant-level cardinality down to gene/class-level -- a far smaller
    # expansion than exploding at the variant level would have been.
    # DictExpression.items() -> array of (key, value) tuples, indexable
    # via [0]/[1] -- used deliberately over hl.array(dict_expr) (which
    # instead produces struct{key, value} elements, a different access
    # pattern) to keep this unambiguous.
    cols_ht = mt.cols()
    cols_ht = cols_ht.annotate(_kv=cols_ht._gene_class_counts.items())
    cols_ht = cols_ht.explode("_kv")
    result = cols_ht.select(
        gene_symbol=cols_ht._kv[0][0],
        variant_class=cols_ht._kv[0][1],
        n_variants=cols_ht._kv[1],
    )
    result = result.key_by("gene_symbol", "variant_class", "s")

    # Write FIRST, then read the written table back for the preview/count.
    # (Previously the --gene preview called result.show() before
    # result.write() -- two separate executions of the whole upstream
    # pipeline, VEP join and all, for exactly the same rows. Writing once
    # and reading the materialized .ht back is free by comparison: a
    # written Hail Table stores each partition's row count in its
    # metadata, so .count() on it is a metadata lookup, not a rescan.)
    result.write(out_path, overwrite=True)
    written = hl.read_table(out_path)
    n_result_rows = written.count()

    if n_result_rows == 0:
        # An empty output is indistinguishable from a silent bug unless
        # it's called out explicitly -- e.g. a gene with no qualifying
        # canonical-transcript variants in this cohort/version, or an
        # overly narrow interval/filter combination.
        scope = f"gene {gene}" if gene else f"chrom {chrom}" if chrom else "the requested scope"
        print(
            f"NOTE: 0 rows written to {out_path} -- no (gene, variant_class, sample) "
            f"combinations found for {scope} ({gnomad_version}). This is a real "
            "result (no qualifying variants survived filtering), not a write failure "
            "-- verify the gene/interval and QC flags (--release-only, "
            "--high-quality-only, --skip-filter-pass) are what you intended."
        )
    else:
        print(f"Wrote {n_result_rows} per-gene-per-class-per-individual variant count rows to {out_path}")
        if gene:
            print(f"\n(variant_class, sample) -> counts for {gene}:")
            written.show(25)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mt-path", default=None,
        help="Explicit MatrixTable/VDS path. If omitted: v4 loads via "
             "get_gnomad_v4_vds() (recommended, applies QC filtering); "
             "v2 falls back to the exomes hardcalls MT path.",
    )
    parser.add_argument("--out-path", required=True, help="Output Hail Table path")
    parser.add_argument(
        "--gnomad-version", choices=["v2", "v4"], default="v4",
        help="Controls loading path/method and canonical-field typing (int in v2, bool in v4)",
    )
    parser.add_argument(
        "--interval-path", default=None,
        help="Optional BED/interval file to restrict to specific genes/regions before processing",
    )
    parser.add_argument(
        "--chrom", default=None,
        help="Restrict to a single chromosome, e.g. '19' or 'chr19' (either "
             "form works for both versions -- normalized internally: v2 is "
             "GRCh37/'19', v4 is GRCh38/'chr19'). Applied before splitting/"
             "reading full data, for efficiency.",
    )
    parser.add_argument(
        "--release-only", action="store_true",
        help="(v4 only, via get_gnomad_v4_vds) Restrict to release samples only",
    )
    parser.add_argument(
        "--high-quality-only", action="store_true",
        help="(v4 only, via get_gnomad_v4_vds) Restrict to high-quality samples only",
    )
    parser.add_argument(
        "--skip-filter-pass", action="store_true",
        help="(v2 only) Skip joining to the release HT and filtering to filters==PASS sites",
    )
    parser.add_argument(
        "--verbose-counts", action="store_true",
        help="(v2 only) Print before/after site counts for the PASS filter. Off by "
             "default because it forces an extra full execution of the join+filter "
             "pipeline (Hail is lazy) -- only enable for debugging/small runs.",
    )
    parser.add_argument(
        "--gene", default=None,
        help="Restrict to a single gene (by canonical-transcript gene_symbol), e.g. BRCA1 -- "
             "for quickly validating the pipeline on a cheap, eyeballable slice instead of a "
             "full chromosome/genome run. Restricts to the gene's locus interval immediately "
             "(before VEP, before entry filtering), then additionally filters to "
             "gene_symbol == --gene after canonical-transcript assignment (interval overlap "
             "alone can pull in a neighboring gene too). Requires --gene-interval unless the "
             "gene is in KNOWN_GENE_INTERVALS (currently just BRCA1, PCSK9). Prints a preview "
             "of the result before writing.",
    )
    parser.add_argument(
        "--gene-interval", default=None,
        help="Locus interval for --gene, e.g. chr17:43044295-43125364 (v4/GRCh38) or "
             "17:41196312-41277500 (v2/GRCh37). Required for --gene unless the gene is in "
             "KNOWN_GENE_INTERVALS.",
    )
    parser.add_argument(
        "--gcp-project", default=None,
        help="GCP project ID to bill for reads from the requester-pays gs://gnomad and "
             "gs://gnomad_v2 buckets (e.g. your project ID from `gcloud config get-value "
             "project`). Required -- reads will fail with a 400 error without it.",
    )
    parser.add_argument(
        "--skip-v4-qc-wrapper", action="store_true",
        help="(v4 only, ignored if --mt-path is set) Bypass get_gnomad_v4_vds() and read the "
             "raw VDS directly. Skips ALL of that function's QC steps, not just the "
             "chr19:5787204 multiallelic-site drop -- also skips duplicate/withdrawn UKB "
             "sample removal and hard-filtered sample removal, since they're bundled into "
             "the same function with no separate toggle. Fine for a quick --gene smoke test; "
             "reconsider for a real production run.",
    )
    args = parser.parse_args()
    main(
        args.mt_path, args.out_path, args.gnomad_version, args.interval_path, args.chrom,
        args.release_only, args.high_quality_only, args.skip_filter_pass,
        args.verbose_counts, args.gene, args.gene_interval, args.gcp_project,
        args.skip_v4_qc_wrapper,
    )

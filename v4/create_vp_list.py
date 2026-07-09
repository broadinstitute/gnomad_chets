
import argparse
import logging
import os
import tempfile
from typing import List, Optional, Tuple

import hail as hl
from gnomad.resources.grch38.reference_data import gencode
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.vep import CSQ_ORDER, filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds

from gnomad_chets.v4.resources import (
    CLINVAR_CATEGORIES,
    CLINVAR_CATEGORY_FIELD_FMT,
    CLINVAR_CATEGORY_SOURCE_TAG,
    CLINVAR_VERSION,
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_LEAST_CONSEQUENCE,
    DEFAULT_MAX_FREQ,
    DEFAULT_EXON_DOWNSTREAM_PADDING,
    DEFAULT_EXON_UPSTREAM_PADDING,
    DEFAULT_IN_TRANS_OE_MAX_AF,
    DEFAULT_MIN_PANGOLIN,
    DEFAULT_MIN_SPLICE_AI,
    DEFAULT_TMP_DIR,
    IN_TRANS_OE_CANDIDATE_SOURCES,
    SITES_FIELD_CLINVAR,
    SITES_FIELD_PANGOLIN,
    SITES_FIELD_SPLICEAI,
    SOURCE_IN_TRANS_OE_CANDIDATE,
    SOURCE_IN_TRANS_OE_INTRONIC_PADDING,
    TEST_INTERVALS,
    get_variant_pair_resources,
)
from gnomad_chets.v4.utils import (
    AN_CUTOFFS,
    clinvar_category_match_expr,
    clinvar_review_flags_expr,
    filter_for_testing,
    get_an_percent_expr,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_vp_matrix")
logger.setLevel(logging.INFO)


def assemble_sites_ht(
    filter_ht: hl.Table,
    freq_ht: hl.Table,
    vep_ht: hl.Table,
    release_filter_ht: Optional[hl.Table] = None,
    an_ht: Optional[hl.Table] = None,
    spliceai_ht: Optional[hl.Table] = None,
    pangolin_ht: Optional[hl.Table] = None,
    clinvar_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """Join every per-variant annotation source needed downstream.

    Filters out variants that fail a *hard* QC filter — ``AS_VQSR``, or
    ``InbreedingCoeff`` when *both* filter sources agree on it — while
    retaining ``AC0`` failures and ``InbreedingCoeff`` *disagreements*, and
    joins the union of fields the variant-pair pipeline reads from
    per-variant annotations. Optional sources are omitted from the schema
    entirely when their HT argument is ``None`` —
    :func:`create_variant_filter_ht` validates that the sources needed by
    each ``include_*`` flag are present.

    **Two filter sets are carried side by side.** ``filter_ht``
    (``only_filters.ht``) and ``release_filter_ht`` (``final_filter.ht``)
    disagree *only* on the ``InbreedingCoeff`` token: they were built from
    different freq HTs, so InbreedingCoeff was recomputed — a genuine value
    change (bidirectional, ~16k variants), not float noise. ``AC0`` and
    ``AS_VQSR`` are byte-identical between them. ``final_filter.ht``'s
    ``filters`` match the v4.1.1 public release exactly (the release is
    populated from it); ``only_filters.ht`` instead matches gnomAD's current
    freq HT. Rather than commit to one InbreedingCoeff definition here, both
    are stored (``filters`` and ``filters_release``) and the
    ``InbreedingCoeff`` filter is left for downstream consumers to apply. The
    QC gate below drops a variant only when *both* sources agree it
    hard-fails: variants both flag for ``InbreedingCoeff`` are removed as
    unambiguous failures, while variants only one source flags are retained
    (the ambiguous set) for downstream resolution.

    Row-set anchor is ``filter_ht`` (``only_filters.ht``, genome-wide over
    all called variants). ``vep`` / ``freq`` / ``an`` / ``splice`` /
    ``clinvar`` become nullable left-join annotations; ``ac`` / ``af`` /
    ``an`` are null for variants absent from ``freq_ht`` (i.e. carried
    only by non-release samples, so AC=0 in the release-scoped freq).

    Always-present fields:

    * ``ac`` / ``af`` / ``an`` — global adj callstats (``freq_ht.freq[0]``
      ``.AC`` / ``.AF`` / ``.AN``); nullable.
    * ``ac_raw`` / ``af_raw`` / ``an_raw`` — global raw callstats
      (``freq_ht.freq[1]``, i.e. pre-adj); nullable.
    * ``vep`` — the full VEP struct; nullable.

    Optional fields (present iff the corresponding source HT was provided):

    * ``filters_release`` — the ``filters`` set from ``release_filter_ht``
      (``final_filter.ht``, == the v4.1.1 release ``filters``). Nullable for
      variants absent from it (the ~16M ``only_filters``-only / unreleased
      variants). Differs from ``filters`` only by ``InbreedingCoeff``.
    * ``an_pct`` — per-locus AN_percent (X/Y-aware), from
      :func:`get_an_percent_expr` on ``an_ht``.
    * ``spliceai_ds_max`` — max SpliceAI Δ.
    * ``pangolin_largest_ds`` — largest Pangolin Δ.
    * ``clinvar`` — struct with one ``is_<category>`` boolean per
      :data:`CLINVAR_CATEGORIES` (membership precomputed via
      :func:`clinvar_category_match_expr` with **relaxed**
      ``remove_no_assertion=False`` / ``remove_conflicting=False`` — 0-star
      and conflicting records are kept), a ``clinvar_review`` flag set
      (:func:`clinvar_review_flags_expr`: ``{"no_assertion"}`` /
      ``{"conflicting"}``, empty when clean — the strict set is
      ``is_<cat> & clinvar_review`` empty), plus ``GENEINFO`` for later
      VEP-symbol cross-referencing in :func:`_get_clinvar_gene_id_expr`.
      Missing struct for variants not in ClinVar.

    :param filter_ht: gnomAD final-filter Table (``only_filters.ht``, the
        all-variants variant.)
    :param freq_ht: gnomAD frequency Table.
    :param vep_ht: gnomAD VEP Table.
    :param release_filter_ht: gnomAD ``final_filter.ht``, the
        release-consistent filter Table. When provided, its ``filters`` are
        stored as ``filters_release`` so the ``InbreedingCoeff`` filter can
        be resolved downstream (see above).
    :param an_ht: gnomAD all-sites AN Table (per-locus). When provided,
        ``an_pct`` is computed via :func:`get_an_percent_expr`.
    :param spliceai_ht: SpliceAI predictor Table.
    :param pangolin_ht: Pangolin predictor Table.
    :param clinvar_ht: Unfiltered ClinVar HT.
    :return: Sites HT keyed by ``(locus, alleles)``, restricted to
        QC-PASS-or-AC0-only variants.
    """
    ht = filter_ht.select("filters")
    freq_arr = freq_ht[ht.locus, ht.alleles].freq
    freq_expr = freq_arr[0]
    raw_freq_expr = freq_arr[1]
    ann_expr = {
        "vep": vep_ht[ht.locus, ht.alleles].vep,
        "ac": freq_expr.AC,
        "af": freq_expr.AF,
        "an": freq_expr.AN,
        "ac_raw": raw_freq_expr.AC,
        "af_raw": raw_freq_expr.AF,
        "an_raw": raw_freq_expr.AN,
    }
    if release_filter_ht is not None:
        ann_expr["filters_release"] = release_filter_ht[ht.locus, ht.alleles].filters
    if an_ht is not None:
        ann_expr["an_pct"] = get_an_percent_expr(an_ht, ht.locus)
    if spliceai_ht is not None:
        ann_expr[SITES_FIELD_SPLICEAI] = spliceai_ht[ht.locus, ht.alleles].spliceai_ds_max
    if pangolin_ht is not None:
        ann_expr[SITES_FIELD_PANGOLIN] = pangolin_ht[ht.locus, ht.alleles].pangolin_largest_ds
    if clinvar_ht is not None:
        cv_row = clinvar_ht[ht.locus, ht.alleles]
        ann_expr[SITES_FIELD_CLINVAR] = hl.or_missing(
            hl.is_defined(cv_row),
            hl.struct(
                # Relaxed membership (keep 0-star + conflicting records); the
                # `clinvar_review` flag below labels why a record is borderline
                # so downstream can drop them. Strict set == is_<cat> &
                # (clinvar_review is empty).
                **{
                    CLINVAR_CATEGORY_FIELD_FMT.format(
                        category=c
                    ): clinvar_category_match_expr(
                        clnsig=cv_row.info.CLNSIG,
                        category=c,
                        remove_no_assertion=False,
                        remove_conflicting=False,
                    )
                    for c in CLINVAR_CATEGORIES
                },
                clinvar_review=clinvar_review_flags_expr(
                    cv_row.info.CLNREVSTAT, cv_row.info.CLNSIGCONF
                ),
                GENEINFO=cv_row.info.GENEINFO,
            ),
        )

    ht = ht.annotate(**ann_expr)
    
    # Keep a variant if it is PASS-or-``AC0``-only under *either* filter
    # source, dropping it only when both sources agree it hard-fails. AC0-only
    # variants are real, QC-clean, just absent from the release-cohort
    # carriers (see gnomad_qc final_filter.py, where AC0 joins the ``filters``
    # set alongside AS_VQSR / InbreedingCoeff). AS_VQSR is byte-identical
    # across the two sources, so it drops the same variants either way; the
    # sources disagree only on InbreedingCoeff. Variants both sources flag for
    # InbreedingCoeff are dropped here as unambiguous IC failures, while
    # variants only one source flags are retained (the ambiguous set) for
    # downstream resolution via ``filters`` vs ``filters_release``. When
    # ``release_filter_ht`` was not provided, gate on only_filters alone. AC /
    # AF / VEP / splice may be null and downstream consumers handle nulls.
    if "filters_release" in ht.row:
        of_pass = ht.filters.difference(hl.set(["AC0"])).length() == 0
        rel_pass = hl.is_defined(ht.filters_release) & (
            ht.filters_release.difference(hl.set(["AC0"])).length() == 0
        )
        ht = ht.filter(of_pass | rel_pass)
    else:
        ht = ht.filter(ht.filters.difference(hl.set(["AC0"])).length() == 0)

    return ht

def _validate_clinvar_categories(
    include_clinvar_categories: Optional[List[str]],
    sites_ht: hl.Table,
) -> None:
    """
    Raise if ``include_clinvar_categories`` is set but ``sites_ht`` lacks ClinVar.

    Validates each requested category against :data:`CLINVAR_CATEGORIES` and
    confirms the sites HT carries the ``clinvar`` annotation built by
    :func:`assemble_sites_ht`.

    :param include_clinvar_categories: ClinVar significance categories
        requested by the caller, or ``None`` / empty if disabled.
    :param sites_ht: Pre-assembled sites HT.
    :raises ValueError: If the ``clinvar`` field is missing from
        ``sites_ht`` or any category is unknown.
    """
    if not include_clinvar_categories:
        return
    if SITES_FIELD_CLINVAR not in sites_ht.row:
        raise ValueError(
            f"sites_ht is missing the {SITES_FIELD_CLINVAR!r} field; rebuild "
            "it with assemble_sites_ht(..., clinvar_ht=...) before requesting "
            "include_clinvar_categories."
        )
    miss_clinvar_categories = [
        c for c in include_clinvar_categories if c not in CLINVAR_CATEGORIES
    ]
    if miss_clinvar_categories:
        raise ValueError(
            f"Unknown ClinVar categories: {miss_clinvar_categories!r}. Valid: "
            f"{CLINVAR_CATEGORIES}."
        )


def _validate_pathogenic_splice(
    include_pathogenic_splice: bool,
    sites_ht: hl.Table,
) -> None:
    """
    Raise if ``include_pathogenic_splice`` is set but ``sites_ht`` lacks splice scores.

    :param include_pathogenic_splice: Whether the pathogenic-splice source
        is enabled.
    :param sites_ht: Pre-assembled sites HT.
    :raises ValueError: If either ``spliceai_ds_max`` or
        ``pangolin_largest_ds`` is missing from ``sites_ht``.
    """
    if not include_pathogenic_splice:
        return
    missing_fields = [
        f for f in (SITES_FIELD_SPLICEAI, SITES_FIELD_PANGOLIN)
        if f not in sites_ht.row
    ]
    if missing_fields:
        raise ValueError(
            f"sites_ht is missing {missing_fields!r}; rebuild it with "
            "assemble_sites_ht(..., spliceai_ht=..., pangolin_ht=...) before "
            "setting include_pathogenic_splice."
        )


def _validate_least_consequence(least_consequence: str) -> None:
    """
    Raise if ``least_consequence`` is not a known VEP consequence.

    :param least_consequence: VEP consequence string to validate against
        :data:`gnomad.utils.vep.CSQ_ORDER` (e.g. ``"3_prime_UTR_variant"``).
    :raises ValueError: If the value is not present in ``CSQ_ORDER``.
    """
    if least_consequence not in CSQ_ORDER:
        raise ValueError(
            f"least_consequence '{least_consequence}' not in CSQ_ORDER"
        )


def _get_vep_gene_id_expr(
    csq_expr: hl.expr.ArrayExpression, least_consequence: str
) -> hl.expr.ArrayExpression:
    """
    Per-variant gene_id array from PC-Ensembl transcripts at-or-above ``least_consequence``.

    Filters ``csq_expr`` (a protein-coding Ensembl
    transcript-consequence array) to transcripts carrying at least one
    ``consequence_term`` no less severe than ``least_consequence``, then
    returns the deduplicated gene_id array.

    :param csq_expr: Protein-coding Ensembl transcript_consequences array.
    :param least_consequence: Lowest-severity VEP consequence to keep, in
        the canonical :data:`gnomad.utils.vep.CSQ_ORDER` order.
    :return: ArrayExpression of ``gene_id`` strings per row.
    """
    allowed_csqs = hl.literal(
        set(CSQ_ORDER[0 : CSQ_ORDER.index(least_consequence) + 1])
    )
    return hl.array(
        hl.set(
            csq_expr.filter(
                lambda tc: tc.consequence_terms.any(lambda c: allowed_csqs.contains(c))
            ).map(lambda csq: csq.gene_id)
        )
    )


def _build_intronic_padding_interval_ht(
    gencode_ht: hl.Table,
    acceptor_padding: int,
    donor_padding: int,
    use_cache: bool = False,
    gencode_version: str = "v39",
) -> hl.Table:
    """
    Build a GENCODE exon-flanking interval Table (optionally cached).

    For each protein-coding exon in the requested GENCODE release, builds
    up to two intervals representing the strand-aware intronic flanking
    zones beyond VEP's built-in splice region. With default padding (3
    acceptor / 8 donor) these exactly overlap VEP's
    ``splice_region_variant`` range; expanding e.g. donor to 12 captures
    intron +9..+12 positions VEP doesn't tag. Padding is clipped to the
    gene body to exclude upstream/downstream regions.

    Caching is opt-in (``use_cache=True``): when enabled the result is
    read from / written to
    ``{hl.tmp_dir()}/region_intervals/region_intervals_gencode_{ver}_a{ap}_d{dp}.ht``.
    With the default ``use_cache=False`` the cache is neither read nor
    written — guard against stale or accidental reuse.

    :param gencode_ht: GENCODE annotation Table (with ``feature`` /
        ``transcript_type`` / ``interval`` / ``gene_id`` / ``strand`` fields)
        to derive protein-coding exons from.
    :param acceptor_padding: Bp to pad on the acceptor (intron) side.
    :param donor_padding: Bp to pad on the donor (intron) side.
    :param use_cache: Whether to consult / refresh the cached interval
        HT at the canonical path under :func:`hl.tmp_dir`. Default
        ``False`` always rebuilds from scratch.
    :param gencode_version: Label used only in the ``use_cache`` cache
        path; the GENCODE data itself comes from ``gencode_ht``.
    :return: Table keyed by interval with ``gene_id`` field.
    """
    cached_path = (
        f"{hl.tmp_dir()}/region_intervals/region_intervals_gencode_"
        f"{gencode_version}_a{acceptor_padding}_d{donor_padding}.ht"
    )
    if use_cache:
        try:
            ht = hl.read_table(cached_path)
            logger.info("Reusing cached region intervals from %s", cached_path)
            return ht
        except Exception:
            pass

    exons = gencode_ht.filter(
        (gencode_ht.feature == "exon")
        & (gencode_ht.transcript_type == "protein_coding")
    )

    # Per-gene genomic extent of the exon body — used to skip padding on
    # the genomic-leftmost / -rightmost exon (padding past those edges
    # would be upstream/downstream of the gene, not intronic). This is
    # strand-unaware on purpose: strand only determines whether each
    # genomic side is acceptor or donor (see ``start_pad`` / ``end_pad``
    # below), not where the gene body ends.
    gene_ranges = exons.group_by(exons.gene_id).aggregate(
        gene_body_min=hl.agg.min(exons.interval.start.position),
        gene_body_max=hl.agg.max(exons.interval.end.position),
    )
    exons = exons.annotate(**gene_ranges[exons.gene_id])

    is_plus = exons.strand == "+"
    start_pad = hl.if_else(is_plus, acceptor_padding, donor_padding)
    end_pad = hl.if_else(is_plus, donor_padding, acceptor_padding)

    contig = exons.interval.start.contig
    exon_start = exons.interval.start.position
    exon_end = exons.interval.end.position
    before_pad_start = hl.max(1, exon_start - start_pad)

    # Build at most two flanking intervals per exon. ``hl.or_missing``
    # returns missing for clipped cases (genomic-edge exon, zero padding,
    # or padding clipped to gene boundary); we drop missing entries before
    # exploding. GENCODE intervals are 1-based fully-closed, so
    # ``exon_start`` is the first exon base and ``exon_end`` is the last.
    before_interval = hl.or_missing(
        (start_pad > 0)
        & (exon_start > exons.gene_body_min)
        & (before_pad_start < exon_start),
        hl.struct(
            interval=hl.interval(
                hl.locus(contig, before_pad_start, "GRCh38"),
                hl.locus(contig, exon_start, "GRCh38"),
            ),
            gene_id=exons.gene_id,
        ),
    )
    after_interval = hl.or_missing(
        (end_pad > 0) & (exon_end < exons.gene_body_max),
        hl.struct(
            interval=hl.interval(
                hl.locus(contig, exon_end, "GRCh38"),
                hl.locus(contig, exon_end + end_pad, "GRCh38"),
                includes_end=True,
            ),
            gene_id=exons.gene_id,
        ),
    )

    all_interval_ht = exons.select(
        _intervals=hl.array([before_interval, after_interval]).filter(
            lambda x: hl.is_defined(x)
        )
    )
    all_interval_ht = all_interval_ht.explode("_intervals").key_by()
    all_interval_ht = all_interval_ht.select(
        interval=all_interval_ht._intervals.interval,
        gene_id=all_interval_ht._intervals.gene_id,
    )
    all_interval_ht = all_interval_ht.key_by("interval")

    if use_cache:
        all_interval_ht = all_interval_ht.checkpoint(
            cached_path, overwrite=True,
        )
        logger.info(
            "Cached region intervals (acceptor: %d bp, donor: %d bp) to %s",
            acceptor_padding, donor_padding, cached_path,
        )
    return all_interval_ht


def _get_intronic_padding_gene_id_expr(
    ht: hl.Table,
    gencode_ht: hl.Table,
    acceptor_padding: int = DEFAULT_EXON_UPSTREAM_PADDING,
    donor_padding: int = DEFAULT_EXON_DOWNSTREAM_PADDING,
    use_cache: bool = False,
    gencode_version: str = "v39",
) -> hl.expr.ArrayExpression:
    """
    Per-variant gene_id array for variants in GENCODE exon-flanking intronic zones.

    VEP already handles the exon body and its built-in splice region (±8 bp
    on the intron side). This function builds the GENCODE exon-flanking
    interval HT for the requested padding (optionally cached; see
    ``use_cache``), then returns the deduplicated array of overlapping
    ``gene_id`` values per row of ``ht``. With the default padding
    (acceptor=3, donor=8) the intervals exactly overlap VEP's
    ``splice_region_variant`` range and add zero new variants; expanding
    e.g. ``donor_padding=12`` captures positions +9..+12 that VEP misses.
    Padding is clipped to the gene body so upstream/downstream regions
    are excluded. Variants outside every interval get an empty array.

    :param ht: Table with a ``locus`` field; the returned expression is
        bound to ``ht``'s row context.
    :param gencode_ht: GENCODE annotation Table to derive exons from
        (passed through to :func:`_build_intronic_padding_interval_ht`).
    :param acceptor_padding: Bp to pad on the acceptor (intron) side of
        each exon boundary.
    :param donor_padding: Bp to pad on the donor (intron) side of each
        exon boundary.
    :param use_cache: Whether to consult / refresh the cached interval
        HT under :func:`hl.tmp_dir`. Default ``False`` always rebuilds
        from scratch — opt in to reuse the GCS cache.
    :param gencode_version: GENCODE release for the exon source.
    :return: ArrayExpression of overlapping ``gene_id`` strings per row.
    """
    all_interval_ht = _build_intronic_padding_interval_ht(
        gencode_ht,
        acceptor_padding=acceptor_padding,
        donor_padding=donor_padding,
        use_cache=use_cache,
        gencode_version=gencode_version,
    )
    return hl.array(
        hl.set(all_interval_ht.index(ht.locus, all_matches=True).gene_id)
    )


def _get_clinvar_gene_id_expr(
    clinvar_expr: hl.expr.StructExpression,
    csq_expr: hl.expr.ArrayExpression,
    category: str,
) -> Tuple[hl.expr.ArrayExpression, hl.expr.BooleanExpression]:
    """
    Per-variant ClinVar gene_id array (+ fallback flag) for one category.

    Returns a pair ``(gene_ids, used_fallback)``:

    * ``gene_ids`` — gene_ids assigned to the variant for ``category``.
      Empty when the variant isn't in ClinVar for this category. When
      the variant *is* in ClinVar for the category, gene_ids are
      restricted to the intersection of VEP protein-coding Ensembl
      transcript gene_symbols and the symbols listed in ClinVar's
      ``GENEINFO``. If that intersection is empty (ClinVar symbol drift,
      retired gene names, etc.) we fall back to the full VEP gene_id
      list so the variant isn't silently dropped.
    * ``used_fallback`` — ``True`` only when the fallback path above
      was taken (i.e. variant is in ClinVar for the category but no VEP
      transcript matched ``GENEINFO`` symbols).

    Category membership is read off the precomputed
    ``clinvar_expr.is_<category>`` boolean (set during
    :func:`assemble_sites_ht`); no per-call category gating runs here.
    Unlike :func:`_get_vep_gene_id_expr`, no consequence filter is
    applied — ClinVar curation status is the signal; gene assignment
    just follows VEP.

    :param clinvar_expr: ``clinvar`` struct on the sites HT (see
        :func:`assemble_sites_ht`). Carries the per-category booleans
        and ``GENEINFO``; missing for variants not in ClinVar.
    :param csq_expr: Protein-coding Ensembl transcript_consequences array.
    :param category: A single value from :data:`CLINVAR_CATEGORIES`.
    :return: ``(gene_ids, used_fallback)`` ArrayExpression /
        BooleanExpression pair.
    """
    vep_gene_ids = hl.array(hl.set(csq_expr.map(lambda csq: csq.gene_id)))

    # GENEINFO is a single string "SYMBOL1:Entrez1|SYMBOL2:Entrez2"; missing
    # -> empty set. Drop empty splits defensively. Note: when ``clinvar_expr``
    # itself is missing the propagation makes ``clinvar_expr.GENEINFO`` missing
    # too, which the coalesce flattens to empty before the split.
    cv_symbols = hl.set(
        hl.coalesce(clinvar_expr.GENEINFO, "")
        .split(r"\|")
        .map(lambda entry: entry.split(":")[0])
        .filter(lambda sym: sym != "")
    )
    matched_gene_ids = hl.array(hl.set(
        csq_expr.filter(
            lambda tc: cv_symbols.contains(tc.gene_symbol)
        ).map(lambda tc: tc.gene_id)
    ))
    is_in_category = hl.coalesce(
        clinvar_expr[CLINVAR_CATEGORY_FIELD_FMT.format(category=category)], False
    )
    used_fallback = is_in_category & (hl.len(matched_gene_ids) == 0)
    gene_ids = hl.if_else(
        is_in_category,
        hl.if_else(used_fallback, vep_gene_ids, matched_gene_ids),
        hl.empty_array(hl.tstr),
    )
    return gene_ids, used_fallback


def _get_hc_lof_gene_id_expr(
    csq_expr: hl.expr.ArrayExpression
) -> hl.expr.ArrayExpression:
    """
    Per-variant gene_id array from HC-LoF Ensembl protein-coding transcripts.

    Returns the deduplicated array of ``gene_id`` values for transcripts
    where ``lof == "HC"``. The ``lof_filter`` flag is intentionally not
    enforced — variants with LOFTEE flags (``NAGNAG_SITE``,
    ``single_exon``, etc.) still count. Gene IDs come only from the
    qualifying HC-LoF transcripts so a variant is gene-assigned only to
    genes where it actually disrupts function. Variants without a
    qualifying transcript get an empty array.

    :param csq_expr: Protein-coding Ensembl transcript_consequences array.
    :return: ArrayExpression of ``gene_id`` strings per row.
    """
    return hl.array(
        hl.set(
            csq_expr.filter(lambda tc: tc.lof == "HC").map(
                lambda csq: csq.gene_id
            )
        )
    )


def create_variant_filter_ht(
    ht: hl.Table,
    least_consequence: str = DEFAULT_LEAST_CONSEQUENCE,
    max_freq: float = DEFAULT_MAX_FREQ,
    include_extra_padding: bool = False,
    acceptor_padding: int = DEFAULT_EXON_UPSTREAM_PADDING,
    donor_padding: int = DEFAULT_EXON_DOWNSTREAM_PADDING,
    include_clinvar_categories: Optional[List[str]] = None,
    include_pathogenic_splice: bool = False,
    include_hc_lof: bool = False,
    min_splice_ai: float = DEFAULT_MIN_SPLICE_AI,
    min_pangolin: float = DEFAULT_MIN_PANGOLIN,
    include_in_trans_oe_candidates: bool = False,
    include_in_trans_oe_intronic_padding: bool = False,
    in_trans_oe_max_af: float = DEFAULT_IN_TRANS_OE_MAX_AF,
    in_trans_oe_acceptor_padding: int = 50,
    in_trans_oe_donor_padding: int = 15,
    gencode_ht: Optional[hl.Table] = None,
    use_cache: bool = False,
) -> hl.Table:
    """
    Create a filter Table for variant pair determination.

    Reads the pre-assembled sites HT produced by :func:`assemble_sites_ht`
    and filters variants to those that pass variant QC, have a consequence
    at least as severe as ``least_consequence``, have a global AF in
    ``(0, max_freq]``, and have an associated gene ID.

    Optionally unions additional variant sets on top of the base VEP filter,
    each gated by its own ``include_*`` flag. Each source has a single
    well-scoped tag emitted into the ``source`` set:

    - ``include_extra_padding`` -> ``gencode_extra_padding``: GENCODE exon
      flanking intronic zones beyond VEP's built-in splice region (±8 bp).
    - ``include_clinvar_categories`` -> one tag per category, e.g.
      ``clinvar_plp`` / ``clinvar_blb`` / ``clinvar_vus`` (AF-capped at
      ``in_trans_oe_max_af`` rather than ``max_freq``, to retain
      common-but-pathogenic alleles). Requires the sites HT to carry the
      ``clinvar`` annotation.
    - ``include_pathogenic_splice`` -> ``splice_path``: SpliceAI Δ >
      ``min_splice_ai`` OR Pangolin Δ > ``min_pangolin``. Requires the
      sites HT to carry ``spliceai_ds_max`` and ``pangolin_largest_ds``.
    - ``include_hc_lof`` -> ``hc_lof``: VEP ``lof == "HC"`` on at least one
      protein-coding Ensembl transcript.
    - ``include_in_trans_oe_candidates`` -> ``in_trans_oe_candidate``:
      higher-AF (up to ``in_trans_oe_max_af``) candidates for the in-trans
      observed-vs-expected analysis.
    - ``include_in_trans_oe_intronic_padding`` ->
      ``in_trans_oe_intronic_padding``: deep intronic windows (branch-point
      region + cryptic 5' splice signals).

    When multiple sources are included, gene IDs are unioned across sources.
    Each variant is annotated with a ``source`` set indicating which
    filter(s) included it. ``af`` (global adj AF) and ``an_pct`` are
    propagated to the output so an AF threshold can be re-applied downstream
    (the source tags already encode the ``max_freq`` / ``in_trans_oe_max_af``
    AF caps at include time). The
    variant's ``filters`` (only_filters) and ``filters_release``
    (release/final_filter, when present on the input) are carried through so
    InbreedingCoeff / QC filtering can be applied at the very end of the
    pipeline rather than gated here. The build parameters are recorded in the
    ``variant_filter_params`` global for provenance.

    :param ht: Pre-assembled sites HT produced by
        :func:`assemble_sites_ht`.
    :param least_consequence: Lowest-severity VEP consequence kept by the
        base filter (must be in :data:`gnomad.utils.vep.CSQ_ORDER`).
    :param max_freq: Upper AF bound (inclusive) for the base filter.
    :param include_extra_padding: Include GENCODE exon-flanking variants
        (tag ``gencode_extra_padding``).
    :param acceptor_padding: Bp upstream of acceptor splice site for the
        extra-padding source.
    :param donor_padding: Bp downstream of donor splice site for the
        extra-padding source.
    :param include_clinvar_categories: ClinVar significance categories to
        union in (any combination of ``"plp"``, ``"blb"``, ``"vus"``).
        Each emits its own ``clinvar_<cat>`` source tag.
    :param include_pathogenic_splice: Include SpliceAI/Pangolin-flagged
        variants (tag ``splice_path``).
    :param include_hc_lof: Include HC-LoF variants (tag ``hc_lof``).
    :param min_splice_ai: Minimum SpliceAI delta score for
        ``include_pathogenic_splice``.
    :param min_pangolin: Minimum Pangolin delta score for
        ``include_pathogenic_splice``.
    :param include_in_trans_oe_candidates: Include in-trans-OE higher-AF
        candidates (tag ``in_trans_oe_candidate``).
    :param include_in_trans_oe_intronic_padding: Include in-trans-OE
        intronic padding (tag ``in_trans_oe_intronic_padding``).
    :param in_trans_oe_max_af: Upper AF bound for in-trans-OE candidates
        and intronic padding (inclusive).
    :param in_trans_oe_acceptor_padding: Bp upstream of acceptor splice
        site for the in-trans-OE intronic padding source (typically 50,
        covering the branch-point region).
    :param in_trans_oe_donor_padding: Bp downstream of donor splice site
        for the in-trans-OE intronic padding source (typically 15,
        covering cryptic 5' splice signals just past VEP's +8).
    :param gencode_ht: GENCODE annotation Table, required when
        ``include_extra_padding`` or ``include_in_trans_oe_intronic_padding``
        is set (used to derive exon-flanking intronic intervals). May be
        ``None`` otherwise.
    :param use_cache: Whether to consult / refresh the cached GENCODE
        interval HTs under :func:`hl.tmp_dir`. Default ``False`` rebuilds
        from scratch each call (guard against stale cache reuse); set
        ``True`` to opt into the cache.
    :return: Table keyed by ``(locus, alleles)`` with ``gene_id``
        (array<str>) and ``source`` (set<str>) fields, plus ``an_pct``
        (int32) and the ``an_cutoffs`` global when ``ht`` carries ``an_pct``.
    """
    _validate_clinvar_categories(include_clinvar_categories, ht)
    _validate_pathogenic_splice(include_pathogenic_splice, ht)
    _validate_least_consequence(least_consequence)
    if (
        include_extra_padding or include_in_trans_oe_intronic_padding
    ) and gencode_ht is None:
        raise ValueError(
            "gencode_ht is required when include_extra_padding or "
            "include_in_trans_oe_intronic_padding is set."
        )

    # NULL af (variant absent from the release-scoped freq HT — e.g. AC0-in-
    # release variants carried only by non-release samples) is treated as
    # passing the AF cap. Otherwise the NA propagates through the per-tag
    # ``hl.if_else`` in add_filters_expr, NA-poisons the unioned ``source``
    # set, and the variant is silently dropped even when it qualifies for a
    # source (ClinVar P/LP, HC-LoF, padding, ...). Inclusion-safe: keep it.
    af_filter_expr = hl.is_missing(ht.af) | (ht.af <= max_freq)
    in_trans_oe_af_expr = hl.is_missing(ht.af) | (ht.af <= in_trans_oe_max_af)

    # Filter VEP transcripts to protein-coding Ensembl once and share
    # across every per-source gene_id helper.
    csq_expr = filter_vep_transcript_csqs_expr(
        ht.vep.transcript_consequences,
        protein_coding=True,
        ensembl_only=True,
    )
    vep_gene_id_expr = _get_vep_gene_id_expr(csq_expr, least_consequence)
    has_vep_gene = hl.len(vep_gene_id_expr) > 0
    gene_id_set_expr = hl.set(vep_gene_id_expr)
    source_tag_expr = {"vep_csq": af_filter_expr & has_vep_gene}

    if include_in_trans_oe_candidates:
        logger.info(
            "Including in-trans-OE candidates (AF in (0, %g], consequence ≥ %s)...",
            in_trans_oe_max_af, least_consequence,
        )
        # Same gene_id derivation as vep_csq; differs only in AF threshold.
        source_tag_expr[SOURCE_IN_TRANS_OE_CANDIDATE] = (
            in_trans_oe_af_expr & has_vep_gene
        )
    if include_extra_padding:
        logger.info(
            "Including GENCODE extra padding variants (acceptor=%d, donor=%d)...",
            acceptor_padding, donor_padding,
        )
        _gene_id_expr = _get_intronic_padding_gene_id_expr(
            ht, gencode_ht, acceptor_padding, donor_padding, use_cache=use_cache
        )
        source_tag_expr["gencode_extra_padding"] = (
            af_filter_expr & (hl.len(_gene_id_expr) > 0)
        )
        gene_id_set_expr = gene_id_set_expr.union(hl.set(_gene_id_expr))
    if include_in_trans_oe_intronic_padding:
        logger.info(
            "Including in-trans-OE intronic padding (acceptor=%dbp, donor=%dbp)...",
            in_trans_oe_acceptor_padding, in_trans_oe_donor_padding,
        )
        _gene_id_expr = _get_intronic_padding_gene_id_expr(
            ht,
            gencode_ht,
            in_trans_oe_acceptor_padding,
            in_trans_oe_donor_padding,
            use_cache=use_cache,
        )
        source_tag_expr[SOURCE_IN_TRANS_OE_INTRONIC_PADDING] = (
            in_trans_oe_af_expr & (hl.len(_gene_id_expr) > 0)
        )
        gene_id_set_expr = gene_id_set_expr.union(hl.set(_gene_id_expr))
    clinvar_fallback_by_tag = {}
    if include_clinvar_categories:
        clinvar_expr = ht[SITES_FIELD_CLINVAR]
        for category in include_clinvar_categories:
            tag = CLINVAR_CATEGORY_SOURCE_TAG[category]
            logger.info(
                "Including ClinVar %s variants (tagged %s)...", category.upper(), tag
            )
            _gene_id_expr, _fallback_expr = _get_clinvar_gene_id_expr(
                clinvar_expr, csq_expr, category
            )
            # ClinVar tags use the higher in-trans-OE AF cap (not max_freq):
            # keep clinically-important common-but-pathogenic alleles, bounded
            # at in_trans_oe_max_af so pairs don't explode.
            source_tag_expr[tag] = in_trans_oe_af_expr & (hl.len(_gene_id_expr) > 0)
            gene_id_set_expr = gene_id_set_expr.union(hl.set(_gene_id_expr))
            clinvar_fallback_by_tag[tag] = _fallback_expr
    if include_pathogenic_splice:
        logger.info(
            "Including pathogenic-splice variants (tagged splice_path)..."
        )
        spliceai_expr = ht[SITES_FIELD_SPLICEAI]
        pangolin_expr = ht[SITES_FIELD_PANGOLIN]
        is_splice_path = (
            (hl.is_defined(spliceai_expr) & (spliceai_expr > min_splice_ai))
            | (hl.is_defined(pangolin_expr) & (pangolin_expr > min_pangolin))
        )
        source_tag_expr["splice_path"] = af_filter_expr & is_splice_path
        # Splice-path gene_id is the PC-Ensembl gene_id set (no consequence
        # filter — deep-intronic variants qualify when the predictors flag
        # them). Contribute when the flag is True.
        pc_gene_ids = hl.set(csq_expr.map(lambda c: c.gene_id))
        gene_id_set_expr = gene_id_set_expr.union(
            hl.if_else(is_splice_path, pc_gene_ids, hl.empty_set(hl.tstr))
        )
    if include_hc_lof:
        logger.info("Including HC LoF variants (tagged hc_lof)...")
        _gene_id_expr = _get_hc_lof_gene_id_expr(csq_expr)
        source_tag_expr["hc_lof"] = af_filter_expr & (hl.len(_gene_id_expr) > 0)
        gene_id_set_expr = gene_id_set_expr.union(hl.set(_gene_id_expr))

    # Collapse per-category clinvar gene-match-fallback bools into a set
    # of the source tags where fallback was triggered (empty when nothing
    # fell back). Only emitted when ClinVar is in the picture.
    clinvar_fallback_expr = hl.empty_set(hl.tstr)
    for tag, fallback_expr in clinvar_fallback_by_tag.items():
        clinvar_fallback_expr = clinvar_fallback_expr.union(
            hl.if_else(
                fallback_expr,
                hl.set([tag]),
                hl.empty_set(hl.tstr),
            )
        )
    ht = ht.select(
        gene_id=hl.array(gene_id_set_expr),
        source=add_filters_expr(source_tag_expr),
        af=ht.af,
        filters=ht.filters,
        **(
            {"filters_release": ht.filters_release}
            if "filters_release" in ht.row else {}
        ),
        **({"an_pct": ht.an_pct} if "an_pct" in ht.row else {}),
        **(
            {
                # Carry the ClinVar quality flags so relaxed (0-star /
                # conflicting) records are labeled and droppable downstream;
                # empty for clean or non-ClinVar variants.
                "clinvar_review": hl.coalesce(
                    ht[SITES_FIELD_CLINVAR].clinvar_review, hl.empty_set(hl.tstr)
                ),
                "clinvar_gene_match_fallback": clinvar_fallback_expr,
            }
            if include_clinvar_categories else {}
        ),
    )
    ht = ht.filter(hl.len(ht.source) > 0)
    if "an_pct" in ht.row:
        ht = ht.annotate_globals(an_cutoffs=hl.literal(AN_CUTOFFS))

    # Record the parameters this filter was built with (provenance).
    ht = ht.annotate_globals(
        variant_filter_params=hl.struct(
            least_consequence=least_consequence,
            max_freq=max_freq,
            include_extra_padding=include_extra_padding,
            acceptor_padding=acceptor_padding,
            donor_padding=donor_padding,
            include_clinvar_categories=hl.literal(
                include_clinvar_categories or [], hl.tarray(hl.tstr)
            ),
            include_pathogenic_splice=include_pathogenic_splice,
            include_hc_lof=include_hc_lof,
            min_splice_ai=min_splice_ai,
            min_pangolin=min_pangolin,
            include_in_trans_oe_candidates=include_in_trans_oe_candidates,
            include_in_trans_oe_intronic_padding=include_in_trans_oe_intronic_padding,
            in_trans_oe_max_af=in_trans_oe_max_af,
            in_trans_oe_acceptor_padding=in_trans_oe_acceptor_padding,
            in_trans_oe_donor_padding=in_trans_oe_donor_padding,
        )
    )
    return ht


def _get_ordered_vp_struct(
    v1: hl.expr.StructExpression, v2: hl.expr.StructExpression
) -> hl.expr.StructExpression:
    """
    Create an ordered variant pair struct ensuring consistent ordering.

    Orders variants by position first, then by alt allele if positions are equal.
    This ensures that (v1, v2) and (v2, v1) are treated as the same pair.

    :param v1: First variant struct with fields 'locus' and 'alleles'.
    :param v2: Second variant struct with fields 'locus' and 'alleles'.
    :return: Struct with fields 'v1' and 'v2' in canonical order.
    """
    return hl.if_else(
        v1 <= v2,
        hl.struct(v1=v1, v2=v2),
        hl.struct(v1=v2, v2=v1),
    )


def filter_pair_ht_to_in_trans_oe_pairs(
    pair_ht: hl.Table,
    filter_ht: hl.Table,
    candidate_sources: Optional[List[str]] = None,
) -> hl.Table:
    """
    Drop pairs where BOTH sides are OE-candidate-only.

    A side is "OE-candidate-only" if every source tag on it is in
    ``candidate_sources`` (default :data:`IN_TRANS_OE_CANDIDATE_SOURCES`);
    i.e., the variant is in the filter HT *only* because of the
    in-trans-OE candidate extensions (higher AF range, intronic padding)
    and has no other annotation justifying inclusion.

    The candidate × candidate space is the explosion to avoid: pairs of
    two OE-only variants only exist because of the OE expansion and would
    inflate the pair list combinatorially. Pairs where AT LEAST ONE side
    has a non-candidate source are retained — this includes all baseline
    (``vep_csq``) pairs and partner-side pairs regardless of the other
    side (including OE-candidate × partner pairs, which are exactly what
    the in-trans-OE feature needs to test higher-AF candidates against).

    :param pair_ht: Output of :func:`create_variant_pair_ht`, keyed by
        ``(locus1, alleles1, locus2, alleles2)``.
    :param filter_ht: Variant filter HT with ``source`` field (set<str>),
        keyed by ``(locus, alleles)``.
    :param candidate_sources: Source tags considered "OE-candidate-only".
        Defaults to :data:`IN_TRANS_OE_CANDIDATE_SOURCES`.
    :return: Filtered pair Table.
    """
    candidate_set = hl.set(candidate_sources or IN_TRANS_OE_CANDIDATE_SOURCES)
    pair_ht = pair_ht.annotate(
        _v1_source=filter_ht[pair_ht.locus1, pair_ht.alleles1].source,
        _v2_source=filter_ht[pair_ht.locus2, pair_ht.alleles2].source,
    )
    pair_ht = pair_ht.filter(
        (hl.len(pair_ht._v1_source.difference(candidate_set)) > 0)
        | (hl.len(pair_ht._v2_source.difference(candidate_set)) > 0)
    )
    return pair_ht.drop("_v1_source", "_v2_source")


def create_variant_pair_ht(
    mt: hl.MatrixTable,
    filter_ht: hl.Table,
    *,
    drop_oe_only_pairs: bool = False,
) -> hl.Table:
    """
    Create a Hail Table of unique ordered variant pairs per sample per gene.

    Pair rows are annotated with per-side ``an_pct1`` / ``an_pct2`` and the
    standard cutoff list is copied into globals as ``an_cutoffs``. When
    ``drop_oe_only_pairs`` is True (typically set with the
    ``--include-in-trans-oe-candidates`` flag), pairs where ≥1 side is
    OE-candidate-only are dropped via
    :func:`gnomad_chets.v4.in_trans_oe.filter_pair_ht_to_in_trans_oe_pairs`
    to avoid candidate × candidate explosion while preserving baseline
    pairs.

    :param mt: MatrixTable with filtered variant data.
    :param filter_ht: Variant filter Table (output of
        :func:`create_variant_filter_ht` + AN annotation). Must be keyed by
        ``(locus, alleles)`` with ``gene_id``, ``an_pct``, and (when
        ``drop_oe_only_pairs`` is True) ``source`` fields.
    :param drop_oe_only_pairs: Drop pairs where both sides are
        OE-candidate-only.
    :return: Hail Table keyed by ``(locus1, alleles1, locus2, alleles2)``
        with one row per unique variant pair plus ``gene_id``, per-side
        ``an_pct{1,2}`` annotations, and ``an_cutoffs`` global.
    """
    n_partitions = filter_ht.n_partitions()
    mt = mt.add_row_index("variant_idx")
    variant_index_ht = mt.rows().key_by("variant_idx").cache()

    mt = mt.annotate_rows(gene_id=filter_ht[mt.locus, mt.alleles].gene_id)

    # Note: We do it this way because a row grouping by gene_id results in an 
    # aggregation with one gene per partition, and this can lead to memory issues.
    # Convert to entries table and explode on gene_id so each variant-gene combination
    # is a separate row.
    ht = mt.select_cols().select_rows("variant_idx", "gene_id").entries()
    ht = ht.filter(ht.GT.is_non_ref())
    ht = ht.explode("gene_id")

    # Group by gene and sample, collecting unique variants per gene/sample.
    # Using collect_as_set ensures each variant appears only once per gene/sample.
    ht = ht.group_by("gene_id", "s").aggregate(
        variants=hl.array(hl.agg.collect_as_set(ht.variant_idx))
    )

    # Filter to samples with at least 2 variants (needed to form pairs).
    ht = ht.filter(ht.variants.length() >= 2)
    ht = ht.checkpoint(
        hl.utils.new_temp_file("create_variant_pair_ht.gene_sample_grouped", "ht")
    )

    # Generate all ordered pairs of variants within each gene/sample.
    # The nested flatmap/map creates all combinations (i, j) where i < j, ensuring
    # each pair is created exactly once.
    ht = ht.annotate(
        pairs=(
            hl.range(0, hl.len(ht.variants)).flatmap(
                lambda i1: (
                    hl.range(i1 + 1, hl.len(ht.variants)).map(
                        lambda i2: _get_ordered_vp_struct(
                            ht.variants[i1], ht.variants[i2]
                        )
                    )
                )
            )
        )
    )

    # Explode pairs.
    ht = ht.explode("pairs")

    # Key by variant pair and select distinct pairs.
    # Use new shuffle method for apply models to prevent shuffle errors.
    hl._set_flags(use_new_shuffle="1")
    ht = ht.group_by(v1=ht.pairs.v1, v2=ht.pairs.v2).aggregate(
        gene_id=hl.agg.collect_as_set(ht.gene_id)
    )
    # Restore partition count; group_by shuffle often coalesces to few partitions
    # (e.g. spark.sql.shuffle.partitions=24), which would carry through to the
    # written variant pair table and downstream steps.
    ht = ht.repartition(n_partitions, shuffle=True).cache()
    hl._set_flags(use_new_shuffle=None)

    # Add a unique index id to each variant pair.
    # This is used later so both variants can be annotated with genotype
    # info separately and then joined together by the common index. This helps with
    # performance issues observed when trying to annotate both variants with genotype
    # info simultaneously.
    ht = ht.add_index("vp_ht_idx").key_by("vp_ht_idx")

    variant_index_keyed_v1 = variant_index_ht[ht.v1]
    variant_index_keyed_v2 = variant_index_ht[ht.v2]
    ht = ht.select(
        "gene_id",
        locus1=variant_index_keyed_v1.locus,
        alleles1=variant_index_keyed_v1.alleles,
        locus2=variant_index_keyed_v2.locus,
        alleles2=variant_index_keyed_v2.alleles,
    )

    if drop_oe_only_pairs:
        logger.info(
            "Filtering variant pair list to drop pairs where both sides are "
            "OE-candidate-only (avoids candidate × candidate explosion while "
            "preserving baseline pairs)."
        )
        ht = filter_pair_ht_to_in_trans_oe_pairs(ht, filter_ht)

    # Propagate per-side AN_percent from the filter HT and record the
    # standard cutoff list in globals so downstream code can threshold.
    ht = ht.annotate(
        an_pct1=filter_ht[ht.locus1, ht.alleles1].an_pct,
        an_pct2=filter_ht[ht.locus2, ht.alleles2].an_pct,
    )
    ht = ht.annotate_globals(an_cutoffs=hl.literal(AN_CUTOFFS))

    return ht


def main(args):
    """Create variant pair matrix from gnomAD v4 VDS."""

    # ------------------------------------------------------------
    # Initialize arguments.
    # ------------------------------------------------------------
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix
    data_type = args.data_type
    least_consequence = args.least_consequence
    max_freq = args.max_freq
    min_an_pct = args.min_an_pct
    
    # ------------------------------------------------------------
    # Initialize test intervals.
    # ------------------------------------------------------------
    scope_flags = [bool(args.gene), bool(args.test_genes), bool(args.interval)]
    if sum(scope_flags) > 1:
        raise ValueError(
            "--gene, --test-genes, and --interval are mutually exclusive."
        )
    
    test = args.test or any(scope_flags)
    if args.gene:
        test_intervals = {args.gene: TEST_INTERVALS[args.gene]}
    elif args.test_genes:
        wanted = [g.strip() for g in args.test_genes.split(",") if g.strip()]
        missing = [g for g in wanted if g not in TEST_INTERVALS]
        if missing:
            raise ValueError(
                f"--test-genes contains unknown TEST_INTERVALS keys: {missing}. "
                f"Available: {sorted(TEST_INTERVALS)}"
            )
        test_intervals = {g: TEST_INTERVALS[g] for g in wanted}
    elif args.interval:
        test_intervals = {args.interval: args.interval}
    else:
        test_intervals = TEST_INTERVALS
    
    # Normalize --test-chrom (e.g. "5" -> "chr5"); None when not set.
    test_chrom = args.test_chrom
    if test_chrom and not test_chrom.startswith("chr"):
        test_chrom = f"chr{test_chrom}"

    hl.init(
        log=os.path.join(tempfile.gettempdir(), "create_vp_matrix.log"),
        tmp_dir=tmp_dir,
    )

    logger.info(
        f"""
        Running script with the following parameters:

            Data type: {data_type}
            Test: {test}
            Gene: {args.gene or 'all test intervals'}
            Output postfix: {output_postfix}
            Overwrite: {overwrite}
            Tmp dir: {tmp_dir}
        """
    )

    # ------------------------------------------------------------
    # Initialize resources.
    # ------------------------------------------------------------
    resources = get_variant_pair_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir if test else None,
        output_postfix=output_postfix,
        overwrite=overwrite,
    )
    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )

    # ------------------------------------------------------------
    # Run pipeline steps.
    # ------------------------------------------------------------
    if args.preprocess_sites_ht:
        logger.info("Assembling per-variant sites HT...")
        res = resources.preprocess_sites_ht
        res.check_resource_existence()

        filter_ht = res.filter_ht.ht()
        release_filter_ht = res.release_filter_ht.ht()
        freq_ht = res.freq_ht.ht()
        vep_ht = res.vep_ht.ht()
        an_ht = res.an_ht.ht()
        spliceai_ht = res.spliceai_ht.ht()
        pangolin_ht = res.pangolin_ht.ht()
        clinvar_ht = res.clinvar_ht.ht()

        # Filter all input HTs to the test interval in one place.
        if test:
            logger.info("Filtering input HTs to test interval...")
            filter_ht = filter_for_testing(filter_ht, test_intervals)
            release_filter_ht = filter_for_testing(release_filter_ht, test_intervals)
            freq_ht = filter_for_testing(freq_ht, test_intervals)
            vep_ht = filter_for_testing(vep_ht, test_intervals)
            an_ht = filter_for_testing(an_ht, test_intervals)
            spliceai_ht = filter_for_testing(spliceai_ht, test_intervals)
            pangolin_ht = filter_for_testing(pangolin_ht, test_intervals)
            clinvar_ht = filter_for_testing(clinvar_ht, test_intervals)

        sites_ht = assemble_sites_ht(
            filter_ht=filter_ht,
            freq_ht=freq_ht,
            vep_ht=vep_ht,
            release_filter_ht=release_filter_ht,
            an_ht=an_ht,
            spliceai_ht=spliceai_ht,
            pangolin_ht=pangolin_ht,
            clinvar_ht=clinvar_ht,
        )

        # Stamp source provenance onto the sites HT globals: the resolved
        # path of every input (each path encodes its gnomAD / annotation
        # version) plus the pinned ClinVar release version.
        sites_ht = sites_ht.annotate_globals(
            source_paths=hl.struct(
                filter_ht=res.filter_ht.path,
                release_filter_ht=res.release_filter_ht.path,
                freq_ht=res.freq_ht.path,
                vep_ht=res.vep_ht.path,
                an_ht=res.an_ht.path,
                spliceai_ht=res.spliceai_ht.path,
                pangolin_ht=res.pangolin_ht.path,
                clinvar_ht=res.clinvar_ht.path,
            ),
            clinvar_version=CLINVAR_VERSION,
        )

        # Coalesce down before the write so we don't produce tens of
        # thousands of tiny part files from VEP's input partitioning.
        sites_ht = sites_ht.naive_coalesce(5000)
        sites_ht = sites_ht.checkpoint(res.sites_ht.path, overwrite=overwrite)
        logger.info("Number of variants in the sites HT: %d", sites_ht.count())

    if args.create_variant_filter_ht:
        logger.info("Creating variant filter Table...")
        res = resources.create_variant_filter_ht
        res.check_resource_existence()

        sites_ht = res.sites_ht.ht()

        # Auto-enable extra padding when non-default padding is specified.
        include_extra_padding = args.include_extra_padding or (
            args.exon_acceptor_padding != DEFAULT_EXON_UPSTREAM_PADDING
            or args.exon_donor_padding != DEFAULT_EXON_DOWNSTREAM_PADDING
        )

        logger.info(
            f"""
            create-variant-filter-ht parameters:

                Least consequence: {least_consequence}
                Max freq: {max_freq}
                Extra padding: {include_extra_padding} (acceptor={args.exon_acceptor_padding}, donor={args.exon_donor_padding})
                ClinVar categories: {args.include_clinvar_categories}
                Pathogenic splice: {args.include_pathogenic_splice} (min_splice_ai={args.min_splice_ai}, min_pangolin={args.min_pangolin})
                HC LoF: {args.include_hc_lof}
                In-trans-OE candidates: {args.include_in_trans_oe_candidates} (max_af={args.in_trans_oe_max_af})
                In-trans-OE intronic padding: {args.include_in_trans_oe_intronic_padding} (acceptor={args.in_trans_oe_acceptor_padding}, donor={args.in_trans_oe_donor_padding})
            """
        )

        # GENCODE is only needed for the exon-flanking intronic padding sources.
        gencode_ht = None
        if include_extra_padding or args.include_in_trans_oe_intronic_padding:
            gencode_ht = gencode.versions["v39"].ht()

        ht = create_variant_filter_ht(
            sites_ht,
            least_consequence=least_consequence,
            max_freq=max_freq,
            include_extra_padding=include_extra_padding,
            acceptor_padding=args.exon_acceptor_padding,
            donor_padding=args.exon_donor_padding,
            include_clinvar_categories=args.include_clinvar_categories,
            include_pathogenic_splice=args.include_pathogenic_splice,
            include_hc_lof=args.include_hc_lof,
            min_splice_ai=args.min_splice_ai,
            min_pangolin=args.min_pangolin,
            include_in_trans_oe_candidates=args.include_in_trans_oe_candidates,
            include_in_trans_oe_intronic_padding=args.include_in_trans_oe_intronic_padding,
            in_trans_oe_max_af=args.in_trans_oe_max_af,
            in_trans_oe_acceptor_padding=args.in_trans_oe_acceptor_padding,
            in_trans_oe_donor_padding=args.in_trans_oe_donor_padding,
            gencode_ht=gencode_ht,
            use_cache=args.use_region_interval_cache,
        )

        ht = ht.checkpoint(res.variant_filter_ht.path, overwrite=overwrite)
        logger.info("Number of variants in the variant filter Table: %d", ht.count())

    if args.filter_vmt:
        logger.info(f"Filtering gnomAD v4 {data_type} variant data MatrixTable...")
        res = resources.filter_vmt

        if test_chrom:
            logger.info(
                "Single-chromosome test mode: restricting --filter-vmt to %s.",
                test_chrom,
            )
            chrom_interval = hl.parse_locus_interval(
                test_chrom, reference_genome="GRCh38"
            )
            filter_intervals = [chrom_interval]
            out_path = (
                f"{DEFAULT_TMP_DIR}/exomes.filtered_vmt.{test_chrom}_test.mt"
            )
        elif test:
            filter_intervals = list(test_intervals.values())
            out_path = res.filtered_vmt.path
            res.check_resource_existence()
        else:
            filter_intervals = None
            out_path = res.filtered_vmt.path
            res.check_resource_existence()

        vp_release_only = args.vp_release_only
        logger.info("filter-vmt parameters: vp_release_only=%s", vp_release_only)
        vds = get_vds_func(
            release_only=vp_release_only,
            high_quality_only=not vp_release_only,
            split=True,
            filter_intervals=filter_intervals,
            filter_variant_ht=res.variant_filter_ht.ht(),
            entries_to_keep=["GT"],
            split_reference_blocks=False,
        )
        vds.variant_data.write(out_path, overwrite=overwrite)
        logger.info("The filtered VDS has been written to %s", out_path)

    if args.create_variant_pair_list_ht:
        logger.info("Creating variant pair list Table...")
        logger.info(
            "create-variant-pair-list-ht parameters: drop_oe_only_pairs=%s",
            args.include_in_trans_oe_candidates,
        )
        res = resources.create_variant_pair_list_ht

        if test_chrom:
            logger.info(
                "Single-chromosome test mode: reading chrom-filtered VMT "
                "+ restricting --create-variant-pair-list-ht to %s.",
                test_chrom,
            )
            interval = hl.parse_locus_interval(
                test_chrom, reference_genome="GRCh38"
            )
            chrom_vmt_path = (
                f"{DEFAULT_TMP_DIR}/exomes.filtered_vmt.{test_chrom}_test.mt"
            )
            mt = hl.read_matrix_table(chrom_vmt_path)
            filter_ht = hl.filter_intervals(
                res.variant_filter_ht.ht(), [interval]
            )
            out_path = (
                f"{DEFAULT_TMP_DIR}/exomes.variant_pairs.{test_chrom}_test.ht"
            )
        else:
            res.check_resource_existence()
            mt = res.filtered_vmt.mt()
            filter_ht = res.variant_filter_ht.ht()
            out_path = res.vp_list_ht.path

        ht = create_variant_pair_ht(
            mt, filter_ht,
            drop_oe_only_pairs=args.include_in_trans_oe_candidates,
        )

        ht = ht.checkpoint(out_path, overwrite=overwrite)
        logger.info(
            "The variant pair list Table has been written to %s.\n"
            "The number of unique variant pairs is %d",
            out_path, ht.count(),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tmp-dir",
        default=DEFAULT_TMP_DIR,
        help="Temporary directory for intermediate files.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to test intervals (all genes in TEST_INTERVALS) for testing.",
    )
    parser.add_argument(
        "--gene",
        choices=list(TEST_INTERVALS),
        help=(
            "Run on a single gene; uses that gene's interval from TEST_INTERVALS and "
            "implies --test."
        ),
    )
    parser.add_argument(
        "--test-genes",
        help=(
            "Comma-separated subset of TEST_INTERVALS keys (e.g. "
            "'SGCA,CAPN3,HFE'). Restricts the test intervals to just those "
            "genes; implies --test. Mutually exclusive with --gene."
        ),
    )
    parser.add_argument(
        "--interval",
        help=(
            "Explicit locus interval (e.g. 'chr19:1-58617616') for chromosome-"
            "scale scans. Bypasses the TEST_INTERVALS lookup; implies --test. "
            "Mutually exclusive with --gene and --test-genes."
        ),
    )
    parser.add_argument(
        "--output-postfix",
        help=(
            'Postfix to append to output file names (e.g., "pcnt_test" for files like '
            "exomes.vp_list.pcnt_test.ht)."
        ),
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Whether to overwrite existing files."
    )
    parser.add_argument(
        "--use-region-interval-cache",
        action="store_true",
        help=(
            "Opt in to caching the GENCODE intronic-padding interval HT "
            f"under {{hl.tmp_dir()}}/region_intervals/ (consulted by "
            "--create-variant-filter-ht). Off by default to guard against "
            "stale-cache reuse; rebuild takes ~2 min."
        ),
    )
    parser.add_argument(
        "--overwrite-cache",
        action="store_true",
        help=(
            "Force recomputation of cached intermediates inside "
            "genotype_count_intermediates.{postfix}/ (gt_info, var_idx, "
            "vp_map_flat, etc.). Without this flag, those intermediates are "
            "reused if already present, which is wrong when the upstream "
            "dense MT or pair list has been regenerated."
        ),
    )
    parser.add_argument(
        "--data-type",
        default=DEFAULT_DATA_TYPE,
        choices=DATA_TYPE_CHOICES,
        help=(
            f'Data type to use. Must be one of {", ".join(DATA_TYPE_CHOICES)}. Default '
            f"is {DEFAULT_DATA_TYPE}.",
        ),
    )
    parser.add_argument(
        "--preprocess-sites-ht",
        action="store_true",
        help=(
            "Assemble the per-variant sites HT (joins filter / freq / "
            "vep / all-sites-AN / SpliceAI / Pangolin / ClinVar into one "
            "wide Table). Required upstream of --create-variant-filter-ht."
        ),
    )
    parser.add_argument(
        "--create-variant-filter-ht",
        action="store_true",
        help=(
            "Create the variant filter Table from the assembled sites HT. "
            "Run --preprocess-sites-ht first."
        ),
    )
    parser.add_argument(
        "--least-consequence",
        default=DEFAULT_LEAST_CONSEQUENCE,
        choices=CSQ_ORDER,
        help=(
            "Lowest-severity consequence to keep. Default "
            f"is {DEFAULT_LEAST_CONSEQUENCE}."
        ),
    )
    parser.add_argument(
        "--max-freq",
        type=float,
        default=DEFAULT_MAX_FREQ,
        help=f"Maximum global AF to keep (inclusive). Default is {DEFAULT_MAX_FREQ}.",
    )
    parser.add_argument(
        "--min-an-pct",
        type=int,
        default=0,
        help=(
            "Exclusive AN_percent floor required on both sides of a pair, "
            "applied when a pair list is consumed (dense-MT / encode / "
            "count steps). Default 0 drops only uncallable (an_pct==0) "
            "endpoints; raise to e.g. 10 for a stricter floor. Pass a "
            "negative value to keep every pair. Must be passed consistently "
            "across the steps that read the pair list."
        ),
    )
    parser.add_argument(
        "--include-extra-padding",
        action="store_true",
        help=(
            "Include variants in GENCODE exon flanking intronic zones beyond "
            "VEP's built-in splice region (±8bp). Auto-enabled when "
            "--exon-acceptor-padding or --exon-donor-padding exceed defaults. "
            "Used with --create-variant-filter-ht."
        ),
    )
    parser.add_argument(
        "--exon-acceptor-padding",
        type=int,
        default=DEFAULT_EXON_UPSTREAM_PADDING,
        help=(
            "Bp to pad on the acceptor (intron) side of exon boundaries. "
            "Strand-aware: applied before exon start on + strand, after exon "
            f"end on - strand. Default is {DEFAULT_EXON_UPSTREAM_PADDING}."
        ),
    )
    parser.add_argument(
        "--exon-donor-padding",
        type=int,
        default=DEFAULT_EXON_DOWNSTREAM_PADDING,
        help=(
            "Bp to pad on the donor (intron) side of exon boundaries. "
            "Strand-aware: applied after exon end on + strand, before exon "
            f"start on - strand. Default is {DEFAULT_EXON_DOWNSTREAM_PADDING}."
        ),
    )
    parser.add_argument(
        "--include-clinvar-categories",
        type=lambda s: [x.strip().lower() for x in s.split(",") if x.strip()],
        default=[],
        help=(
            "Comma-separated ClinVar significance categories to union into "
            "the variant filter HT. Choices: plp (Pathogenic/Likely "
            "pathogenic), blb (Benign/Likely benign), vus (Uncertain "
            "significance). Emits per-category source tags (clinvar_plp, "
            "clinvar_blb, clinvar_vus). Each category is also a valid "
            "--partner-set value in run_in_trans_oe.py. Example: "
            "--include-clinvar-categories plp,blb,vus."
        ),
    )
    parser.add_argument(
        "--include-pathogenic-splice",
        action="store_true",
        help=(
            "Include variants with SpliceAI Δ > --min-splice-ai OR Pangolin "
            "Δ > --min-pangolin in the variant filter HT (tag splice_path)."
        ),
    )
    parser.add_argument(
        "--include-hc-lof",
        action="store_true",
        help=(
            "Include variants with VEP lof == 'HC' on a protein-coding "
            "Ensembl transcript in the variant filter HT (tag hc_lof)."
        ),
    )
    parser.add_argument(
        "--min-splice-ai",
        type=float,
        default=DEFAULT_MIN_SPLICE_AI,
        help=(
            f"Minimum spliceAI delta score for --include-noncoding-pathogenic. "
            f"Default is {DEFAULT_MIN_SPLICE_AI}."
        ),
    )
    parser.add_argument(
        "--min-pangolin",
        type=float,
        default=DEFAULT_MIN_PANGOLIN,
        help=(
            f"Minimum pangolin delta score for --include-noncoding-pathogenic. "
            f"Default is {DEFAULT_MIN_PANGOLIN}."
        ),
    )
    parser.add_argument(
        "--include-in-trans-oe-candidates",
        action="store_true",
        help=(
            "Union higher-AF (up to --in-trans-oe-max-af) candidates into "
            "the variant filter HT (tag in_trans_oe_candidate). When set, "
            "the variant pair list is also filtered to keep only pairs "
            "where ≥1 side has a non-OE-candidate source tag (avoids "
            "candidate × candidate explosion). Partner sets are not "
            "auto-included anymore — add them explicitly via "
            "--include-clinvar-categories / --include-pathogenic-splice / "
            "--include-hc-lof."
        ),
    )
    parser.add_argument(
        "--in-trans-oe-max-af",
        type=float,
        default=DEFAULT_IN_TRANS_OE_MAX_AF,
        help=(
            "Upper AF bound (inclusive) for in-trans-OE candidates. Default "
            f"{DEFAULT_IN_TRANS_OE_MAX_AF}; deliberately above the standard "
            "pipeline's 5%% cap so common-but-suspect variants enter the "
            "analysis when paired with a P/LP."
        ),
    )
    parser.add_argument(
        "--include-in-trans-oe-intronic-padding",
        action="store_true",
        help=(
            "Union intronic-padding variants (branch-point region + cryptic "
            "5' splice signals, see --in-trans-oe-acceptor-padding / "
            "--in-trans-oe-donor-padding) into the variant filter HT (tag "
            "in_trans_oe_intronic_padding). Variants here are usually "
            "VEP-classified intron_variant and would otherwise be invisible "
            "to consequence-based filters."
        ),
    )
    parser.add_argument(
        "--in-trans-oe-acceptor-padding",
        type=int,
        default=50,
        help=(
            "Bp upstream of the acceptor splice site to include for "
            "--in-trans-oe-include-intronic-padding. Default 50 (covers the "
            "branch-point region, typically -18 to -40)."
        ),
    )
    parser.add_argument(
        "--in-trans-oe-donor-padding",
        type=int,
        default=15,
        help=(
            "Bp downstream of the donor splice site to include for "
            "--in-trans-oe-include-intronic-padding. Default 15 (covers "
            "cryptic 5' splice signals just past VEP's +8 splice region)."
        ),
    )
    parser.add_argument(
        "--filter-vmt",
        action="store_true",
        help="Filter the MatrixTable for determining variant pairs.",
    )
    parser.add_argument(
        "--vp-release-only",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Use release-only samples for variant pair discovery (--filter-vmt). "
            "When False (default), uses high-quality samples instead."
        ),
    )
    parser.add_argument(
        "--create-variant-pair-list-ht",
        action="store_true",
        help="first create just the list of possible variant pairs.",
    )
    parser.add_argument(
        "--test-chrom",
        nargs="?",
        const="chr19",
        default=None,
        help=(
            "Restrict --filter-vmt and --create-variant-pair-list-ht to "
            "a single chromosome (useful for full-genome runtime / cost "
            "estimation). Pass without a value to default to chr19, or "
            "specify e.g. '--test-chrom 5'. Each step reads its inputs "
            "from production paths and writes its output to a "
            "chrom-postfixed path under DEFAULT_TMP_DIR — production "
            "locations are untouched."
        ),
    )
 
    args = parser.parse_args()
    main(args)

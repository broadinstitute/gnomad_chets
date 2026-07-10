"""
Compute observed-vs-expected in-trans co-occurrence for candidate variants
against curated partner sets (ClinVar P/LP, predicted-damaging) within the
same gene.

For a candidate variant ``c`` in gene ``g`` with partner set ``S``, this module
produces, for each candidate × gene × partner_set:

* ``E_total`` = sum over partners ``i`` in ``S`` of
  ``E_i = 2 · N_pair · AF_c · (1 - AF_c) · AF_i · (1 - AF_i)`` (expected
  count of in-trans compound heterozygotes under HWE + independence between
  loci). ``N_pair = sum(gt_counts)`` is the per-pair joint callable count
  emitted by the gt_counts pipeline (the encoder tracks ``n_callable`` per
  variant so the AABB cell correctly excludes samples with no GT entry).
  The ``(1 - AF)`` factors restrict the count to the "true in-trans"
  configuration (one haplotype carries the candidate, the other carries the
  partner, neither carries both); they are ≈ 1 in the rare-variant regime
  but matter for higher-AF candidates entering via the in-trans-OE candidate
  sources (AF range up to 0.5).
* ``O_total`` = sum over partners ``i`` of ``O_i = p_chet_i · double_carriers_i``
  where ``p_chet_i`` is the EM-derived Ptrans probability for the pair and
  ``double_carriers_i = gt_counts_adj[4]`` (the AaBb cell of the 9-element
  count array, ordered ``[AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb,
  aabb]``).
* ``poisson_lower_tail_p`` = ``P(X <= O_total | mean = E_total)`` — one-sided
  depletion p-value.

Pipeline status: skeleton — function signatures, schema, and per-pair math
are in place; per-candidate aggregation and the Poisson p-value step are
stubs (``NotImplementedError``) pending the GraphQL contract being locked
on the browser side.

Reuses, does not reimplement:

* :func:`gnomad_chets.v4.create_vp_list.create_variant_filter_ht` for
  partner-set membership (variants tagged with ``clinvar_<cat>`` /
  ``splice_path`` / ``hc_lof`` per the corresponding ``--include-*`` flag).
* :func:`filter_clinvar_by_category` (local helper, TODO upstream) for the
  ClinVar significance categorical filter.
* ``hl.experimental.haplotype_freq_em`` for the EM (mirrors
  :func:`gnomad_chets.v2.phasing.get_em_expr`).
* :func:`gnomad_chets.v4.resources.get_variant_pair_genotype_counts_ht` for
  the input pair-counts Table.
"""

import logging
from typing import Optional

import hail as hl
from scipy.stats import poisson
from gnomad.utils.filtering import filter_to_clinvar_pathogenic
from gnomad.utils.vep import filter_vep_transcript_csqs_expr

from gnomad_chets.v4.resources import (
    DEFAULT_MAX_FREQ,
    DEFAULT_MIN_PANGOLIN,
    DEFAULT_MIN_SPLICE_AI,
)

logger = logging.getLogger("in_trans_oe")
logger.setLevel(logging.INFO)


########################################################################################
### Partner-set string constants (mirror the GraphQL InTransPartnerSet enum)
########################################################################################

PARTNER_SET_CLINVAR_PLP = "CLINVAR_PLP"
PARTNER_SET_CLINVAR_BLB = "CLINVAR_BLB"
PARTNER_SET_CLINVAR_VUS = "CLINVAR_VUS"
PARTNER_SET_PREDICTED_DAMAGING = "PREDICTED_DAMAGING"
PARTNER_SET_CLINVAR_PLP_OR_PREDICTED = "CLINVAR_PLP_OR_PREDICTED"

PARTNER_SETS = (
    PARTNER_SET_CLINVAR_PLP,
    PARTNER_SET_CLINVAR_BLB,
    PARTNER_SET_CLINVAR_VUS,
    PARTNER_SET_PREDICTED_DAMAGING,
    PARTNER_SET_CLINVAR_PLP_OR_PREDICTED,
)

NULL_MODEL_POISSON_HWE = (
    "Poisson(E_total) under HWE + independence between candidate and partner "
    "loci; E_i = 2 * N_pair * AF_c * (1 - AF_c) * AF_i * (1 - AF_i), summed "
    "over partners. N_pair is the per-pair joint callable sample count "
    "(samples with non-missing GT at both variants), which corrects for "
    "callable-coverage differences between low-coverage positions (e.g. 3' "
    "UTR) and the rest of the gene."
)

NULL_MODEL_POISSON_HWE_LD_ADJUSTED = (
    "E_ld = N * (1 - Prod_i (1 - 2 * MAF_A * MAF_Bi * (1 - D*_i))), "
    "where D*_i = 1 if |D'_i| >= tau else 0 (tau = ld_threshold), "
    "D_i = MAF_{A,Bi} - MAF_A * MAF_Bi from the gt_counts haplotype-count "
    "estimator, and D'_i = D_i / D_max sign-branched. "
    "poisson_lower_tail_p_ld_adjusted = P(X <= O | X ~ Poisson(E_ld)). "
    "Unobserved partners default to D* = 0 (no LD correction)."
)


########################################################################################
### Partner-set construction (Phase 1 = ClinVar P/LP only)
########################################################################################

def create_clinvar_plp_partner_ht(
    vep_ht: hl.Table,
    freq_ht: hl.Table,
    filter_ht: hl.Table,
    clinvar_path_ht: hl.Table,
    max_freq: float = DEFAULT_MAX_FREQ,
) -> hl.Table:
    """
    Build the ClinVar P/LP partner-set Table for in-trans depletion.

    A variant is a ClinVar P/LP partner candidate when:

    * it appears in ``clinvar_path_ht`` (already filtered to P/LP via
      :func:`gnomad.utils.filtering.filter_to_clinvar_pathogenic`),
    * it passes variant QC (``filter_ht.filters`` empty),
    * it has a global AF in ``(0, max_freq]``,
    * it has at least one protein-coding Ensembl VEP transcript consequence
      (gene IDs are taken from those transcripts).

    Output is a Table keyed by ``(locus, alleles)`` with fields ``gene_id``
    (array of ENSG strings), ``af``, ``an`` (callable allele count, used to
    estimate per-pair sample count for cross-join expectations), and
    ``partner_set`` (the literal :data:`PARTNER_SET_CLINVAR_PLP`).

    :param vep_ht: VEP Table for gene_id assignment.
    :param freq_ht: Frequency Table.
    :param filter_ht: Final filter Table for QC pass.
    :param clinvar_path_ht: ClinVar Table already filtered to P/LP.
    :param max_freq: Maximum global AF to keep (inclusive).
    :return: Partner-set Table keyed by ``(locus, alleles)``.
    """
    ht = vep_ht.annotate(
        filters=filter_ht[vep_ht.locus, vep_ht.alleles].filters,
        af=freq_ht[vep_ht.locus, vep_ht.alleles].freq[0].AF,
        an=freq_ht[vep_ht.locus, vep_ht.alleles].freq[0].AN,
        gene_id=hl.array(hl.set(filter_vep_transcript_csqs_expr(
            vep_ht.vep.transcript_consequences,
            protein_coding=True,
            ensembl_only=True,
        ).map(lambda csq: csq.gene_id))),
        _is_clinvar_path=hl.is_defined(
            clinvar_path_ht[vep_ht.locus, vep_ht.alleles]
        ),
    )

    ht = ht.filter(
        (ht.filters.length() == 0)
        & (ht.af > 0)
        & (ht.af <= max_freq)
        & hl.is_defined(ht.gene_id)
        & (hl.len(ht.gene_id) > 0)
        & ht._is_clinvar_path
    )

    return ht.select(
        gene_id=ht.gene_id,
        af=ht.af,
        an=ht.an,
        partner_set=hl.literal(PARTNER_SET_CLINVAR_PLP),
    )


def create_predicted_damaging_partner_ht(
    vep_ht: hl.Table,
    freq_ht: hl.Table,
    filter_ht: hl.Table,
    spliceai_ht: hl.Table,
    pangolin_ht: hl.Table,
    max_freq: float = DEFAULT_MAX_FREQ,
    min_splice_ai: float = DEFAULT_MIN_SPLICE_AI,
    min_pangolin: float = DEFAULT_MIN_PANGOLIN,
    min_revel: Optional[float] = None,
) -> hl.Table:
    """
    Build the predicted-damaging partner-set Table (Phase 2).

    Includes high-confidence LoF + REVEL/SpliceAI/Pangolin thresholds.
    Schema matches :func:`create_clinvar_plp_partner_ht`; ``partner_set``
    is :data:`PARTNER_SET_PREDICTED_DAMAGING`.

    :param vep_ht: VEP Table.
    :param freq_ht: Frequency Table.
    :param filter_ht: Final filter Table.
    :param spliceai_ht: SpliceAI in silico predictor Table.
    :param pangolin_ht: Pangolin in silico predictor Table.
    :param max_freq: Maximum global AF.
    :param min_splice_ai: Minimum SpliceAI delta score.
    :param min_pangolin: Minimum Pangolin delta score.
    :param min_revel: Optional minimum REVEL score for missense partners.
    :return: Partner-set Table keyed by ``(locus, alleles)``.
    """
    raise NotImplementedError("Skeleton — Phase 2.")


########################################################################################
### Per-pair O/E annotations (pure transform on the existing pair-counts Table)
########################################################################################

def annotate_pair_oe_terms(
    vp_gt_counts_ht: hl.Table,
    freq_ht: hl.Table,
    n_samples: int,
    use_adj: bool = True,
    phased_ht: Optional[hl.Table] = None,
    include_ld_adjusted: bool = False,
    ld_threshold: float = 0.5,
) -> hl.Table:
    """
    Annotate the variant-pair genotype-counts Table with per-pair O and E.

    For each pair (v1, v2), adds:

    * ``af1``, ``af2``: global allele frequencies pulled from ``freq_ht``.
    * ``double_carriers``: ``gt_counts[4]`` (the AaBb cell — number of
      individuals heterozygous for both v1 and v2).
    * ``p_chet``: Ptrans, the EM-derived probability that an observed
      double-het is in trans. Computed as
      ``(h01 * h10) / (h00 * h11 + h01 * h10)`` from
      ``hl.experimental.haplotype_freq_em(gt_counts)`` — the same expression
      used by :func:`gnomad_chets.v2.phasing.get_em_expr`. When ``phased_ht``
      is supplied, ``p_chet`` is read from
      ``phased_ht.em.{adj|raw}.p_chet`` instead (skipping the EM recompute).
    * ``n_pair``: per-pair callable sample count, ``sum(gt_counts)`` -- the
      number of individuals with non-missing genotype for both v1 and v2.
      The gt_counts pipeline tracks per-variant ``n_callable`` in the encoder
      so the AABB cell correctly excludes missing-GT samples (samples with
      no entry at the variant); ``sum(gt_counts)`` therefore equals the joint
      callable count, not the global cohort size. Falls back to ``n_samples``
      if ``sum(gt_counts) == 0``.
    * ``e_pair``: ``2 * n_pair * af1 * (1 - af1) * af2 * (1 - af2)``
      (expected in-trans count under HWE + independence -- the "true in-trans"
      configuration where one haplotype carries v1 only and the other carries
      v2 only). Using the per-pair callable ``n_pair`` corrects for
      callable-coverage differences at low-coverage positions (e.g. 3' UTR)
      vs the rest of the gene. The ``(1 - AF)`` factors round to 1 in the
      rare-variant regime but matter for higher-AF candidates from the
      in-trans-OE candidate sources (AF up to 0.5).
    * ``o_pair``: ``p_chet * double_carriers`` (observed in-trans count).

    The ``populations`` array on the input is preserved if present; per-pop
    O/E annotations will be added by a follow-up function once the v4
    genotype-counts table carries ancestry-stratified counts (Phase 2.5).

    :param vp_gt_counts_ht: Output of
        :func:`gnomad_chets.v4.compute_vp_counts.create_variant_pair_genotype_counts`.
        Keyed by ``(locus1, alleles1, locus2, alleles2)`` with
        ``gt_counts_raw`` and ``gt_counts_adj`` 9-element arrays.
    :param freq_ht: Frequency Table for AF lookup. Same data type
        (exomes/genomes) as the pair-counts Table.
    :param n_samples: Fallback sample count used to scale ``e_pair`` when
        ``sum(gt_counts) == 0`` (rare; would only happen for malformed pair
        rows). The per-pair ``n_pair`` from ``sum(gt_counts)`` is preferred.
        Should match the sample count behind ``gt_counts_adj`` (or ``raw`` if
        ``use_adj=False``).
    :param use_adj: If ``True`` (default), use ``gt_counts_adj``; otherwise
        ``gt_counts_raw``.
    :param phased_ht: Optional pre-computed phased Table (output of
        :func:`gnomad_chets.v4.phase_gnomad.get_phased_gnomad_ht`), keyed by
        the same 4-tuple ``(locus1, alleles1, locus2, alleles2)`` with
        ``em.{adj,raw}.p_chet``. When supplied, EM is skipped and
        ``p_chet`` is looked up from ``phased_ht`` (partition-local indexed
        read). Same result, avoids redundant EM in per-gene OE calls.
    :param include_ld_adjusted: If ``True``, adds LD-adjusted expected-in-trans
        annotations following Rachel Unger's formulation: ``maf1``, ``maf2``,
        ``n_hap_ab``, ``n_chrom``, ``maf_ab_hap``, ``d_ab``, ``d_max``,
        ``d_prime``, ``d_star``, ``ld_term``, ``log_ld_term``. These are
        consumed by :func:`aggregate_oe_per_candidate` to compute
        ``total_expected_in_trans_ld_adjusted`` via the product-over-partners
        form ``N * (1 - exp(sum(log(ld_term))))``. Default ``False`` (schema
        unchanged when off).
    :param ld_threshold: |D'| threshold (tau) above which a candidate/partner
        pair is treated as fully in-cis (D* = 1 -> partner contributes 0 to
        expected-in-trans). Only used when ``include_ld_adjusted=True``.
        Must be in [0, 1].
    :return: Annotated pair-counts Table.
    """
    if not (0 <= ld_threshold <= 1):
        raise ValueError(
            f"ld_threshold must be in [0, 1], got {ld_threshold}"
        )
    gt_counts = (
        vp_gt_counts_ht.gt_counts_adj if use_adj else vp_gt_counts_ht.gt_counts_raw
    )
    # The per-sample-grouping pipeline emits int64 counts; haplotype_freq_em
    # requires int32. Cast (counts are individual genotypes, well below 2^31).
    gt_counts = gt_counts.map(hl.int32)
    # freq[0] = adj subset, freq[1] = raw subset (per gnomad_qc convention).
    freq_idx = 0 if use_adj else 1
    v1_freq = freq_ht[vp_gt_counts_ht.locus1, vp_gt_counts_ht.alleles1].freq[freq_idx]
    v2_freq = freq_ht[vp_gt_counts_ht.locus2, vp_gt_counts_ht.alleles2].freq[freq_idx]
    af1 = v1_freq.AF
    af2 = v2_freq.AF

    if phased_ht is not None:
        # Read p_chet from the pre-computed phased HT (same 4-tuple key as
        # gt-counts). Partition-local indexed lookup, no shuffle; avoids
        # re-running EM in every per-gene OE call.
        em_field = phased_ht[
            vp_gt_counts_ht.locus1, vp_gt_counts_ht.alleles1,
            vp_gt_counts_ht.locus2, vp_gt_counts_ht.alleles2,
        ].em
        p_chet = em_field.adj.p_chet if use_adj else em_field.raw.p_chet
    else:
        hap_counts = hl.experimental.haplotype_freq_em(gt_counts)
        p_chet_raw = (hap_counts[1] * hap_counts[2]) / (
            hap_counts[0] * hap_counts[3] + hap_counts[1] * hap_counts[2]
        )
        # The EM denominator can collapse to 0 (all zero double-hets, or extreme
        # haplotype configurations), producing NaN p_chet. Pairs with no double
        # carriers can't contribute to observed in-trans regardless, so coerce
        # NaN to missing -- the o_pair guard below sets the contribution to 0.
        p_chet = hl.if_else(hl.is_nan(p_chet_raw), hl.missing(hl.tfloat64), p_chet_raw)

    double_carriers = gt_counts[4]
    # Per-pair callable sample count: ``sum(gt_counts)`` = samples
    # adj-PASS at both variants. The encoder tracks ``entries_set`` /
    # ``n_callable`` per variant so the AABB cell correctly excludes
    # samples with no entry at either side; ``sum(gt_counts)`` therefore
    # equals the joint adj-PASS-callable count. Falls back to
    # ``n_samples`` only if ``sum(gt_counts) == 0`` (malformed pair rows).
    gt_sum = hl.float64(hl.sum(gt_counts))
    n_pair = hl.if_else(gt_sum > 0, gt_sum, hl.float64(n_samples))

    # Exact in-trans expectation under HWE + independence between loci:
    # 2 * N_pair * AF_c * (1 - AF_c) * AF_p * (1 - AF_p), counting only the
    # "true in-trans" configuration (one haplotype carries the candidate, the
    # other carries the partner, neither carries both). The (1 - AF) factors
    # absorb candidate-hom and cis double-carrier individuals, who can't
    # contribute in-trans pairs. They round to 1 for the rare-variant pipeline
    # default (AF cap 5%) but matter for higher-AF candidates flagged via the
    # in-trans-OE candidate sources (AF range up to 0.5).
    e_pair = 2 * n_pair * af1 * (1 - af1) * af2 * (1 - af2)
    # Pairs with no double-het carriers contribute 0 to observed in-trans
    # regardless of (possibly missing) p_chet. Guarding here prevents NaN
    # propagation through the per-candidate sum aggregation downstream.
    o_pair = hl.if_else(
        (double_carriers > 0) & hl.is_defined(p_chet),
        p_chet * double_carriers,
        0.0,
    )

    annotations = dict(
        af1=af1,
        af2=af2,
        n_pair=n_pair,
        double_carriers=double_carriers,
        p_chet=p_chet,
        e_pair=e_pair,
        o_pair=o_pair,
    )

    if include_ld_adjusted:
        # LD-adjusted expected-in-trans following Rachel Unger's formulation.
        # Adds per-pair fields consumed by aggregate_oe_per_candidate to
        # compute the product-form E_ld = N * (1 - Prod_i (1 - 2 * MAF_A *
        # MAF_Bi * (1 - D*_i))) via a log-sum trick.
        # Defensive against ALT-AF > 0.5 (shouldn't happen after max_freq cap
        # but cheap; matches Rachel's MAF-convention formulas).
        maf1 = hl.min(af1, 1 - af1)
        maf2 = hl.min(af2, 1 - af2)

        # Joint haplotype MAF from gt_counts (orientation-invariant: indices
        # 4/5/7 are symmetric under v1<->v2 swap, and index 8 is trivially
        # symmetric). n_{A,B} = 2*gt[8] + gt[7] + gt[5] + 0.5*gt[4] --
        # moment-of-methods estimator that doesn't need EM; uses uniform
        # 0.5 prior for the double-het's two configurations.
        n_hap_ab = (
            hl.float64(2 * gt_counts[8])
            + hl.float64(gt_counts[7])
            + hl.float64(gt_counts[5])
            + 0.5 * hl.float64(gt_counts[4])
        )
        n_chrom = 2.0 * n_pair
        maf_ab_hap = hl.if_else(
            n_chrom > 0, n_hap_ab / n_chrom, hl.missing(hl.tfloat64)
        )

        # Composite D and sign-branched D_max.
        d_ab = maf_ab_hap - maf1 * maf2
        d_max_pos = hl.min(maf1 * (1 - maf2), (1 - maf1) * maf2)
        d_max_neg = hl.min(maf1 * maf2, (1 - maf1) * (1 - maf2))
        d_max = hl.if_else(d_ab >= 0, d_max_pos, d_max_neg)

        # D' clipped to [-1, 1]; missing if D_max <= 0 (either MAF is 0 or 1).
        d_prime_raw = hl.if_else(
            d_max > 0, d_ab / d_max, hl.missing(hl.tfloat64)
        )
        d_prime = hl.if_else(
            hl.is_defined(d_prime_raw),
            hl.max(-1.0, hl.min(1.0, d_prime_raw)),
            hl.missing(hl.tfloat64),
        )
        # D* = 1 iff |D'| >= tau; 0 otherwise (including missing D_prime).
        d_star = hl.if_else(
            hl.is_defined(d_prime) & (hl.abs(d_prime) >= ld_threshold),
            1.0, 0.0,
        )

        # Per-pair LD-adjusted expected-in-trans term. Guard the log:
        # cap to [eps, 1]. raw_term hits 1 exactly when MAF=0 or D*=1 (then
        # log(1)=0 -> partner contributes 0). raw_term < 0 requires MAF_A *
        # MAF_B > 0.5 which is impossible under the max_freq <= 0.5 cap;
        # clamp anyway.
        ld_term_eps = 1e-300
        raw_term = 1.0 - 2.0 * maf1 * maf2 * (1.0 - d_star)
        ld_term = hl.max(hl.min(raw_term, 1.0), ld_term_eps)
        log_ld_term = hl.log(ld_term)

        annotations.update(
            maf1=maf1,
            maf2=maf2,
            n_hap_ab=n_hap_ab,
            n_chrom=n_chrom,
            maf_ab_hap=maf_ab_hap,
            d_ab=d_ab,
            d_max=d_max,
            d_prime=d_prime,
            d_star=d_star,
            ld_term=ld_term,
            log_ld_term=log_ld_term,
        )

    return vp_gt_counts_ht.annotate(**annotations)


########################################################################################
### Aggregation: per-pair → per-(candidate, gene, partner_set)
########################################################################################

def aggregate_oe_per_candidate(
    annotated_pair_ht: hl.Table,
    partner_ht: hl.Table,
    candidate_ht: hl.Table,
    freq_ht: hl.Table,
    partner_set: str,
    n_samples: int,
    include_ld_adjusted: bool = False,
    ld_threshold: float = 0.5,
) -> hl.Table:
    """
    Roll up per-pair O/E to per-(candidate, gene, partner_set), including
    cross-joined unobserved (candidate, partner) pairs.

    For each candidate × gene combination, enumerates every P/LP partner in
    that gene (from ``partner_ht``) and joins against the observed per-pair
    O/E table. Pairs without a matching observed entry — i.e. no individual
    is a co-carrier of both variants — get ``O = 0`` and a synthetic
    ``E = 2 * N_pair_est * AF_c * (1 - AF_c) * AF_p * (1 - AF_p)``, where
    ``N_pair_est = min(AN_c, AN_p) / 2`` (the upper bound on jointly-callable
    samples; falls back to ``n_samples`` if AN is missing on either side).

    This corrects the conceptual gap in the prior version, which only summed
    E over partners with observed pairs and therefore underestimated total E
    (and inflated p-values) for low-power candidates with rare partners.

    Self-pairs (``v1 == v2``) are excluded.

    Output schema (one row per (candidate locus, candidate alleles, gene_id,
    partner_set)):

    .. code-block::

        key: (locus, alleles, gene_id, partner_set)
        candidate_af: float
        n_partners: int                             # all P/LP partners in the gene
        total_expected_in_trans: float              # sum(e_pair) over all partners
        total_observed_in_trans: float              # sum(o_pair) over all partners
        partners: array<struct{
            locus, alleles, af, double_carriers, p_chet,
            expected_in_trans, observed_in_trans
        }>
        populations: array<struct{...}>             # empty in Phase 2.0
        warnings: array<str>

    :param annotated_pair_ht: Output of :func:`annotate_pair_oe_terms`.
    :param partner_ht: Partner-set Table (e.g. from
        :func:`create_clinvar_plp_partner_ht`). Must have ``af`` and ``an``
        fields for the per-pair E synthesis.
    :param candidate_ht: Table of candidate variants. For Phase 2.0 this is
        the output of
        :func:`gnomad_chets.v4.create_vp_list.create_variant_filter_ht`.
    :param freq_ht: Frequency Table for candidate AF/AN lookup.
    :param partner_set: One of :data:`PARTNER_SETS`.
    :param n_samples: Fallback sample count if AN is unavailable on either
        side of an unobserved pair. Should match the cohort behind the gt
        counts (e.g. 730947 for v4.1 exomes).
    :param include_ld_adjusted: If ``True``, also emits
        ``total_expected_in_trans_ld_adjusted`` (product form over partners
        with a D' threshold), ``n_ld_adjusted_used`` (count of partners with
        D* = 1 -- LD-linked and thus zeroed out), and adds ``d_prime``,
        ``d_star``, ``ld_term`` to each ``partners`` struct. ``ld_threshold``
        is stamped as a global. Requires that ``annotated_pair_ht`` was
        produced by :func:`annotate_pair_oe_terms` with
        ``include_ld_adjusted=True``. Default ``False``.
    :param ld_threshold: |D'| threshold (tau) stamped as a global on the
        output when ``include_ld_adjusted=True``. Must match the value
        passed to :func:`annotate_pair_oe_terms`.
    :return: Aggregated per-candidate Table with the schema above.
    """
    if partner_set not in PARTNER_SETS:
        raise ValueError(f"partner_set {partner_set!r} not in {PARTNER_SETS}")
    if not (0 <= ld_threshold <= 1):
        raise ValueError(
            f"ld_threshold must be in [0, 1], got {ld_threshold}"
        )

    # ----- Step 1: Build the per-(candidate, partner, gene) OBSERVED view. -----
    # The pair-counts Table is keyed by (v1, v2) with v1 <= v2, so each pair
    # gets two views (candidate=v1 and candidate=v2). The union doubles the
    # row count temporarily but lets us use a single (candidate, partner)
    # lookup downstream.
    ht = annotated_pair_ht.key_by()
    ld_view_extras_a = {}
    ld_view_extras_b = {}
    if include_ld_adjusted:
        # d_prime / d_star / ld_term / log_ld_term are orientation-invariant
        # by construction (indices 4/5/7 in gt_counts are symmetric under
        # v1<->v2 swap, and MAF_A * MAF_B is commutative), so the same values
        # go into both views.
        ld_view_extras_a = dict(
            d_prime=ht.d_prime,
            d_star=ht.d_star,
            ld_term=ht.ld_term,
            log_ld_term=ht.log_ld_term,
        )
        ld_view_extras_b = dict(
            d_prime=ht.d_prime,
            d_star=ht.d_star,
            ld_term=ht.ld_term,
            log_ld_term=ht.log_ld_term,
        )
    view_a = ht.select(
        candidate_locus=ht.locus1,
        candidate_alleles=ht.alleles1,
        partner_locus=ht.locus2,
        partner_alleles=ht.alleles2,
        candidate_af=ht.af1,
        partner_af=ht.af2,
        double_carriers=ht.double_carriers,
        p_chet=ht.p_chet,
        e_pair=ht.e_pair,
        o_pair=ht.o_pair,
        **ld_view_extras_a,
    )
    view_b = ht.select(
        candidate_locus=ht.locus2,
        candidate_alleles=ht.alleles2,
        partner_locus=ht.locus1,
        partner_alleles=ht.alleles1,
        candidate_af=ht.af2,
        partner_af=ht.af1,
        double_carriers=ht.double_carriers,
        p_chet=ht.p_chet,
        e_pair=ht.e_pair,
        o_pair=ht.o_pair,
        **ld_view_extras_b,
    )
    observed_pairs = view_a.union(view_b)

    # Filter to (candidate in candidate_ht) and (partner in partner_ht), then
    # explode by shared gene_id.
    observed_pairs = observed_pairs.annotate(
        partner_gene_ids=partner_ht[observed_pairs.partner_locus, observed_pairs.partner_alleles].gene_id,
        candidate_gene_ids=candidate_ht[observed_pairs.candidate_locus, observed_pairs.candidate_alleles].gene_id,
    )
    observed_pairs = observed_pairs.filter(
        hl.is_defined(observed_pairs.partner_gene_ids)
        & hl.is_defined(observed_pairs.candidate_gene_ids)
    )
    observed_pairs = observed_pairs.annotate(
        shared_gene_ids=hl.array(
            hl.set(observed_pairs.candidate_gene_ids)
            .intersection(hl.set(observed_pairs.partner_gene_ids))
        )
    )
    observed_pairs = observed_pairs.filter(hl.len(observed_pairs.shared_gene_ids) > 0)
    observed_pairs = observed_pairs.filter(
        ~(
            (observed_pairs.candidate_locus == observed_pairs.partner_locus)
            & (observed_pairs.candidate_alleles == observed_pairs.partner_alleles)
        )
    )
    observed_pairs = (
        observed_pairs.explode("shared_gene_ids")
        .rename({"shared_gene_ids": "gene_id"})
    )
    observed_pairs = observed_pairs.key_by(
        "candidate_locus", "candidate_alleles",
        "partner_locus", "partner_alleles", "gene_id",
    )
    observed_pairs = observed_pairs.cache()

    # ----- Step 2: Build the (candidate × partner × gene) cross-product. -----
    # Index candidates by (variant, gene), enrich with AF/AN from freq_ht.
    cand_with_freq = candidate_ht.annotate(
        candidate_af=freq_ht[candidate_ht.locus, candidate_ht.alleles].freq[0].AF,
        candidate_an=freq_ht[candidate_ht.locus, candidate_ht.alleles].freq[0].AN,
    )
    cand_exp = cand_with_freq.key_by()
    cand_exp = cand_exp.select(
        candidate_locus=cand_exp.locus,
        candidate_alleles=cand_exp.alleles,
        candidate_af=cand_exp.candidate_af,
        candidate_an=cand_exp.candidate_an,
        gene_id=cand_exp.gene_id,
    )
    cand_exp = cand_exp.explode("gene_id")
    cand_exp = cand_exp.filter(hl.is_defined(cand_exp.candidate_af))

    # Group partners by gene_id so each gene gets an array of its partners.
    part_exp = partner_ht.key_by()
    part_exp = part_exp.select(
        partner_locus=part_exp.locus,
        partner_alleles=part_exp.alleles,
        partner_af=part_exp.af,
        partner_an=part_exp.an,
        gene_id=part_exp.gene_id,
    )
    part_exp = part_exp.explode("gene_id")
    partners_by_gene = part_exp.group_by("gene_id").aggregate(
        partner_list=hl.agg.collect(
            hl.struct(
                partner_locus=part_exp.partner_locus,
                partner_alleles=part_exp.partner_alleles,
                partner_af=part_exp.partner_af,
                partner_an=part_exp.partner_an,
            )
        )
    )

    # Attach the partner list to each (candidate, gene); explode to one row
    # per (candidate, partner, gene); exclude self-pairs.
    cand_exp = cand_exp.annotate(
        partner_list=partners_by_gene[cand_exp.gene_id].partner_list
    )
    cand_exp = cand_exp.filter(
        hl.is_defined(cand_exp.partner_list) & (hl.len(cand_exp.partner_list) > 0)
    )
    all_pairs = cand_exp.explode("partner_list")
    all_pairs = all_pairs.transmute(
        partner_locus=all_pairs.partner_list.partner_locus,
        partner_alleles=all_pairs.partner_list.partner_alleles,
        partner_af=all_pairs.partner_list.partner_af,
        partner_an=all_pairs.partner_list.partner_an,
    )
    all_pairs = all_pairs.filter(
        ~(
            (all_pairs.candidate_locus == all_pairs.partner_locus)
            & (all_pairs.candidate_alleles == all_pairs.partner_alleles)
        )
    )

    # ----- Step 3: Left-join with observed pairs; synthesize unobserved. -----
    obs = observed_pairs[
        all_pairs.candidate_locus,
        all_pairs.candidate_alleles,
        all_pairs.partner_locus,
        all_pairs.partner_alleles,
        all_pairs.gene_id,
    ]
    n_pair_synth = hl.coalesce(
        hl.min(all_pairs.candidate_an, all_pairs.partner_an) / 2,
        hl.float64(n_samples),
    )
    n_pair_synth = hl.if_else(n_pair_synth > 0, n_pair_synth, hl.float64(n_samples))
    af_c = all_pairs.candidate_af
    af_p = all_pairs.partner_af
    synthetic_e = 2 * n_pair_synth * af_c * (1 - af_c) * af_p * (1 - af_p)

    ld_pair_extras = {}
    if include_ld_adjusted:
        # Unobserved partners default to D* = 0 (no LD correction available
        # since we have no gt_counts row to estimate MAF_{A,B} from).
        # ld_term_synth = 1 - 2 * MAF_c * MAF_p; clamp to [eps, 1].
        maf_c = hl.min(all_pairs.candidate_af, 1 - all_pairs.candidate_af)
        maf_p = hl.min(all_pairs.partner_af, 1 - all_pairs.partner_af)
        raw_term_synth = 1.0 - 2.0 * maf_c * maf_p
        ld_term_synth = hl.max(hl.min(raw_term_synth, 1.0), 1e-300)
        log_ld_term_synth = hl.log(ld_term_synth)
        ld_pair_extras = dict(
            ld_term=hl.coalesce(obs.ld_term, ld_term_synth),
            log_ld_term=hl.coalesce(obs.log_ld_term, log_ld_term_synth),
            d_prime=obs.d_prime,
            d_star=hl.coalesce(obs.d_star, 0.0),
        )

    pairs = all_pairs.annotate(
        e_pair=hl.coalesce(obs.e_pair, synthetic_e),
        o_pair=hl.coalesce(obs.o_pair, 0.0),
        double_carriers=hl.coalesce(obs.double_carriers, 0),
        p_chet=obs.p_chet,
        **ld_pair_extras,
    )
    pairs = pairs.cache()

    ld_aggregate_extras = {}
    ld_partner_struct_extras = {}
    if include_ld_adjusted:
        ld_aggregate_extras = dict(
            sum_log_ld_term=hl.agg.sum(pairs.log_ld_term),
            n_ld_adjusted_used=hl.int32(
                hl.agg.count_where(pairs.d_star == 1.0)
            ),
            _candidate_an_take=hl.agg.take(pairs.candidate_an, 1),
        )
        ld_partner_struct_extras = dict(
            d_prime=pairs.d_prime,
            d_star=pairs.d_star,
            ld_term=pairs.ld_term,
        )

    result = pairs.group_by(
        pairs.candidate_locus,
        pairs.candidate_alleles,
        pairs.gene_id,
    ).aggregate(
        candidate_af=hl.agg.take(pairs.candidate_af, 1)[0],
        n_partners=hl.int32(hl.agg.count()),
        total_expected_in_trans=hl.agg.sum(pairs.e_pair),
        total_observed_in_trans=hl.agg.sum(pairs.o_pair),
        partners=hl.agg.collect(
            hl.struct(
                locus=pairs.partner_locus,
                alleles=pairs.partner_alleles,
                af=pairs.partner_af,
                double_carriers=pairs.double_carriers,
                p_chet=pairs.p_chet,
                expected_in_trans=pairs.e_pair,
                observed_in_trans=pairs.o_pair,
                **ld_partner_struct_extras,
            )
        ),
        **ld_aggregate_extras,
    )

    if include_ld_adjusted:
        # Rachel's outer N: use the candidate's own callable count
        # (candidate_an / 2), which upper-bounds jointly-callable samples for
        # the candidate. Falls back to n_samples if AN is missing.
        cand_an_taken = result._candidate_an_take
        cand_an_val = hl.if_else(
            hl.len(cand_an_taken) > 0,
            cand_an_taken[0],
            hl.missing(hl.tint64),
        )
        n_ld = hl.coalesce(
            hl.float64(cand_an_val) / 2.0, hl.float64(n_samples)
        )
        n_ld = hl.if_else(n_ld > 0, n_ld, hl.float64(n_samples))
        total_expected_in_trans_ld_adjusted = (
            n_ld * (1.0 - hl.exp(result.sum_log_ld_term))
        )
        result = result.annotate(
            total_expected_in_trans_ld_adjusted=total_expected_in_trans_ld_adjusted,
        ).drop("sum_log_ld_term", "_candidate_an_take")

    populations_type = hl.tstruct(
        id=hl.tstr,
        n_partners=hl.tint32,
        total_expected_in_trans=hl.tfloat64,
        total_observed_in_trans=hl.tfloat64,
        poisson_lower_tail_p=hl.tfloat64,
    )

    warning_exprs = [
        hl.if_else(
            result.candidate_af < 0.005,
            hl.literal("low_candidate_af"),
            hl.missing(hl.tstr),
        ),
        hl.if_else(
            result.total_expected_in_trans < 1,
            hl.literal("low_total_expected"),
            hl.missing(hl.tstr),
        ),
    ]
    if include_ld_adjusted:
        warning_exprs.append(
            hl.if_else(
                result.total_expected_in_trans_ld_adjusted < 1,
                hl.literal("low_total_expected_ld_adjusted"),
                hl.missing(hl.tstr),
            )
        )

    result = result.annotate(
        partner_set=hl.literal(partner_set),
        populations=hl.empty_array(populations_type),
        warnings=hl.array(warning_exprs).filter(lambda w: hl.is_defined(w)),
    )

    result = result.rename(
        {"candidate_locus": "locus", "candidate_alleles": "alleles"}
    ).key_by("locus", "alleles", "gene_id", "partner_set")

    if include_ld_adjusted:
        result = result.annotate_globals(
            ld_threshold=hl.float64(ld_threshold),
            ld_null_model=hl.literal(NULL_MODEL_POISSON_HWE_LD_ADJUSTED),
        )

    return result


########################################################################################
### Poisson p-value (depletion test)
########################################################################################

def compute_poisson_p_ht(
    ht: hl.Table,
    include_ld_adjusted: bool = False,
) -> hl.Table:
    """
    Annotate ``poisson_lower_tail_p`` and ``null_model`` on an aggregated
    Table.

    Given fields ``total_expected_in_trans`` (E) and ``total_observed_in_trans``
    (O), computes the one-sided depletion p-value
    ``P(X <= O | X ~ Poisson(E))`` via ``scipy.stats.poisson.cdf``.

    ``O`` is fractional (``p_chet * double_carriers``); ``poisson.cdf`` floors
    it to an integer internally, giving the conservative
    ``P(X <= floor(O) | E)``.

    The aggregated Table is typically 1k–5k rows (one row per candidate ×
    gene × partner_set), so a collect / scipy / re-parallelize round-trip
    is the cleanest implementation. ``poisson_lower_tail_p`` is set to
    missing when ``E <= 0`` or either E/O is missing — the upstream
    ``low_total_expected`` warning (added by :func:`aggregate_oe_per_candidate`
    when ``E < 1``) flags underpowered tests for downstream consumers.

    Also annotates ``null_model = NULL_MODEL_POISSON_HWE``.

    When ``include_ld_adjusted=True``, additionally emits
    ``poisson_lower_tail_p_ld_adjusted`` computed against
    ``total_expected_in_trans_ld_adjusted`` (same collect / scipy /
    re-parallelize round-trip; single round-trip, not two).

    :param ht: Table from :func:`aggregate_oe_per_candidate`.
    :param include_ld_adjusted: If ``True``, also emits
        ``poisson_lower_tail_p_ld_adjusted`` from
        ``total_expected_in_trans_ld_adjusted``.
    :return: Same Table with ``poisson_lower_tail_p`` and ``null_model``.
    """
    select_kwargs = dict(
        _e=ht.total_expected_in_trans,
        _o=ht.total_observed_in_trans,
    )
    if include_ld_adjusted:
        select_kwargs["_e_ld"] = ht.total_expected_in_trans_ld_adjusted
    collected = ht.select(**select_kwargs).collect()

    p_rows = []
    for r in collected:
        if r._e is None or r._o is None or r._e <= 0:
            p = None
        else:
            p = float(poisson.cdf(r._o, r._e))
        row_kwargs = dict(
            locus=r.locus,
            alleles=r.alleles,
            gene_id=r.gene_id,
            partner_set=r.partner_set,
            poisson_lower_tail_p=p,
        )
        if include_ld_adjusted:
            e_ld = r._e_ld
            if e_ld is None or r._o is None or e_ld <= 0:
                p_ld = None
            else:
                p_ld = float(poisson.cdf(r._o, e_ld))
            row_kwargs["poisson_lower_tail_p_ld_adjusted"] = p_ld
        p_rows.append(hl.Struct(**row_kwargs))

    schema_fields = dict(
        locus=ht.locus.dtype,
        alleles=ht.alleles.dtype,
        gene_id=ht.gene_id.dtype,
        partner_set=ht.partner_set.dtype,
        poisson_lower_tail_p=hl.tfloat64,
    )
    if include_ld_adjusted:
        schema_fields["poisson_lower_tail_p_ld_adjusted"] = hl.tfloat64

    p_ht = hl.Table.parallelize(
        p_rows,
        schema=hl.tstruct(**schema_fields),
        key=list(ht.key),
    )

    annotate_kwargs = dict(
        poisson_lower_tail_p=p_ht[ht.key].poisson_lower_tail_p,
        null_model=hl.literal(NULL_MODEL_POISSON_HWE),
    )
    if include_ld_adjusted:
        annotate_kwargs["poisson_lower_tail_p_ld_adjusted"] = (
            p_ht[ht.key].poisson_lower_tail_p_ld_adjusted
        )

    return ht.annotate(**annotate_kwargs)


########################################################################################
### Top-level orchestration (pure transform; CLI / file I/O lives in main script)
########################################################################################

def in_trans_oe_pipeline(
    vp_gt_counts_ht: hl.Table,
    freq_ht: hl.Table,
    partner_ht: hl.Table,
    candidate_ht: hl.Table,
    n_samples: int,
    partner_set: str,
    use_adj: bool = True,
    include_ld_adjusted: bool = False,
    ld_threshold: float = 0.5,
) -> hl.Table:
    """
    Run the in-trans observed-vs-expected pipeline end-to-end as a pure
    transform: pair-counts + AFs + partner set + candidates → per-candidate
    depletion Table.

    Equivalent to::

        annotate_pair_oe_terms → aggregate_oe_per_candidate → compute_poisson_p_ht

    Callers are responsible for reading the inputs and writing the output;
    this function does no I/O.

    :param vp_gt_counts_ht: Pair-counts Table. See
        :func:`annotate_pair_oe_terms`.
    :param freq_ht: Frequency Table.
    :param partner_ht: Partner-set Table. See
        :func:`create_clinvar_plp_partner_ht`.
    :param candidate_ht: Candidate-variant Table.
    :param n_samples: Sample count behind the genotype counts.
    :param partner_set: One of :data:`PARTNER_SETS`.
    :param use_adj: Whether to use the adj or raw counts.
    :param include_ld_adjusted: If ``True``, also compute Rachel Unger's
        LD-adjusted expected-in-trans (product form) and the corresponding
        Poisson p-value, threaded through all three stages.
    :param ld_threshold: |D'| threshold used when
        ``include_ld_adjusted=True``.
    :return: Per-candidate depletion Table.
    """
    annotated = annotate_pair_oe_terms(
        vp_gt_counts_ht, freq_ht, n_samples=n_samples, use_adj=use_adj,
        include_ld_adjusted=include_ld_adjusted,
        ld_threshold=ld_threshold,
    )
    aggregated = aggregate_oe_per_candidate(
        annotated, partner_ht, candidate_ht, freq_ht,
        partner_set=partner_set,
        n_samples=n_samples,
        include_ld_adjusted=include_ld_adjusted,
        ld_threshold=ld_threshold,
    )
    return compute_poisson_p_ht(
        aggregated, include_ld_adjusted=include_ld_adjusted
    )

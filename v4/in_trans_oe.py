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

* :func:`gnomad_chets.v4.create_vp_matrix.create_variant_filter_ht` for
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
      used by :func:`gnomad_chets.v2.phasing.get_em_expr`.
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
        :func:`gnomad_chets.v4.create_vp_matrix.create_variant_pair_genotype_counts`.
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
    :return: Annotated pair-counts Table.
    """
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

    return vp_gt_counts_ht.annotate(
        af1=af1,
        af2=af2,
        n_pair=n_pair,
        double_carriers=double_carriers,
        p_chet=p_chet,
        e_pair=e_pair,
        o_pair=o_pair,
    )


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
        :func:`gnomad_chets.v4.create_vp_matrix.create_variant_filter_ht`.
    :param freq_ht: Frequency Table for candidate AF/AN lookup.
    :param partner_set: One of :data:`PARTNER_SETS`.
    :param n_samples: Fallback sample count if AN is unavailable on either
        side of an unobserved pair. Should match the cohort behind the gt
        counts (e.g. 730947 for v4.1 exomes).
    :return: Aggregated per-candidate Table with the schema above.
    """
    if partner_set not in PARTNER_SETS:
        raise ValueError(f"partner_set {partner_set!r} not in {PARTNER_SETS}")

    # ----- Step 1: Build the per-(candidate, partner, gene) OBSERVED view. -----
    # The pair-counts Table is keyed by (v1, v2) with v1 <= v2, so each pair
    # gets two views (candidate=v1 and candidate=v2). The union doubles the
    # row count temporarily but lets us use a single (candidate, partner)
    # lookup downstream.
    ht = annotated_pair_ht.key_by()
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

    pairs = all_pairs.annotate(
        e_pair=hl.coalesce(obs.e_pair, synthetic_e),
        o_pair=hl.coalesce(obs.o_pair, 0.0),
        double_carriers=hl.coalesce(obs.double_carriers, 0),
        p_chet=obs.p_chet,
    )
    pairs = pairs.cache()

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
            )
        ),
    )

    populations_type = hl.tstruct(
        id=hl.tstr,
        n_partners=hl.tint32,
        total_expected_in_trans=hl.tfloat64,
        total_observed_in_trans=hl.tfloat64,
        poisson_lower_tail_p=hl.tfloat64,
    )

    result = result.annotate(
        partner_set=hl.literal(partner_set),
        populations=hl.empty_array(populations_type),
        warnings=hl.array(
            [
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
        ).filter(lambda w: hl.is_defined(w)),
    )

    return result.rename(
        {"candidate_locus": "locus", "candidate_alleles": "alleles"}
    ).key_by("locus", "alleles", "gene_id", "partner_set")


########################################################################################
### Poisson p-value (depletion test)
########################################################################################

def compute_poisson_p_ht(ht: hl.Table) -> hl.Table:
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

    :param ht: Table from :func:`aggregate_oe_per_candidate`.
    :return: Same Table with ``poisson_lower_tail_p`` and ``null_model``.
    """
    collected = ht.select(
        _e=ht.total_expected_in_trans,
        _o=ht.total_observed_in_trans,
    ).collect()

    p_rows = []
    for r in collected:
        if r._e is None or r._o is None or r._e <= 0:
            p = None
        else:
            p = float(poisson.cdf(r._o, r._e))
        p_rows.append(hl.Struct(
            locus=r.locus,
            alleles=r.alleles,
            gene_id=r.gene_id,
            partner_set=r.partner_set,
            poisson_lower_tail_p=p,
        ))

    p_ht = hl.Table.parallelize(
        p_rows,
        schema=hl.tstruct(
            locus=ht.locus.dtype,
            alleles=ht.alleles.dtype,
            gene_id=ht.gene_id.dtype,
            partner_set=ht.partner_set.dtype,
            poisson_lower_tail_p=hl.tfloat64,
        ),
        key=list(ht.key),
    )

    return ht.annotate(
        poisson_lower_tail_p=p_ht[ht.key].poisson_lower_tail_p,
        null_model=hl.literal(NULL_MODEL_POISSON_HWE),
    )


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
    :return: Per-candidate depletion Table.
    """
    annotated = annotate_pair_oe_terms(
        vp_gt_counts_ht, freq_ht, n_samples=n_samples, use_adj=use_adj
    )
    aggregated = aggregate_oe_per_candidate(
        annotated, partner_ht, candidate_ht, freq_ht,
        partner_set=partner_set,
        n_samples=n_samples,
    )
    return compute_poisson_p_ht(aggregated)

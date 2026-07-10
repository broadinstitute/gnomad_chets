"""
Run the in-trans observed-vs-expected aggregation against the chets pipeline
output, optionally exporting a JSON blob for the gnomad-browser demo panel.

Inputs (must already exist; produced by ``v4/compute_vp_counts.py``):

* ``variant_pair_genotype_counts.ht`` — output of step 4
* ``variant_filter.ht`` — output of step 1 (the candidate set)

Plus gnomad_qc resources read directly:

* ``filter_ht``, ``freq_ht``, ``vep_ht`` — for partner-set construction
* ``clinvar.ht()`` filtered through ``filter_to_clinvar_pathogenic`` — partner set

Outputs:

* ``in_trans_oe.{partner_set}.{output_postfix}.ht`` — per-(candidate, gene)
  depletion table (full HT for re-use).
* ``--output-json`` (optional) — single-candidate row enriched with HGVSp /
  consequence / ClinVar significance / gene symbol, formatted to match the
  ``DepletionResult`` TypeScript type used by the browser panel.

Submit on Dataproc after the pipeline finishes::

    hailctl dataproc submit test-chets-highmem gnomad_chets/v4/run_in_trans_oe.py \\
        --gene CAPN3 \\
        --output-postfix capn3_demo \\
        --candidate-variant-id 15-42403721-C-G \\
        --output-json gs://gnomad-tmp-30day/in_trans_demo.capn3.json \\
        --pyfiles gnomad_chets

Then ``gsutil cp gs://gnomad-tmp-30day/in_trans_demo.capn3.json
browser/src/VariantPage/in_trans_demo_data.json`` on the local repo.
"""

import argparse
import json
import logging
import math
from typing import Any, Optional

import hail as hl
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.vep import filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.annotations import get_freq, get_vep
from gnomad_qc.v4.resources.variant_qc import final_filter

from gnomad_chets.v4.in_trans_oe import (
    NULL_MODEL_POISSON_HWE,
    PARTNER_SET_CLINVAR_BLB,
    PARTNER_SET_CLINVAR_PLP,
    PARTNER_SET_CLINVAR_VUS,
    aggregate_oe_per_candidate,
    annotate_pair_oe_terms,
    compute_poisson_p_ht,
)
from gnomad_chets.v4.resources import (
    CLINVAR_CATEGORIES,
    CLINVAR_CATEGORY_PARTNER_SET,
    DEFAULT_DATA_TYPE,
    DEFAULT_MAX_FREQ,
    DEFAULT_TMP_DIR,
    TEST_INTERVALS,
    VARIANT_COOCCURRENCE_ROOT,
    _get_output_postfix,
    get_excluded_genes_ht,
    get_phase,
    get_variant_filter_ht,
    get_variant_pair_genotype_counts_ht,
)
from gnomad_chets.v4.utils import (
    annotate_filters_af,
    filter_clinvar_by_category,
    filter_for_testing,
)

# Maps the --partner-set CLI value to the ClinVar category string used by
# create_clinvar_category_partner_ht / filter_clinvar_by_category.
_PARTNER_SET_TO_CLINVAR_CATEGORY = {
    PARTNER_SET_CLINVAR_PLP: "plp",
    PARTNER_SET_CLINVAR_BLB: "blb",
    PARTNER_SET_CLINVAR_VUS: "vus",
}

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_in_trans_oe")
logger.setLevel(logging.INFO)


def create_clinvar_category_partner_ht(
    vep_ht: hl.Table,
    clinvar_ht: hl.Table,
    category: str,
    max_freq: float = DEFAULT_MAX_FREQ,
) -> hl.Table:
    """Build a ClinVar partner-set Table for one significance category.

    Expects ``vep_ht`` to already have ``filters``, ``af``, and ``an``
    annotated (via :func:`annotate_filters_af` with ``include_an=True``)
    and be pre-filtered to QC-PASS + AF > 0. Takes the *unfiltered*
    ClinVar HT and filters internally via :func:`filter_clinvar_by_category`,
    keeping the category gating in one place. The output adds a
    ``partner_set`` literal (e.g. ``"CLINVAR_PLP"``) so the same HT can
    be consumed downstream as a partner set.

    :param vep_ht: VEP Table with ``filters`` / ``af`` / ``an`` pre-annotated.
    :param clinvar_ht: Unfiltered ClinVar HT.
    :param category: One of :data:`CLINVAR_CATEGORIES`.
    :param max_freq: Upper AF bound (inclusive).
    :return: Partner-set Table keyed by ``(locus, alleles)`` with
        ``gene_id``, ``af``, ``an``, ``partner_set``.
    :raises ValueError: If ``category`` is not in :data:`CLINVAR_CATEGORIES`.
    """
    if category not in CLINVAR_CATEGORIES:
        raise ValueError(
            f"Unknown ClinVar category: {category!r}. Valid: "
            f"{CLINVAR_CATEGORIES}."
        )
    clinvar_filtered = filter_clinvar_by_category(clinvar_ht, category)
    ht = vep_ht.annotate(
        gene_id=hl.array(hl.set(filter_vep_transcript_csqs_expr(
            vep_ht.vep.transcript_consequences,
            protein_coding=True,
            ensembl_only=True,
        ).map(lambda csq: csq.gene_id))),
        _is_in_clinvar_category=hl.is_defined(
            clinvar_filtered[vep_ht.locus, vep_ht.alleles]
        ),
    )
    ht = ht.filter(
        (ht.af <= max_freq)
        & hl.is_defined(ht.gene_id)
        & (hl.len(ht.gene_id) > 0)
        & ht._is_in_clinvar_category
    )
    return ht.select(
        gene_id=ht.gene_id,
        af=ht.af,
        an=ht.an,
        partner_set=hl.literal(CLINVAR_CATEGORY_PARTNER_SET[category]),
    )


def _poisson_lower_tail_p(o: float, e: float) -> Optional[float]:
    """One-sided depletion p: P(X <= round(O) | X ~ Poisson(E)). Forward recursion."""
    if e is None or e <= 0:
        return None
    k = max(int(round(o)), 0)
    term = math.exp(-e)
    s = term
    for i in range(1, k + 1):
        term *= e / i
        s += term
    return min(s, 1.0)


def _enrich_partner(
    partner: hl.struct, vep_ht: hl.Table, clinvar_ht: hl.Table, gene_id: hl.expr.StringExpression
) -> hl.expr.StructExpression:
    """Add HGVSp / consequence / ClinVar significance to a partner struct."""
    vep_row = vep_ht[partner.locus, partner.alleles]
    clinvar_row = clinvar_ht[partner.locus, partner.alleles]
    matching_csq = (
        vep_row.vep.transcript_consequences
        .filter(lambda csq: csq.gene_id == gene_id)
        .filter(
            lambda csq: csq.canonical == 1
        )
    )
    fallback_csq = vep_row.vep.transcript_consequences.filter(
        lambda csq: csq.gene_id == gene_id
    )
    chosen_csq = hl.if_else(
        hl.len(matching_csq) > 0, matching_csq[0], fallback_csq[0]
    )
    return partner.annotate(
        hgvsp=chosen_csq.hgvsp,
        consequence=chosen_csq.consequence_terms[0],
        clinvar_clinical_significance=hl.delimit(
            clinvar_row.info.CLNSIG, ", "
        ),
    )


def _resolve_output_path(
    resource_name: str,
    extension: str,
    test: bool,
    output_postfix: Optional[str],
) -> str:
    output_dir = DEFAULT_TMP_DIR if test else VARIANT_COOCCURRENCE_ROOT
    postfix = _get_output_postfix(output_postfix, test)
    return f"{output_dir}/exomes.{resource_name}{postfix}{extension}"


def _row_to_demo_json(
    row: Any,
    candidate_variant_id: str,
    dataset: str,
    partner_set: str,
    n_total_partners_in_gene: int,
    partner_limit: int,
) -> dict:
    """Convert a single aggregation row into the DEMO_DATA shape for the browser."""
    sorted_partners = sorted(row.partners, key=lambda p: -(p.af or 0))
    truncated = len(sorted_partners) > partner_limit
    selected = sorted_partners[:partner_limit]

    partners_out = []
    for p in selected:
        chrom_str = p.locus.contig.replace("chr", "")
        ref, alt = p.alleles
        partners_out.append({
            "variant_id": f"{chrom_str}-{p.locus.position}-{ref}-{alt}",
            "hgvsp": p.hgvsp if "hgvsp" in p else None,
            "consequence": p.consequence if "consequence" in p else None,
            "clinvar_clinical_significance": (
                p.clinvar_clinical_significance if "clinvar_clinical_significance" in p else None
            ),
            "allele_frequency": p.af,
            "double_carriers": int(p.double_carriers),
            "p_compound_heterozygous": p.p_chet if not (p.p_chet is None or math.isnan(p.p_chet)) else None,
            "expected_in_trans": p.expected_in_trans,
            "observed_in_trans": p.observed_in_trans,
        })

    e_total = row.total_expected_in_trans
    o_total = row.total_observed_in_trans
    poisson_p = _poisson_lower_tail_p(o_total, e_total)

    warnings = list(row.warnings) if row.warnings else []
    if truncated:
        warnings.append(
            f"Top {partner_limit} of {n_total_partners_in_gene} partners shown (sorted by AF)."
        )

    power_note = None
    if e_total < 1:
        power_note = (
            f"Total expected in-trans count is {e_total:.2f} — below the threshold "
            f"(E < 1) where the depletion test is meaningful."
        )

    return {
        "candidate_variant_id": candidate_variant_id,
        "gene_id": row.gene_id,
        "gene_symbol": row.gene_symbol if "gene_symbol" in row else None,
        "dataset": dataset,
        "partner_set": partner_set,
        "null_model": NULL_MODEL_POISSON_HWE,
        "n_partners": int(row.n_partners),
        "total_expected_in_trans": e_total,
        "total_observed_in_trans": o_total,
        "poisson_lower_tail_p": poisson_p,
        "log_bayes_factor_hypomorphic": None,
        "power_note": power_note,
        "partners": partners_out,
        "partners_truncated": truncated,
        "populations": [],
        "warnings": warnings,
    }


def main(args):
    if args.gene and args.interval:
        raise ValueError("--gene and --interval are mutually exclusive.")
    if not (0 <= args.ld_threshold <= 1):
        raise ValueError(
            f"--ld-threshold must be in [0, 1], got {args.ld_threshold}"
        )
    if args.ld_threshold != 0.5 and not args.with_ld_adjustment:
        logger.warning(
            "--ld-threshold set but --with-ld-adjustment is off; ignoring."
        )
    test = args.test or bool(args.gene) or bool(args.interval)
    if args.interval:
        # Use a synthetic label; filter_for_testing only cares about the value.
        test_intervals = {args.interval: args.interval}
    elif args.gene:
        test_intervals = {args.gene: TEST_INTERVALS[args.gene]}
    else:
        test_intervals = TEST_INTERVALS

    # Explicitly init Hail so cache()/checkpoint() spill to GCS, not the
    # executors' local disk. On small clusters (jg3 = 2×n1-standard-8), the
    # chr19-scale aggregation blows past the local /tmp allotment and the
    # executors are killed by YARN before Spark can retry cleanly.
    hl.init(tmp_dir=DEFAULT_TMP_DIR, log="/tmp/run_in_trans_oe.log")

    # Step 1: Read pipeline outputs. --gt-counts-ht-path / --variant-filter-ht-path
    # let a per-gene run reuse a shared upstream bundle while writing a
    # per-gene OE output HT (matches the --phased-ht-path pattern).
    logger.info("Reading variant pair genotype counts and variant filter HTs...")
    if args.gt_counts_ht_path:
        logger.info("Reading gt-counts HT from override path: %s", args.gt_counts_ht_path)
        vp_gt_counts_ht = hl.read_table(args.gt_counts_ht_path)
    else:
        vp_gt_counts_ht = get_variant_pair_genotype_counts_ht(
            data_type=args.data_type,
            test=test,
            output_postfix=args.output_postfix,
        ).ht()
    if args.variant_filter_ht_path:
        logger.info("Reading variant_filter HT from override path: %s", args.variant_filter_ht_path)
        candidate_ht = hl.read_table(args.variant_filter_ht_path)
    else:
        candidate_ht = get_variant_filter_ht(
            data_type=args.data_type,
            test=test,
            output_postfix=args.output_postfix,
        ).ht()

    excluded_gene_set = None
    if args.skip_excluded_genes:
        if args.excluded_genes_ht_path:
            excluded_path = args.excluded_genes_ht_path
        else:
            excluded_res = get_excluded_genes_ht(
                data_type=args.data_type,
                test=test,
                output_postfix=args.output_postfix,
            )
            excluded_path = excluded_res.path
        try:
            excluded_ht = hl.read_table(excluded_path)
            excluded_gene_set = excluded_ht.aggregate(
                hl.agg.collect_as_set(excluded_ht.gene_id)
            )
            logger.info(
                "Loaded excluded-genes set from %s (n=%d) — will drop "
                "candidates and partners whose gene_id intersects it.",
                excluded_path, len(excluded_gene_set),
            )
        except Exception as e:
            logger.warning(
                "Could not load excluded-genes HT at %s (%s). Continuing "
                "without the excluded-gene filter; pass "
                "--no-skip-excluded-genes to suppress this warning.",
                excluded_path, e,
            )
            excluded_gene_set = None

    # Step 2: Read gnomad_qc resources for partner-set construction.
    logger.info("Reading gnomad_qc filter / freq / vep HTs and ClinVar...")
    filter_ht = final_filter(data_type=args.data_type).ht()
    freq_ht = get_freq(data_type=args.data_type).ht()
    vep_ht = get_vep(data_type=args.data_type).ht()
    clinvar_full_ht = clinvar.ht()

    if test:
        logger.info("Filtering input HTs to test intervals: %s", list(test_intervals))
        filter_ht = filter_for_testing(filter_ht, test_intervals)
        freq_ht = filter_for_testing(freq_ht, test_intervals)
        vep_ht = filter_for_testing(vep_ht, test_intervals)
        clinvar_full_ht = filter_for_testing(clinvar_full_ht, test_intervals)

    # Step 3: Build partner set for the requested ClinVar significance
    # category (default CLINVAR_PLP; also supports CLINVAR_BLB / CLINVAR_VUS).
    # Pre-annotate vep_ht with filters / af / an and drop non-PASS or AF=0
    # rows once, matching the pattern used in create_variant_filter_ht.
    annotated_vep_ht = annotate_filters_af(
        vep_ht, filter_ht, freq_ht, include_an=True
    )
    annotated_vep_ht = annotated_vep_ht.filter(
        (annotated_vep_ht.filters.length() == 0) & (annotated_vep_ht.af > 0)
    )
    clinvar_category = _PARTNER_SET_TO_CLINVAR_CATEGORY[args.partner_set]
    logger.info("Building %s partner set...", args.partner_set)
    partner_ht = create_clinvar_category_partner_ht(
        vep_ht=annotated_vep_ht,
        clinvar_ht=clinvar_full_ht,
        category=clinvar_category,
        max_freq=args.max_freq,
    )

    # AF filter on candidates — partners are already ≤ max_freq via
    # create_clinvar_category_partner_ht above, but candidate_ht (from the
    # variant_filter) may include in_trans_oe_candidate entries up to AF 0.5.
    # Rare-only default gives a cleaner statistic — common-AF candidates
    # dominate the top hits with mostly LD-driven signal.
    candidate_af_expr = freq_ht[candidate_ht.locus, candidate_ht.alleles].freq[0].AF
    n_cand_pre = candidate_ht.count()
    candidate_ht = candidate_ht.annotate(_af=candidate_af_expr)
    candidate_ht = candidate_ht.filter(
        hl.is_defined(candidate_ht._af)
        & (candidate_ht._af > 0)
        & (candidate_ht._af <= args.max_freq)
    ).drop("_af")
    n_cand_post = candidate_ht.count()
    logger.info(
        "Candidate AF filter (≤ %s): %d → %d rows (%.1f%% kept).",
        args.max_freq, n_cand_pre, n_cand_post,
        100.0 * n_cand_post / max(n_cand_pre, 1),
    )

    # Excluded-gene filter on candidates and partners: drop entries whose
    # gene_id array intersects the pipeline's --exclude-gene-ids set. Prevents
    # false-positive top hits with O=0 for MUC16/RYR1/FBN3/ABCA7 candidates
    # whose pairs were dropped from gt-counts upstream.
    if excluded_gene_set:
        excl_lit = hl.literal(excluded_gene_set)
        n_pre = candidate_ht.count()
        candidate_ht = candidate_ht.filter(
            hl.len(hl.set(candidate_ht.gene_id).intersection(excl_lit)) == 0
        )
        logger.info(
            "Excluded-gene filter on candidates: %d → %d rows.",
            n_pre, candidate_ht.count(),
        )
        n_pre = partner_ht.count()
        partner_ht = partner_ht.filter(
            hl.len(hl.set(partner_ht.gene_id).intersection(excl_lit)) == 0
        )
        logger.info(
            "Excluded-gene filter on partners: %d → %d rows.",
            n_pre, partner_ht.count(),
        )

    # Step 4: Per-pair O/E annotations and per-candidate aggregation.
    # Reuse the pre-computed phased HT for p_chet if it's available — the
    # EM step already ran when phase_gnomad.py wrote it, so re-running
    # haplotype_freq_em here is wasted work (~1-2 min per per-gene call).
    phased_ht = None
    if args.use_phased:
        if args.phased_ht_path:
            phased_path = args.phased_ht_path
        else:
            phased_path = get_phase(
                data_type=args.data_type,
                test=test,
                output_postfix=args.output_postfix,
            ).path
        try:
            phased_ht = hl.read_table(phased_path)
            logger.info(
                "Loaded phased HT from %s — will read p_chet from "
                "em.%s.p_chet instead of re-running EM.",
                phased_path, "adj" if args.use_adj else "raw",
            )
        except Exception as e:
            logger.warning(
                "Could not load phased HT at %s (%s). Falling back to "
                "inline EM.", phased_path, e,
            )
            phased_ht = None

    logger.info(
        "Annotating pair O/E and aggregating per candidate × gene "
        "(LD-adjusted=%s, ld_threshold=%s)...",
        args.with_ld_adjustment, args.ld_threshold,
    )
    annotated = annotate_pair_oe_terms(
        vp_gt_counts_ht, freq_ht, n_samples=args.n_samples, use_adj=args.use_adj,
        phased_ht=phased_ht,
        include_ld_adjusted=args.with_ld_adjustment,
        ld_threshold=args.ld_threshold,
    )
    result = aggregate_oe_per_candidate(
        annotated, partner_ht, candidate_ht, freq_ht,
        partner_set=args.partner_set,
        n_samples=args.n_samples,
        include_ld_adjusted=args.with_ld_adjustment,
        ld_threshold=args.ld_threshold,
    )
    result = compute_poisson_p_ht(
        result, include_ld_adjusted=args.with_ld_adjustment
    )

    # Powered-rows filter: drop rows below the min-E floor. Below ~1
    # expected pair, depletion is unresolvable regardless of the observed
    # count. This shrinks the output by ~99% on chr19 (most PLP candidate
    # rows have E << 1) but concentrates the Poisson-tail statistic on the
    # small fraction of rows where it can actually distinguish signal.
    if args.min_expected_in_trans > 0:
        n_pre = result.count()
        result = result.filter(
            result.total_expected_in_trans >= args.min_expected_in_trans
        )
        n_post = result.count()
        logger.info(
            "Powered-rows filter (E ≥ %s): %d → %d rows (%.2f%% kept).",
            args.min_expected_in_trans, n_pre, n_post,
            100.0 * n_post / max(n_pre, 1),
        )

    # Step 5: Annotate gene_symbol (from VEP) on the result; enrich partners.
    logger.info("Enriching with gene_symbol, hgvsp, consequence, and ClinVar significance...")
    # gene_symbol from any tx_csq matching gene_id (canonical preferred).
    candidate_vep = vep_ht[result.locus, result.alleles].vep
    gene_symbol_match = candidate_vep.transcript_consequences.filter(
        lambda csq: csq.gene_id == result.gene_id
    )
    result = result.annotate(
        gene_symbol=hl.if_else(
            hl.len(gene_symbol_match) > 0,
            gene_symbol_match[0].gene_symbol,
            hl.missing(hl.tstr),
        ),
    )
    # NOTE: partner enrichment (hgvsp / consequence / clinvar significance)
    # used to happen here via ``result.partners.map(_enrich_partner)``, but
    # that idiom — array-of-struct map with cross-table joins on vep_ht and
    # clinvar_full_ht — produces an IR Hail can't render at write time
    # (``KeyError: '__uid_*'``). For the demo path we only need enrichment
    # for the one candidate's partners in step 7, which is small enough
    # to do in Python after collecting. The full HT therefore stores
    # bare partner structs; downstream consumers wanting enriched data can
    # call ``_enrich_partner_python`` (defined below) themselves.

    # Step 6: Write the full HT.
    full_ht_path = _resolve_output_path(
        f"in_trans_oe.{args.partner_set.lower()}",
        ".ht",
        test=test,
        output_postfix=args.output_postfix,
    )
    logger.info("Writing aggregated HT → %s", full_ht_path)
    result = result.checkpoint(full_ht_path, overwrite=args.overwrite)
    n_rows = result.count()
    logger.info("Aggregated HT has %d rows.", n_rows)

    # Step 7: If a specific candidate was requested, dump JSON for the browser demo.
    if args.candidate_variant_id and args.output_json:
        chrom, pos, ref, alt = args.candidate_variant_id.split("-")
        chrom_normalized = f"chr{chrom}" if not chrom.startswith("chr") else chrom
        candidate_locus = hl.locus(chrom_normalized, int(pos), "GRCh38")
        candidate_alleles_py = [ref, alt]
        candidate_alleles = hl.literal(candidate_alleles_py)

        # 7a: Pick the candidate's gene. Prefer a gene that already has an
        # aggregation row (i.e. ≥1 observed-pair P/LP partner); otherwise fall
        # back to the first gene_id from the candidate's variant_filter row.
        sub = result.filter(
            (result.locus == candidate_locus) & (result.alleles == candidate_alleles)
        )
        observed_rows = sub.collect()

        if observed_rows:
            chosen_observed_row = observed_rows[0]
            if len(observed_rows) > 1:
                logger.info(
                    "Candidate is in %d genes with observed P/LP pairs; picking %s.",
                    len(observed_rows), chosen_observed_row.gene_id,
                )
            chosen_gene_id = chosen_observed_row.gene_id
            candidate_af_value = chosen_observed_row.candidate_af
            existing_warnings = list(chosen_observed_row.warnings) if chosen_observed_row.warnings else []
            observed_partners = list(chosen_observed_row.partners)
        else:
            # No co-carrier pairs with any P/LP partner. Resolve gene_id via
            # candidate_ht when present, otherwise via VEP. Resolve AF via freq_ht.
            chosen_observed_row = None
            cand_filter_rows = candidate_ht.filter(
                (candidate_ht.locus == candidate_locus)
                & (candidate_ht.alleles == candidate_alleles)
            ).collect()
            cand_gene_ids = (
                list(cand_filter_rows[0].gene_id) if cand_filter_rows else []
            )
            if not cand_gene_ids:
                # Fall back to VEP for gene_id assignment.
                vep_lookup_ht = vep_ht.filter(
                    (vep_ht.locus == candidate_locus)
                    & (vep_ht.alleles == candidate_alleles)
                )
                vep_lookup_rows = vep_lookup_ht.select(
                    transcript_consequences=vep_lookup_ht.vep.transcript_consequences,
                ).collect()
                if vep_lookup_rows:
                    cand_gene_ids = list({
                        tc.gene_id for tc in vep_lookup_rows[0].transcript_consequences
                        if tc.gene_id is not None
                    })
            if not cand_gene_ids:
                logger.warning(
                    "Candidate %s has no gene_id in candidate_ht or VEP; "
                    "cannot build demo data.",
                    args.candidate_variant_id,
                )
                return
            # Pick the gene_id with the most P/LP partners — this naturally
            # favors the disease gene over neighboring annotations (overlapping
            # ncRNAs, antisense transcripts) that have no ClinVar variants.
            partner_counts_by_gene = {
                g: partner_ht.filter(partner_ht.gene_id.contains(g)).count()
                for g in cand_gene_ids
            }
            chosen_gene_id = max(partner_counts_by_gene, key=partner_counts_by_gene.get)
            logger.info(
                "Candidate %s gene-id candidates: %s; picked %s (%d P/LP partners).",
                args.candidate_variant_id,
                partner_counts_by_gene,
                chosen_gene_id,
                partner_counts_by_gene[chosen_gene_id],
            )
            cand_freq_ht = freq_ht.filter(
                (freq_ht.locus == candidate_locus)
                & (freq_ht.alleles == candidate_alleles)
            )
            cand_freq_rows = cand_freq_ht.select(
                af=cand_freq_ht.freq[0].AF,
            ).collect()
            if not cand_freq_rows:
                logger.warning(
                    "Candidate %s has no freq_ht entry; cannot build demo data.",
                    args.candidate_variant_id,
                )
                return
            candidate_af_value = cand_freq_rows[0].af
            existing_warnings = [
                "No co-carrier pair list entries found for this candidate; "
                "observed counts default to 0 across all partners. The "
                "depletion signal here may reflect missing upstream pair-list "
                "coverage rather than true biological depletion."
            ]
            observed_partners = []
            logger.info(
                "Candidate %s has no observed P/LP pair partners — building "
                "expected-only demo from cross-join with the gene's P/LP set.",
                args.candidate_variant_id,
            )

        # 7a-bis: Look up candidate AN (used by the cross-join for unobserved
        # partners). Fetch unconditionally so both the observed and the
        # expected-only branches have it.
        cand_an_ht = freq_ht.filter(
            (freq_ht.locus == candidate_locus)
            & (freq_ht.alleles == candidate_alleles)
        )
        cand_an_rows = cand_an_ht.select(an=cand_an_ht.freq[0].AN).collect()
        candidate_an_value = cand_an_rows[0].an if cand_an_rows else None

        # 7b: Build the full P/LP partner set for this candidate's gene.
        # ``partner_ht.gene_id`` is an ``array<str>``; semi-join via filter.
        partners_in_gene_ht = partner_ht.filter(
            partner_ht.gene_id.contains(chosen_gene_id)
        )
        # Exclude the candidate itself (a variant can't be its own in-trans partner).
        partners_in_gene_ht = partners_in_gene_ht.filter(
            ~((partners_in_gene_ht.locus == candidate_locus)
              & (partners_in_gene_ht.alleles == candidate_alleles))
        )
        partner_rows = partners_in_gene_ht.collect()
        n_total_in_gene = len(partner_rows)

        # Index observed partners by (chrom, pos, alleles) for O() lookup.
        observed_by_key = {
            (p.locus.contig, p.locus.position, tuple(p.alleles)): p
            for p in observed_partners
        }

        # 7c: Cross-join — every P/LP partner in the gene gets a row. If the
        # pair is observed (in the pair list), use its O / double_carriers /
        # p_chet; otherwise default to 0 (no co-carriers seen) and compute E
        # from marginal AFs.
        cross_join_partners = []
        for partner in partner_rows:
            key = (partner.locus.contig, partner.locus.position, tuple(partner.alleles))
            obs = observed_by_key.get(key)
            if obs is not None:
                cross_join_partners.append(obs)
            else:
                # Match annotate_pair_oe_terms:
                #   2 * N_pair * AF_c * (1-AF_c) * AF_p * (1-AF_p)
                # For unobserved pairs we don't have a per-pair gt_count, so
                # estimate N_pair as min(AN_candidate, AN_partner) / 2 -- the
                # number of jointly-callable samples (upper-bounded by the
                # less-callable variant). Falls back to args.n_samples if AN
                # isn't available for either side.
                if (candidate_an_value is not None and partner.an is not None):
                    n_pair_est = min(candidate_an_value, partner.an) / 2
                else:
                    n_pair_est = args.n_samples
                e_pair = (
                    2 * n_pair_est
                    * candidate_af_value * (1 - candidate_af_value)
                    * partner.af * (1 - partner.af)
                )
                cross_join_partners.append(
                    hl.Struct(
                        locus=partner.locus,
                        alleles=list(partner.alleles),
                        af=partner.af,
                        double_carriers=0,
                        p_chet=None,
                        expected_in_trans=e_pair,
                        observed_in_trans=0.0,
                    )
                )

        # 7d: Enrich each partner with hgvsp / consequence / clinvar_clinical_significance.
        # Doing this in Python (after a semi_join + collect of a tiny VEP / ClinVar
        # slice) avoids the Hail IR-rendering bug we'd hit doing it via partners.map.
        partner_keys_ht = hl.Table.parallelize(
            [
                hl.struct(locus=p.locus, alleles=list(p.alleles))
                for p in cross_join_partners
            ],
            schema=hl.tstruct(locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr)),
            key=["locus", "alleles"],
        )
        vep_partners_ht = vep_ht.semi_join(partner_keys_ht)
        vep_rows = vep_partners_ht.select(
            transcript_consequences=vep_partners_ht.vep.transcript_consequences,
        ).collect()
        clinvar_partners_ht = clinvar_full_ht.semi_join(partner_keys_ht)
        clinvar_rows = clinvar_partners_ht.select(
            clnsig=hl.delimit(clinvar_partners_ht.info.CLNSIG, ", "),
        ).collect()

        vep_lookup = {
            (r.locus.contig, r.locus.position, tuple(r.alleles)): r
            for r in vep_rows
        }
        clinvar_lookup = {
            (r.locus.contig, r.locus.position, tuple(r.alleles)): r.clnsig
            for r in clinvar_rows
        }

        def _pick_tx_csq(tx_csqs, gene_id):
            canonical = [c for c in tx_csqs if c.gene_id == gene_id and c.canonical == 1]
            if canonical:
                return canonical[0]
            any_match = [c for c in tx_csqs if c.gene_id == gene_id]
            return any_match[0] if any_match else None

        enriched_partners = []
        for p in cross_join_partners:
            key = (p.locus.contig, p.locus.position, tuple(p.alleles))
            v = vep_lookup.get(key)
            csq = _pick_tx_csq(v.transcript_consequences, chosen_gene_id) if v else None
            enriched_partners.append(
                hl.Struct(
                    **{f: p[f] for f in p},
                    hgvsp=csq.hgvsp if csq else None,
                    consequence=csq.consequence_terms[0] if csq else None,
                    clinvar_clinical_significance=clinvar_lookup.get(key),
                )
            )

        # 7e: Recompute totals over the full cross-joined partner set.
        total_e = sum(p.expected_in_trans for p in enriched_partners)
        total_o = sum(
            p.observed_in_trans for p in enriched_partners
            if p.observed_in_trans is not None
        )

        # Look up gene_symbol for the candidate from VEP.
        cand_vep_ht = vep_ht.filter(
            (vep_ht.locus == candidate_locus) & (vep_ht.alleles == candidate_alleles)
        )
        candidate_vep_rows = cand_vep_ht.select(
            transcript_consequences=cand_vep_ht.vep.transcript_consequences,
        ).collect()
        candidate_vep_csqs = (
            list(candidate_vep_rows[0].transcript_consequences)
            if candidate_vep_rows else []
        )
        gene_symbol_csq = _pick_tx_csq(candidate_vep_csqs, chosen_gene_id)
        gene_symbol = gene_symbol_csq.gene_symbol if gene_symbol_csq else None

        # Build a synthetic chosen_row for _row_to_demo_json.
        chosen_row = hl.Struct(
            locus=candidate_locus,
            alleles=candidate_alleles_py,
            gene_id=chosen_gene_id,
            gene_symbol=gene_symbol,
            candidate_af=candidate_af_value,
            n_partners=len(enriched_partners),
            total_expected_in_trans=total_e,
            total_observed_in_trans=total_o,
            partners=enriched_partners,
            warnings=existing_warnings,
        )

        demo_dict = _row_to_demo_json(
            chosen_row,
            candidate_variant_id=args.candidate_variant_id,
            dataset=args.dataset_label,
            partner_set=PARTNER_SET_CLINVAR_PLP,
            n_total_partners_in_gene=n_total_in_gene,
            partner_limit=args.partner_limit,
        )

        # Note this is real chets-derived data — no need for the demo-mode warning.
        # The browser-side panel will treat it as if it came from the resolver.
        with hl.hadoop_open(args.output_json, "w") as f:
            json.dump(demo_dict, f, indent=2)
        logger.info("Wrote demo JSON → %s", args.output_json)


def get_argparser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-type", choices=["exomes", "genomes"], default=DEFAULT_DATA_TYPE)
    p.add_argument("--test", action="store_true",
                   help="Filter inputs to TEST_INTERVALS. Implied by --gene.")
    p.add_argument("--gene", choices=list(TEST_INTERVALS),
                   help="Run only this single gene (uses its TEST_INTERVALS interval).")
    p.add_argument("--interval",
                   help="Explicit locus interval (e.g. 'chr19:1-58617616'). "
                        "Bypasses the TEST_INTERVALS lookup for chromosome-scale "
                        "or custom scans. Mutually exclusive with --gene; the "
                        "supplied string is passed to hl.parse_locus_interval "
                        "for annotation filtering, and --output-postfix must "
                        "match a pipeline run that emitted HTs at that scope.")
    p.add_argument("--output-postfix",
                   help="Postfix used by the compute_vp_counts run that produced the HTs.")
    p.add_argument("--max-freq", type=float, default=DEFAULT_MAX_FREQ,
                   help="Upper AF bound (inclusive) applied to BOTH partners "
                        "and candidates. Default: %(default)s. Raise this to "
                        "keep the OE-candidate high-AF expansion; but note "
                        "common-AF candidates dominate top hits with mostly "
                        "LD-driven signal, so the rare-only default gives a "
                        "cleaner statistic for biologically relevant "
                        "recessive-disease-risk cases.")
    p.add_argument("--min-expected-in-trans", type=float, default=1.0,
                   help="Drop rows where total_expected_in_trans is below "
                        "this floor before writing the aggregated HT. Below "
                        "~1 expected pair, depletion is unresolvable "
                        "regardless of the observed count, and these rows "
                        "just pollute the QQ tail. Default: %(default)s. "
                        "Set to 0 to keep everything.")
    p.add_argument("--skip-excluded-genes", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Drop candidates and partners whose gene_id intersects "
                        "the pipeline's --exclude-gene-ids set (loaded via "
                        "get_excluded_genes_ht for the same --output-postfix, "
                        "or from --excluded-genes-ht-path if provided). "
                        "Prevents false-positive top hits with O=0 caused by "
                        "excluded-gene pairs being missing from gt-counts "
                        "upstream. Default on; pass --no-skip-excluded-genes "
                        "to disable (e.g., when no excluded_genes HT exists "
                        "for the postfix).")
    p.add_argument("--excluded-genes-ht-path",
                   help="Explicit path to an excluded_genes HT (keyed by "
                        "gene_id). Overrides the postfix-based lookup; useful "
                        "when --output-postfix has a suffix (e.g. "
                        "'chr19_test.tier3_fixed') for which no excluded HT "
                        "was written but the base postfix's ('chr19_test') "
                        "excluded HT still applies.")
    p.add_argument("--gt-counts-ht-path",
                   help="Explicit path to the variant-pair gt-counts HT. "
                        "Overrides the postfix-based lookup. Use together "
                        "with a per-gene --output-postfix to reuse a shared "
                        "upstream bundle while writing per-gene OE outputs.")
    p.add_argument("--variant-filter-ht-path",
                   help="Explicit path to the variant-filter (candidate) HT. "
                        "Overrides the postfix-based lookup; pairs naturally "
                        "with --gt-counts-ht-path.")
    p.add_argument("--n-samples", type=int, default=730947,
                   help="Total exome sample count behind gt_counts_adj. Default: 730947 (v4.1).")
    p.add_argument("--use-adj", action="store_true", default=True)
    p.add_argument("--no-use-adj", dest="use_adj", action="store_false",
                   help="Use gt_counts_raw instead of gt_counts_adj.")
    p.add_argument(
        "--partner-set",
        choices=["CLINVAR_PLP", "CLINVAR_BLB", "CLINVAR_VUS"],
        default="CLINVAR_PLP",
        help="ClinVar significance category to use as the partner set. "
             "Output HT path includes the lower-case partner-set name.",
    )
    p.add_argument("--phased-ht-path",
                   help="Explicit GCS path to a phased HT (output of "
                        "phase_gnomad.py). When set (or the default "
                        "postfix-based resource exists), p_chet is read "
                        "from em.{adj|raw}.p_chet on this HT instead of "
                        "re-running haplotype_freq_em. Pass an empty string "
                        "or --no-use-phased to force inline EM.")
    p.add_argument("--use-phased", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Whether to reuse the pre-computed phased HT for "
                        "p_chet (default: True). --no-use-phased forces "
                        "re-running EM inline.")
    p.add_argument("--with-ld-adjustment", action=argparse.BooleanOptionalAction,
                   default=False,
                   help="Also compute Rachel Unger's LD-adjusted "
                        "expected-in-trans statistic (product form over "
                        "partners with a D' threshold). Adds "
                        "total_expected_in_trans_ld_adjusted and "
                        "poisson_lower_tail_p_ld_adjusted alongside the "
                        "existing sum-form fields; existing outputs are "
                        "unchanged. Default off (opt-in).")
    p.add_argument("--ld-threshold", type=float, default=0.5,
                   help="|D'| threshold above which a candidate/partner pair "
                        "is treated as fully in-cis (D* = 1 -> partner "
                        "contributes 0 to expected-in-trans). Only used with "
                        "--with-ld-adjustment. Common choices in the "
                        "LD-pruning literature: 0.5 (default, permissive), "
                        "0.8 (strict). Must be in [0, 1].")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing output HT.")
    p.add_argument("--candidate-variant-id",
                   help="e.g. 15-42403721-C-G; if set together with --output-json, "
                   "writes a single-candidate JSON in DepletionResult shape.")
    p.add_argument("--output-json",
                   help="GCS or local path for single-candidate JSON output.")
    p.add_argument("--partner-limit", type=int, default=25,
                   help="Cap partners in the output JSON (top N by AF).")
    p.add_argument("--dataset-label", default="gnomad_r4",
                   help="DatasetId string baked into the JSON ('gnomad_r4' for v4.1).")
    return p


if __name__ == "__main__":
    main(get_argparser().parse_args())

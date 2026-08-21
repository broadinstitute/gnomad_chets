"""
Build a candidate x ClinVar-partner variant pair list for a gene.

Every qualifying variant in a gene is paired against every ClinVar variant of a chosen
significance category in the same gene (P/LP by default), excluding self-pairs. The
result is the input format ``count_supplied_variant_pairs.py`` expects.

This is a different question from the one ``create_vp_list.py`` answers. That builds the
*co-occurrence* list -- pairs actually observed together in at least one sample, which
is what the released co-occurrence tables and the in-trans-OE feature are built on.
This one enumerates pairs by annotation alone, so it includes pairs no sample carries
together. For FKRP that is 146,584 pairs against 10,190, from a slightly smaller variant
set: the co-occurrence list averages 8.5 partners per variant, this averages 146.

Ported from Rachel Ungar's hypomorph_stat_v4 notebook, which built these lists by
re-deriving MAF, a canonical-transcript gene call and ClinVar category straight from the
public v4 sites table. All three are already on the pipeline's sites HT
(:func:`gnomad_chets.v4.create_vp_list.assemble_sites_ht`), so this reads them instead:

* ``af`` -- global adj AF, same source as the notebook's ``freq[0].AF``.
* ``clinvar.is_plp`` / ``is_blb`` / ``is_vus`` -- precomputed by
  :func:`gnomad_chets.v4.utils.clinvar_category_match_expr`, which is stricter than the
  notebook's CLNSIG substring test: it can drop zero-star and conflicting records, and
  its ``blb`` category excludes records that also carry a pathogenic assertion.
  Zero-star and conflicting records are kept in the field but labelled in
  ``clinvar.clinvar_review``; they are excluded here unless
  ``--keep-no-assertion`` / ``--keep-conflicting`` are passed.
* ``vep`` -- for the gene call.

Example
-------
.. code-block:: bash

    hailctl dataproc submit CLUSTER \\
      --pyfiles=/abs/path/to/gnomad_chets \\
      /abs/path/to/gnomad_chets/v4/build_candidate_partner_pairs.py -- \\
      --gene-symbol FKRP \\
      --output gs://my-bucket/FKRP_variant_pairs.ht
"""
import argparse
import logging
from typing import Optional

import hail as hl

from gnomad_chets.v4.create_vp_list import _get_ordered_vp_struct
from gnomad_chets.v4.resources import (
    CLINVAR_CATEGORIES,
    CLINVAR_CATEGORY_FIELD_FMT,
    CLINVAR_CATEGORY_PLP,
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    SITES_FIELD_CLINVAR,
    get_sites_ht,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("build_candidate_partner_pairs")
logger.setLevel(logging.INFO)

DEFAULT_MAX_AF = 0.05

CLINVAR_REVIEW_NO_ASSERTION = "no_assertion"
CLINVAR_REVIEW_CONFLICTING = "conflicting"
"""Flags emitted into the sites HT's ``clinvar.clinvar_review`` set by
:func:`gnomad_chets.v4.utils.clinvar_review_flags_expr`; an empty set means a clean,
reviewed, non-conflicting record."""


def canonical_transcript_gene_expr(
    vep_expr: hl.expr.StructExpression,
) -> hl.expr.StructExpression:
    """
    Pick one gene per variant, from the canonical transcript for its worst consequence.

    Deliberately *not* :func:`gnomad_chets.v4.create_vp_list._get_vep_gene_id_expr`,
    which returns every gene id whose transcript consequence is at least as severe as a
    threshold. That is the right answer for the co-occurrence pipeline, where a variant
    genuinely belongs to each gene it affects and pairs are formed per gene. Here each
    variant needs exactly one gene so the cross is well defined and the output carries a
    scalar ``gene``/``gene_id``, matching the notebook this replaces.

    Preference order: canonical transcript whose consequence terms include the variant's
    worst consequence, then any transcript with that consequence, then the first
    transcript.

    :param vep_expr: The full VEP struct.
    :return: Struct with ``gene_id``, ``gene_symbol`` and ``consequence``.
    """
    tcs = vep_expr.transcript_consequences
    worst = tcs.filter(
        lambda tc: tc.consequence_terms.contains(vep_expr.most_severe_consequence)
    )
    canonical = worst.filter(lambda tc: tc.canonical == 1)
    best = (
        hl.case()
        .when(hl.len(canonical) > 0, canonical[0])
        .when(hl.len(worst) > 0, worst[0])
        .default(tcs[0])
    )

    return hl.struct(
        gene_id=best.gene_id,
        gene_symbol=best.gene_symbol,
        consequence=vep_expr.most_severe_consequence,
    )


def annotate_sites_for_pairing(
    sites_ht: hl.Table,
    category: str = CLINVAR_CATEGORY_PLP,
    max_af: Optional[float] = DEFAULT_MAX_AF,
    *,
    remove_no_assertion: bool = True,
    remove_conflicting: bool = True,
) -> hl.Table:
    """
    Reduce the sites HT to the fields the cross needs, with an ``is_partner`` flag.

    :param sites_ht: Sites Table from :func:`get_sites_ht`.
    :param category: ClinVar category the partner side is drawn from.
    :param max_af: Keep variants with ``0 < af <= max_af``; None keeps all.
    :param remove_no_assertion: Drop zero-star ClinVar records from the partner set.
    :param remove_conflicting: Drop conflicting-interpretation records from the
        partner set.
    :return: Table keyed by (locus, alleles) with gene_id / gene_symbol /
        consequence / af / is_partner.
    """
    ht = sites_ht
    if max_af is not None:
        ht = ht.filter(hl.is_defined(ht.af) & (ht.af > 0) & (ht.af <= max_af))

    ht = ht.filter(hl.is_defined(ht.vep) & (hl.len(ht.vep.transcript_consequences) > 0))
    ht = ht.annotate(_g=canonical_transcript_gene_expr(ht.vep))

    # The sites HT stores relaxed ClinVar membership -- zero-star and conflicting
    # records are kept in is_<category>, with clinvar_review labelling why a record is
    # borderline (empty set == clean). The strict set is therefore
    # ``is_<category> & clinvar_review is empty``; the raw CLNSIG/CLNREVSTAT fields are
    # not carried on the sites HT, so there is nothing to recompute from here.
    cv = ht[SITES_FIELD_CLINVAR]
    is_category = hl.or_else(
        cv[CLINVAR_CATEGORY_FIELD_FMT.format(category=category)], False
    )
    drop_flags = {
        flag
        for flag, drop in (
            (CLINVAR_REVIEW_NO_ASSERTION, remove_no_assertion),
            (CLINVAR_REVIEW_CONFLICTING, remove_conflicting),
        )
        if drop
    }
    if drop_flags:
        review = hl.or_else(cv.clinvar_review, hl.empty_set(hl.tstr))
        is_partner = is_category & (
            hl.len(review.intersection(hl.literal(drop_flags))) == 0
        )
    else:
        is_partner = is_category

    ht = ht.select(
        gene_id=ht._g.gene_id,
        gene_symbol=ht._g.gene_symbol,
        consequence=ht._g.consequence,
        af=ht.af,
        is_partner=is_partner,
    )

    return ht.filter(hl.is_defined(ht.gene_id))


def create_candidate_partner_pair_ht(
    ht: hl.Table, canonical_order: bool = False
) -> hl.Table:
    """
    Cross every variant against every partner in the same gene, dropping self-pairs.

    With ``n`` variants in a gene of which ``p`` are partners, this emits
    ``n * p - p`` rows: each partner pairs with every other variant, and each
    non-partner pairs with every partner. Partner-partner combinations therefore appear
    twice, once in each orientation.

    :param ht: Annotated sites Table from :func:`annotate_sites_for_pairing`.
    :param canonical_order: Reorder each pair so v1 <= v2 by locus position then alt
        allele, matching the released co-occurrence tables, and drop the resulting
        duplicate partner-partner rows. Off by default, which preserves the
        candidate-first orientation (``locus1`` is always the variant under test).
    :return: Table keyed by (locus1, alleles1, locus2, alleles2) with gene / gene_id.
    """
    partners = ht.filter(ht.is_partner)
    by_gene = partners.group_by(partners.gene_id).aggregate(
        partners=hl.agg.collect(
            hl.struct(locus=partners.locus, alleles=partners.alleles)
        )
    )

    pairs = ht.annotate(_partners=by_gene[ht.gene_id].partners)
    pairs = pairs.filter(hl.is_defined(pairs._partners) & (hl.len(pairs._partners) > 0))
    pairs = pairs.explode("_partners")
    # Drop self-pairs: a partner variant does not pair with itself.
    pairs = pairs.filter(
        (pairs.locus != pairs._partners.locus)
        | (pairs.alleles != pairs._partners.alleles)
    )

    pairs = pairs.key_by()
    if canonical_order:
        ordered = _get_ordered_vp_struct(
            hl.struct(locus=pairs.locus, alleles=pairs.alleles),
            hl.struct(locus=pairs._partners.locus, alleles=pairs._partners.alleles),
        )
        pairs = pairs.select(
            locus1=ordered.v1.locus,
            alleles1=ordered.v1.alleles,
            locus2=ordered.v2.locus,
            alleles2=ordered.v2.alleles,
            gene=pairs.gene_symbol,
            gene_id=pairs.gene_id,
        )
        # Partner-partner combinations are emitted from both sides; canonicalising
        # collapses them onto the same key, so dedupe.
        return pairs.key_by("locus1", "alleles1", "locus2", "alleles2").distinct()

    pairs = pairs.select(
        locus1=pairs.locus,
        alleles1=pairs.alleles,
        locus2=pairs._partners.locus,
        alleles2=pairs._partners.alleles,
        gene=pairs.gene_symbol,
        gene_id=pairs.gene_id,
    )

    return pairs.key_by("locus1", "alleles1", "locus2", "alleles2")


def main(args):
    """Build a candidate x ClinVar-partner pair list."""
    hl.init(
        log="/tmp/build_candidate_partner_pairs.log",
        tmp_dir=args.tmp_dir,
        default_reference="GRCh38",
    )

    sites_ht = (
        hl.read_table(args.sites_ht)
        if args.sites_ht
        else get_sites_ht(data_type=args.data_type).ht()
    )
    if args.interval:
        sites_ht = hl.filter_intervals(
            sites_ht,
            [hl.parse_locus_interval(args.interval, reference_genome="GRCh38")],
        )

    ht = annotate_sites_for_pairing(
        sites_ht,
        category=args.partner_category,
        max_af=args.max_af,
        remove_no_assertion=not args.keep_no_assertion,
        remove_conflicting=not args.keep_conflicting,
    )
    if args.gene_symbol:
        symbols = hl.literal(set(args.gene_symbol.split(",")))
        ht = ht.filter(symbols.contains(ht.gene_symbol))
    ht = ht.checkpoint(hl.utils.new_temp_file("sites_for_pairing", "ht"))

    n_variants, n_partners = ht.aggregate(
        (hl.agg.count(), hl.agg.count_where(ht.is_partner))
    )
    logger.info(
        "%d variants, of which %d are ClinVar %s partners.",
        n_variants, n_partners, args.partner_category,
    )
    if n_partners == 0:
        raise ValueError(
            f"No ClinVar {args.partner_category} partners found -- with no partner set "
            "the cross is empty. Check --gene-symbol / --interval, or relax with "
            "--keep-no-assertion / --keep-conflicting."
        )

    pairs = create_candidate_partner_pair_ht(ht, canonical_order=args.canonical_order)
    pairs = pairs.checkpoint(args.output, overwrite=args.overwrite)
    logger.info("Wrote %d variant pairs to %s", pairs.count(), args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--output", required=True, help="Path to write the pair list to."
    )
    parser.add_argument(
        "--tmp-dir",
        required=True,
        help="Hail scratch directory. Must be a GCS path when running on Dataproc.",
    )
    parser.add_argument(
        "--sites-ht",
        help=(
            "Explicit sites HT path. Defaults to the pipeline's sites resource for "
            "--data-type."
        ),
    )
    parser.add_argument(
        "--data-type",
        default=DEFAULT_DATA_TYPE,
        choices=DATA_TYPE_CHOICES,
        help=f"Data type. Default is {DEFAULT_DATA_TYPE}.",
    )
    parser.add_argument(
        "--gene-symbol",
        help=(
            "Comma-separated gene symbol(s) to restrict to, e.g. 'FKRP' or "
            "'FKRP,CAPN3'."
        ),
    )
    parser.add_argument(
        "--interval",
        help=(
            "Locus interval to restrict to before annotating, e.g. "
            "'chr19:46746046-46776988'. Combine with --gene-symbol to keep the scan "
            "bounded on a large sites HT."
        ),
    )
    parser.add_argument(
        "--partner-category",
        default=CLINVAR_CATEGORY_PLP,
        choices=list(CLINVAR_CATEGORIES),
        help=(
            f"ClinVar category the partner side is drawn from. Default is "
            f"{CLINVAR_CATEGORY_PLP}."
        ),
    )
    parser.add_argument(
        "--max-af",
        type=float,
        default=DEFAULT_MAX_AF,
        help=(
            f"Keep variants with 0 < adj AF <= this value. Default {DEFAULT_MAX_AF}. "
            "Pass a negative value to keep every variant."
        ),
    )
    parser.add_argument(
        "--canonical-order",
        action="store_true",
        help=(
            "Order each pair by locus position (tie-broken on alt allele) to match the "
            "released co-occurrence tables, deduping partner-partner rows. Off by "
            "default, which keeps locus1 as the variant under test."
        ),
    )
    parser.add_argument(
        "--keep-no-assertion",
        action="store_true",
        help="Keep zero-star ClinVar records in the partner set.",
    )
    parser.add_argument(
        "--keep-conflicting",
        action="store_true",
        help="Keep conflicting-interpretation ClinVar records in the partner set.",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite the output if it exists."
    )

    parsed = parser.parse_args()
    if parsed.max_af is not None and parsed.max_af < 0:
        parsed.max_af = None
    main(parsed)

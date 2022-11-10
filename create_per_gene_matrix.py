import argparse
from collections import defaultdict
import hail as hl
from itertools import combinations_with_replacement
import logging
from typing import Dict, List, Tuple, Union

from gnomad.resources.grch37.gnomad import public_release
from gnomad_qc.v2.resources import get_gnomad_meta, get_gnomad_data, annotations_ht_path

from resources import (
    full_mt_path,
    get_revel_annotations_path,
    het_hom_per_gene_path,
    phased_vp_count_ht_path,
    vp_ann_ht_path,
    vp_per_gene_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gene_matrix")
logger.setLevel(logging.INFO)

CHET_THRESHOLD = 0.55
SAME_HAP_THRESHOLD = 0.1
ALLELE_FREQUENCY_CUTOFFS = {
    0.05,
    0.02,
    0.015,
    0.01,
    0.005,
    0.001,
    0.0005,
    0.0001,
    0.00005,
    0.00001,
}
STRONG_REVEL_CUTOFF = 0.932
MODERATE_REVEL_CUTOFF = 0.773
SUPPORTING_REVEL_CUTOFF = 0.644
CSQ_CODES = [
    "lof",
    "strong_revel_missense",
    "moderate_to_strong_revel_missense",
    "supporting_to_strong_revel_missense",
    "missense",
    "synonymous",
]
LEN_CSQ_CODES = len(CSQ_CODES)
CSQ_IDX = list(range(LEN_CSQ_CODES))
CSQ_COMBOS = list(combinations_with_replacement(CSQ_CODES, 2))
CSQ_IDX_COMBOS = list(combinations_with_replacement(range(LEN_CSQ_CODES), 2))
LEN_CSQ_COMBOS = len(CSQ_COMBOS)
PHASE_GROUPS = ["chet", "same_hap", "unphased"]


def filter_to_test(
    tables: List[Union[hl.Table, hl.MatrixTable]]
) -> List[Union[hl.Table, hl.MatrixTable]]:
    """
    Filter any Tables in the list to chr20, and any MatrixTables in the list to the first 20 partitions on chr20.

    :param tables: List of Tables and/or MatrixTables to filter for testing
    :return: List of filtered Tables and/or MatrixTables
    """
    test_tables = []
    for t in tables:
        t = hl.filter_intervals(t, [hl.parse_locus_interval("20")])
        if isinstance(t, hl.MatrixTable):
            t = t._filter_partitions(range(20))
        test_tables.append(t)

    return test_tables


def get_csq_pair_combo_map() -> Tuple[List[str], hl.expr.DictExpression]:
    """
    Get list of possible CSQ pairs and a mapping of the CSQ pair code index to all cumulative CSQs that it belongs to.

    For example, 0 (lof_lof) is included in the following cumulative consequences:
        - lof_strong_revel_missense_or_worse
        - lof_moderate_revel_missense_or_worse
        - lof_supporting_revel_missense_or_worse
        - lof_missense_or_worse
        - lof_synonymous_or_worse
        - strong_revel_missense_or_worse_strong_revel_missense_or_worse
        - strong_revel_missense_or_worse_moderate_revel_missense_or_worse
        - strong_revel_missense_or_worse_supporting_revel_missense_or_worse
        - strong_revel_missense_or_worse_missense_or_worse
        - strong_revel_missense_or_worse_synonymous_or_worse
        - moderate_revel_missense_or_worse_moderate_revel_missense_or_worse
        - moderate_revel_missense_or_worse_supporting_revel_missense_or_worse
        - moderate_revel_missense_or_worse_missense_or_worse
        - moderate_revel_missense_or_worse_synonymous_or_worse
        - supporting_revel_missense_or_worse_supporting_revel_missense_or_worse
        - supporting_revel_missense_or_worse_missense_or_worse
        - supporting_revel_missense_or_worse_synonymous_or_worse
        - missense_or_worse_missense_or_worse
        - missense_or_worse_synonymous_or_worse
        - synonymous_or_worse_synonymous_or_worse

    :return: Tuple with the list of possible CSQ pairs and a Dictionary expression mapping a CSQ code index to all
        cumulative CSQs that it belongs to.
    """
    csq_pair_codes = []
    csq_pair_map = defaultdict(list)
    for i, (csq_idx_combo, csq_combo) in enumerate(zip(CSQ_IDX_COMBOS, CSQ_COMBOS)):
        # Add CSQ pair combo to list of all possible pair combinations
        csq_pair_codes.append(
            "_".join(
                [
                    csq.replace("_to_strong", "") + "_or_worse" if csq != "lof" else csq
                    for csq in csq_combo
                ]
            )
        )

        # Loop through all pair combos already added to the dict and determine if the index of the new combo
        # should be included in the old combo index list
        for past_csq_idx_combo, cum_idx in csq_pair_map.items():
            # If worst CSQ in combo is equal to the previously added combo or the best CSQ in the new combo
            # is the same as or better than the previous add index of the new CSQ combo in the'csq_pair_codes' list
            # to the previously added combo index list
            if (past_csq_idx_combo[0] == csq_idx_combo[0]) or (
                past_csq_idx_combo[1] <= csq_idx_combo[1]
            ):
                cum_idx.append(i)
        # All combos include their own index in their index list
        csq_pair_map[csq_idx_combo].append(i)

    for csq_idx1, csq_idx2 in combinations_with_replacement(range(1, LEN_CSQ_CODES), 2):
        csq_idx_combo = (csq_idx1, csq_idx2)
        if csq_idx1 == csq_idx2:
            csq_pair_codes.append("_".join([CSQ_CODES[csq_idx1], CSQ_CODES[csq_idx1]]))
            i = len(csq_pair_codes) - 1
            if csq_idx1 == (LEN_CSQ_CODES - 1):
                csq_pair_map[csq_idx_combo].append(i)
        csq_pair_map[csq_idx_combo].extend(
            list(range(i, LEN_CSQ_COMBOS + LEN_CSQ_CODES - 2))
        )
        i += 1

    csq_pair_map = {k: {csq_pair_codes[i] for i in v} for k, v in csq_pair_map.items()}

    return csq_pair_codes, hl.literal(csq_pair_map)


def get_worst_gene_csq_code_expr_revel(
    vep_expr: hl.expr.StructExpression, revel_expr: hl.expr.StringExpression
) -> hl.expr.DictExpression:
    """
    Filter VEP transcript consequences to canonical and protein coding and annotate with the worst consequence.

    Consequence order worst to least:
        - lof
        - strong_revel_missense
        - moderate_to_strong_revel_missense
        - supporting_to_strong_revel_missense
        - missense
        - synonymous

    :param vep_expr: Expression containing VEP information for the variant.
    :param revel_expr: Expression containing the Revel score for the variant.
    :return: Dictionary expression mapping gene ID to the worst consequence.
    """
    worst_gene_csq_expr = vep_expr.transcript_consequences.filter(
        lambda tc: (tc.canonical == 1) & (tc.biotype == "protein_coding")
    ).map(
        lambda ts: ts.select(
            "gene_id",
            "gene_symbol",
            csq=(
                hl.case(missing_false=True)
                .when(ts.lof == "HC", CSQ_CODES.index("lof"))
                .when(
                    (ts.consequence_terms.all(lambda x: x == "missense_variant") & (revel_expr >= STRONG_REVEL_CUTOFF)),
                    CSQ_CODES.index("strong_revel_missense"),
                )
                .when(
                    (ts.consequence_terms.all(lambda x: x == "missense_variant") & (revel_expr >= MODERATE_REVEL_CUTOFF)),
                    CSQ_CODES.index("moderate_to_strong_revel_missense"),
                )
                .when(
                    (ts.consequence_terms.all(lambda x: x == "missense_variant") & (revel_expr >= SUPPORTING_REVEL_CUTOFF)),
                    CSQ_CODES.index("supporting_to_strong_revel_missense"),
                )
                .when(
                    ts.consequence_terms.all(lambda x: x == "missense_variant"),
                    CSQ_CODES.index("missense"),
                )
                .when(
                    ts.consequence_terms.all(lambda x: x == "synonymous_variant"),
                    CSQ_CODES.index("synonymous"),
                )
                .or_missing()
            ),
        )
    )

    worst_gene_csq_expr = worst_gene_csq_expr.filter(lambda x: hl.is_defined(x.csq))
    worst_gene_csq_expr = worst_gene_csq_expr.group_by(lambda x: x.gene_id)
    worst_gene_csq_expr = worst_gene_csq_expr.map_values(
        lambda x: hl.sorted(x, key=lambda y: y.csq)[0]
    )

    return worst_gene_csq_expr


def compute_from_vp_mt(test: bool, overwrite: bool) -> None:
    """
    Compute chet, same_hap, unphased, and all het_het sample counts by functional csq from variant pair MatrixTable.
    :param test: Whether to filter the variant pair MatrixTable to the first 20 partitions for testing.
    :param overwrite: Whether to overwrite the final output.
    :return: None
    """
    meta = get_gnomad_meta("exomes")
    vp_mt = hl.read_matrix_table(full_mt_path("exomes"))
    revel_ht = hl.read_table(get_revel_annotations_path("exomes"))
    vp_mt = vp_mt.filter_cols(meta[vp_mt.col_key].release)
    ann_ht = hl.read_table(vp_ann_ht_path("exomes"), _n_partitions=20000)
    phase_ht = hl.read_table(phased_vp_count_ht_path("exomes"))
    af_groups = hl.literal(ALLELE_FREQUENCY_CUTOFFS)
    csq_pair_codes, csq_pair_map = get_csq_pair_combo_map()

    if test:
        logger.info(
            "Filtering variant pair MatrixTable to the first 20 partitions of chromosome 20, and the annotation Table "
            "and phase Table to chr20..."
        )
        vp_mt, ann_ht, phase_ht = filter_to_test([vp_mt, ann_ht, phase_ht])

    logger.info(
        "Getting the expression (one for each variant in the pair) for the worst gene consequence of the canonical "
        "transcript per gene using Revel scores for missense variants..."
    )
    vep1_expr = get_worst_gene_csq_code_expr_revel(
        ann_ht.vep1, revel_ht[ann_ht.locus1, ann_ht.alleles1].revel.revel_score
    )
    vep2_expr = get_worst_gene_csq_code_expr_revel(
        ann_ht.vep2, revel_ht[ann_ht.locus2, ann_ht.alleles2].revel.revel_score
    )

    logger.info(
        "Reducing annotation Table to necessary annotations and adding an annotation indicating all relevant "
        "consequence groupings for the variant pair..."
    )

    ann_ht = ann_ht.select(
        "snv1",
        "snv2",
        is_singleton_vp=(ann_ht.freq1["all"].AC < 2) & (ann_ht.freq2["all"].AC < 2),
        popmax_or_global_af=hl.max(ann_ht.popmax1.AF, ann_ht.popmax2.AF, ann_ht.freq1["all"].AF, ann_ht.freq2["all"].AF, filter_missing=True),
        bottlenecked_af=hl.max(ann_ht.freq1["fin"].AF, ann_ht.freq1["asj"].AF, ann_ht.freq1["oth"].AF, ann_ht.freq2["fin"].AF, ann_ht.freq2["asj"].AF, ann_ht.freq2["oth"].AF, filter_missing=True),
        filtered=(hl.len(ann_ht.filters1) > 0) | (hl.len(ann_ht.filters2) > 0),
        vep=vep1_expr.keys()
        .filter(lambda k: vep2_expr.contains(k))
        .map(
            lambda k: vep1_expr[k].annotate(
                csq=csq_pair_map.get(
                    hl.tuple(
                        [
                            hl.min(vep1_expr[k].csq, vep2_expr[k].csq),
                            hl.max(vep1_expr[k].csq, vep2_expr[k].csq),
                        ]
                    )
                )
            )
        ),
    )

    logger.info(
        "Annotating variant pair MatrixTable with the annotation Table and phase information..."
    )
    vp_mt = vp_mt.annotate_rows(
        **ann_ht[vp_mt.row_key], phase_info=phase_ht[vp_mt.row_key].phase_info
    )

    logger.info(
        "Filtering variant pair MatrixTable rows to variant pairs with both PASS variants and "
        "both AF <= 0.05 in bottlenecked populations (fin, asj, oth)..."
    )
    vp_mt = vp_mt.filter_rows(~vp_mt.filtered)
    vp_mt = vp_mt.filter_rows((vp_mt.bottlenecked_af <= 0.05) | hl.is_missing(vp_mt.bottlenecked_af))

    logger.info(
        "Filtering variant pair MatrixTable entries keeping only entries where both variants have a het GT and pass "
        "adj filtering..."
    )
    vp_mt = vp_mt.filter_entries(
        vp_mt.GT1.is_het() & vp_mt.GT2.is_het() & vp_mt.adj1 & vp_mt.adj2
    )
    vp_mt = vp_mt.select_entries(x=True)
    vp_mt = vp_mt.filter_rows(hl.agg.any(vp_mt.x))

    logger.info(
        "Annotate variant pair MatrixTable with all relevant allele frequencies (max variant pair popmax or global AF <= AF cutoff "
        "for the following cutoffs: %s) and the co-occurrence grouping of the variant pair "
        "(chet cutoff: em.adj.p_chet >= %f , same hap cutoff: em.adj.p_chet <= %f)...",
        ALLELE_FREQUENCY_CUTOFFS,
        CHET_THRESHOLD,
        SAME_HAP_THRESHOLD,
    )
    vp_mt = vp_mt.annotate_rows(
        af_cutoff=af_groups.filter(lambda af: vp_mt.popmax_or_global_af <= af),
        chet_group=hl.case()
        .when(
            ~vp_mt.is_singleton_vp
            & (vp_mt.phase_info["all"].em.adj.p_chet >= CHET_THRESHOLD),
            "chet",
        )
        .when(
            ~vp_mt.is_singleton_vp
            & (vp_mt.phase_info["all"].em.adj.p_chet <= SAME_HAP_THRESHOLD),
            "same_hap",
        )
        .default("unphased"),
    )

    logger.info(
        "Exploding rows by VEP annotation (list of structs with 'gene_id', 'gene_symbol', 'csq' annotations)..."
    )
    vp_mt = vp_mt.explode_rows(vp_mt.vep)
    vp_mt = vp_mt.transmute_rows(**vp_mt.vep)
    vp_mt = vp_mt.explode_rows(vp_mt.csq)

    logger.info(
        "Aggregating variant pairs by gene_id, gene_symbol, csq, and af_cutoff and for each sample adding entry "
        "annotations for 'chet', 'same_hap', 'unphased', and 'het_het' if the sample contains any "
        "variant pair within those co-occurrence groupings..."
    )
    vp_mt = (
        vp_mt.group_rows_by("gene_id", "gene_symbol", "csq", "af_cutoff")
        .aggregate(
            **{
                chet_group: hl.agg.filter(
                    vp_mt.chet_group == chet_group, hl.agg.any(vp_mt.x)
                )
                for chet_group in PHASE_GROUPS
            },
            het_het=hl.agg.any(vp_mt.x),
        )
        .key_rows_by("gene_id", "gene_symbol", "csq")
    )
    vp_mt = vp_mt.explode_rows(vp_mt.af_cutoff)

    vp_mt = vp_mt.group_rows_by("gene_id", "gene_symbol", "csq", "af_cutoff").aggregate(
        **{chet_group: hl.agg.any(vp_mt[chet_group]) for chet_group in PHASE_GROUPS},
        het_het=hl.agg.any(vp_mt.het_het),
    )

    logger.info(
        "Writing MatrixTable aggregated at the gene-level for input into combine_and_count_two_rare_var..."
    )
    vp_mt.write('gs://gnomad-sarah/tmp/vp_two_rare_variant.mt')


def compute_from_full_mt(test: bool, overwrite: bool) -> None:
    """
    Compute hom sample counts by functional csq from full variant MatrixTable.
    :param test: Whether to filter the full gnomAD MatrixTable to the first 20 partitions for testing.
    :param overwrite: Whether to overwrite the final output.
    :return: None
    """
    mt = get_gnomad_data("exomes", adj=True, release_samples=True)
    freq_ht = public_release("exomes").ht().select("freq", "popmax")
    revel_ht = hl.read_table(get_revel_annotations_path("exomes"))
    vep_ht = public_release("exomes").ht().select("vep")
    rf_ht = hl.read_table(annotations_ht_path("exomes", "rf"))
    af_groups = hl.literal(ALLELE_FREQUENCY_CUTOFFS)
    csq_pair_codes, csq_pair_map = get_csq_pair_combo_map()

    if test:
        logger.info(
            "Filtering full MatrixTable to the first 20 partitions of chromosome 20, and the frequency Table, "
            "VEP Table, and random forest variant QC Table to chr20..."
        )
        mt, freq_ht, vep_ht, rf_ht = filter_to_test([mt, freq_ht, vep_ht, rf_ht])
    
    logger.info(
        "Getting the expression for the worst gene consequence of the canonical "
        "transcript per gene using Revel scores for missense variants..."
    )
    vep_expr = get_worst_gene_csq_code_expr_revel(
        vep_ht[mt.row_key].vep, revel_ht[mt.row_key].revel.revel_score
    )

    logger.info(
        "Annotating full MatrixTable with VEP, filter information, annotation indicating all relevant consequence "
        "groupings for the variant (matching the annotation used in compute_from_vp_mt), "
        "bottlenecked population AF (fin, asj, oth), and all relevant allele frequencies for filtering "
        "(variant popmax or global AF <= AF cutoff for the following cutoffs: %s)...",
        ALLELE_FREQUENCY_CUTOFFS,
    )

    mt = mt.select_rows(
        vep=vep_expr.values()
        .map(
            lambda k: k.annotate(
                csq=csq_pair_map.get(
                    hl.tuple(
                        [
                            k.csq,
                            k.csq,
                        ]
                    )
                )
            )
        ),
        filters=rf_ht[mt.row_key].filters,
        af_cutoff=af_groups.filter(lambda af: hl.max(freq_ht[mt.row_key].popmax[0].AF,freq_ht[mt.row_key].freq[0].AF, filter_missing=True) <= af),
        bottlenecked_af=hl.max(freq_ht[mt.row_key].freq[3].AF, freq_ht[mt.row_key].freq[4].AF, freq_ht[mt.row_key].freq[5].AF, filter_missing=True)
    )

    logger.info(
        "Filtering MatrixTable to PASS variants with VEP information, bottlenecked population AF <= 0.05,"
        "and popmax or global AF <= %f)...",
        max(ALLELE_FREQUENCY_CUTOFFS),
    )
    mt = mt.filter_rows(
        (hl.len(mt.af_cutoff) > 0) & 
        ((mt.bottlenecked_af <= 0.05) | hl.is_missing(mt.bottlenecked_af)) &
        (hl.len(mt.vep) > 0) & 
        (hl.len(mt.filters) == 0)
    )

    logger.info(
        "Selecting MatrixTable entries with non-reference alleles and annotating "
        "het GT entries ('is_het')..."
    )
    mt = mt.select_entries(mt.GT)

    logger.info(
        "Exploding rows by VEP annotation (list of structs with 'gene_id' and 'gene_symbol' annotations)..."
    )
    mt = mt.explode_rows(mt.vep)
    mt = mt.transmute_rows(**mt.vep)
    mt = mt.explode_rows(mt.csq)

    logger.info(
        "Grouping rows by gene_id, gene_symbol, csq and the set of possible af_cutoffs and adding 'het' and 'hom' "
        "entry annotations per sample..."
    )
    mt = (
        mt.group_rows_by("gene_id", "gene_symbol", "csq", "af_cutoff")
        .aggregate(
            het=hl.agg.any(mt.GT.is_het()),
            hom=hl.agg.any(mt.GT.is_hom_var()),
        )
        .key_rows_by("gene_id", "gene_symbol", "csq")
    )
    mt = mt.explode_rows(mt.af_cutoff)
    mt = mt.group_rows_by("gene_id", "gene_symbol", "csq", "af_cutoff").aggregate(
        het=hl.agg.any(mt.het),
        hom=hl.agg.any(mt.hom),
    )

    logger.info(
        "Writing MatrixTable aggregated at the gene-level for input into combine_and_count_two_rare_var..."
    )
    mt = mt.drop('meta')
    mt.write('gs://gnomad-sarah/tmp/hom_two_rare_variant.mt')


def combine_and_count_two_rare_var(test: bool, overwrite: bool) -> None:
    """
    Combine the variant pair count MatrixTable with the homozygous count Matrixtable generated by compute_from_vp_mt
    and compute_from_full_mt to count individuals with two rare variants by phase.
    :param test: Whether to load MatrixTables from compute_from_vp_mt and compute_from_full_mt filtered to the first 20 partitions for testing.
    :param overwrite: Whether to overwrite the final output.
    :return: None
    """
    
    if test:
        logger.info(
            "Loading gene-level MatrixTables from analysis filtered to the first 20 partitions of chromosome 20..."
        )
        vp_mt = hl.read_matrix_table('gs://gnomad-sarah/tmp/vp_two_rare_variant.mt')
        hom_mt = hl.read_matrix_table('gs://gnomad-sarah/tmp/hom_two_rare_variant.mt')

    logger.info(
        "Merging the entries of the variant pair MatrixTable and the full MatrixTable"
        "to give 'chet', 'same_hap', 'unphased', 'het_het', and 'hom' annotations per sample..."
    )

    merge_mt = hom_mt.annotate_entries(
        chet = vp_mt.index_entries(hom_mt.row_key, hom_mt.col_key).chet,
        same_hap = vp_mt.index_entries(hom_mt.row_key, hom_mt.col_key).same_hap,
        unphased = vp_mt.index_entries(hom_mt.row_key, hom_mt.col_key).unphased,
        het_het = vp_mt.index_entries(hom_mt.row_key, hom_mt.col_key).het_het
    )

    logger.info(
        "Performing final aggregation to get sample counts for co-occurrence groupings ('chet', 'same_hap', "
        "'unphased', 'het_het', 'hom', 'chet_or_hom', and 'any_two') by gene_id, gene_symbol, "
        "csq, and af_cutoff..."
    )

    gene_ht = merge_mt.annotate_rows(
        **{
            f"n_{chet_group}": hl.agg.count_where(merge_mt[chet_group])
            for chet_group in PHASE_GROUPS
        },
        n_het_het=hl.agg.count_where(merge_mt.het_het),
        n_hom=hl.agg.count_where(merge_mt.hom),
        n_chet_or_hom=hl.agg.count_where((merge_mt.hom) | (merge_mt.chet)),
        n_any_two=hl.agg.count_where((merge_mt.het_het) | (merge_mt.hom))
    ).rows()

    #gene_ht = gene_ht.checkpoint(
        #vp_per_gene_path("exomes", test=test), overwrite=overwrite
    #)
    #gene_ht.flatten().export(vp_per_gene_path("exomes", test=test, extension="tsv.gz"))
    if test:
        gene_ht.flatten().export("gs://gnomad-sarah/tmp/two_rare_variants_counts.tsv.gz")


def main(args):
    hl.init(
        log="gs://gnomad-sarah/tmp/gene_matrix_update.log",
        default_reference="GRCh37",
        tmp_dir="gs://gnomad-sarah/tmp/",
    )

    if args.compute_from_vp_mt:
        compute_from_vp_mt(args.test, args.overwrite)
    if args.compute_from_full_mt:
        compute_from_full_mt(args.test, args.overwrite)
    if args.combine_and_count_two_rare_var:
        combine_and_count_two_rare_var(args.test, args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compute-from-vp-mt",
        help="Computes chet, same_hap, and unphased sample counts by functional csq from variant pair MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--compute-from-full-mt",
        help="Computes het and hom sample counts by functional csq from full variant MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--combine-and-count-two-rare-var",
        help="Computes chet, same_hap, unphased, het_het, any_two, hom, and hom_or_chet sample counts by functional csq from variant pair and full MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Computes on first 20 partitions of chr20 only",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)

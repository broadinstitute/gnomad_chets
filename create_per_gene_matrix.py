import argparse
from collections import defaultdict
from itertools import combinations_with_replacement
import logging
from typing import Dict, List, Tuple, Union
from resources import (
    get_revel_annotations_path,
    full_mt_path,
    phased_vp_count_ht_path,
    vp_ann_ht_path,
)
from gnomad.resources.grch37.gnomad import public_release
from gnomad_qc.v2.resources import get_gnomad_meta, get_gnomad_data, annotations_ht_path
import hail as hl


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


def filter_to_test(
    tables: List[Union[hl.Table, hl.MatrixTable]]
) -> List[Union[hl.Table, hl.MatrixTable]]:
    test_tables = []
    for t in tables:
        t = hl.filter_intervals(t, [hl.parse_locus_interval('20')])
        if isinstance(t, hl.MatrixTable):
            t = t._filter_partitions(range(20))
        test_tables.append(t)

    return test_tables


def get_csq_pair_combo_map() -> Tuple[List[str], hl.expr.DictExpression]:
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

    csq_pair_map = {k: [csq_pair_codes[i] for i in v] for k, v in csq_pair_map.items()}

    return csq_pair_codes, hl.literal(csq_pair_map)

def get_cum_csq_map() -> hl.expr.DictExpression:
    # Adds cumulative CSQ codes to indivudual CSQ codes
    cum_csq_codes = CSQ_CODES[0:len(CSQ_CODES)]
    for csq in CSQ_CODES[1:len(CSQ_CODES)]:
        cum_csq_codes.append(
            csq.replace("_to_strong", "") + "_or_worse" if csq != "lof" else csq
        )
    len_cum_csq_codes = len(cum_csq_codes)
    cum_csq_idx = list(range(len_cum_csq_codes))
    cum_csq_map = defaultdict(list)
    # Adds cumulative CSQ index from 'cum_csq_codes' to the cum_csq_map
    for csq in CSQ_IDX:
        for cum_csq in cum_csq_idx:
            if (cum_csq == csq):
                    cum_csq_map[csq].append(cum_csq_idx[cum_csq])
            if (cum_csq < LEN_CSQ_CODES-1) and (cum_csq > csq) and csq != 0:
                    cum_csq_map[csq].append(cum_csq_idx[cum_csq])
            if (cum_csq >= LEN_CSQ_CODES) and (csq < LEN_CSQ_CODES) and (csq <= cum_csq-(LEN_CSQ_CODES-1)):
                    cum_csq_map[csq].append(cum_csq_idx[cum_csq])
                    
    cum_csq_map = {k: [cum_csq_codes[i] for i in v] for k, v in cum_csq_map.items()}

    return hl.literal(cum_csq_map)


def get_worst_gene_csq_code_expr_revel(
    vep_expr: hl.expr.StructExpression, revel_expr: hl.expr.StringExpression
) -> hl.expr.DictExpression:
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
                    (revel_expr >= STRONG_REVEL_CUTOFF),
                    CSQ_CODES.index("strong_revel_missense"),
                )
                .when(
                    (revel_expr >= MODERATE_REVEL_CUTOFF),
                    CSQ_CODES.index("moderate_to_strong_revel_missense"),
                )
                .when(
                    (revel_expr >= SUPPORTING_REVEL_CUTOFF),
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


def compute_from_vp_mt(test: bool, overwrite: bool):
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
        popmax_af=hl.max(ann_ht.popmax1.AF, ann_ht.popmax2.AF, filter_missing=False),
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

    logger.info("Checkpointing annotation Table...")
    ann_ht = ann_ht.checkpoint(
        f"gs://gnomad-tmp/compound_hets/chet_per_gene.annotation{'.test' if test else ''}.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    logger.info("Annotating variant pair MatrixTable with the annotation Table and phase information...")
    vp_mt = vp_mt.annotate_rows(
        **ann_ht[vp_mt.row_key], phase_info=phase_ht[vp_mt.row_key].phase_info
    )

    logger.info("Filtering variant pair MatrixTable to variant pairs with both PASS variants...")
    vp_mt = vp_mt.filter_rows(~vp_mt.filtered)

    logger.info(
        "Filtering variant pair MatrixTable entries keeping only entries where both variants have a het GT and pass "
        "adj filtering..."
    )
    vp_mt = vp_mt.filter_entries(
        vp_mt.GT1.is_het() & vp_mt.GT2.is_het() & vp_mt.adj1 & vp_mt.adj2
    )
    vp_mt = vp_mt.select_entries(x=True)
    vp_mt = vp_mt.filter_rows(hl.agg.any(vp_mt.x))
    vp_mt = vp_mt.checkpoint(
        f"gs://gnomad-tmp/compound_hets/chet_per_gene.filter_entries{'.test' if test else ''}.mt",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )
    
    logger.info(
        "Annotate variant pair MatrixTable with all relevant allele frequencies (max variant pair popmax AF <= AF cutoff "
        "for the following cutoffs: %s) and the co-occurrence grouping of the variant pair "
        "(chet cutoff: em.adj.p_chet >= %f , same hap cutoff: em.adj.p_chet <= %f)...",
        ALLELE_FREQUENCY_CUTOFFS,
        CHET_THRESHOLD,
        SAME_HAP_THRESHOLD
    )
    vp_mt = vp_mt.annotate_rows(
        af_cutoff=af_groups.filter(lambda af: vp_mt.popmax_af <= af),
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
        "Exploding rows by VEP annotation (list of structs with 'gene_id', 'gene_symbol', 'csq' annotations), allele "
        "frequency cutoffs..."
    )
    vp_mt = vp_mt.explode_rows("vep")
    vp_mt = vp_mt.transmute_rows(**vp_mt.vep)
    vp_mt = vp_mt.explode_rows(vp_mt.af_cutoff)
    vp_mt = vp_mt.explode_rows(vp_mt.csq)

    logger.info("Checkpointing exploded variant pair MatrixTable...")
    vp_mt = vp_mt.checkpoint(
        f"gs://gnomad-tmp/compound_hets/chet_per_gene.rows_exploded{'.test' if test else ''}.mt",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    logger.info(
        "Aggregating variant pairs by gene_id, gene_symbol, csq, and af_cutoff and for each sample adding entry "
        "annotations indicating if the sample contains any variant pair (het_het) and annotations for 'chet', "
        "'same_hap', and 'unphased' if the sample contains any variant pair within those co-occurrence groupings..."
    )
    vp_mt = vp_mt.group_rows_by("gene_id", "gene_symbol", "csq", "af_cutoff").aggregate(
        **{
            chet_group: hl.agg.filter(
                vp_mt.chet_group == chet_group, hl.agg.any(vp_mt.x)
            )
            for chet_group in ["chet", "same_hap", "unphased"]
        },
        het_het=hl.agg.any(vp_mt.x),
    )
    
    # NOTE: first part runs with autoscaling and only completes by restarting from this checkpoint and using all workers to finish
    # may want to split the function

    logger.info(
        "Performing final aggregation to get sample counts for co-occurrence groupings ('chet', 'same_hap', and "
        "'unphased') and 'het_het' by gene_id, gene_symbol, csq, and af_cutoff"
    )
    gene_ht = vp_mt.annotate_rows(
        **{
            f"n_{chet_group}": hl.agg.count_where(vp_mt[chet_group])
            for chet_group in ["chet", "same_hap", "unphased"]
        },
        n_het_het=hl.agg.count_where(vp_mt.het_het),
    ).rows()
    gene_ht = gene_ht.annotate(pop="all")

    gene_ht = gene_ht.checkpoint(
        f"gs://gnomad-tmp/compound_hets/chet_per_gene{'.test' if test else ''}.ht",
        overwrite=overwrite,
    )
    gene_ht.flatten().export(
        f"gs://gnomad-sarah/compound_hets/chet_per_gene{'.test' if test else ''}.tsv.gz"
    )

def compute_from_full_mt(test: bool, overwrite: bool):
    mt = get_gnomad_data('exomes', adj=True, release_samples=True)
    freq_ht = public_release('exomes').ht().select('freq','popmax')
    revel_ht = hl.read_table(get_revel_annotations_path('exomes'))
    vep_ht = public_release('exomes').ht().select('vep')
    rf_ht = hl.read_table(annotations_ht_path('exomes', 'rf'))
    af_groups = hl.literal(ALLELE_FREQUENCY_CUTOFFS)
    cum_csq_map = get_cum_csq_map()

    if test:
        logger.info(
            "Filtering variant MatrixTable to the first 20 partitions of chromosome 20, and the frequency Table, "
            "vep Table, and rf Table to chr20..."
        )
        mt, freq_ht, vep_ht, rf_ht = filter_to_test([mt, freq_ht, vep_ht, rf_ht])

    logger.info(
        "Getting the expression for the worst gene consequence of the canonical "
        "transcript per gene using Revel scores for missense variants and adding an annotation "
        "indicating all relevant consequence groupings for the variant..."
    )
    vep_ht = vep_ht.annotate(
        vep=vep_ht.vep.annotate(revel_score=revel_ht[vep_ht.key].revel.revel_score)
    )
    vep_expr = get_worst_gene_csq_code_expr_revel(
        vep_ht.vep, vep_ht.vep.revel_score
    )
    vep_ht = vep_ht.annotate(
        vep=vep_expr.keys()
        .map(
            lambda k: vep_expr[k].annotate(
                csq=cum_csq_map.get(
                            vep_expr[k].csq
                )
            )
        ),
    )

    logger.info("Annotating MatrixTable with the popmax AF, VEP, and filters information...")
    mt = mt.annotate_rows(
        popmax_af=freq_ht[mt.row_key].popmax[0].AF,
        vep=vep_ht[mt.row_key].vep,
        filters=rf_ht[mt.row_key].filters
    )

    logger.info("Filtering MatrixTable to PASS variants <= popmax AF threshold of 0.05...")
    mt = mt.filter_rows(
        (mt.popmax_af <= 0.05) &
        (hl.len(mt.vep) > 0) &
        (hl.len(mt.filters) == 0)
    )

    logger.info(
        "Filtering MatrixTable entries to non-reference alleles and annotating "
        "het GT entries ('is_het')..."
    )
    mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = mt.select_entries(
        is_het=mt.GT.is_het()
    )

    logger.info(
        "Annotating variant MatrixTable with all relevant allele frequencies (variant popmax AF <= AF cutoff "
        "for the following cutoffs: %s)...",
        ALLELE_FREQUENCY_CUTOFFS
    )
    mt = mt.annotate_rows(
        af_cutoff=af_groups.filter(lambda af: mt.popmax_af <= af),
    )

    logger.info(
        "Exploding rows by VEP annotation (list of structs with 'gene_id' and 'gene_symbol' annotations) and allele "
        "frequency cutoffs..."
    )
    mt = mt.explode_rows("vep")
    mt = mt.transmute_rows(**mt.vep)
    mt = mt.explode_rows(mt.csq)
    mt = mt.explode_rows(mt.af_cutoff)

    logger.info("Checkpointing exploded variant MatrixTable...")
    mt = mt.checkpoint(
        f"gs://gnomad-tmp/compound_hets/het_and_hom_per_gene.rows_exploded{'.test' if test else ''}.mt",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    # NOTE: first part runs with autoscaling and only completes by restarting from this checkpoint and using all workers to finish
    # may want to split the function

    logger.info(
        "Performing final aggregation to get sample counts for GT groupings 'het' and 'hom' "
        "by gene_id, gene_symbol, csq, and af_cutoff..."
    )
    mt = mt.group_rows_by("gene_id", "gene_symbol", "csq", "af_cutoff").aggregate(
                het=hl.agg.any(mt.is_het),
                hom=hl.agg.any(~mt.is_het),
    )
    gene_ht = mt.annotate_rows(
        n_het=hl.agg.count_where(mt.het),
        n_hom=hl.agg.count_where(mt.hom),
    ).rows()
    gene_ht = gene_ht.annotate(pop="all")
    
    gene_ht = gene_ht.checkpoint(
        'gs://gnomad-tmp/compound_hets/het_and_hom_per_gene{}.ht'.format(
            '.test' if test else ''
        ),
        overwrite=overwrite
    )
    
    gene_ht.flatten().export(
        'gs://gnomad-sarah/compound_hets/het_and_hom_per_gene{}.tsv.gz'.format(
            '.test' if test else ''
        )
    )

def main(args):
    hl.init(
        log="/tmp/gene_matrix.log",
        default_reference="GRCh37",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    if args.compute_from_vp_mt:
        compute_from_vp_mt(args.test, args.overwrite)
    if args.compute_from_full_mt:
        compute_from_full_mt(args.test, args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compute_from_vp_mt",
        help="Computes chet, same_hap, and unphased sample counts by functional csq from variant pair MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--compute_from_full_mt", 
        help="Computes het and hom sample counts by functional csq from full variant MatrixTable", 
        action="store_true")
    parser.add_argument(
        "--test", 
        help="Computes on first 20 partitions of chr20 only", 
        action="store_true"
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
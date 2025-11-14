from gnomad_hail import *
from gnomad_hail.utils.sample_qc import add_filters_expr

TX_ANN_HT_PATH = 'gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.tx_annotations.ht'

def annotate_tx(mt: hl.MatrixTable, overwrite: bool):
    from tx_annotation import tx_annotate_mt, gtex_v7_tx_summary_mt_path, all_coding_csqs
    mt = mt.select_rows()
    vep = hl.read_table(TX_ANN_HT_PATH)

    mt = mt.annotate_rows(
        vep=vep[mt.row_key].vep
    )

    mt = tx_annotate_mt(
        mt,
        gtex=hl.read_table(gtex_v7_tx_summary_mt_path),
        filter_to_csqs=all_coding_csqs,
        tx_annotation_type="proportion"
    )

    ht = mt.rows()
    ht = ht.annotate(
        vep=ht.vep.drop('colocated_variants')
    )

    ht.write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.tx_annotations.ht', overwrite=overwrite)


def add_vep_missense_metrics_to_tx_annotations(tx_annotation_root: hl.expr.ArrayExpression, vep_root: hl.expr.ArrayExpression) -> hl.expr.ArrayExpression:
    return tx_annotation_root.map(
        lambda tx: hl.cond(
            tx.csq == 'missense_variant',
            tx.annotate(
                **hl.sorted(
                    vep_root.transcript_consequences.filter(
                        lambda vep: (vep.gene_symbol == tx.symbol) & (vep.gene_id == tx.ensg) & (vep.consequence_terms.contains('missense_variant'))
                    ),
                    key=lambda x: x.polyphen_score,
                    reverse=True
                )[0].select('polyphen_prediction', 'sift_prediction')
            ),
            tx.annotate(
                polyphen_prediction=hl.null(hl.tstr),
                sift_prediction=hl.null(hl.tstr)
            )
        )
    )


def main(args):
    mt = hl.read_matrix_table(args.input_mt)
    meta = hl.read_table(args.meta_ht)

    if args.tx_annotations:
        annotate_tx(mt, args.overwrite)

    mt = mt.annotate_cols(
        **meta[mt.col_key]
    )
    mt = mt.filter_cols(hl.len(mt.sample_filters) == 0)
    mt = mt.filter_entries(get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD, haploid_adj_dp=5))
    mt = mt.select_rows('filters')

    filters_expr = hl.bind(
        lambda case_call_rate, control_call_rate: hl.struct(
            case_call_rate=case_call_rate,
            control_call_rate=control_call_rate,
            filters=add_filters_expr(
                filters={
                    "call_rate": (case_call_rate < args.min_call_rate) | (control_call_rate < args.min_call_rate),
                    'AC0': ~hl.agg.any(mt.GT.is_non_ref())
                },
                current_filters=mt.filters
            )
        ),
        hl.agg.filter(mt.is_case, hl.agg.fraction(hl.is_defined(mt.GT))),
        hl.agg.filter(~mt.is_case, hl.agg.fraction(hl.is_defined(mt.GT)))
    )

    gnomad_exomes = get_gnomad_public_data('exomes')[mt.row_key]
    gnomad_genomes = get_gnomad_public_data('genomes')[mt.row_key]

    tx_annotations = hl.read_table(TX_ANN_HT_PATH)
    tx_annotations = tx_annotations.annotate(
        tx_annotation=add_vep_missense_metrics_to_tx_annotations(tx_annotations.tx_annotation, tx_annotations.vep)
    )

    mt = mt.annotate_rows(
        **tx_annotations[mt.row_key],
        **filters_expr,
        gnomad_exomes_popmax=gnomad_exomes.popmax[0],
        gnomad_exomes_freq=gnomad_exomes.freq[0],
        gnomad_genomes_popmax=gnomad_genomes.popmax[0],
        gnomad_genomes_freq=gnomad_genomes.freq[0],
        gnomad_exomes_filters=gnomad_exomes.filters,
        gnomad_genomes_filters=gnomad_genomes.filters,
    )

    mt.rows().write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.annotations.ht', overwrite=args.overwrite)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_mt', help='Location of the input MT', default="gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.mt")
    parser.add_argument('--meta_ht', help='Location of samples meta HT.', default='gs://gnomad/projects/compound_hets/myoseq/sample_qc/MacArthur_LGMD_Callset_Jan2019.full_meta.ht')
    parser.add_argument('--tx_annotations', help='Generates pext/vep annotations', action='store_true')
    parser.add_argument('--min_call_rate', help='Minimum call rate to keep variant', default=0.9, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)



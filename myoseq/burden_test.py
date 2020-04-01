import hail as hl
import argparse
from gnomad.utils.slack import try_slack
from gnomad.utils.annotations import get_adj_expr


def main(args):
    if args.create_gene_sample_mt:
        mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.mt')
        meta = hl.read_table('gs://gnomad/projects/compound_hets/myoseq/sample_qc/MacArthur_LGMD_Callset_Jan2019.full_meta.ht')
        pop_distance = hl.read_table('gs://gnomad-lfran/compound_hets/myoseq/sample_qc/myoseq_pop_distance_to_max_kde.ht')
        variant_annotations_ht = hl.read_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.annotations.ht')
        variant_annotations_ht.drop('was_split', 'a_index')
        mt = mt.annotate_cols(
            **meta[mt.col_key],
            **pop_distance[mt.col_key],
        )
        mt = mt.annotate_rows(
            **variant_annotations_ht[mt.row_key]
        )

        # Filter samples failing QC
        mt = mt.filter_cols(
            (hl.len(mt.sample_filters) == 0) &
            (mt.distance < args.pop_distance) # NFE pop-distance away from densest point in KDE in pc-space (selects only NFEs)
        )
        counts = mt.aggregate_cols(hl.agg.counter(mt.is_case))
        print(f'Found {counts[True]} cases and {counts[False]} controls for gene aggregation.')


        # Filter sites failing QC, without any tx_annotation (i.e. without a protein-coding variant) or too common
        mt = mt.filter_rows(
            (hl.len(mt.filters)==0) &
            hl.is_defined(mt.tx_annotation) &
            (hl.or_else(mt.gnomad_exomes_popmax.AF, hl.or_else(mt.gnomad_genomes_popmax.AF, 0.0)) < args.max_gnomad_af)
        )

        # Keep non-ref entries only
        entries_filter_expr = mt.GT.is_non_ref()
        if not args.raw:
            entries_filter_expr = mt.GT.is_non_ref() & get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD, haploid_adj_dp=5)
        mt = mt.filter_entries(entries_filter_expr)

        # Annotate genes and
        mt = mt.annotate_rows(
            gene=hl.set(mt.tx_annotation.map(
                lambda x: hl.struct(gene_symbol=x.symbol, gene_id=x.ensg)
            ))
        )

        # Aggregate by gene
        mt = mt.explode_rows(mt.gene)
        mt = mt.annotate_rows(tx_annotation=mt.tx_annotation.filter(lambda x: (x.symbol == mt.gene.gene_symbol) & (x.ensg == mt.gene.gene_id)))
        # mt.write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_filtered_gene_exploded.mt', overwrite=True)

        # TODO: Add pext to missense counts

        # mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_filtered_gene_exploded.mt')
        mt = mt.group_rows_by(**mt.gene).aggregate(
            locus_interval=hl.locus_interval(hl.agg.take(mt.locus, 1)[0].contig, hl.agg.min(mt.locus.position), hl.agg.max(mt.locus.position),includes_end=True),
            n_het_lof=hl.agg.count_where(mt.GT.is_het() & mt.tx_annotation.any(lambda x: x.lof == 'HC')),
            n_hom_lof=hl.agg.count_where(mt.GT.is_hom_var() & mt.tx_annotation.any(lambda x: x.lof == 'HC')),
            n_het_lof_pext=hl.agg.count_where(mt.GT.is_het() & mt.tx_annotation.any(lambda x: (x.lof == 'HC') & (x.Muscle_Skeletal >= args.pext_cutoff))),
            n_hom_lof_pext=hl.agg.count_where(mt.GT.is_hom_var() & mt.tx_annotation.any(lambda x: (x.lof == 'HC') & (x.Muscle_Skeletal >= args.pext_cutoff))),
            n_het_missense=hl.agg.count_where(mt.GT.is_het() & mt.tx_annotation.any(lambda x: x.csq == 'missense_variant')),
            n_hom_missense=hl.agg.count_where(mt.GT.is_hom_var() & mt.tx_annotation.any(lambda x: x.csq == 'missense_variant')),
            n_het_damaging_missense=hl.agg.count_where(mt.GT.is_het() & mt.tx_annotation.any(
                lambda x: (x.polyphen_prediction == 'probably damaging') | (x.sift_prediction == 'deleterious')
            )),
            n_hom_damaging_missense=hl.agg.count_where(mt.GT.is_hom_var() & mt.tx_annotation.any(
                lambda x: (x.polyphen_prediction == 'probably damaging') | (x.sift_prediction == 'deleterious')
            )),
            n_het_synonymous=hl.agg.count_where(mt.GT.is_het() & mt.tx_annotation.any(lambda x: x.csq == 'synonymous_variant')),
            n_hom_synonymous=hl.agg.count_where(mt.GT.is_hom_var() & mt.tx_annotation.any(lambda x: x.csq == 'synonymous_variant'))
        ).write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_gene_burden.mt', overwrite=args.overwrite)

    if args.run_burden_tests:
        mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_gene_burden.mt')

        def fet_expr(het_count_exp: hl.expr.Int64Expression, hom_count_expr: hl.expr.Int64Expression):
            return hl.bind(
                lambda x: hl.struct(
                    counts=x,
                    dominant=
                        hl.fisher_exact_test(
                        x[0][0],
                        x[0][1] + x[0][2],
                        x[1][0],
                        x[1][1] + x[1][2]
                    ),
                    recessive=hl.fisher_exact_test(
                        x[0][0] + x[0][1],
                        x[0][2],
                        x[1][0] + x[1][1],
                        x[1][2]
                    )
                ),
                hl.bind(
                    lambda x: [
                        [hl.int32(hl.cond(x.contains(False), x[False].get(0, 0), 0)), hl.int32(hl.cond(x.contains(False), x[False].get(1, 0), 0)), hl.int32(hl.cond(x.contains(False), x[False].get(2, 0), 0))],
                        [hl.int32(hl.cond(x.contains(True), x[True].get(0, 0), 0)), hl.int32(hl.cond(x.contains(True), x[True].get(1, 0), 0)), hl.int32(hl.cond(x.contains(True), x[True].get(2, 0), 0))],
                    ],
                    hl.agg.group_by(mt.is_case, hl.agg.counter(hl.min(2, het_count_exp + 2 * hom_count_expr)))
                )
            )

        mt = mt.annotate_rows(
            **{
                'lof': fet_expr(mt.n_het_lof, mt.n_hom_lof),
                'lof_pext': fet_expr(mt.n_het_lof_pext, mt.n_hom_lof_pext),
                'lof_missense': fet_expr(mt.n_het_lof + mt.n_het_missense, mt.n_het_lof + mt.n_hom_missense),
                'lof_damaging_missense': fet_expr(mt.n_het_lof + mt.n_het_damaging_missense, mt.n_het_lof + mt.n_hom_damaging_missense),
                'synonymous': fet_expr(mt.n_het_synonymous, mt.n_hom_synonymous)
            }
        )

        mt.write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_gene_burden_tests.mt', overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--create_gene_sample_mt', help='Creates a gene x sample MT with counts', action='store_true')
    parser.add_argument('--run_burden_tests', help='Runs burden test on the gene x sample MT', action='store_true')
    parser.add_argument('--pop_distance', help='Population distance in PCA space for NFEs', default=0.07, type=float)
    parser.add_argument('--max_gnomad_af', help='Maximum gnmomAD AF (popmax) to consider.', default=0.01, type=float)
    parser.add_argument('--pext_cutoff', help='Minimum pext to consider.', default=0.1, type=float)
    parser.add_argument('--raw', help='If set, all genotypes are used regardless of their quality, otherwise only adj are used.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


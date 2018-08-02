from gnomad_hail import *
# from gnomad_hail.resources.phasing import *
import hail as hl
# from gnomad_hail.utils.constants import CSQ_CODING_HIGH_IMPACT, CSQ_CODING_MEDIUM_IMPACT, POLYPHEN_PREDICTIONS, SIFT_PREDICTIONS


def variant_pair_mt2(mt, row_groups):

    # Convert entries to minimal representation
    mt = mt.select_entries(
        GT=hl.or_missing(
            hl.is_missing(mt.GT) | mt.GT.is_non_ref(),
            hl.struct(
                is_missing=hl.or_missing(hl.is_missing(mt.GT), hl.struct()),
                is_het=hl.or_missing(mt.GT.is_het(), hl.struct()),
                adj=hl.or_missing(mt.adj, hl.struct())
            )
        )
    )

    #Try to join table to itself?

    et = mt.select_cols().select_rows(*row_groups).entries()
    et.describe()
    et = et.group_by(*row_groups, *mt.col_key).aggregate( #Make sample ID instead of s // possibly for the gene too
        vgt=(hl.agg.collect(hl.struct(locus=et.locus, alleles=et.alleles, GT=et.GT)))
    )

    et = et.annotate(
        vgt=hl.range(0, hl.len(et.vgt))
                      .flatmap(lambda i1: hl.range(i1+1, hl.len(et.vgt))
                               .map(lambda i2: hl.struct(v1=et.vgt[i1], v2=et.vgt[i2])))
    )

    et.describe()
    et = et.explode(et.vgt)
    et = et.key_by(*mt.col_key, locus1=et.vgt.v1.locus, alleles1=et.vgt.v1.alleles, locus2=et.vgt.v2.locus, alleles2=et.vgt.v2.alleles)
    et = et.select(gt1=et.vgt.v1.GT, gt2=et.vgt.v2.GT)

    vp_mt = et.to_matrix_table(row_key=['locus1','alleles1','locus2','alleles2'], col_key=['s'])
    return vp_mt


def variant_pair_mt(mt, row_groups):

    def variant_pairs_ht(mt, row_groups):
        mt = mt.add_row_index(name='_row_idx')
        mt = mt.add_col_index(name='_col_idx')

        mt = mt.group_rows_by(*row_groups).aggregate(
            variant_pairs_entry=hl.set(
                hl.agg.collect(
                    hl.tuple([mt._col_idx, mt._row_idx])
                )
                    # [col_idx, [row_idx]]
                    .group_by(lambda x: x[0])
                    # [[col_idx, row_idx), (col_idx, row_idx), ...],  ...]
                    .values()
                    # [[row_idx, row_idx, ...],  ...]
                    .map(lambda x: x.map(lambda y: y[0]))
                    .flatmap(lambda x: hl.range(0, hl.len(x))
                               .flatmap(lambda i1: hl.range(i1 + 1, hl.len(x))
                                        .map(lambda i2: hl.tuple([x[i1], x[i2]]))))
            )
        )

        mt.describe()

        ht = mt.annotate_rows(
            variant_pairs=hl.agg.take(mt.variant_pairs_entry, 1)[0]
        ).rows()

        ht = ht.explode('variant_pairs')
        ht = ht.key_by(v1_idx=ht.variant_pairs[0], v2_idx=ht.variant_pairs[1])
        return ht.select()

    pass




def per_gene_bruden(mt):

    #v1, v2 => gt1, gt2
    mt = mt.select_entries(

    )




def continue_work(args):
    data_type = 'exomes' if args.exomes else 'genomes'
    agg_rows_name: str = 'agg_rows'
    agg_cols_name: str = 'agg_cols'
    row_groups = ['gene']
    tmp_dir = 'gs://gnomad/projects/compound_hets/grouped_sparse_genotype_matrix_by_gene'

    mt = hl.read_matrix_table(f'{tmp_dir}/cols_agg.mt')

    mt = mt.group_rows_by(*row_groups).aggregate(
        # [row_idx, [(col_idx, entry)]]
        non_refs=hl.agg.collect(hl.tuple([mt._row_idx, mt.non_refs.map(lambda x: hl.tuple([x[0], hl.struct(
            is_missing=hl.or_missing(hl.is_missing(x[1].GT), hl.struct()),
            is_het=hl.or_missing(x[1].GT.is_het(), hl.struct()),
            adj=hl.or_missing(x[1].adj, hl.struct())
        )]))]))
    )

    mt.write(f'{tmp_dir}/rows_agg.mt', overwrite=True)
    mt = hl.read_matrix_table(f'{tmp_dir}/rows_agg.mt')

    rows_ann = hl.read_table(f'{tmp_dir}/rows_ann.ht')

    mt = mt.annotate_rows(
        **rows_ann[mt.row_key]
    )

    mt = mt.select_entries(
        **{
            # {row_idx: [(col_idx, call)}
            agg_rows_name + "_non_ref": hl.dict(mt.non_refs),
            # {col_idx: [(row_idx, call)}
            agg_cols_name + "_non_ref": (
                mt.non_refs
                    .flatmap(
                    lambda x: x[1].map(
                        lambda y: hl.tuple([
                            y[0],
                            hl.tuple([x[0], y[1]])
                        ])
                    )
                )
                    .group_by(lambda x: x[0])
                    .map_values(lambda x: x.map(lambda y: y[1]))
            )
        }
    )

    write_temp_gcs(mt, grouped_sparse_genotype_matrix_path(data_type), overwrite=args.overwrite)


def main(args):

    data_type = 'exomes' if args.exomes else 'genomes'
    
    mt = get_gnomad_data(data_type, release_samples=True)

    if args.chrom20:
        print("Selecting chrom20")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval('20')])

    mt = mt.select_cols(pop=mt.meta.pop)
    mt = mt.select_rows('a_index')

    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    #cadd = hl.read_table(cadd_ht_path)
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))

    consequences_to_keep = hl.literal(set(["inframe_insertion", "inframe_deletion", "missense_variant"]))
    # consequences_index = hl.literal({pred: i for i, pred in enumerate(CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT)})
    # polyphen_pred_index = hl.literal({pred: i for i, pred in enumerate(POLYPHEN_PREDICTIONS)})
    # sift_pred_index = hl.literal({pred: i for i, pred in enumerate(SIFT_PREDICTIONS)})

    mt = mt.annotate_globals(
        vep_consequences_to_keep=consequences_to_keep
        # polyphen_pred=POLYPHEN_PREDICTIONS,
        # sift_pred=SIFT_PREDICTIONS
    )

    mt = mt.select_rows(
        vep=(
            vep_ht[mt.row_key].vep.transcript_consequences
            .filter(
                lambda tc: (tc.biotype == 'protein_coding') &
                           (tc.allele_num == mt.a_index) &
                           ((hl.is_defined(tc.lof) & (tc.lof == "HC")) | tc.consequence_terms.any(lambda c: consequences_to_keep.contains(c)))
            )
            .map(
                lambda x: x.select(
                    'gene_symbol',
                    'transcript_id',
                    lof=hl.or_missing(x.lof == "HC", hl.struct()),
                    #lof=hl.cond(hl.is_defined(x.lof), hl.struct(), hl.null(hl.tstruct())),
                    # polyphen_prediction=polyphen_pred_index[x.polyphen_prediction],
                    # sift_prediction=sift_pred_index[x.sift_prediction],
                    canonical=hl.or_missing(x.canonical == 1, hl.struct())
                )
            )
            .group_by(lambda x: x.gene_symbol)
            .values()
        ),
        #cadd=hl.int(cadd[mt.row_key].phred_score),
        af=hl.float32(hl.find(lambda x: x.meta.get('group') == 'adj', freq[mt.row_key].freq).AF[1])
    )

    mt = mt.filter_rows(hl.is_defined(mt.vep) & (hl.len(mt.vep) > 0) & (mt.af > 0) & (mt.af <= 0.01))
    mt = mt.explode_rows(mt.vep)

    mt = mt.annotate_rows(
        gene=mt.vep[0].gene_symbol #,
        # vep=mt.vep.map(lambda x:
        #                x.select(
        #                    'transcript_id', 'lof', 'canonical'
        #                )
        #                )
    )

    #mt = mt.transmute_rows(**mt.vep)
    mt.describe()
    mt = variant_pair_mt2(mt, ['gene'])
    mt.write('gs://gnomad/projects/compound_hets/gnomad_rare_vp{}.mt'.format('_chrom20' if args.chrom20 else ''), overwrite=args.overwrite)
    #mt = create_grouped_sparse_genotypes_matrix_table(mt, ['gene'], ['pop'], tmp_dir='gs://gnomad/projects/compound_hets/grouped_sparse_genotype_matrix_by_gene')
    #write_temp_gcs(mt, grouped_sparse_genotype_matrix_path(data_type), overwrite=args.overwrite)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chrom20', help='Only run chrom20', action='store_true')


    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
        #continue_work(args)




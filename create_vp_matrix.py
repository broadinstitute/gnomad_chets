from gnomad_hail import *
import hail as hl

def create_variant_pair_mt(mt: hl.MatrixTable, row_groups: List[str]):
    """
    Create a variant-pair MatrixTable containing all variant-pairs that appear both in the same individual within the row group.
    E.g., if the row group is a gene, this creates a variant-pair MT with all variant-pairs in a given gene.

    :param MatrixTable mt: Input MatrixTable
    :param list of str row_groups: Row annotations for delimiting variant-pairs territory
    :return: Variant-pair MT
    :rtype: MatrixTable
    """

    mt = mt.filter_entries(mt.GT.is_non_ref())
    et = mt.select_cols().select_rows(*row_groups).entries()
    et = et.group_by(*row_groups, *mt.col_key).aggregate(
        vgt=(hl.agg.collect(hl.struct(locus=et.locus, alleles=et.alleles, gt=hl.struct(GT=et.GT, adj=et.adj))))
    )

    et = et.annotate(
        vgt=hl.range(0, hl.len(et.vgt))
                      .flatmap(lambda i1: hl.range(i1+1, hl.len(et.vgt))
                               .map(lambda i2: hl.struct(v1=et.vgt[i1], v2=et.vgt[i2])))
    )

    et = et.explode(et.vgt)
    et = et.key_by(*mt.col_key, locus1=et.vgt.v1.locus, alleles1=et.vgt.v1.alleles, locus2=et.vgt.v2.locus, alleles2=et.vgt.v2.alleles)
    et = et.select(gt1=et.vgt.v1.gt.GT, adj1=et.vgt.v1.gt.adj, gt2=et.vgt.v2.gt.GT, adj2=et.vgt.v2.gt.adj)
    et = et.distinct()

    vp_mt = et.to_matrix_table(row_key=['locus2','alleles2','locus1','alleles1'], col_key=['s'])
    return vp_mt


def filter_freq_and_csq(mt: hl.MatrixTable, data_type: str, max_freq: float, least_consequence: str):
    """
    Filters MatrixTable to include variants that:
    1. Have a global AF <= `max_freq`
    2. Have a consequence at least as severe as `least_consequence` (based on ordering from CSQ_ORDER)

    :param MatrixTable mt: Input MT
    :param str data_type: One of 'exomes' or 'genomes'
    :param float max_freq: Max. AF to keep
    :param str least_consequence: Least consequence to keep.
    :return: Filtered MT
    :rtype: MatrixTable
    """

    # mt = mt.select_cols(pop=mt.meta.pop)
    mt = mt.select_rows('a_index')

    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    vep_consequences = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1]))

    mt = mt.select_rows(
        vep=(
            hl.set(
                vep_ht[mt.row_key].vep.transcript_consequences
                    .filter(
                    lambda tc: (tc.biotype == 'protein_coding') &
                               (tc.allele_num == mt.a_index) &
                               (tc.consequence_terms.any(lambda c: vep_consequences.contains(c)))
                )
                    .map(lambda x: x.gene_id)
            )
        ),
        af=hl.float32(freq[mt.row_key].freq[0].AF)
    )

    mt = mt.filter_rows(hl.is_defined(mt.vep) & (hl.len(mt.vep) > 0) & (mt.af > 0) & (mt.af <= max_freq))
    mt = mt.explode_rows(mt.vep)
    mt = mt.rename({'vep': 'gene_id'})
    return mt


def main(args):

    def get_out_path(stage: str = ''):
        return 'gs://gnomad{}/compound_hets/{}{}{}_{}_{}_vp{}.mt'.format(
            '-tmp/' if stage == 'mini_mt' else '/projects',
            data_type,
            '_pbt' if args.pbt else '',
            f'_{stage}' if stage else '',
            args.max_freq,
            args.least_consequence,
            f'_chrom{args.chrom}' if args.chrom else '')

    data_type = 'exomes' if args.exomes else 'genomes'

    if args.create_mini_mt:
        if args.pbt:
            mt = hl.read_matrix_table(pbt_phased_trios_mt_path('genomes'))
            mt = mt.filter_entries(mt.PBT_GT.is_non_ref())
            mt = mt.annotate_entries(GT=mt.PBT_GT)
        else:
            mt = get_gnomad_data(data_type, non_refs_only=True)
            mt = mt.filter_cols(mt.meta.high_quality)
            mt = mt.select_entries('GT', 'adj', 'is_missing')

        if args.chrom:
            print(f"Selecting chrom {args.chrom}")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.chrom)])

        mt = filter_freq_and_csq(mt, data_type, args.max_freq, args.least_consequence)
        mt.write(get_out_path('mini_mt'), overwrite=args.overwrite)

    if args.create_non_missing_vp:
        mt = hl.read_matrix_table(get_out_path('mini_mt'))
        mt = create_variant_pair_mt(mt, ['gene_id'])
        mt.write(get_out_path('non_missing'), overwrite=args.overwrite)

    if args.create_full_vp:
        non_missing_mt = hl.read_matrix_table(get_out_path('non_missing'))
        missing_mt = hl.read_matrix_table(get_out_path('mini_mt'))
        non_missing_mt = non_missing_mt.key_rows_by('locus2', 'alleles2')
        non_missing_mt = non_missing_mt.select_entries(**non_missing_mt.entry, is_missing2=missing_mt[non_missing_mt.row_key, non_missing_mt.col_key].is_missing)
        non_missing_mt = non_missing_mt.key_rows_by('locus1', 'alleles1')
        non_missing_mt = non_missing_mt.select_entries(**non_missing_mt.entry, is_missing1=missing_mt[non_missing_mt.row_key, non_missing_mt.col_key].is_missing)
        non_missing_mt.write(get_out_path(), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--pbt', help='Runs on PBT-phased data instead of the entire gnomAD. Note that the PBT_GT will be renamed as GT',
                        action='store_true')
    parser.add_argument('--create_mini_mt', help='Creates a filtered, minimal MT that is then used to create the VP MT.',
                        action='store_true')
    parser.add_argument('--create_non_missing_vp', help='Creates the VP MT containing non-missing GTs only.', action='store_true')
    parser.add_argument('--create_full_vp', help='Creates the VP MT.', action='store_true')
    parser.add_argument('--least_consequence', help='Includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: 3_prime_UTR_variant)',
                        default='3_prime_UTR_variant')
    parser.add_argument('--max_freq', help='If specified, maximum global adj AF for genotypes table to emit. (default: 0.02)', default=0.02, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chrom', help='Only run on given chromosome')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)




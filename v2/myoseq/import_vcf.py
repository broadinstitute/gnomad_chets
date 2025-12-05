from gnomad_hail import *
import hail as hl


def main(args):
    output_prefix = args.vcf.split(".vcf")[0]

    if args.import_vcf:
        mt = hl.import_vcf(args.vcf, force_bgz=True)
        mt.write(f'{output_prefix}.unsplit.mt', overwrite=args.overwrite)
        mt = hl.read_matrix_table(f'{output_prefix}.unsplit.mt')
        mt = hl.split_multi_hts(mt)
        mt.write(f'{output_prefix}.mt', overwrite=args.overwrite)

    # if args.write_hardcalls:
    #     mt = hl.read_matrix_table(f'{output_prefix}.unsplit.mt')
    #     mt = annotate_adj(mt.select_cols(sex=ht[hl.literal(data_type), mt.s].sex))
    #     mt = mt.select_entries(GT=hl.case(missing_false=True).when(hl.call(mt.PGT[0], mt.PGT[1]) == mt.GT, mt.PGT).default(mt.GT),
    #                            PID=mt.PID, adj=mt.adj)
    #     mt = adjust_sex_ploidy(mt, mt.sex)
    #     mt = mt.select_cols().naive_coalesce(10000)
    #     mt.write(f'{output_prefix}.hardcalls.unsplit.mt', args.overwrite)
    #
    #     mt = get_gnomad_data(f'{output_prefix}.hardcalls.unsplit.mt', split=False, meta_root=None)
    #     mt = hl.split_multi_hts(mt)
    #     mt.write(get_gnomad_data_path(f'{output_prefix}.hardcalls.mt', hardcalls=True, split=True), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', help='Location of the VCF to import', default="gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.gz")
    parser.add_argument('--import_vcf', help='Imports the input VCF.', action='store_true')
    parser.add_argument('--write_hardcalls', help='Generates hardcalls data from the original import. NOT IMPLEMENTED/USED', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
from gnomad_hail import *
import hail as hl


def main(args):
    if args.import_vcf:
        mt = hl.import_vcf('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.gz', force_bgz=True)
        mt.write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.mt', overwrite=args.overwrite)
        mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.mt')
        mt = hl.split_multi_hts(mt)
        mt.write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.split.mt', overwrite=args.overwrite)

    if args.write_hardcalls:
        mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.mt')
        mt = annotate_adj(mt.select_cols(sex=ht[hl.literal(data_type), mt.s].sex))
        mt = mt.select_entries(GT=hl.case(missing_false=True).when(hl.call(mt.PGT[0], mt.PGT[1]) == mt.GT, mt.PGT).default(mt.GT),
                               PID=mt.PID, adj=mt.adj)
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(10000)
        mt.write(get_gnomad_data_path(data_type, hardcalls=True, split=False), args.overwrite)

    if args.split_hardcalls:
        mt = get_gnomad_data(data_type, split=False, meta_root=None)
        mt = hl.split_multi_hts(mt)
        mt.write(get_gnomad_data_path(data_type, hardcalls=True, split=True), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--import_vcf', help='Imports the MYOSEQ data from  VCF.', action='store_true')
    parser.add_argument('--write_hardcalls', help='Generates hardcalls data.', action='store_true')
    parser.add_argument('--split_hardcalls', help='Generates split hardcalls data.', action='store_true')
    parser.add_argument('--split', help='Generates hardcalls.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
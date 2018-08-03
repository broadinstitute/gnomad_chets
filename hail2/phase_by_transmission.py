from gnomad_hail import *
import hail as hl


def main(args):
    data_types = []
    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:
        mt = get_gnomad_data(data_type, split=False)
        ped = hl.Pedigree.read(fam_path(data_type), delimiter='\\t')

        ped_samples = hl.literal(set([s for trio in ped.complete_trios() for s in [trio.s, trio.pat_id, trio.mat_id]]))

        mt = mt.filter_cols(ped_samples.contains(mt.s) & mt.meta.high_quality)
        mt = mt.select_cols().select_rows()
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

        tm = hl.trio_matrix(mt, ped)
        tm = hl.experimental.phase_trio_matrix_by_transmission(tm)

        pmt = explode_trio_matrix(tm)
        pmt.write(pbt_phased_trios_mt_path(data_type, split=False), stage_locally=True, overwrite=args.overwrite)

        pmt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type, split=False))
        pmt = split_multi_dynamic(pmt)
        pmt.write(pbt_phased_trios_mt_path(data_type), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)




import argparse
from

#create 
def main(args):
    tmp_dir=args.tmp_dir
    data_type = 'exomes' if args.exomes else 'genomes'

    hl.init(log="{tmp_dir}/hail_vp_matrix.log")


if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument(
        '--exomes', 
        action='store_true',
        help='Run on exomes. One and only one of --exomes or --genomes is required.')
    data_grp.add_argument('--genomes',
        action='store_true',
        help='Run on genomes. One and only one of --exomes or --genomes is required.')
    parser.add_argument(
        '--testing',
        action='store_true',
        help='if you are testing the pipeline, for developers only')
    parser.add_argument(
        '--tmp-dir',
        required=True,
        help='Temporary directory for intermediate files.'
    )

    args = parser.parse_args()
    main(args)

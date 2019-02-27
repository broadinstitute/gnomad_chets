from gnomad_hail import *
from resources import *
import pandas as pd

def main(args):

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, True, args.least_consequence, args.max_freq]

    pbt = hl.read_table(pbt_phase_count_ht_path(*path_args))
    # Apply relevant filters
    if not args.keep_sex_chrom:
        autosomes = hl.parse_locus_interval('1-22')
        pbt = pbt.filter(autosomes.contains(pbt.locus1))  # locus1 and locus2 are always on the same contig

    if not args.keep_filtered:
        pbt = pbt.filter((hl.len(pbt.filters1) == 0) & (hl.len(pbt.filters2) == 0))
        logger.debug(f'Rows remaining after keeping non-filtered variants: {pbt.count()}')

    if not args.keep_raw:
        pbt = pbt.filter(pbt.adj1 & pbt.adj2)
        logger.debug(f'Rows remaining after keeping adj-only: {pbt.count()}')

    if not args.keep_other_pop:
        pbt = pbt.filter(pbt.is_defined(pbt.pop) & (pbt.pop != 'oth'))
        logger.debug(f'Rows remaining after removing oth samples: {pbt.count()}')




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                          action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                          action='store_true')
    parser.add_argument('--least_consequence', help=f'Least consequence for the input (just to get the path right). (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--max_freq', help=f'Maximum global adj AF for the input (just to get the path right). (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', help='Path to write output tsv file')

    parser.add_argument('--keep_filtered', help="When set, keeps filtered variants. (by default if either of the two variants was filtered it is excluded.", action='store_true')
    parser.add_argument('--keep_raw', help="When set, keeps non-ADJ GTs (by default, all genotypes in the trio need to be ADJ to keep).", action='store_true')
    parser.add_argument('--keep_other_pop', help="When set, keeps trios that don't have an assigned population (either 'oth' or mixed).", action='store_true')
    parser.add_argument('--keep_sex_chrom', help="When set, keeps sex chromosomes.", action='store_true')
    parser.add_argument('--include_discordant_vps', help="When set, includes variant-pairs for which different trios suggest different phasing within the same population.", action='store_true')
    parser.add_argument('--max_pop_freq', help="Maximum population AF to include.", default=0.01, type=float)

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
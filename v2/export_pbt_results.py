from gnomad.utils.slack import try_slack
from resources import *
import hail as hl
import logging
import argparse

logger = logging.getLogger("export_pbt_results")

def main(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, True, args.least_consequence, args.max_freq]
    ht = hl.read_table(pbt_trio_et_path(*path_args))

    # Apply relevant filters
    if not args.export_sex_chrom:
        autosomes = hl.parse_locus_interval('1-22')
        ht = ht.filter(autosomes.contains(ht.locus1)) # locus1 and locus2 are always on the same contig

    if not args.export_filtered:
        ht = ht.filter((hl.len(ht.filters1) == 0) & (hl.len(ht.filters2) == 0))
        logger.debug(f'Rows remaining after keeping non-filtered variants: {ht.count()}')

    if not args.export_raw:
        ht = ht.filter(ht.adj1 & ht.adj2)
        logger.debug(f'Rows remaining after keeping adj-only: {ht.count()}')

    if not args.export_other_pop:
        ht = ht.filter(hl.is_defined(ht.pop) & (ht.pop != 'oth'))
        logger.debug(f'Rows remaining after removing oth samples: {ht.count()}')

    pbt_vp_summary = hl.read_table(pbt_phase_count_ht_path(*path_args))
    pbt_vp_summary = pbt_vp_summary.filter(pbt_vp_summary.adj.n_same_hap + pbt_vp_summary.adj.n_chet > 0)
    indexed_pbt_vp_summary = pbt_vp_summary[ht.key]
    discordant_expr = (indexed_pbt_vp_summary.adj.n_same_hap > 0) & (indexed_pbt_vp_summary.adj.n_chet > 0)
    if args.exclude_discordant_vps:
        ht = ht.filter(discordant_expr)
    else:
        ht = ht.annotate(
            trio_phase_discordant=discordant_expr
        )

    ht = ht.filter((ht.pop_freq1.af <= args.max_pop_freq) & (ht.pop_freq2.af <= args.max_pop_freq))
    logger.debug(f'Rows remiaining after removing sites with freq > {args.max_pop_freq}: {ht.count()}')

    # Annotate phase from gnomAD
    vp_ht = hl.read_table(phased_vp_count_ht_path(*path_args))
    vp_ht = vp_ht.select(
        'em',
        'em_plus_one',
        'likelihood_model',
        'singlet_het_ratio'
    )
    ht = ht.annotate(
        **vp_ht[ht.key]
    )

    ht = ht.flatten()
    ht.export(args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                          action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                          action='store_true')
    parser.add_argument('--export_filtered', help="When set, also exports filtered variants. (by default if either of the two variants was filtered it is excluded.", action='store_true')
    parser.add_argument('--export_raw', help="When set, also exports non-ADJ GTs (by default, all genotypes in the trio need to be ADJ to export).", action='store_true')
    parser.add_argument('--export_other_pop', help="When set, also exports trios that don't have an assigned population (either 'oth' or mixed).", action='store_true')
    parser.add_argument('--export_sex_chrom', help="When set, also exports sex chromosomes.", action='store_true')
    parser.add_argument('--exclude_discordant_vps', help="When set, excludes variant-pairs for which different trios suggest different phasing within the same population.", action='store_true')
    parser.add_argument('--max_pop_freq', help="Maximum population AF to include.", default=0.01, type=float)
    parser.add_argument('--least_consequence', help=f'Least consequence for the input (just to get the path right). (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--max_freq', help=f'Maximum global adj AF for the input (just to get the path right). (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', help='Path to write output tsv file')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
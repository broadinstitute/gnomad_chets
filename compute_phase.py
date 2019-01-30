from gnomad_hail import *
from resources import *
from phasing import *

def main(args):

    hl.init(log="/tmp/phasing.hail.log")

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]

    ht = hl.read_table(vp_count_ht_path(*path_args))

    if not args.no_em:

        def get_em_expr(gt_counts):
            hap_counts = hl.experimental.haplotype_freq_em(gt_counts)
            return hl.bind(
                lambda x:  hl.struct(
                    hap_counts=x,
                    p_chet=(x[1] * x[2]) / (x[0] * x[3] + x[1] * x[2])
                ),
                hap_counts
            )

        ht = ht.annotate(
            em=hl.struct(
                raw=get_em_expr(ht.gt_counts.raw),
                adj=get_em_expr(ht.gt_counts.adj),
            ),
            em_plus_one=hl.struct(
                raw=get_em_expr(ht.gt_counts.raw + [0,0,0,0,1,0,0,0,0]),
                adj=get_em_expr(ht.gt_counts.adj + [0,0,0,0,1,0,0,0,0]),
            )
        )

    if not args.no_lr:
        def get_lr_annotation(gt_counts):
            same_hap_likelihood = same_hap_likelihood_expr(gt_counts)
            chet_likelihood = chet_likelihood_expr(gt_counts)
            return hl.bind(
                lambda x, y: hl.struct(
                    same_hap_like=x,
                    chet_like=y,
                    lr_chet=y - x
                ),
                same_hap_likelihood,
                chet_likelihood
            )

        ht = ht.annotate(
            likelihood_model=hl.struct(
                raw=get_lr_annotation(ht.gt_counts.raw),
                adj=get_lr_annotation(ht.gt_counts.adj)
            )
        )

    if not args.no_shr:
        def get_single_het_expr(gt_counts):
            return hl.bind(
                lambda x:
                hl.cond(x[1] > x[3],
                        (x[1] + x[2]) / hl.sum(hl.range(1, 9).filter(lambda x: x % 3 > 0).map(lambda y: x[y])),
                        (x[3] + x[6]) / hl.sum(x[3:])
                        ),
                gt_counts
            )

        ht = ht.annotate(
            singlet_het_ratio=hl.struct(
                raw=get_single_het_expr(ht.gt_counts.raw),
                adj=get_single_het_expr(ht.gt_counts.adj)
            )
        )

    ht.write(phased_vp_count_ht_path(*path_args), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--pbt', help='If --pbt is specified, then only sites present in PBT samples are used and counts exclude PBT samples.',
                        action='store_true')
    parser.add_argument('--least_consequence', help=f'Least consequence for the input (just to get the path right). (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--no_em', help=f'Do not compute EM phase', action='store_true')
    parser.add_argument('--no_lr', help=f'Do not compute likelihood-ratio phase', action='store_true')
    parser.add_argument('--no_shr', help=f'Do not compute single het ratio phase', action='store_true')
    parser.add_argument('--max_freq', help=f'Maximum global adj AF for the input (just to get the path right). (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chrom', help='Only run on given chromosome')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
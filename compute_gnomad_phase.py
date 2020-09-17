from resources import *
from phasing import *
import argparse


def main(args):
    hl.init(log="/tmp/phasing.hail.log")

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]

    ht = hl.read_table(vp_count_ht_path(*path_args))
    expr_fun = []

    if not args.no_em:
        expr_fun.append(get_em_expressions)

    if not args.no_lr:
        expr_fun.append(get_lr_expressions)

    if not args.no_shr:
        expr_fun.append(get_single_het_expressions)

    if not expr_fun:
        raise (Exception("No expressions to annotate"))

    # Support for both exploded or dict versions of gt_counts
    # dict
    if isinstance(ht.gt_counts, hl.expr.DictExpression):
        ht = ht.select(
            phase_info=ht.gt_counts.map_values(
                lambda pop_count: hl.bind(
                    lambda x: hl.struct(
                        gt_counts=x,
                        **{
                            k: v for f in expr_fun for k, v in f(x).items()
                        }
                    ),
                    hl.struct(
                        raw=pop_count.raw.map(lambda y: hl.int32(y)),
                        adj=pop_count.adj.map(lambda z: hl.int32(z))
                    )
                )
            )
        )
    # exploded
    else:
        ht = ht.annotate(
            **{
                k: v for f in expr_fun for k, v in f(ht.gt_counts).items()
            }
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
    parser.add_argument('--chrom', help='Only run on given chromosome')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    main(args)

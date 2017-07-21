import argparse
from utils import *
from hail import *



def main(args):

    hc = HailContext()

    if args.exomes:
        vds = hc.read(full_exome_hardcalls_adj_vds_path)
        vds = vds.filter_intervals(Interval.parse("20"))
        vds = vds.filter_samples_expr('sa.meta.drop_status == "keep"')

        subsets = { "n{}".format(size): size for size in args.subset }

        logger.info("Creating the following subsets: {}".format(",".join(["{}:{}".format(name, size) for name,size in subsets.iteritems()])))

        vds = create_sample_subsets(vds, subsets)

        vds = vds.annotate_variants_expr([
            "va.calldata.{0} = gs.filter(g => sa.{0}).callStats(g => v)".format(name) for name in subsets.keys()
        ] + [
            "va.calldata.n{} = gs.callStats(g => v)".format(vds.num_samples)
        ])

        nsamples = vds.num_samples

        vds = vds.drop_samples()

        vds = vds.filter_alleles("va.calldata.n{}.AC[aIndex] == 0".format(nsamples), keep=False)

        vds.write(args.output, overwrite=args.overwrite)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input VDS is exomes',
                        action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes',
                        action='store_true')
    parser.add_argument('--subset', help='Subset size(s)', type=int, nargs="+", required=True)
    parser.add_argument('--output', '-o', help='Output VDS', required=True)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.genomes:
        sys.exit("Not implemented yet -- need ADJ hardcalls ;)") #TODO: fix that

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)



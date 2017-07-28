import argparse
from utils import *
from hail import *


def main(args):

    hc = HailContext()

    if args.exomes:
        vds_path = full_exome_vds_path if args.exomes else full_genome_vds_path
        vds = hc.read(vds_path)

        vds = vds.filter_samples_expr('sa.meta.drop_status == "keep"')

        subsets = { "n{}".format(size): size for size in args.subset }

        logger.info("Creating the following subsets: {}".format(",".join(["{}:{}".format(name, size) for name,size in subsets.iteritems()])))

        vds = vds.annotate_variants_expr("va.calldata.full = gs.callStats(g => v)")
        vds = vds.filter_alleles("va.calldata.full.AC[aIndex] == 0", keep=False)

        vds = create_sample_subsets(vds, subsets)

        vds = vds.annotate_variants_expr([
            "va.calldata.raw.{0} = gs.filter(g => sa.{0}).callStats(g => v)".format(name) for name in subsets.keys()
        ] + [
            "va.calldata.adj.{0} = gs.filter(g => sa.{0} && {1}).callStats(g => v)".format(name, ADJ_CRITERIA) for name in subsets.keys()
        ] + [
            "va.calldata.raw.n{0} = gs.callStats(g => v)".format(vds.num_samples),
            "va.calldata.adj.n{0} = gs.filter(g => {1}).callStats(g => v)".format(vds.num_samples, ADJ_CRITERIA)
        ])

        vds = vds.drop_samples()

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

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


# Note 2017-07-27, how I actually got this to work given memory issues:
# python ~/gnomad_qc/hail/pyhail.py --cluster exomes --script create_subsample_ACs.py --exomes --subset 55000 60000 65000 70000 75000 80000 85000 90000 95000 100000 105000 110000 115000 120000 --output gs://gnomad-exomes/subsets/random_subsamples/gnomad.exomes.subsamples.sites_55k_123k.vds --overwrite --spark_conf spark.memory.fraction=0.33,spark.executor.memory=15g,spark.yarn.executor.memoryOverhead=25g --jar gs://hail-common/hail-hail-is-master-all-spark2.0.2-ff26e57.jar --pyhail gs://hail-common/pyhail-hail-is-master-ff26e57.zip
# python ~/gnomad_qc/hail/pyhail.py --cluster exomes2 --script create_subsample_ACs.py --exomes --subset 10 20 50 100 200 500 1000 2000 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000 --output gs://gnomad-exomes/subsets/random_subsamples/gnomad.exomes.subsamples.sites_10_50k.vds --overwrite --spark_conf spark.memory.fraction=0.33,spark.executor.memory=15g,spark.yarn.executor.memoryOverhead=25g --jar gs://hail-common/hail-hail-is-master-all-spark2.0.2-ff26e57.jar --pyhail gs://hail-common/pyhail-hail-is-master-ff26e57.zip
# Then, later:
# vds1 = hc.read('gs://gnomad-exomes/subsets/random_subsamples/gnomad.exomes.subsamples.sites_10_50k.vds')
# vds2 = hc.read('gs://gnomad-exomes/subsets/random_subsamples/gnomad.exomes.subsamples.sites_55k_123k.vds')
# vds1 = vds1.annotate_variants_expr(['va.calldata.raw.n{0} = select(va.calldata.raw.n{0}, AC, AN), va.calldata.adj.n{0} = select(va.calldata.adj.n{0}, AC, AN)'.format(x) for x in ['10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '15000', '20000', '25000', '30000', '35000', '40000', '45000', '50000']])
# vds2 = vds2.annotate_variants_expr(['va.calldata.raw.n{0} = select(va.calldata.raw.n{0}, AC, AN), va.calldata.adj.n{0} = select(va.calldata.adj.n{0}, AC, AN)'.format(x) for x in ['55000', '60000', '65000', '70000', '75000', '80000', '85000', '90000', '95000', '100000', '105000', '110000', '115000', '120000', '123136']])
# vds = vds1.annotate_variants_vds(vds2, 'va.calldata.raw = merge(drop(va.calldata.raw, n123136), vds.calldata.raw), va.calldata.adj = merge(drop(va.calldata.adj, n123136), vds.calldata.adj)')
# vds.write('gs://gnomad-exomes/subsets/random_subsamples/gnomad.exomes.subsamples.sites.vds', overwrite=True)

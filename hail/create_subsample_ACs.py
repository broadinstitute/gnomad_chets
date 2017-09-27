import argparse
from utils import *
from hail import *
from collections import Counter


def create_sample_subsets(vds, subsets, pop_location, sa_root="sa.downsampled_subset", global_root="global.downsampled_subset", seed=42):
    """
    Creates sample subsets annotations by taking random subsets of the desired size(s).

    :param VariantDataset vds: Input VDS
    :param dict of str:int subsets: A dict containing the desired subset names and sizes
    :param str sa_root: Root where sample annotations are written. Should start with `sa`
    :param str global_root: Root where global annotations are written. Should start with `global`
    :return: Annotated VDS
    :rtype: VariantDataset
    """
    random.seed(seed)
    sample_pops = vds.query_samples('samples.map(s => [s, {}]).collect()'.format(pop_location))
    pop_output = {}
    expr = []
    for name, n_samples in subsets.iteritems():
        random.shuffle(sample_pops)
        samples_this_subset = set(sample_pops[:n_samples])
        vds = vds.annotate_global("{}.{}".format(global_root, name), [x[0] for x in samples_this_subset], TSet(TString()))
        expr.append('{0}.{1} = {2}.{1}.contains(s)'.format(sa_root, name, global_root))

        pop_output[name] = Counter([x[1] for x in samples_this_subset])
    vds = vds.annotate_samples_expr(expr)

    return vds, pop_output


def main(args):

    hc = HailContext()

    vds_path = full_exome_vds_path if args.exomes else full_genome_vds_path
    pops = EXOME_POPS if args.exomes else GENOME_POPS
    pop_location = 'sa.meta.final_pop' if args.genomes else 'sa.meta.population'

    vds = hc.read(vds_path)

    keep_text = 'sa.meta.drop_status == "keep"' if args.exomes else 'sa.meta.keep'
    vds = vds.filter_samples_expr(keep_text)

    if args.pcr_free_only:
        if args.exomes:
            logger.warning('PCR-free not available for exomes')
        else:
            vds = vds.filter_samples_expr('sa.meta.pcr_free')

    subsets = {"n{}".format(size): size for size in args.subset}

    logger.info("Creating the following subsets: {}".format(",".join(["{}:{}".format(name, size) for name, size in subsets.iteritems()])))

    vds = vds.annotate_variants_expr("va.calldata.full = gs.callStats(g => v)")
    vds = vds.filter_alleles("va.calldata.full.AC[aIndex] == 0", keep=False)

    vds, pop_counts = create_sample_subsets(vds, subsets, pop_location)

    vds = vds.annotate_variants_expr([
        "va.calldata.raw.{0} = gs.filter(g => sa.{0}).callStats(g => v)".format(name) for name in subsets.keys()
    ] + [
        "va.calldata.adj.{0} = gs.filter(g => sa.{0} && {1}).callStats(g => v)".format(name, ADJ_CRITERIA) for name in subsets.keys()
    ] + [
        "va.calldata.raw.n{0} = gs.callStats(g => v)".format(vds.num_samples),
        "va.calldata.adj.n{0} = gs.filter(g => {1}).callStats(g => v)".format(vds.num_samples, ADJ_CRITERIA)
    ])

    if args.populations:
        total_pops = vds.query_samples('samples.map(s => {}).counter()'.format(pop_location))
        vds = vds.annotate_variants_expr([
            "va.calldata.pop_raw.{0}.n{1} = gs.filter(g => {2} == '{0}' && sa.{3}).callStats(g => v)".format(pop, pop_counts[name][pop.lower()], pop_location, name) for name in subsets.keys() for pop in pops
        ] + [
            "va.calldata.pop_adj.{0}.n{1} = gs.filter(g => {2} == '{0}' && sa.{3} && {4}).callStats(g => v)".format(pop, pop_counts[name][pop.lower()], pop_location, name, ADJ_CRITERIA) for name in subsets.keys() for pop in pops
        ] + [
            "va.calldata.pop_raw.{0}.n{1} = gs.filter(g => {2} == '{0}').callStats(g => v)".format(pop, total_pops[pop.lower()], pop_location) for pop in pops
        ] + [
            "va.calldata.pop_adj.{0}.n{1} = gs.filter(g => {2} == '{0}' && {3}).callStats(g => v)".format(pop, total_pops[pop.lower()], pop_location, ADJ_CRITERIA) for pop in pops
        ])

    vds = vds.drop_samples()

    vds.write(args.output, overwrite=args.overwrite)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input VDS is exomes', action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes', action='store_true')
    parser.add_argument('--populations', help='Also calculate population-specific metrics', action='store_true')
    parser.add_argument('--pcr_free_only', help='Use only PCR-free genomes (only applicable for genomes)', action='store_true')
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

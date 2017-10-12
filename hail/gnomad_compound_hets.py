import argparse
from compound_hets_utils import *
from hail import *

def get_vds(hc, args):
    data_type = "exomes" if args.exomes else "genomes"
    vds = get_gnomad_data(hc,
                          data_type,
                          hardcalls="adj" if args.adj else "raw",
                          split=True,
                          release_samples=True)

    if args.chrom20:
        vds = vds.filter_intervals(Interval.parse("20:1-10000000"))

    vds = filter_low_conf_regions(vds,
                                  high_conf_regions=[exomes_high_conf_regions_intervals_path] if "exomes" else None)
    vds = vds.annotate_variants_vds(get_gnomad_public_data(hc, data_type, split=True), root='va.release')
    vds = annotate_gene_impact(vds, vep_root='va.vep')

    # Add methylated CpG annotation
    vds = vds.annotate_variants_vds(hc.read(context_vds_path), expr='va.methylated_cpg = vds.methylation.value >= 0.25,'
                                                                    'va.coverage = vds.coverage')

    vds = vds.filter_variants_expr(
        'isDefined(va.gene) && isDefined(va.impact) && va.release.info.AF_POPMAX <= {} && va.release.filters.isEmpty()'.format(
            args.max_af))
    vds = vds.annotate_variants_expr(
        'va.release.info = select(va.release.info, AC, AN, AF, POPMAX, AC_POPMAX, AN_POPMAX)')
    vds = vds.annotate_variants_expr(['va.release = select(va.release, info, filters)',
                                      'va.vep = drop(va.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences)'])
    vds = vds.annotate_variants_expr(['va = drop(va, stats, pass, calldata, tdt, rsid)'])
    vds = vds.annotate_samples_expr('sa = {pop: sa.meta.population}')
    pprint(vds.variant_schema)

    return vds.persist()


def main(args):

    hc = HailContext(log='/gnomad_compound_hets.log')
    vds = None

    if args.write_pairs_kt:

        vds = get_vds(hc, args)

        n_partitions = 50 if args.chrom20 else 4000

        kt = vds.phase_em(['va.gene'], n_partitions, sa_keys='sa.pop', by_sample=args.by_sample)

        kt.write(args.output + ".variant_pairs.kt", overwrite = args.overwrite)

    if args.write_single_kt:
        if vds is None:
            vds = get_vds(hc, args)

        kt = vds.genotypes_table()
        kt = kt.filter('g.isHomVar')
        kt.write(args.output + ".single_variants.kt", overwrite = args.overwrite)


    if args.export_pairs:
        kt = hc.read_table(args.output + ".variant_pairs.kt")
        kt, count_cols = flatten_counts(kt,gt_anns=['g1', 'g2'])

        variant_cols = {'chrom': 'v{}.contig',
                        'pos': 'v{}.start',
                        'ref': 'v{}.ref',
                        'alt': 'v{}.alt',
                        'cpg': 'va{}.methylated_cpg',
                        'pass': 'va{}.release.filters.isEmpty',
                        'impact': 'va{}.impact',
                        'ac_raw': 'va{0}.info.AC[va{0}.aIndex -1]',
                        'an_raw': 'va{}.info.AN',
                        'ac_release': 'va{}.release.info.AC',
                        'an_release': 'va{}.release.info.AN',
                        'exome_coverage': 'va{}.coverage.exome',
                        'genome_coverage': 'va{}.coverage.genome',
                        'wasSplit': 'va{}.wasSplit',
                        'alleleType': 'va{}.alleleType'}

        (
            # kt.filter('prob_same_haplotype < 0.5 || single_het_ratio > 0.5')
            kt.filter('va1.impact == "high" && va2.impact == "high"')
                .annotate(['gene = `va.gene`',
                           'pop = `sa.pop`',
                           'sample = s',
                           'distance = (v1.start - v2.start).abs()',
                           'single_het_ratio = if(AaBB > AABb) (AABb + AAbb) / (AABb + AAbb + AaBb + aaBb + Aabb + aabb) else (AaBB + aaBB) / (AaBB + aaBB + AaBb + aaBb + Aabb + aabb)'
                           ] +
                          ['{0}{1} = {2}'.format(name, n, expr.format(n)) for n in ["1", "2"] for name, expr in
                           variant_cols.iteritems()]
                          )
                .select(
                ['gene', 'sample', 'pop', 'prob_same_haplotype', 'distance', 'single_het_ratio'] +
                count_cols +
                [col + n for n in ["1", "2"] for col in variant_cols.keys()])
                .export(args.output + ".variant_pairs.txt.bgz")
        )

    if args.export_single:
        kt = hc.read_table(args.output + ".single_variants.kt")
        (
            kt.annotate(['pop = sa.pop', 'gene = va.gene', 'impact  = va.impact',
                         'sample = s', 'alleleType = va.alleleType',
                         'ac_raw = va.info.AC[va.aIndex -1]', 'ac = va.release.info.AC',
                         'pass = va.release.filters.isEmpty', 'cpg = va.methylated_cpg',
                         'wasSplit = va.wasSplit',
                         'ref = v.ref', 'alt = v.alt',
                         'chrom = v.contig',
                         'pos = v.start',
                         'gq = g.gq',
                         'dp = g.dp',
                         'ad0 = g.ad[0]',
                         'ad1 = g.ad[1]'
                         ])
                .select(['gene', 'chrom', 'pos', 'ref', 'alt', 'cpg', 'wasSplit',
                         'impact', 'alleleType', 'ac_raw', 'ac', 'pass', 'sample', 'pop'])
                .export(
                args.output + ".single_variants.txt.bgz")
        )

    if args.slack_channel:
        send_message(args.slack_channel, 'gnomAD compound hets {} is done processing!'.format(args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Use exomes', required=False, action='store_true')
    parser.add_argument('--genomes', help='Use genomes', required=False, action='store_true')
    parser.add_argument('--write_pairs_kt', help='Writes variant pairs KeyTable.', required=False, action='store_true')
    parser.add_argument('--write_single_kt', help='Writes single variants KeyTable.', required=False, action='store_true')
    parser.add_argument('--export_pairs', help='Writes variant pairs tsv.', required=False, action='store_true')
    parser.add_argument('--export_single', help='Writes single variants tsv.', required=False, action='store_true')
    parser.add_argument('--max_af', help='Maximum AF for a site to be retained (default 0.01).', required=False,
                        type=float, default=0.01)
    parser.add_argument('--by_sample', help='Compute compound het per sample rather than per site', required=False, action='store_true')
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrites existing results.', required=False,
                        action='store_true')
    parser.add_argument('--adj', help='Use ADJ genotypes only', required=False, action='store_true')
    parser.add_argument('--chrom20', help='Process chrom 20 only', required=False, action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.exomes and args.genomes:
        sys.exit("Only one of --exomes, --genomes can be specified.")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
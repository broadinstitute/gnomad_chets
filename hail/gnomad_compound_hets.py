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

    vds = vds.annotate_variants_vds(get_gnomad_public_data(hc, data_type="exomes", split=True),
                                    root='va')
    vds = vds.filter_variants_expr('va.filters.isEmpty')
    #vds = annotate_gene_impact(vds, vep_root='va.vep')

    # Add methylated CpG annotation
    vds = vds.annotate_variants_vds(hc.read(context_vds_path), expr='va.methylated_cpg = vds.methylation.value >= 0.25,'
                                                                    'va.coverage = vds.coverage')
    vds = filter_vep_to_canonical_transcripts(vds)
    vds = process_consequences(vds)
    vds = vds.annotate_variants_expr(
        'va.vep.transcript_consequences = va.vep.transcript_consequences.filter(csq => {} csq.lof == "HC")'.format(
            '' if args.lof_only else 'csq.most_severe_consequence == "missense_variant" ||'
        ))

    vds = vds.filter_variants_expr('va.info.AF <= {} && !va.vep.transcript_consequences.isEmpty()'.format(args.max_af))

    # vds = vds.filter_variants_expr(
    #     'isDefined(va.gene) && isDefined(va.impact) && va.release.info.AF <= {} && va.release.filters.isEmpty()'.format(
    #         args.max_af))
    vds = vds.annotate_variants_expr([
        'va.info = select(va.info, AC, AN, AF, POPMAX, AC_POPMAX, AN_POPMAX, AC_raw, AN_raw)',
        'va.vep = drop(va.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences)',
        'va.gene = va.vep.transcript_consequences[0].gene_symbol'
    ])
    vds = vds.annotate_samples_expr('sa = {pop: sa.meta.population}')
    pprint(vds.variant_schema)

    return vds.persist()


def main(args):

    hc = HailContext(log='/gnomad_compound_hets.log')
    vds = None

    if args.write_pairs_kt:

        vds = get_vds(hc, args)

        n_partitions = 50 if args.chrom20 else 1000 if args.lof_only else 4000

        kt = vds.phase_em(['va.gene'], n_partitions, sa_keys='sa.pop', by_sample=args.by_sample)

        kt.write(args.output + ".variant_pairs.kt", overwrite = args.overwrite)

    if args.write_single_kt:
        if vds is None:
            vds = get_vds(hc, args)

        kt = vds.genotypes_table()
        kt = kt.filter('g.isHomVar')
        kt.write(args.output + ".single_variants.kt", overwrite = args.overwrite)

    variant_cols = {'chrom': 'v{}.contig',
                    'pos': 'v{}.start',
                    'ref': 'v{}.ref',
                    'alt': 'v{}.alt',
                    'cpg': 'va{}.methylated_cpg',
                    # 'pass': 'va{}.filters.isEmpty',
                    # 'impact': 'va{}.impact',
                    # 'ac_raw': 'va{0}.info.AC_raw[va{0}.aIndex -1]',
                    # 'an_raw': 'va{}.info.AN_raw',
                    'ac_release': 'va{}.info.AC',
                    'an_release': 'va{}.info.AN',
                    # 'exome_coverage': 'va{}.coverage.exome',
                    # 'genome_coverage': 'va{}.coverage.genome',
                    'wasSplit': 'va{}.wasSplit',
                    'snv': 'v{}.altAllele.isSNP()',
                    'indel': 'v{}.altAllele.isIndel()'}

    if args.export_pairs:
        kt = hc.read_table(args.output + ".variant_pairs.kt")
        #Order variant by position
        kt = conditional_column_swap(kt,
                                                 swap_expr='v1.start > v2.start',
                                                 columns=[
                                                     ('{}1'.format(col),
                                                      '{}2'.format(col)) for col in [
                                                         'v',
                                                         'va',
                                                     ]
                                                 ],
                                     gt_counts_col="genotype_counts",
                                     hc_counts_col="haplotype_counts")

        kt, hc_count_cols = flatten_haplotype_counts(kt)

        (
            # kt.filter('prob_same_haplotype < 0.5 || single_het_ratio > 0.5')
            kt
             #   .filter('va1.impact == "high" && va2.impact == "high"')
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
                hc_count_cols +
                [col + n for n in ["1", "2"] for col in variant_cols.keys()])
                .export(args.output + ".variant_pairs.txt.bgz")
        )

    if args.export_single:

        kt = hc.read_table(args.output + ".single_variants.kt")
        (
            kt.annotate(['gene = va.gene',
                         'pop = sa.pop',
                         'sample = s'] +
                        ['{0} = {1}'.format(name, expr.format("")) for name, expr in variant_cols.iteritems()]
                        )
                .select(['gene', 'sample', 'pop'] +
                        variant_cols.keys())
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
    parser.add_argument('--lof_only', help='Only output LoFs', required=False,
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.exomes and args.genomes:
        sys.exit("Only one of --exomes, --genomes can be specified.")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
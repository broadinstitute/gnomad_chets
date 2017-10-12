__author__ = 'konradk'
import argparse
from utils import *


def main():
    hc = HailContext()

    vds = get_gnomad_data(hc, data_type="exomes", hardcalls="adj", split=True, release_samples=True)
    vds = vds.filter_intervals(Interval.parse('1-22'))
    vds = filter_low_conf_regions(vds,
                                  high_conf_regions=exomes_high_conf_regions_intervals_path)
    vds = vds.annotate_variants_vds(get_gnomad_public_data(hc, data_type="exomes", split=True),
                                    root='va')

    vds = vds.filter_variants_expr('va.filters.isEmpty')

    # Filter consequences
    vds = filter_vep_to_canonical_transcripts(vds)
    vds = process_consequences(vds)
    vds = vds.annotate_variants_expr(
        'va.vep.transcript_consequences = va.vep.transcript_consequences.filter(csq => {} csq.lof == "HC")'.format(
            '' if args.lof_only else 'csq.most_severe_consequence == "missense_variant" ||'
        ))

    vds = vds.filter_variants_expr('va.info.AF <= {} && !va.vep.transcript_consequences.isEmpty()'.format(args.max_af))

    # Create keytable, filtering early
    annotations = {'ac': 'va.info.AC',
                   'an': 'va.info.AN',
                   'ac_popmax': 'va.info.AC_POPMAX',
                   'an_popmax': 'va.info.AN_POPMAX'}

    exome_kt = (vds.genotypes_table().filter('g.isHet()')
                .annotate(['csqs = va.vep.transcript_consequences'] +
                            ['{} = [}'.format(k,v) for k,v in annotations.iteritems()])
                .select(['v', 'g', 's','csqs'] + annotations.keys())
                .explode('csqs'))

    kt = exome_kt.aggregate_by_key(['sample = s', 'gene = csqs.gene_symbol'],
                                   ['variants = v.collect()'] +
                                   ['{0} = {0}.collect()'.format(ann) for ann in annotations.keys()])

    kt = kt.filter('variants.length > 1')
    kt = kt.persist()
    max_n_variants = kt.query('variants.map(x => x.kength).max()')

    kt_vp = kt.annotate(['v1 = variants[0]',
                         'n_variants = variants.length',
                          'variants = variants[1:]'])
    kt_vp = kt_vp.explode('variants')
    kt_vp = kt_vp.annotate(['v2 = variants'] +
                           ['{0}{1} = {0}[{2}]'.format(ann, n, n-1) for n in [1,2] for ann in annotations.keys()])

    for i in range(1,max_n_variants-1):
        temp_kt = (
            kt.filter('all_variants.length > i')
            .annotate(['v1 = variants[{}]'.format(i),
                          'variants = variants[{}:]'.format(i + 1)])
        )
        temp_kt = temp_kt.explode('variants')
        temp_kt = temp_kt.annotate(['v2 = variants'] +
                               ['{0}{1} = {0}[{2}]'.format(ann, n,  i + n -1) for n in [1, 2] for ann in
                                annotations.keys()])
        kt_vp = kt_vp.union(temp_kt)

    kt_vp = kt_vp.drop(['variants'])
    kt_vp.write(args.output + ".kt")
    kt_vp = hc.read_table(args.output + ".kt")
    kt_vp.export(args.output + ".kt.tsv.bgz")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max_af', help='Maximum AF for a site to be retained (default 0.01).', required=False,
                        type=float, default=0.01)
    parser.add_argument('--lof_only', help='Only output LoFs', required=False,
                        action='store_true')
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
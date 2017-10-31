from utils import *
import hail
import argparse
import utils


def get_nNonRefSamples(vds, vep_vds, af_bins = [0.01, 0.05]):

    def get_af_bin_expr(expr, af_bins):
        if af_bins:
            if not expr:
                expr = "%s"
            if len(af_bins) >1:
                expr = get_af_bin_expr(expr % ('if(va.AF <= %s) "%s" else %%s' % (af_bins[0], af_bins[0])), af_bins[1:])
            else:
                expr = expr % ('if(va.AF <= %s) "%s" else ">%s"' % (af_bins[0], af_bins[0], af_bins[0]))

        return expr

    vds = vds.annotate_variants_expr(['va.nNonRef = gs.filter(g => g.isCalledNonRef).map(g => g.oneHotAlleles(v).map(x => min([x, 1]))).sum()',
                                      'va.AC  = gs.filter(g => g.isCalled).count() * 2'])
    vds = vds.drop_samples().split_multi().annotate_variants_expr('va.nNonRef = va.nNonRef[va.aIndex], va.AF = va.nNonRef[va.aIndex] / va.AC')
    vds = vds.persist()
    vds = vds.annotate_variants_vds(vep_vds, expr = 'va.vep = vds.vep')
    kt = vds.variants_keytable()
    af_bin_expr = get_af_bin_expr("", af_bins)
    kt_annotations = ['most_severe_consequence = va.vep.most_severe_consequence']
    agg_keys = ['most_severe_consequence = most_severe_consequence']
    if af_bin_expr:
        kt_annotations.append('ac_bin = ' + af_bin_expr)
        agg_keys.append('ac_bin = ac_bin')

    res = kt.annotate(kt_annotations).aggregate_by_key(agg_keys, 'nNonRefGTs = va.map(x => x.nNonRef.toLong).sum()')
    res.persist()
    print(res.collect())
    res.to_dataframe().orderBy(['nNonRefGTs']).show(truncate=False, n=5000)

        #print(vds.filter_variants_expr(filter_criteria).query_genotypes('gs.filter(g => g.isCalledNonRef).count()'))


    #

def main(args):

    hc = hail.HailContext(log='/test.log')
    logger.setLevel(logging.DEBUG)

    if args.get_n_nonref_samples:
        get_nNonRefSamples(hc.read(full_genome_vds_path), hc.read(full_genomes_vep_split_vds_path))
        get_nNonRefSamples(hc.read(full_exome_vds_path), hc.read(full_genomes_vep_split_vds_path))

    elif args.export_variants:
        if args.exomes:
            release = hc.read(final_exome_split_vds_path)
            vep = hc.read(full_exomes_vep_split_vds_path)
            output_prefix = args.output + 'exomes'
        else:
            release = hc.read(final_genome_split_vds_path)
            vep = hc.read(full_genomes_vep_split_vds_path)
            output_prefix = args.output + 'genomes'

        #Load release variant annotations to export
        annotations_to_ignore = ['DB', 'GQ_HIST_ALL', 'DP_HIST_ALL', 'AB_HIST_ALL', 'GQ_HIST_ALT', 'DP_HIST_ALT',
                                 'AB_HIST_ALT', 'GC.*', '.*_Male', '.*_Female', 'CSQ']

        annotations, a_annotations, g_annotations, dot_annotations = get_numbered_annotations(release, 'va.info')
        a_annotations = filter_annotations_regex(a_annotations, annotations_to_ignore)
        annotations = filter_annotations_regex(annotations, annotations_to_ignore)
        release = release.annotate_variants_expr([
            'va.RF=va.info.AS_FilterStatus.contains("RF")',
            'va.AC0=va.info.AS_FilterStatus.contains("AC0")',
            'va.segdup = va.filters.contains("SEGDUP")',
            'va.lcr = va.filters.contains("LCR")',
            'va.filtered = !va.info.AS_FilterStatus.isEmpty || va.filters.contains("LCR") || va.filters.contains("SEGDUP") || va.info.InbreedingCoeff < -0.3',
            'va.info.AS_FilterStatus = va.info.AS_FilterStatus.mkString(",")'
        ])

        #Use VEP-annotated full sites to get (1) VEP annotations and (2) Original annotations (useful for sites not in release)
        vep = vep.annotate_variants_expr(['va.info.AC_orig = va.info.AC[va.aIndex-1]',
                                          'va.info.AN_orig = va.info.AN',
                                          'va.info.AF_orig = va.info.AF[va.aIndex-1]',
                                          'va.info = drop(va.info, AC, AN, AF, DB, END, MLEAC, MLEAF, HaplotypeScore)'])
        orig_annotations, _, _, _ = get_numbered_annotations(vep, 'va.info')
        orig_annotations = [a.name for a in orig_annotations]
        vep = vep.annotate_variants_vds(release,
                                        expr = ",".join(['va.info.%s = vds.info.%s' % (a.name,a.name) for a in annotations + a_annotations if a.name not in orig_annotations]) +
                                        ", va.inRelease = isDefined(vds.info), va.RF = vds.RF, va.AC0 = vds.AC0, va.segdup = vds.segdup, va.lcr = vds.lcr, va.filtered = vds.filtered")

        annotations_to_export = ['chrom=v.contig',
             'pos=v.start',
             'ref=v.ref',
             'alt=v.alt',
             'was_split = va.wasSplit',
             'RF_filtered=va.RF',
             'AC0_filtered=va.AC0',
             'segdup_filtered=va.segdup',
             'lcr_filtered=va.lcr',
             'inbreedingCoeff_filtered=va.info.InbreedingCoeff < -0.3',
             'filtered = va.filtered',
             'most_severe_consequence = va.vep.most_severe_consequence',
             'in_release = va.inRelease'
             ]
        annotations_to_export.extend(['{0} = va.info.{0}'.format(a.name) for a in annotations + a_annotations ])
        annotations_to_export.extend(['{0} = va.info.{0}'.format(a) for a in orig_annotations])

        vep = vep.persist()
        #vep.export_variants(output_prefix + '.variants.tsv.gz', ','.join(annotations_to_export), parallel=True)

        #Export transcripts table
        kt = vep.filter_variants_expr('isDefined(va.vep.transcript_consequences) && !va.vep.transcript_consequences.isEmpty').variants_table()
        kt = kt.annotate(['chrom = v.contig',
                          'pos = v.start',
                          'ref = v.ref',
                          'alt = v.alt',
                          'transcript_consequences = va.vep.transcript_consequences',
                          'most_severe_consequence = va.vep.most_severe_consequence'])
        kt = kt.select(['chrom', 'pos', 'ref', 'alt', 'transcript_consequences', 'most_severe_consequence'])
        kt = kt.explode(['transcript_consequences'])
        kt = kt.flatten()
        kt = kt.rename({ x: x.replace(".","_") for x in kt.columns if x.startswith("transcript_consequences.")})
        kt = kt.explode(['transcript_consequences_consequence_terms'])
        kt = kt.annotate(['transcript_consequences_domains = transcript_consequences_domains.map(x => x.db + ": " + x.name).mkString(",")'])
        kt.export(output_prefix + '.transcript_consequences.tsv.gz', parallel=True)

    elif args.export_genotypes:
        if args.exomes:
            vds = add_exomes_sa(hc.read(full_exome_vds_path))
            #vds = exomes_sites_vcf.preprocess_vds(vds, hc.read(vqsr_vds_path), vds_pops=exomes_sites_vcf.pops)
            release = hc.read(final_exome_split_vds_path)
            vep = hc.read(full_exomes_vep_split_vds_path)
            output_prefix = args.output + 'exomes'
        else:
            vds = add_genomes_sa(hc.read(full_genome_vds_path))
            #vds = genomes_sites_vcf.preprocess_vds(vds, vds_pops=genomes_sites_vcf.pops)
            release = hc.read(final_genome_split_vds_path)
            vep = hc.read(full_genomes_vep_split_vds_path)
            output_prefix = args.output + 'genomes'

        vds = vds.filter_genotypes('g.isCalledNonRef()')

        high_impact = ",".join(["transcript_ablation",
                       "splice_acceptor_variant",
                       "splice_donor_variant",
                       "stop_gained" ,
                       "frameshift_variant",
                       "stop_lost",
                       "start_lost",
                       "transcript_amplification"])

        vds = vds.split_multi(propagate_gq=True)
        vds = vds.annotate_variants_vds(release, expr='va.info.AC = vds.info.AC,'
                                                     'va.info.AN = vds.info.AN')
        vds = vds.annotate_variants_vds(vep, expr = 'va.vep = vds.vep')
        vds = vds.filter_variants_expr('va.info.AN == 0 || va.info.AC / va.info.AN <= 0.01 || [%s].toSet.contains(va.vep.most_severe_consequence)' % high_impact)
        vds = vds.persist()

        cols = [
            'chrom=v.contig',
            'pos=v.start',
            'ref=v.ref',
            'alt=v.alt',
            'RF=va.RF',
            'sample=s',
            'gt_isHet=g.isHet',
            'gt_isHQ=%s' % ADJ_CRITERIA,
            'gt_fakeRef=g.fakeRef',
            'gt_dp=g.dp', #Need prefix because va.info.DP exists
            'gt_ad0=g.ad[0]',
            'gt_ad1=g.ad[1]',
            'gt_gq=g.gq',
            'gt_pl0=g.pl[0]',
            'gt_pl1=g.pl[1]',
            'gt_pl2=g.pl[2]'
        ]

        if args.genomes:
            vds = vds.annotate_global('global.non_coding_csq', set(CSQ_NON_CODING), TSet(TString()))
            coding = vds.filter_variants_expr('isDefined(va.vep.most_severe_consequence) && !global.non_coding_csq.contains(va.vep.most_severe_consequence)')
            coding.export_genotypes(output_prefix + '.coding.tsv.gz', ",".join(cols), parallel=True)
            non_coding = vds.filter_variants_expr('!isDefined(va.vep.most_severe_consequence) || global.non_coding_csq.contains(va.vep.most_severe_consequence)')
            non_coding.export_genotypes(output_prefix + '.non_coding.tsv.gz', ",".join(cols), parallel=True)
        else:
            vds.export_genotypes(output_prefix + '.tsv.gz', ",".join(cols), parallel=True)

    send_message(channel='@laurent', message='BigQuery export is done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Export exomes genotypes', action='store_true')
    parser.add_argument('--genomes', help='Export genomes genotypes', action='store_true')
    parser.add_argument('--export_genotypes', help='Export genotypes', action='store_true')
    parser.add_argument('--export_variants', help='Export variants', action='store_true')
    parser.add_argument('--get_n_nonref_samples', help='Counts the number of non-ref genotypes by AC and functional annotations', action='store_true')
    parser.add_argument('--output', '-o', help='Output directory')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if (args.export_genotypes or args.export_variants) and not args.output:
        sys.exit('Error: --output not specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


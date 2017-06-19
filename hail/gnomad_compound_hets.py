import argparse
from compound_hets_utils import *
from hail import *


def main(args):

    hc = HailContext(log='/gnomad_compound_hets.log')

    if not args.vds:
        #Select exomes / genomes files
        release = final_exome_split_vds if args.exomes else final_genome_split_vds
        vep = full_exomes_vep_split_vds if args.exomes else full_genomes_vep_split_vds
        high_conf_regions = exomes_high_conf_regions_path if args.exomes else None

        vds = add_exomes_sa(hc.read(full_exome_vds)) if args.exomes else add_genomes_sa(hc.read(full_genome_vds))
        vds = vds.filter_intervals([Interval.parse(x) for x in ['2:27424706-27424707','9:27429756-27429757']])
        logger.info(vds.count_variants())
        if args.genomes:
            vds = vds.annotate_samples_expr(['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop'])
        vds = filter_low_conf_regions(vds, high_conf_regions=high_conf_regions)
        logger.info(vds.count_variants())
        vds = vds.split_multi()
        logger.info(vds.count_variants())
        vds = vds.annotate_variants_vds(hc.read(release), root='va.release')
        vds = vds.annotate_variants_vds(hc.read(vep), expr='va.vep = vds.vep')
        vds = annotate_gene_impact(vds)
        vds = annotate_methylation(vds)
        #af_expr = 'va.calldata.all_samples_raw' if args.raw_AF else 'va.release.info.AF'
        vds = vds.filter_variants_expr('isDefined(va.gene) && isDefined(va.impact) && va.impact != "low" && va.release.info.AF_POPMAX <= {}'.format(args.max_af))
        logger.info(vds.count_variants())
        vds = filter_to_adj(vds)
        logger.info(vds.count_variants())
       # vds.hardcalls().write(args.output + '.vds', overwrite = args.overwrite)
       # vds = hc.read(args.output + '.vds')
        sys.exit("boom")
    else:
        vds = hc.read(args.vds)

    #TODO: Move up if regenerating data
    if args.exomes:
        vds = vds.filter_samples_expr('sa.meta.drop_status == "keep"')
    else:
        vds = vds.filter_samples_expr('sa.meta.keep')

    if args.chrom20:
        vds = vds.filter_intervals(Interval.parse("20:1-10000000"))

    if args.write_pairs_kt:
        n_partitions = 50 if args.chrom20 else 7000
        vds = vds.annotate_variants_expr('va.release.info = select(va.release.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)')
        vds = vds.annotate_variants_expr(['va.release = select(va.release, info, filters)',
                                         'va.vep = drop(va.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences) '])
        logger.info(str(vds.variant_schema))
        #TODO REMOVE AFTER TESTING!!!!
        n_partitions = 10
        #vds = vds.filter_samples_list(['10C113973', '10C113780', '10C113781'])
        vds = vds.filter_intervals([Interval.parse(x) for x in ['2:27424706-27424707', '2:27429756-27429757']])
        vds = vds.annotate_variants_expr('va = select(va, impact, gene)')
        vds = vds.annotate_samples_expr('sa = {pop: sa.meta.population}')
        kt = vds.phase_em(['va.gene'], n_partitions)
        kt = kt.persist()
        kt.to_dataframe().show()
        kt.export("gs://gnomad-lfran/tmp/ch_test.txt")
        sys.exit("boom")
        #TODO END
        kt.write(args.output + ".variant_pairs.kt", overwrite = args.overwrite)
        kt = hc.read_table(args.output + ".variant_pairs.kt")
        (
            kt.filter('prob_same_haplotype < 0.5')
                .annotate(['gene = `va.gene`', 'pop = sa.meta.population', 'sample = s',
                           'impact1 = va1.impact', 'impact2 = va2.impact',
                           'alleleType1 = va1.alleleType', 'alleleType2 = va2.alleleType',
                           'ac_raw1 = va1.info.AC[va1.aIndex -1]', 'ac_raw2 = va2.info.AC[va2.aIndex -1]',
                           'ac1 = va1.release.info.AC', 'ac2 = va2.release.info.AC',
                           'pass1 = va1.release.filters.isEmpty', 'pass2 = va2.release.filters.isEmpty',
                           'distance = (v1.start - v2.start).abs()',
                           'wasSplit1 = va1.wasSplit', 'wasSplit2 = va2.wasSplit',
                           'cpg1 = va1.methylated_cpg', 'cpg2 = va2.methylated_cpg',
                           'ref1 = v1.ref', 'alt1 = v1.alt',
                           'ref2 = v2.ref', 'alt2 = v2.alt',
                           'chrom1 = v1.contig', 'chrom2 = v2.contig',
                           'pos1 = v1.start', 'pos2 = v2.start'
                           ])
                .select(['gene', 'chrom1','pos1','ref1','alt1','cpg1', 'pass1', 'impact1', 'alleleType1', 'ac1', 'ac_raw1',
                         'chrom2','pos2','ref2','alt2','cpg2', 'pass2', 'impact2', 'alleleType2', 'ac2', 'ac_raw2',
                         'sample', 'pop', 'prob_same_haplotype','wasSplit1','wasSplit2',
                         'distance'])
                .export(args.output + ".variant_pairs.txt.bgz")
        )

    if args.write_single_kt:
        kt = vds.genotypes_table()
        kt = kt.filter('g.isHomVar')
        #kt = kt.aggregate_by_key(['gene = va.gene', 'impact = va.impact', 'alleleType = va.alleleType' , 'sample = s'],
        #                         ['nHomVar = g.count()', 'sa = sa.take(1)[0]'])
        kt.write(args.output + ".single_variants.kt", overwrite = args.overwrite)
        kt = hc.read_table(args.output + ".single_variants.kt")
        (
            kt.annotate(['pop = sa.meta.population','gene = va.gene','impact  = va.impact',
                         'sample = s', 'alleleType = va.alleleType',
                         'ac_raw = va.info.AC[va.aIndex -1]', 'ac = va.release.info.AC',
                         'pass = va.release.filters.isEmpty', 'cpg = va.methylated_cpg',
                         'wasSplit = va.wasSplit',
                         'cpg = va.methylated_cpg',
                         'ref = v.ref', 'alt = v.alt',
                         'chrom = v.contig',
                         'pos = v.start'
                         ])
                .select(['gene','chrom','pos','ref','alt','cpg','wasSplit',
                         'impact','alleleType','ac_raw','ac','pass','sample', 'pop'])
                .export(
        'gs://gnomad/compound_hets/gnomad_exomes.single_variants.txt.bgz')
        )

    #kt = kt.aggregate_by_key(['gene = va.gene', 'impact = va.impact', 'pop = sa.population'],
    #                         ['totHomVar = nHomVar.sum()', 'nHasHomVar = nHomVar.count()'])

    if args.slack_channel:
        send_message(args.slack_channel, 'gnomAD {} {} is done processing!'.format("exomes" if args.exomes else "genomes",
                                                                                   args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Use exomes and writes intermediate input', required=False, action='store_true')
    parser.add_argument('--genomes', help='Use genomes and writes intermediate input', required=False, action='store_true')
    parser.add_argument('--vds', help='Input VDS file for KT generation', required=False)
    parser.add_argument('--write_pairs_kt', help='Writes variant pairs KeyTable.', required=False, action='store_true')
    parser.add_argument('--write_single_kt', help='Writes single variants KeyTable.', required=False, action='store_true')
    parser.add_argument('--max_af', help='Maximum AF for a site to be retained (default 0.01).', required=False,
                        type=float, default=0.01)
    #parser.add_argument('--raw', help='Use raw rather than release AF for site selection', required=False, action='store_true')
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrites existing results.', required=False,
                        action='store_true')
    parser.add_argument('--chrom20', help='Process chrom 20 only', required=False, action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1: #TODO: Interactions between --exomes, --genomes and --vds are not clean
        sys.exit("One and only one of --exomes or --genomes is required.")

    main(args)
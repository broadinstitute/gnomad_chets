import argparse
from hail import *
from compound_hets_utils import *


def main(args):

    hc = HailContext()

    if args.trios:

        if args.filter_to_adj and (args.overwrite) and False: #TODO: rethink the flow of things!
            trios = add_exomes_sa(hc.read(full_exome_vds)) if args.exomes_trios else add_genomes_sa(hc.read(full_genome_hardcalls_split_vds))
            trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')
            ann = hc.read(args.trios)
            trios = trios.annotate_variants_vds(ann, 'va = vds')
            trios = trios.filter_variants_expr('isDefined(va.gene)')
            trios = filter_to_adj(trios)
            trios.write(args.output + '.adj.vds', overwrite=args.overwrite)
            trios = hc.read(args.output + '.adj.vds')
            logger.info("Count() in trios: {}".format(str(trios.count())))

        else:
            trios = hc.read(args.trios)
        ped = Pedigree.read(args.fam_file)
    else:
        if args.exomes_trios:
            trios = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds))
            trios = trios.annotate_variants_vds(hc.read(final_exome_split_vds), root = 'va.release')
            ped = Pedigree.read(exomes_fam)
        else:
            trios = add_genomes_sa(hc.read(full_genome_hardcalls_split_vds))
            trios = trios.annotate_variants_vds(hc.read(final_genome_split_vds), root = 'va.release')
            ped = Pedigree.read(genomes_fam)

        trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')

        ped = ped.filter_to(trios.sample_ids)
        logger.info("Found {0} trios in VDS.".format(len(ped.complete_trios())/3))

        if args.debug:
            n_variants =  trios.query_variants(['variants.map(x => va.release.filters.isEmpty()).counter()',
                                                'variants.map(v => va.calldata.all_samples_raw.AF <= {0}).counter()'.format(args.max_af)])

            logger.debug("Number of variants PASS: {0}, Number of variants with AF below {1}: {2}".format(str(n_variants[0]),args.max_af, str(n_variants[1])))

        trios = trios.filter_variants_expr(
            'va.release.filters.isEmpty() && va.calldata.all_samples_raw.AF <= {0} && gs.filter(g => g.isCalledNonRef).count() > 0'.format(args.max_af))

        if args.exomes_trios:
            trios = trios.annotate_variants_vds(hc.read(full_exomes_vep_split_vds), expr = 'va.vep = vds.vep')
        else:
            trios = trios.annotate_variants_vds(hc.read(full_genomes_vep_split_vds), expr = 'va.vep = vds.vep')

        trios = filter_low_conf_regions(trios, high_conf_regions=exomes_high_conf_regions_path)

        #Add methylated CpG annotation
        trios = annotate_methylation(trios)

        #Add VEP annotations
        trios = annotate_gene_impact(trios)

        trios.write(args.output + '.vds', overwrite=args.overwrite)
        trios = hc.read(args.output + '.vds')

    if args.chrom20:
        trios = trios.filter_intervals(Interval.parse("20:1-10000000"))

    n_partitions = 50 if args.chrom20 else 7000

    if args.exomes_trios: #TODO: Remove once data is re-generated
        trios = trios.annotate_variants_vds(hc.read(final_exome_split_vds), root='va.release')
    else:
        trios = trios.annotate_variants_vds(hc.read(final_genome_split_vds), root='va.release')

    #Select only fields of interest in va
    trios = trios.annotate_variants_expr([
    'va.release.info = select(va.release.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)'
    ])
    trios = trios.annotate_variants_expr(['va.release = select(va.release, info, filters)',
                                  'va.vep = drop(va.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences)',
                                          ])
    trios = trios.annotate_variants_expr(['va = drop(va, stats, calldata, tdt, methylation)',
                                          ])
    trios = trios.annotate_samples_expr(['sa.meta = select(sa.meta, population)'])

    trans_kt = trios.phase_by_transmission(ped, 'va.gene', n_partitions)
    trans_kt = trans_kt.key_by(['va.gene','v1','v2'])
    trans_kt.write("gs://gnomad-lfran/tmp/trans_kt.kt", overwrite=True) #TODO remove when working
    trans_kt = trans_kt.persist()

    trans_variants = trans_kt.query('v1.collect().extend(v2.collect()).toSet')
    logger.info("Found {} variants in pairs in trios.".format(len(trans_variants)))
    #trans_variants_kt = KeyTable.from_py(hc, [{'v': x} for x in trans_variants], TStruct(['v'], [TVariant()]), ['v'], 50)
    #trans_variants_kt = trans_variants_kt.persist()

    if args.reference:
        reference = hc.read(args.reference)
        reference = reference.filter_intervals(
            [Interval(start=Locus(v.contig, v.start), end=Locus(v.contig, v.start + 1)) for v in trans_variants])

    else:
        if args.filter_to_adj:
            reference = add_exomes_sa(hc.read(full_exome_vds))
        else:
            reference = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds))

        reference = reference.filter_samples_expr('!isDefined(sa.fam.famID) && sa.meta.drop_status == "keep"')
        # Changed to Locus to filter before split_multi
        reference = reference.filter_intervals(
            [Interval(start=Locus(v.contig, v.start), end=Locus(v.contig, v.start + 1)) for v in trans_variants])
        logger.info("Found ~{} of these variants in the reference pre-filtering".format(reference.count_variants()))
        # reference = reference.filter_variants_table(trans_variants_kt)

        if args.filter_to_adj:
            reference = filter_to_adj(reference)

        reference = reference.filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')
        reference = reference.split_multi()
        reference = reference.annotate_variants_expr(['va.AC = gs.callStats(g =>v).AC[1]'])
        reference = reference.annotate_variants_vds(trios, expr='va.gene = vds.gene')
        reference = reference.filter_variants_expr('isDefined(va.gene)')
        reference = reference.annotate_variants_expr('va = select(va, gene, AC)')
        reference = reference.annotate_samples_expr('sa = {pop: sa.meta.population}')
        if args.write_reference:
            reference.write(args.output + ".reference.vds", overwrite=True)
            reference = hc.read(args.output + ".reference.vds")
        else:
            reference = reference.persist()

    logger.info("Found ~{} of these variants in the reference post-filtering".format(reference.count_variants()))

    reference_kt = reference.phase_em(['va.gene'], n_partitions,
                                      sa_keys='sa.pop' if args.split_by_pop else None,
                                      variant_pairs= trans_kt)

    reference_kt.write("gs://gnomad-lfran/tmp/ref_kt.kt", overwrite=True) #TODO: Remove when things are working...
    #reference_kt.export("gs://gnomad-lfran/tmp/ref_kt.txt") #TODO: Remove when things are working...
    # pprint(reference_kt.schema())
    #reference_kt.to_dataframe().show()

    reference_keys = ['v1 = v1','v2 = v2']
    if args.split_by_pop:
        reference_keys.append('pop = sa.pop')

    reference_kt = reference_kt.aggregate_by_key(reference_keys,
                                                 ['haplotype_counts = haplotype_counts.take(1)[0]',
                                                  'genotype_counts = genotype_counts.take(1)[0]',
                                                  'prob_same_haplotype = prob_same_haplotype.take(1)[0]',
                                                  'ac1 = va1.take(1)[0].AC',
                                                  'ac2 = va2.take(1)[0].AC'])

    trans_kt = trans_kt.annotate(['gene = `va.gene`', 'pop = kidSA.meta.population'])
    if args.split_by_pop:
        trans_kt = trans_kt.key_by(['v1','v2','pop'])
    else:
        trans_kt = trans_kt.key_by(['v1', 'v2'])

    phase_trios_kt = trans_kt.join(reference_kt, how="left")
    phase_trios_kt = phase_trios_kt.persist()
    phase_trios_kt.write(args.output + '.kt',overwrite=args.overwrite)

    (
        phase_trios_kt.annotate(['fam = kidSA.fam.famID',
                                 'impact1 = va1.impact', 'impact2 = va2.impact',
                                 'alleleType1 = va1.alleleType', 'alleleType2 = va2.alleleType',
                                 'ac_raw1 = va1.info.AC[va1.aIndex -1]', 'ac_raw2 = va2.info.AC[va2.aIndex -1]',
                                 'pass1 = va1.release.filters.isEmpty', 'pass2 = va2.release.filters.isEmpty',
                                 'same_trio_haplotype = same_haplotype',
                                 'distance = (v1.start - v2.start).abs()',
                                 'wasSplit1 = va1.wasSplit', 'wasSplit2 = va2.wasSplit',
                                 'cpg1 = va1.methylated_cpg', 'cpg2 = va2.methylated_cpg',
                                 'ref1 = v1.ref', 'alt1 = v1.alt',
                                 'ref2 = v2.ref', 'alt2 = v2.alt',
                                 'chrom1 = v1.contig', 'chrom2 = v2.contig',
                                 'pos1 = v1.start','pos2 = v2.start'])
            .select(['gene', 'chrom1','pos1','ref1','alt1','cpg1', 'pass1', 'impact1', 'alleleType1', 'ac1', 'ac_raw1',
                     'chrom2','pos2','ref2','alt2','cpg2', 'pass2', 'impact2', 'alleleType2', 'ac2', 'ac_raw2',
                     'fam', 'pop', 'prob_same_haplotype', 'same_trio_haplotype', 'genotype_counts', 'haplotype_counts','distance',
                     'wasSplit1','wasSplit2'])
            .export(args.output + '.txt.bgz')
    )

    if args.slack_channel:
        send_message(args.slack_channel, 'Trio compound hets %s is done processing!' % args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--trios', help='Trio VDS', required=False)
    parser.add_argument('--fam_file', help='Fam file to use when --trios is given', required=False)
    parser.add_argument('--exomes_trios', help='Writes exomes trios VDS or use exome VDS fif filter_to_adj', required=False, action='store_true') #Todo, this is not clean
    parser.add_argument('--genomes_trios', help='Writes genomes trios VDS', required=False, action='store_true')
    parser.add_argument('--write_reference', help='Writes reference VDS', required=False, action='store_true')
    parser.add_argument('--reference', help='Specifies file to use as reference', required=False)
    parser.add_argument('--max_af', help='Maximum AF for a site to be retained (default 0.01).', required=False, type=float, default=0.01)
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--debug', help='Output debug statements', required=False, action='store_true')
    parser.add_argument('--chrom20', help='Process chrom 20 only', required=False, action='store_true')
    parser.add_argument('--filter_to_adj', help='Use Adj genotypes only', required=False, action='store_true')
    parser.add_argument('--split_by_pop', help='Splits data by population when computing EM', required=False, action='store_true')
    parser.add_argument('--overwrite', help='Overwrites existing results.', required=False,
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if int(args.exomes_trios) + int(args.genomes_trios) != 1:
        sys.exit("One and only one of --trios, --exomes_trios or --genomes_trios is required.")

    if args.trios and not args.fam_file:
        sys.exit("Must specify --fam_file when using --trios.")

    main(args)

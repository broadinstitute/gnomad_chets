import argparse
from hail import *
from compound_hets_utils import *


def create_trio_vds(hc, args):
    if args.exomes_trios:
        trios = hc.read(full_exome_vds_path) if args.filter_to_adj else hc.read(full_exome_hardcalls_split_vds_path)
        trios = add_exomes_sa(trios)
        ped = Pedigree.read(exomes_fam_path)
    else:
        trios = hc.read(full_genome_vds_path) if args.filter_to_adj else hc.read(full_genome_hardcalls_split_vds_path)
        trios = add_genomes_sa(trios)
        ped = Pedigree.read(genomes_fam_path)

    trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')
    trios = filter_low_conf_regions(trios, high_conf_regions=[exomes_high_conf_regions_intervals_path])
    ped = ped.filter_to(trios.sample_ids)
    logger.info("Found {0} trios in VDS.".format(len(ped.complete_trios()) / 3))

    if args.filter_to_adj:
        trios = trios.split_multi()
        trios = filter_to_adj(trios)

    if args.exomes_trios:
        trios = trios.annotate_variants_vds(hc.read(final_exome_split_vds_path), expr='va.filters = vds.filters,'
                                                                                      'va.release.info = select(vds.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)')
        trios = trios.annotate_variants_vds(hc.read(full_exomes_vep_split_vds_path), expr='va.vep = drop(vds.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences)')
    else:
        trios = trios.annotate_variants_vds(hc.read(final_genome_split_vds_path), expr='va.release.filters = vds.filters,'
                                                                                      'va.release.info = select(vds.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)')
        trios = trios.annotate_variants_vds(hc.read(full_genomes_vep_split_vds_path), expr='va.vep = drop(vds.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences)')

    if args.debug:
        trios = trios.persist()
        n_variants = trios.query_variants(['variants.map(x => va.release.filters.isEmpty()).counter()',
                                           'variants.map(v => va.calldata.all_samples_raw.AF <= {0}).counter()'.format(
                                               args.max_af)])

        logger.debug(
            "Number of variants PASS: {0}, Number of variants with AF below {1}: {2}".format(str(n_variants[0]),
                                                                                             args.max_af,
                                                                                             str(n_variants[1])))

    trios = trios.filter_variants_expr(
        'va.release.filters.isEmpty() && va.calldata.all_samples_raw.AF <= {0} && gs.filter(g => g.isCalledNonRef).count() > 0'.format(
            args.max_af))

    # Add methylated CpG annotation
    trios = trios.annotate_variants_vds(hc.read(context_vds_path), expr='va.methylated_cpg = vds.methylation.value >= 0.25,'
                                                                        'va.coverage = vds.coverage')

    # Add VEP annotations
    trios = annotate_gene_impact(trios)

    #Drop unused annotations
    trios = trios.annotate_variants_expr(['va = drop(va, stats, calldata, tdt, pass, info)',
                                          ])

    outpath = args.output + '.adj.vds' if args.filter_to_adj else args.output + '.vds'

    trios.write(outpath, overwrite=args.overwrite)
    return hc.read(outpath), ped


def create_reference(hc, args, trios_vds, trans_variants):
    if args.filter_to_adj:
        reference = add_exomes_sa(hc.read(full_exome_vds_path))
    else:
        reference = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds_path))

    reference = reference.filter_samples_expr('!isDefined(sa.fam.famID) && sa.meta.drop_status == "keep"')
    # Changed to Locus to filter before split_multi
    reference = reference.filter_variants_list(trans_variants)
    reference = reference.persist()

    logger.info("Found ~{} of these variants in the reference pre-filtering".format(reference.count_variants()))
    # reference = reference.filter_variants_table(trans_variants_kt)

    if args.filter_to_adj:
        reference = filter_to_adj(reference)

    reference = reference.annotate_variants_vds(trios_vds, expr='va.gene = vds.gene')
    reference = reference.filter_variants_expr('isDefined(va.gene) && gs.filter(g => g.isCalledNonRef).count() > 0')
    reference = reference.split_multi()
    reference = reference.annotate_variants_expr(['va.AC = gs.callStats(g =>v).AC[1]'])
    reference = reference.annotate_variants_expr('va = select(va, gene, AC)')
    reference = reference.annotate_samples_expr('sa = {pop: sa.meta.population}')
    reference.write(args.output + ".reference.vds", overwrite=True)

    return hc.read(args.output + ".reference.vds")


def main(args):

    hc = HailContext()

    if args.trios:
        trios = hc.read(args.trios)
        ped = Pedigree.read(args.fam_file)
    else:
        trios, ped = create_trio_vds(hc, args)

    if args.chrom20:
        trios = trios.filter_intervals(Interval.parse("20:1-10000000"))

    n_partitions = 50 if args.chrom20 else 7000

    #Select only fields of interest in va
    trios = trios.annotate_samples_expr(['sa.meta = select(sa.meta, population)'])

    trans_kt = trios.phase_by_transmission(ped, 'va.gene', n_partitions)
    trans_kt = trans_kt.key_by(['va.gene','v1','v2'])
    trans_kt = trans_kt.persist()

    trans_variants = trans_kt.query('v1.collectAsSet().union(v2.collectAsSet()).toArray()')
    logger.info("Found {} variants in pairs in trios.".format(len(trans_variants)))

    if args.reference:
        reference = hc.read(args.reference)
        reference = reference.filter_variants_list(trans_variants)
    else:
        reference = create_reference(hc, args, trios, trans_variants)

    logger.info("Found ~{} of these variants in the reference post-filtering".format(reference.count_variants()))

    reference_kt = reference.phase_em(['va.gene'], n_partitions,
                                      sa_keys='sa.pop' if args.split_by_pop else None,
                                      variant_pairs= trans_kt)

    logger.info("{} entries in reference_kt".format(reference_kt.count()))

    join_keys = [('v1','v1'), ('v2','v2')]
    if args.split_by_pop:
        join_keys.append(('pop','`sa.pop`'))

    reference_kt = reference_kt.aggregate_by_key(['{} = {}'.format(k,v) for k,v in join_keys],
                                                 ['haplotype_counts = haplotype_counts.takeBy(x => isMissing(x).toInt,1)[0]',
                                                  'genotype_counts = genotype_counts.takeBy(x => isMissing(x).toInt,1)[0]',
                                                  'prob_same_haplotype = prob_same_haplotype.takeBy(x => isMissing(x).toInt,1)[0]',
                                                  'ac1 = va1.takeBy(x => isMissing(x).toInt,1)[0].AC',
                                                  'ac2 = va2.takeBy(x => isMissing(x).toInt,1)[0].AC'])

    logger.info("{} entries in aggregated reference_ky".format(reference_kt.count()))

    trans_kt = trans_kt.annotate(['gene = `va.gene`', 'pop = kidSA.meta.population'])
    trans_kt = trans_kt.key_by([k for k,v in join_keys])

    phase_trios_kt = trans_kt.join(reference_kt, how="left")
    phase_trios_kt = phase_trios_kt.persist()
    phase_trios_kt.write(args.output + '.kt',overwrite=args.overwrite)

    phase_trios_kt, count_cols = flatten_counts(phase_trios_kt, gt_anns=['kid_v1','kid_v2','mom_v1','mom_v2','dad_v1','dad_v2'])

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
                                 'exome_coverage1 = va1.coverage.exome',
                                 'exome_coverage2 = va2.coverage.exome',
                                 'genome_coverage1 = va1.coverage.genome',
                                 'genome_coverage2 = va2.coverage.genome',
                                 'ref1 = v1.ref', 'alt1 = v1.alt',
                                 'ref2 = v2.ref', 'alt2 = v2.alt',
                                 'chrom1 = v1.contig', 'chrom2 = v2.contig',
                                 'pos1 = v1.start', 'pos2 = v2.start'
                                 ])
            .select(
            ['gene', 'chrom1', 'pos1', 'ref1', 'alt1', 'cpg1', 'pass1', 'impact1', 'alleleType1', 'ac1', 'ac_raw1',
             'chrom2', 'pos2', 'ref2', 'alt2', 'cpg2', 'exome_coverage1','exome_coverage2',
             'genome_coverage1','genome_coverage2',
             'pass2', 'impact2', 'alleleType2', 'ac2', 'ac_raw2',
             'fam', 'pop', 'prob_same_haplotype', 'same_trio_haplotype', 'distance',
             'wasSplit1', 'wasSplit2'] + count_cols)
            .export(args.output + '.txt.bgz')
    )

    if args.slack_channel:
        send_message(args.slack_channel, 'Trio compound hets %s is done processing!' % args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #Creates clean VDS with trios filtered down to correct variants
    parser.add_argument('--trios', help='Trio VDS. If not provided, a trios VDS will be generated', required=False)
    parser.add_argument('--fam_file', help='Fam file to use when --trios is given', required=False)
    parser.add_argument('--exomes_trios', help='Writes exomes trios VDS', required=False, action='store_true')
    parser.add_argument('--genomes_trios', help='Writes genomes trios VDS', required=False, action='store_true')
    parser.add_argument('--reference', help='Specifies file to use as reference, otherwise will write reference', required=False)
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

    if int(args.exomes_trios) + int(args.genomes_trios) + int(args.trios is not None) != 1:
        sys.exit("One and only one of --trios, --exomes_trios or --genomes_trios is required.")

    if args.trios and not args.fam_file:
        sys.exit("Must specify --fam_file when using --trios.")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

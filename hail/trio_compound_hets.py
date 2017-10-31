import argparse
from hail import *
from compound_hets_utils import *

subpop_tsv_path = "gs://gnomad-exomes/sampleqc/draft_subpops_exomes.txt"

def create_trio_vds(hc, args):
    trios = get_gnomad_data(hc,"exomes" if args.exomes_trios else "genomes")
    ped = Pedigree.read(exomes_fam_path if args.exomes_trios else genomes_fam_path)

    trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')
    trios = filter_low_conf_regions(trios, high_conf_regions=[exomes_high_conf_regions_intervals_path])
    ped = ped.filter_to(trios.sample_ids)
    logger.info("Found {0} trios in VDS.".format(len(ped.complete_trios()) / 3))

    trios = trios.split_multi()

    if args.filter_to_adj:
        trios = filter_to_adj(trios)

    if args.exomes_trios:
        trios = trios.annotate_variants_vds(hc.read(public_exomes_vds_path(split=True)), expr='va.filters = vds.filters,'
                                                                                      'va.release = select(vds.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)')
    else:
        trios = trios.annotate_variants_vds(hc.read(public_genomes_vds_path(split=True)), expr='va.filters = vds.filters,'
                                                                                      'va.release = select(vds.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)')

    trios = trios.annotate_variants_vds(
        get_gnomad_data(hc,"exomes" if args.exomes_trios else "genomes",hardcalls=True,split=True).drop_samples(),
        expr='va.vep = drop(vds.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences)'
    )

    if args.debug:
        trios = trios.persist()
        n_variants = trios.query_variants(['variants.map(x => va.filters.isEmpty()).counter()',
                                           'variants.map(v => va.release.AF <= {0}).counter()'.format(
                                               args.max_af)])

        logger.debug(
            "Number of variants PASS: {0}, Number of variants with AF below {1}: {2}".format(str(n_variants[0]),
                                                                                             args.max_af,
                                                                                             str(n_variants[1])))

    trios = trios.filter_variants_expr(
        'va.filters.isEmpty() && va.release.AF <= {0} && gs.filter(g => g.isCalledNonRef).count() > 0'.format(
            args.max_af))

    # Add methylated CpG annotation
    trios = trios.annotate_variants_vds(hc.read(context_vds_path), expr='va.methylated_cpg = vds.methylation.value >= 0.25,'
                                                                        'va.coverage = vds.coverage')

    # Add VEP annotations
    trios = annotate_gene_impact(trios)

    #Drop unused annotations
    trios = trios.annotate_variants_expr(['va = select(va, filters, info, coverage, methylated_cpg, release, vep, gene, impact, aIndex, wasSplit)'])

    outpath = args.output + '.adj.vds' if args.filter_to_adj else args.output + '.vds'

    trios.write(outpath, overwrite=args.overwrite)
    return hc.read(outpath), ped


def create_reference(hc, args, trios_vds, trans_variants):
    reference = get_gnomad_data(hc,
                                "exomes",
                                hardcalls= not args.filter_to_adj,
                                split=True)

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

    if args.write_kt:
        if args.trios:
            trios = hc.read(args.trios)
            ped = Pedigree.read(args.fam_file)
        else:
            trios, ped = create_trio_vds(hc, args)

        if args.split_by_subpop:
            trios = trios.annotate_samples_table(
                hc.import_table(subpop_tsv_path,
                                key="sample"),
                root = 'sa.meta.subpop'
            )

        if args.chrom20:
            trios = trios.filter_intervals(Interval.parse("20:1-10000000"))

        n_partitions = 10 if args.chrom20 else 500

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

        if args.split_by_subpop:
            reference = reference.annotate_samples_table(
                hc.import_table(subpop_tsv_path,
                                key="sample"),
                root='sa.subpop'
            )

        logger.info("Found ~{} of these variants in the reference post-filtering".format(reference.count_variants()))

        reference_kt = reference.phase_em(['va.gene'], n_partitions,
                                          sa_keys='sa.pop' if args.split_by_pop else 'sa.subpop' if args.split_by_subpop else None,
                                          variant_pairs= trans_kt)

        logger.debug("{} entries in reference_kt".format(reference_kt.count()))

        join_keys = [('v1','v1'), ('v2','v2')]
        if args.split_by_pop:
            join_keys.append(('pop','`sa.pop`'))
        elif args.split_by_subpop:
            join_keys.append(('subpop', '`sa.subpop`'))

        reference_kt = reference_kt.aggregate_by_key(['{} = {}'.format(k,v) for k,v in join_keys],
                                                     ['haplotype_counts = haplotype_counts.takeBy(x => isMissing(x).toInt,1)[0]',
                                                      'genotype_counts = genotype_counts.takeBy(x => isMissing(x).toInt,1)[0]',
                                                      'prob_same_haplotype = prob_same_haplotype.takeBy(x => isMissing(x).toInt,1)[0]',
                                                      'ac1 = va1.takeBy(x => isMissing(x).toInt,1)[0].AC',
                                                      'ac2 = va2.takeBy(x => isMissing(x).toInt,1)[0].AC'])

        logger.debug("{} entries in aggregated reference_ky".format(reference_kt.count()))

        trans_kt = trans_kt.annotate(['gene = `va.gene`', 'pop = kidSA.meta.population'])
        trans_kt = trans_kt.key_by([k for k,v in join_keys])

        phase_trios_kt = trans_kt.join(reference_kt, how="left")
        phase_trios_kt = phase_trios_kt.persist()

        #Get the full genotypes from the full trios VDS -- no choice :(
        phase_trios_kt = phase_trios_kt.drop(['kid_v1', 'kid_v2', 'mom_v1', 'mom_v2', 'dad_v1', 'dad_v2'])
        trios = trios.filter_variants_list(trans_variants)
        #trios_gt_kt.persist() -- not sure -- probably better to re-generate

        trio_roles = {
        "kid": [trio.proband for trio in ped.complete_trios()],
        "dad" : [trio.father for trio in ped.complete_trios()],
         "mom": [trio.mother for trio in ped.complete_trios()]
        }

        for name, samples in trio_roles.iteritems():
            trio_role_kt = (
                trios.filter_samples_list(samples)
                    .genotypes_table()
                    .select(['v', 's', 'g'])
            )
            for v in ["v1","v2"]:
                phase_trios_kt = phase_trios_kt.key_by([v,name]).join(
                    trio_role_kt.rename({'v': v, 's': name, 'g': "{}_{}".format(name, v)})
                        .key_by([v, name]),
                    how="left"
                )

        phase_trios_kt.write(args.output + '.kt', overwrite=args.overwrite)

    if args.write_results:
        phase_trios_kt = hc.read_table(args.output + '.kt')

        phase_trios_kt = conditional_column_swap(phase_trios_kt,
                                                 swap_expr='(isMissing(ac1) && isDefined(ac2)) || (isDefined(ac1) && isDefined(ac2) && ac1 > ac2)',
                                                 columns=[
                                                     ('{}1'.format(col),
                                                      '{}2'.format(col)) for col in [
                                                         'v',
                                                         'va',
                                                         'ac',
                                                         'kid_v',
                                                         'mom_v',
                                                         'dad_v'
                                                     ]
                                                 ],
                                     gt_counts_col="genotype_counts",
                                     hc_counts_col="haplotype_counts")

        phase_trios_kt, hc_count_cols = flatten_haplotype_counts(phase_trios_kt)
        phase_trios_kt, gt_cols = flatten_genotypes(phase_trios_kt,['kid_v1', 'kid_v2', 'mom_v1', 'mom_v2', 'dad_v1',
                                                                          'dad_v2'])
        variant_cols = {'chrom': 'v{}.contig',
                        'pos': 'v{}.start',
                        'ref': 'v{}.ref',
                        'alt': 'v{}.alt',
                        'cpg': 'va{}.methylated_cpg',
                        'pass': 'va{}.filters.isEmpty',
                        'impact': 'va{}.impact',
                        'ac_notrios': 'ac{}',
                        'ac_raw': 'va{0}.info.AC[va{0}.aIndex -1]',
                        'an_raw': 'va{}.info.AN',
                        'ac_release': 'va{}.release.AC',
                        'an_release': 'va{}.release.AN',
                        'exome_coverage': 'va{}.coverage.exome',
                        'genome_coverage': 'va{}.coverage.genome',
                        'wasSplit': 'va{}.wasSplit'}

        (
            phase_trios_kt.annotate(['fam = kidSA.fam.famID',
                                     'same_trio_haplotype = same_haplotype',
                                     'distance = (v1.start - v2.start).abs()',
                                     ] +
                                    ['{0}{1} = {2}'.format(name, n, expr.format(n)) for n in ["1","2"] for name, expr in variant_cols.iteritems()])
                .select(
                ['gene', 'fam', 'pop', 'prob_same_haplotype', 'same_trio_haplotype', 'distance'] +
                [col + n for n in ["1","2"] for col in variant_cols.keys()] +
                hc_count_cols + gt_cols)
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
    parser.add_argument('--split_by_subpop', help='Splits data by sub-population when computing EM', required=False,
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrites existing results.', required=False,
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--write_results', help='Writes results tsv.', required=False, action='store_true')
    parser.add_argument('--write_kt', help='Writes trio phased KeyTable.', required=False, action='store_true')
    args = parser.parse_args()

    if args.write_kt:
        if int(args.exomes_trios) + int(args.genomes_trios) + int(args.trios is not None) != 1:
            sys.exit("One and only one of --trios, --exomes_trios or --genomes_trios is required.")

        if args.trios and not args.fam_file:
            sys.exit("Must specify --fam_file when using --trios.")

    if args.write_results and not args.write_kt and (args.exomes_trios or
                                                         args.genomes_trios or
                                                         args.trios or
                                                         args.max_af or
                                                         args.split_by_pop or
                                                         args.split_by_subpop or
                                                         args.fam_file):
        print(
        "WARN: Some of the arguments provided will not have any influence since this is only exporting results from a previously computed KT")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

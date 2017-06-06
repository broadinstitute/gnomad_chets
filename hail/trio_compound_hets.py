import argparse
from hail import *
from compound_hets_utils import *


def main(args):

    hc = HailContext()

    if args.trios:

        if args.filter_to_adj and (args.exomes_trios or args.genomes_trios):
            trios = add_exomes_sa(hc.read(full_exome_vds)) if args.exomes_trios else add_genomes_sa(hc.read(full_genome_hardcalls_split_vds))
            trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')
            ann = hc.read(args.trios)
            trios = trios.annotate_variants_vds(ann, 'va = vds')
            trios = trios.filter_variants_expr('isDefined(va.gene)')
            trios = filter_to_adj(trios)
            trios.write(args.output + '.adj.vds', overwrite=True)
            trios = hc.read(args.output + '.adj.vds')
            logger.info("Count() in trios: {}".format(str(trios.count())))

        else:
            trios = hc.read(args.trios)
        ped = Pedigree.read(args.fam_file)
    else:
        if args.exomes_trios:
            trios = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds))
            trios = trios.annotate_variants_vds(hc.read(final_exome_split_vds), expr = 'va.filters = vds.filters')
            ped = Pedigree.read(exomes_fam)
        else:
            trios = add_genomes_sa(hc.read(full_genome_hardcalls_split_vds))
            trios = trios.annotate_variants_vds(hc.read(final_genome_split_vds), expr='va.filters = vds.filters')
            ped = Pedigree.read(genomes_fam)

        trios = trios.filter_samples_expr('isDefined(sa.fam.famID)')

        ped = ped.filter_to(trios.sample_ids)
        logger.info("Found {0} trios in VDS.".format(len(ped.complete_trios())/3))

        if args.debug:
            n_variants =  trios.query_variants(['variants.map(x => va.filters.isEmpty()).counter()',
                                                'variants.map(v => va.calldata.all_samples_raw.AF <= {0}).counter()'.format(args.max_af)])

            logger.debug("Number of variants PASS: {0}, Number of variants with AF below {1}: {2}".format(str(n_variants[0]),args.max_af, str(n_variants[1])))

        trios = trios.filter_variants_expr(
            'va.filters.isEmpty() && va.calldata.all_samples_raw.AF <= {0} && gs.filter(g => g.isCalledNonRef).count() > 0'.format(args.max_af))

        if args.exomes_trios:
            trios = trios.annotate_variants_vds(hc.read(full_exomes_vep_split_vds), expr = 'va.vep = vds.vep')
        else:
            trios = trios.annotate_variants_vds(hc.read(full_genomes_vep_split_vds), expr = 'va.vep = vds.vep')

        trios = filter_low_conf_regions(trios, high_conf_regions=exomes_high_conf_regions_path)

        #Add methylated CpG annotation
        trios = annotate_methylation(trios)

        #Add VEP annotations
        trios = annotate_gene_impact(trios)

        trios.write(args.output + '.vds')
        trios = hc.read(args.output + '.vds')

    if args.filter_to_adj:
        reference = add_exomes_sa(hc.read(full_genome_vds))
    else:
        reference = add_exomes_sa(hc.read(full_exome_hardcalls_split_vds))

    reference = reference.filter_samples_expr('!isDefined(sa.fam.famID) && sa.meta.drop_status == "keep"')

    if args.chrom20:
        trios = trios.filter_intervals(Interval.parse("20:1-10000000"))
        reference = reference.filter_intervals(Interval.parse("20:1-10000000"))

    if args.filter_to_adj:
        reference = filter_to_adj(reference)

    reference = reference.filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')

    n_partitions = 50 if args.chrom20 else 7000

    trios.phase_trios(
        reference_vds = reference,
        pedigree = ped,
        gene_ann = 'va.gene',
        output = args.output + '.stats.txt.bgz',
        num_partitions = n_partitions,
        va_strat = 'va.impact,va.alleleType,va.methylated_cpg',
        sa_strat = 'sa.fam.famID',
        run_coseg = False,
        run_em = True)

    if args.slack_channel:
        send_message(args.slack_channel, 'Trio compound hets %s is done processing!' % args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--trios', help='Trio VDS', required=False)
    parser.add_argument('--fam_file', help='Fam file to use when --trios is given', required=False)
    parser.add_argument('--exomes_trios', help='Writes exomes trios VDS or use exome VDS fif filter_to_adj', required=False, action='store_true') #Todo, this is not clean
    parser.add_argument('--genomes_trios', help='Writes genomes trios VDS', required=False, action='store_true')
    parser.add_argument('--max_af', help='Maximum AF for a site to be retained (default 0.01).', required=False, type=float, default=0.01)
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--debug', help='Output debug statements', required=False, action='store_true')
    parser.add_argument('--chrom20', help='Process chrom 20 only', required=False, action='store_true')
    parser.add_argument('--filter_to_adj', help='Use Adj genotypes only', required=False, action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    #if int(args.trios is not None) + int(args.exomes_trios) + int(args.genomes_trios) != 1:
    #    sys.exit("One and only one of --trios, --exomes_trios or --genomes_trios is required.")

    if args.trios and not args.fam_file:
        sys.exit("Must specify --fam_file when using --trios.")

    main(args)

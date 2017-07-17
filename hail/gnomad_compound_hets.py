import argparse
from compound_hets_utils import *
from hail import *


def main(args):

    hc = HailContext(log='/gnomad_compound_hets.log')

    if not args.vds:
        if args.exomes:
            release = final_exome_split_vds_path
            vep = full_exomes_vep_split_vds_path
            high_conf_regions = [exomes_high_conf_regions_intervals_path]
            vds = add_exomes_sa(hc.read(full_exome_vds_path))
            vds = vds.filter_samples_expr('sa.meta.drop_status == "keep"')
        else:
            release = final_genome_split_vds_path
            vep = full_genomes_vep_split_vds_path
            high_conf_regions = None
            vds = add_genomes_sa(hc.read(full_genome_vds_path))
            vds = vds.filter_samples_expr('sa.meta.keep')
            vds = vds.annotate_samples_expr(
                ['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop'])

        vds = filter_low_conf_regions(vds, high_conf_regions=high_conf_regions)
        vds = vds.split_multi()
        vds = vds.annotate_variants_vds(hc.read(release), root='va.release')
        vds = vds.annotate_variants_vds(hc.read(vep), expr='va.vep = vds.vep')
        vds = annotate_gene_impact(vds)
        vds = annotate_methylation(vds)
        vds = vds.filter_variants_expr('isDefined(va.gene) && isDefined(va.impact) && va.impact != "low" && va.release.info.AF_POPMAX <= {}'.format(args.max_af))
        if args.adj:
            vds = filter_to_adj(vds)

        vds.hardcalls().write(args.output + '.vds', overwrite = args.overwrite)
        vds = hc.read(args.output + '.vds')
    else:
        vds = hc.read(args.vds)

    if args.chrom20:
        vds = vds.filter_intervals(Interval.parse("20:1-10000000"))

    if args.write_pairs_kt:
        n_partitions = 50 if args.chrom20 else 7000
        vds = vds.annotate_variants_expr('va.release.info = select(va.release.info, AC, AN, AF, AC_NFE, AN_NFE, AF_NFE, AC_EAS, AN_EAS, AF_EAS, POPMAX, AC_POPMAX, AN_POPMAX)')
        vds = vds.annotate_variants_expr(['va.release = select(va.release, info, filters)',
                                         'va.vep = drop(va.vep, colocated_variants, intergenic_consequences, regulatory_feature_consequences, motif_feature_consequences) '])
        logger.info(str(vds.variant_schema))

        vds = vds.annotate_samples_expr('sa = {pop: sa.meta.population}')
        kt = vds.phase_em(['va.gene'], n_partitions, sa_keys='sa.pop', by_sample=args.by_sample)
        kt = kt.persist()

        kt.write(args.output + ".variant_pairs.kt", overwrite = args.overwrite)
        kt = hc.read_table(args.output + ".variant_pairs.kt")
        (
            kt.filter('prob_same_haplotype < 0.5')
                .annotate(['gene = `va.gene`', 'pop = `sa.pop`', 'sample = s',
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
                           'pos1 = v1.start', 'pos2 = v2.start',
                           'AABB = genotype_counts[0]',
                           'AaBB = genotype_counts[1]',
                           'aaBB = genotype_counts[2]',
                           'AABb = genotype_counts[3]',
                           'AaBb = genotype_counts[4]',
                           'aaBb = genotype_counts[5]',
                           'AAbb = genotype_counts[6]',
                           'Aabb = genotype_counts[7]',
                           'aabb = genotype_counts[8]',
                           'AB = haplotype_counts[0]',
                           'Ab = haplotype_counts[1]',
                           'aB = haplotype_counts[2]',
                           'ab = haplotype_counts[3]'
                           ])
                .select(
                ['gene', 'chrom1', 'pos1', 'ref1', 'alt1', 'cpg1', 'pass1', 'impact1', 'alleleType1', 'ac1', 'ac_raw1',
                 'chrom2', 'pos2', 'ref2', 'alt2', 'cpg2', 'pass2', 'impact2', 'alleleType2', 'ac2', 'ac_raw2',
                 'sample', 'pop', 'prob_same_haplotype', 'wasSplit1', 'wasSplit2',
                 'distance', 'AABB', 'AaBB', 'aaBB', 'AABb', 'AaBb', 'aaBb', 'AAbb', 'Aabb', 'aabb', 'AB', 'Ab', 'aB',
                 'ab'])
                .export(args.output + ".variant_pairs.txt.bgz")
        )

    if args.write_single_kt:
        kt = vds.genotypes_table()
        kt = kt.filter('g.isHomVar')
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

    if args.slack_channel:
        send_message(args.slack_channel, 'gnomAD compound hets {} is done processing!'.format(args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Use exomes and writes intermediate input', required=False, action='store_true')
    parser.add_argument('--genomes', help='Use genomes and writes intermediate input', required=False, action='store_true')
    parser.add_argument('--vds', help='Input VDS file for KT generation', required=False)
    parser.add_argument('--write_pairs_kt', help='Writes variant pairs KeyTable.', required=False, action='store_true')
    parser.add_argument('--write_single_kt', help='Writes single variants KeyTable.', required=False, action='store_true')
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

    if int(args.exomes) + int(args.genomes) + int(args.vds is not None) != 1:
        sys.exit("One and only one of --exomes, --genomes or --vds is required.")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
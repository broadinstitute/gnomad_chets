#!/usr/bin/env python

__author__ = 'konrad'

import argparse
from utils import *
import subprocess
import os
from hail import *

DOT_ANN_DICT = {
    'AS_RF_POSITIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_POSITIVE_TRAIN in orMissing(isDefined(oldTrain),'
                            'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                            'orMissing(!newTrain.isEmpty(),newTrain))',
    'AS_RF_NEGATIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_NEGATIVE_TRAIN in orMissing(isDefined(oldTrain),'
                            'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                            'orMissing(!newTrain.isEmpty(),newTrain))'
}
vqsr_vds_path = None
date_time = time.strftime("%Y-%m-%d_%H:%M")

genome_sa_drop = ["releasable", "gvcf_date", "qc_pop", "keep", "qc_sample", "BAM", "CLOUD_GVCF", "ON_PREM_GVCF", "final_pop_old"]

def select_annotations(vds):
    annotations_to_ignore = ['DB', 'GQ_HIST_ALL', 'DP_HIST_ALL', 'AB_HIST_ALL', 'GQ_HIST_ALT', 'DP_HIST_ALT',
                             'AB_HIST_ALT', 'A[CN]_..._.*ale']
    info = get_ann_type('va.info', vds.variant_schema)

    return [x.name for x in filter_annotations_regex(info.fields, annotations_to_ignore)]


def get_pops(vds, pop_path, min_count=10):
    subset_pops = vds.query_samples('samples.map(s => %s).counter()' % pop_path)
    return [pop.upper() for (pop, count) in subset_pops.items() if count >= min_count and pop is not None]


def fix_number_attributes(vds):
    info = get_ann_type('va.info',vds.variant_schema)

    for f in info.fields:
        for k,v in f.attributes.iteritems():
            if k == 'Number' and len(v) > 1:
                vds = vds.set_va_attributes('va.info.%s' % f.name, {'Number': v[0]})

    return vds

def get_subset_vds(hc, args):

    if args.exomes:
        vqsr_vds = hc.read(vqsr_vds_path)
        vds = hc.read(full_exome_vds)
    else:
        vqsr_vds = None
        vds = hc.read(full_genome_vds)
    pop_path = 'sa.meta.population' if args.exomes else 'sa.meta.final_pop'

    vds = preprocess_vds(vds, vqsr_vds, [], release=args.release_only)

    if args.genomes:
        vds = vds.annotate_samples_expr('sa.meta = drop(sa.meta, %s)' % ",".join(genome_sa_drop))

    if args.projects:
        list_data = read_list_data(args.projects)
        id_path = "sa.meta.pid" if args.exomes else "sa.meta.project_or_cohort"
        vds = (vds
               .annotate_global_py('global.projects', list_data, TSet(TString()))
               .filter_samples_expr('global.projects.contains(%s)' % id_path, keep=True))
    elif args.samples:
        list_data = read_list_data(args.samples)
        vds = (vds
               .annotate_global_py('global.samples', list_data, TSet(TString()))
               .filter_samples_expr('global.samples.contains(s)', keep=True))
    elif args.expr:
        vds = vds.filter_samples_expr(args.expr)
    else:
        print "Should not have gotten here. Need to add --projects --samples or --expr"
        sys.exit(1)

    vds = (vds
           .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
           .filter_alleles('va.calldata.raw.AC[aIndex] == 0', keep=False)
           .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
    )
    num_samples = vds.query_samples('samples.count()')
    if num_samples:
        logger.info('Got %s samples', num_samples)
    else:
        logger.critical('No samples found! Check input files')
        sys.exit(1)
    pops = get_pops(vds, pop_path)
    logger.info('Populations found: %s', pops)

    vds = (
        vds.annotate_global_py('global.pops', map(lambda x: x.lower(), pops), TArray(TString()))
        .persist()
    )

    return vds, pops


def main(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hc = HailContext(log='/hail.log')

    vds = None

    pops = None

    # Pre
    if not args.skip_pre_process:

        vds, pops = get_subset_vds(hc, args)

        create_sites_vds_annotations(vds, pops, dbsnp_vcf,
                                     filter_alleles=False,
                                     drop_star=False,
                                     generate_hists=False).write(args.output + ".pre.autosomes.sites.vds",
                                                                 overwrite=args.overwrite)
        create_sites_vds_annotations_X(vds, pops, dbsnp_vcf,
                                       filter_alleles=False,
                                       drop_star=False,
                                       generate_hists=False).write(args.output + ".pre.X.sites.vds",
                                                                   overwrite=args.overwrite)
        if args.exomes: create_sites_vds_annotations_Y(vds, pops, dbsnp_vcf,
                                                       filter_alleles=False,
                                                       drop_star=False).write(args.output + ".pre.Y.sites.vds",
                                                                              overwrite=args.overwrite)

    if not args.skip_merge:
        # Combine VDSes
        vdses = [hc.read(args.output + ".pre.autosomes.sites.vds"), hc.read(args.output + ".pre.X.sites.vds")]
        if args.exomes: vdses.append(hc.read(args.output + ".pre.Y.sites.vds"))
        vdses = merge_schemas(vdses)
        sites_vds = vdses[0].union(vdses[1:])
        sites_vds.write(args.output + '.pre.sites.vds', overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.sites.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.sites.vep.vds", overwrite=args.overwrite)
         )

    if not args.skip_post_process:
        # Post
        sites_vds = hc.read(args.output + ".pre.sites.vep.vds")

        if not args.no_annotations_with_both_release:
            release_dict = {
                'exomes': {'out_root': 'va.info.ge_', 'name': 'gnomAD exomes', 'vds': hc.read(final_exome_vds)},
                'genomes': {'out_root': 'va.info.gg_', 'name': 'gnomAD genomes', 'vds': hc.read(final_genome_vds)}
            }
        else:
            if args.exomes:
                release_dict = {
                    'exomes': {'out_root': 'va.info.ge_', 'name': 'gnomAD exomes', 'vds': hc.read(final_exome_vds)}
                }
            else:
                release_dict = {
                    'genomes': {'out_root': 'va.info.gg_', 'name': 'gnomAD genomes', 'vds': hc.read(final_genome_vds)}
                }

        key = 'exomes' if args.exomes else 'genomes'

        sites_vds = post_process_subset(sites_vds, release_dict, key, DOT_ANN_DICT)
        sites_vds = sites_vds.annotate_variants_expr('va.info = select(va.info, %s)' % ",".join(select_annotations(sites_vds)))

        sites_vds.write(args.output + ".sites.vds", overwrite=args.overwrite)

    if not args.skip_write_vds:
        if vds is None:
            vds, pops = get_subset_vds(hc, args)
        sites_vds = hc.read(args.output + ".sites.vds").min_rep()
        vds.annotate_variants_vds(sites_vds, 'va = vds').write(args.output + ".vds", overwrite=args.overwrite)

    vds = hc.read(args.output + ".vds")

    if not args.skip_sanity_checks:
        if pops is None:
            subset_vds, pops = get_subset_vds(hc, args)
        sites_sanity_check_text = run_sites_sanity_checks(vds, pops, return_string=True, skip_star=True)
        samples_sanity_check_text = run_samples_sanity_checks(vds,
                                                              hc.read(full_exome_vds) if args.exomes else hc.read(full_genome_vds),
                                                              n_samples=10, verbose=True, return_string=True)
        if args.slack_channel:
            send_snippet(args.slack_channel, sites_sanity_check_text, 'sites_sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))
            send_snippet(args.slack_channel, samples_sanity_check_text,
                         'samples_sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))
        else:
            logger.info(sites_sanity_check_text)
            logger.info(samples_sanity_check_text)

    if not args.skip_write_vcf:
        vds = fix_number_attributes(vds)

        as_filter_status_fields = ['va.info.AS_FilterStatus']
        if args.no_annotations_with_both_release:
            if args.exomes:
                as_filter_status_fields.append('va.info.ge_AS_FilterStatus')
            else:
                as_filter_status_fields.append('va.info.gg_AS_FilterStatus')
        else:
            as_filter_status_fields.extend(['va.info.ge_AS_FilterStatus','va.info.gg_AS_FilterStatus'])

        if args.exomes:
            vds = vds.filter_variants_intervals(IntervalTree.read(exome_calling_intervals))

        if args.write_vcf_per_chrom:
            for contig in range(1,23):
                write_vcfs(vds, contig, args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                           as_filter_status_fields=as_filter_status_fields,
                           append_to_header=additional_vcf_header)
            write_vcfs(vds, 'X', args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                       as_filter_status_fields=as_filter_status_fields,
                       append_to_header=additional_vcf_header)
            if args.exomes:
                write_vcfs(vds, 'Y', args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                           as_filter_status_fields=as_filter_status_fields,
                           append_to_header=additional_vcf_header)
        else:
            write_vcfs(vds, '', args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                       as_filter_status_fields=as_filter_status_fields,
                       append_to_header=additional_vcf_header)

        vds.export_samples(args.output + '.sample_meta.txt.bgz', 'sa.meta.*')

    if args.slack_channel:
        send_message(args.slack_channel, 'Subset %s is done processing!' % args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--release_only', help='Whether only releaseables should be included in subset (default: False)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--no_annotations_with_both_release', help='If set, only annotates genomes with release genomes or exomes with release exomes.', action='store_true')
    parser.add_argument('--projects', help='File with projects to subset')
    parser.add_argument('--samples', help='File with samples to subset')
    parser.add_argument('--expr', help='Expression to subset')
    parser.add_argument('--skip_pre_process', help='Skip pre-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_merge', help='Skip merge step (assuming already done)', action='store_true')
    parser.add_argument('--skip_vep', help='Skip VEP (assuming already done)', action='store_true')
    parser.add_argument('--skip_post_process', help='Skip post-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_write_vds', help='Skip writing final VDS (assuming already done)', action='store_true')
    parser.add_argument('--skip_sanity_checks', help='Skip sanity checks', action='store_true')
    parser.add_argument('--skip_write_vcf', help='Skip writing VCF', action='store_true')
    parser.add_argument('--write_vcf_per_chrom', help='If set, generates a VCF for each chromosome. Otherwise, creates a single VCF.', action='store_true')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if int(args.samples is not None) + int(args.projects is not None) + int(args.expr is not None) != 1:
        sys.exit('Error: One and only one of --samples or --projects or --expr must be specified')

    if args.exomes:
        from exomes_sites_vcf import preprocess_vds, vqsr_vds_path, RF_SNV_CUTOFF, RF_INDEL_CUTOFF
    else:
        from genomes_sites_vcf import preprocess_vds, RF_SNV_CUTOFF, RF_INDEL_CUTOFF

    main(args)

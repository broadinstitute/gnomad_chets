#!/usr/bin/env python

__author__ = 'konrad'

import argparse
import gzip
from utils import *
import subprocess
import os
from hail import *

DOT_ANN_DICT  = {
    'AS_RF_POSITIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_POSITIVE_TRAIN in orMissing(isDefined(oldTrain),'
                            'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                            'orMissing(!newTrain.isEmpty(),newTrain))',
    'AS_RF_NEGATIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_NEGATIVE_TRAIN in orMissing(isDefined(oldTrain),'
                            'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                            'orMissing(!newTrain.isEmpty(),newTrain))'
}


def read_list_data(input_file):
    if input_file.startswith('gs://'):
        subprocess.check_output(['gsutil', 'cp', input_file, '.'])
        f = gzip.open(os.path.basename(input_file)) if input_file.endswith('gz') else open(os.path.basename(input_file))
    else:
        f = gzip.open(input_file) if input_file.endswith('gz') else open(input_file)
    projects = set()
    for line in f:
        projects.add(line.strip())
    f.close()
    return projects


def get_pops(vds, pop_path, min_count=10):
    subset_pops = vds.query_samples('samples.map(s => %s).counter()' % pop_path)
    return [pop.upper() for (pop, count) in subset_pops.items() if count >= min_count and pop is not None]


def get_subset_vds(hc, args):

    if args.exomes:
        vqsr_vds = hc.read(vqsr_vds_path)
        vds = hc.read(full_exome_vds)
    else:
        vqsr_vds = None
        vds = hc.read(full_genome_vds)
    pop_path = 'sa.meta.population' if args.exomes else 'sa.meta.final_pop'

    if args.projects:
        data_type = 'projects'
        list_data = read_list_data(args.projects)
        id_path = "sa.meta.pid" if args.exomes else "sa.meta.project_or_cohort"
    else:
        data_type = 'samples'
        list_data = read_list_data(args.samples)
        id_path = "s"

    vds = preprocess_vds(vds, vqsr_vds, [], release=args.release_only)

    vds = (vds
            .annotate_global_py('global.%s' % data_type, list_data, TSet(TString()))
            .filter_samples_expr('global.%s.contains(%s)' % (data_type, id_path), keep=True)
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
        release_dict = {
            'exomes': {'out_root': 'va.info.ge_', 'name': 'gnomAD exomes', 'vds': hc.read(final_exome_vds)},
            'genomes': {'out_root': 'va.info.gg_', 'name': 'gnomAD genomes', 'vds': hc.read(final_genome_vds)}
        }
        key = 'exomes' if args.exomes else 'genomes'

        post_process_subset(sites_vds, release_dict, key, DOT_ANN_DICT).write(args.output + ".sites.vds", overwrite=args.overwrite)

    sites_vds = hc.read(args.output + ".sites.vds")

    if vds is None:
        vds = hc.read(full_exome_vds) if args.exomes else hc.read(full_genome_vds)
        vds, pops = get_subset_vds(hc, args)

    if not args.skip_write_vds:
        vds.annotate_variants_vds(sites_vds, 'va = vds').write(args.output + ".vds", overwrite=args.overwrite)

    vds = hc.read(args.output + ".vds")

    if not args.skip_sanity_checks:
        sanity_check_text = run_sanity_checks(vds, pops, return_string=True, skip_star=True)
        if args.slack_channel:
            send_snippet(args.slack_channel, sanity_check_text, 'sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))
        else:
            logger.info(sanity_check_text)

    if args.exomes:
        vds = vds.filter_variants_intervals(IntervalTree.read(exome_calling_intervals))
        write_vcfs(vds, '', args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   as_filter_status_fields=['va.info.AS_FilterStatus', 'va.info.ge_AS_FilterStatus', 'va.info.gg_AS_FilterStatus'],
                   append_to_header=additional_vcf_header)
    else:
        for contig in range(1,23):
            write_vcfs(vds, contig, args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                       as_filter_status_fields=(
                       'va.info.AS_FilterStatus', 'va.info.ge_AS_FilterStatus', 'va.info.gg_AS_FilterStatus'),
                       append_to_header=additional_vcf_header)
        write_vcfs(vds, 'X', args.output, None, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   as_filter_status_fields=(
                   'va.info.AS_FilterStatus', 'va.info.ge_AS_FilterStatus', 'va.info.gg_AS_FilterStatus'),
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
    parser.add_argument('--projects', help='File with projects to subset')
    parser.add_argument('--samples', help='File with samples to subset')
    parser.add_argument('--skip_pre_process', help='Skip pre-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_merge', help='Skip merge step (assuming already done)', action='store_true')
    parser.add_argument('--skip_vep', help='Skip VEP (assuming already done)', action='store_true')
    parser.add_argument('--skip_post_process', help='Skip post-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_write_vds', help='Skip writing final VDS (assuming already done)', action='store_true')
    parser.add_argument('--skip_sanity_checks', help='Skip sanity checks', action='store_true')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    args = parser.parse_args()

    if args.exomes and args.genomes:
        sys.exit('Error: Only one of --exomes and --genomes can be specified')

    if not args.exomes and not args.genomes:
        sys.exit('Error: One of --exomes or --genomes must be specified')

    if args.samples and args.projects:
        sys.exit('Error: Only one of --samples and --projects can be specified')

    if not args.samples and not args.projects:
        sys.exit('Error: One of --samples or --projects must be specified')

    if args.exomes:
        from exomes_sites_vcf import *
    else:
        from genomes_sites_vcf import *

    main(args)

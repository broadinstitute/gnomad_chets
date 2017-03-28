#!/usr/bin/env python

__author__ = 'konrad'

import argparse
import gzip
from utils import *
import subprocess
import os
from hail import *


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


def main(args):
    if args.projects:
        data_type = 'projects'
        list_data = read_list_data(args.projects)
        id_path = "sa.meta.pid" if args.exomes else "sa.meta.project_or_cohort"
    else:
        data_type = 'samples'
        list_data = read_list_data(args.samples)
        id_path = "s.id"

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hc = HailContext(log='/hail.log')

    pop_path = 'sa.meta.population' if args.exomes else 'sa.meta.final_pop'

    # Pre
    if not args.skip_pre_process:
        if args.exomes:
            vqsr_vds = hc.read(vqsr_vds_path)
            vds = hc.read(full_exome_vds)
        else:
            vqsr_vds = None
            vds = hc.read(full_genome_vds)

        vds = preprocess_vds(vds, vqsr_vds, [], release=args.release_only)

        vds = (vds
               .annotate_global_py('global.%s' % data_type, list_data, TSet(TString()))
               .filter_samples_expr('global.%s.contains(%s)' % (data_type, id_path), keep=True))
        logger.info('Got %s samples', vds.query_samples('samples.count()'))
        pops = get_pops(vds, pop_path)
        logger.info('Populations found: %s', pops)

        vds = vds.annotate_global_py('global.pops', map(lambda x: x.lower(), pops), TArray(TString()))

        create_sites_vds_annotations(vds, pops, dbsnp_vcf, False, False).write(args.output + ".pre.autosomes.vds", overwrite=args.overwrite)
        create_sites_vds_annotations_X(vds, pops, dbsnp_vcf, False, False).write(args.output + ".pre.X.vds", overwrite=args.overwrite)
        if args.exomes: create_sites_vds_annotations_Y(vds, pops, dbsnp_vcf, False, False).write(args.output + ".pre.Y.vds", overwrite=args.overwrite)

        # Combine VDSes
        auto_vds = hc.read(args.output + ".pre.autosomes.vds")
        x_vds = hc.read(args.output + ".pre.X.vds")
        y_vds = hc.read(args.output + ".pre.Y.vds")
        vds = auto_vds.union([x_vds, y_vds])
        vds.write(args.output + 'pre.vds', overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.vep.vds", overwrite=args.overwrite)
         )

    if not args.skip_post_process:
        # Post
        vds = hc.read(args.output + ".pre.vep.vds")
        dot_ann_dict = {
            'AS_RF_POSITIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_POSITIVE_TRAIN in orMissing(isDefined(oldTrain),'
                                    'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                                    'orMissing(!newTrain.isEmpty(),newTrain))',
            'AS_RF_NEGATIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_NEGATIVE_TRAIN in orMissing(isDefined(oldTrain),'
                                    'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                                    'orMissing(!newTrain.isEmpty(),newTrain))'
        }
        release_dict = {
            'exomes': {'out_root': 'va.info.ge_', 'name': 'gnomAD exomes', 'vds': hc.read(final_exome_vds)},
            'genomes': {'out_root': 'va.info.gg_', 'name': 'gnomAD genomes', 'vds': hc.read(final_genome_vds)}
        }
        key = 'exomes' if args.exomes else 'genomes'

        post_process_subset(vds, release_dict, key, dot_ann_dict).write(args.output + ".vds", overwrite=args.overwrite)

    vds = hc.read(args.output + ".vds")
    pops = get_pops(vds, pop_path)
    sanity_check = run_sanity_checks(vds, pops, return_string=True, skip_star=True)
    if args.slack_channel:
        send_snippet(args.slack_channel, sanity_check, 'sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))

    vds = hc.read(args.output + ".vds").filter_variants_intervals(IntervalTree.read(exome_calling_intervals))
    write_vcfs(vds, '', args.output + '.internal', args.output, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
               as_filter_status_fields=('va.info.AS_FilterStatus', 'va.info.ge_AS_FilterStatus', 'va.info.gg_AS_FilterStatus'),
               append_to_header=additional_vcf_header)
    vds.export_samples(args.output + '.sample_meta.txt.bgz', 'sa.meta.*')

    if args.slack_channel:
        send_message(channel=args.slack_channel, message='Subset %s is done processing!' % args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--release_only', help='Whether only releaseables should be included in subset (default: False)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--projects', help='File with projects to subset')
    parser.add_argument('--samples', help='File with samples to subset')
    parser.add_argument('--skip_pre_process', help='Skip pre-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_post_process', help='Skip pre-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_vep', help='Skip pre-processing (assuming already done)', action='store_true')
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

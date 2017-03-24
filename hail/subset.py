#!/usr/bin/env python

__author__ = 'konrad'

import argparse
import gzip
from utils import *
import subprocess
import os
from hail import *


def read_projects(project_file):
    if project_file.startswith('gs://'):
        subprocess.check_output(['gsutil', 'cp', project_file, '.'])
        f = gzip.open(os.path.basename(project_file)) if project_file.endswith('gz') else open(os.path.basename(project_file))
    else:
        f = gzip.open(project_file) if project_file.endswith('gz') else open(project_file)
    projects = set()
    for line in f:
        projects.add(line.strip())
    f.close()
    return projects


def get_pops(vds, pop_path, min_count=10):
    subset_pops = vds.query_samples('samples.map(s => %s).counter()' % pop_path)
    return [pop.upper() for (pop, count) in subset_pops.items() if count >= min_count and pop is not None]


def main(args):
    projects = read_projects(args.projects)

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hc = HailContext(log='/hail.log')

    pop_path = 'sa.meta.population' if args.exomes else 'sa.meta.final_pop'

    # Pre
    if not args.skip_pre_process:
        if args.exomes:
            vqsr_vds = hc.read(vqsr_vds_path)
            vds = preprocess_vds(hc.read(full_exome_vds), vqsr_vds, release=args.release_only)
            pid_path = "sa.meta.pid"

        else:
            vds = preprocess_vds(hc.read(full_genome_vds), vqsr_vds=None, release=args.release_only)
            pid_path = "sa.meta.project_or_cohort"

        vds = (vds
               .annotate_global_py('global.projects', projects, TSet(TString()))
               .filter_samples_expr('global.projects.contains(%s)' % pid_path , keep=True))
        pops = get_pops(vds, pop_path)

        create_sites_vds_annotations(
            vds,
            pops,
            dbsnp_path=dbsnp_vcf,
            drop_star=False,
            drop_samples=False
        ).write(args.output + ".pre.autosomes.vds", overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.autosomes.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.vep.autosomes.vds", overwrite=args.overwrite)
         )

    if not args.skip_post_process:
        # Post
        vds = hc.read(args.output + ".pre.vep.autosomes.vds")
        dot_ann_dict = {
            'AS_RF_POSITIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_POSITIVE_TRAIN in orMissing(isDefined(oldTrain),'
                                    'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                                    'orMissing(!newTrain.isEmpty(),newTrain))',
            'AS_RF_NEGATIVE_TRAIN': '%s = let oldTrain = vds.find(x => isDefined(x)).info.AS_RF_NEGATIVE_TRAIN in orMissing(isDefined(oldTrain),'
                                    'let newTrain = range(aIndices.length).filter(i => oldTrain.toSet.contains(aIndices[i])) in '
                                    'orMissing(!newTrain.isEmpty(),newTrain))'
        }
        release_dict = {
            'exomes': {'out_root': 'va.info.ge_', 'name': 'gnomAD genomes', 'vds': hc.read(final_exome_autosomes)},
            'genomes': {'out_root': 'va.info.gg_', 'name': 'gnomAD genomes', 'vds':hc.read(final_genome_autosomes)}
        }
        key = 'exomes' if args.exomes else 'genomes'

        post_process_subset(vds, release_dict,
                            key,
                            dot_annotations_dict=dot_ann_dict).write(args.output + ".autosomes.vds", overwrite=args.overwrite)

    vds = hc.read(args.output + ".autosomes.vds")
    pops = get_pops(vds, pop_path)
    logger.info('Got %s populations', len(pops))
    logger.info(', '.join(pops))
    sanity_check = run_sanity_checks(vds, pops, return_string=True, skip_star=True)
    if args.slack_channel:
        send_snippet(args.slack_channel, sanity_check, 'autosome_sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))

    vds = (hc.read(args.output + ".autosomes.vds").filter_variants_intervals(IntervalTree.read(autosomes_intervals))
           .filter_variants_intervals(IntervalTree.read(exome_calling_intervals)))
    write_vcfs(vds, '', args.output + '.internal', args.output, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header)
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
        sys.exit('Error: One of --exomes and --genomes must be specified')

    if args.exomes:
        from exomes_sites_vcf import *
    else:
        from genomes_sites_vcf import *

    main(args)

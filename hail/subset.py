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


def main(args):
    projects = read_projects(args.projects)

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hc = HailContext(log='/hail.log')

    if args.exomes:
        vds = hc.read(full_exome_vds)
        vqsr_vds = hc.read(vqsr_vds_path)
        pid_path = "sa.meta.pid"
        pop_path = "sa.meta.population"

    else:
        vds = (
            hc.read(full_genome_vds)
                .annotate_samples_table(genomes_meta, 'Sample', root='sa.meta',
                                        config=hail.TextTableConfig(impute=True))
        )
        vds = vds.filter_variants_intervals(IntervalTree.parse_all(['22']))
        vqsr_vds = None
        pid_path = "sa.meta.project_or_cohort"
        pop_path = "sa.meta.final_pop"

    vds = (vds
           .annotate_global_py('global.projects', projects, TSet(TString()))
           .filter_samples_expr('global.projects.contains(%s)' % pid_path , keep=True))
    subset_pops = vds.query_samples('samples.map(s => %s).counter()' % pop_path)
    pops = [pop.upper() for (pop, count) in subset_pops.items() if count >= 10 and pop is not None]

    # Pre
    if not args.skip_pre_process:
        create_sites_vds_annotations(
            preprocess_vds(vds, vqsr_vds, release=args.release_only),
            pops,
            dbsnp_path=dbsnp_vcf,
            drop_star=False
        ).write(args.output + ".pre.autosomes.vds", overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.autosomes.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.vep.autosomes.vds", overwrite=args.overwrite)
         )

    if not args.skip_post_process:
        # Post
        vds = hc.read(args.output + ".pre.vep.autosomes.vds")
        release_dict = {
            'ge_': hc.read(final_exome_autosomes),
            'gg_': hc.read(final_genome_autosomes)
        }
        key = 'ge_' if args.exomes == 'exomes' else 'gg_'
        as_filter_attr = release_dict[key].get_va_attributes('va.info.AS_FilterStatus')

        post_process_subset(vds, release_dict,
                            'va.info.%sAS_FilterStatus' % key,
                            as_filter_attr).write(args.output + ".autosomes.vds", overwrite=args.overwrite)

        vds = hc.read(args.output + ".autosomes.vds")
        sanity_check = run_sanity_checks(vds, pops, return_string=True)
        if args.slack_channel:
            send_snippet(args.slack_channel, sanity_check, 'autosome_sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))

        vds = hc.read(args.output + ".autosomes.vds").filter_variants_intervals(autosomes_intervals).filter_variants_intervals(exome_calling_intervals)
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

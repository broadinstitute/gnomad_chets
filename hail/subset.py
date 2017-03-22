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

    hc = HailContext(log='/hail.log')

    if args.input == 'exomes':
        vds = hc.read(full_exome_vds)
        vqsr_vds = hc.read(vqsr_vds_path)
    else:
        vds = hc.read(full_genome_vds)
        vqsr_vds = None

    vds = (vds
           .annotate_global_py('global.projects', projects, TSet(TString()))
           .filter_samples_expr('global.projects.contains(sa.meta.pid)', keep=True))
    subset_pops = vds.query_samples('samples.map(s => sa.meta.population).counter()')
    pops = [pop for (pop, count) in subset_pops.items() if count >= 10 and pop is not None]

    # Pre
    if not args.skip_pre_process:
        create_sites_vds_annotations(
            preprocess_vds(vds, vqsr_vds, release=args.release_only),
            pops,
            dbsnp_path=dbsnp_vcf,
            drop_star=False
        ).write(args.output + ".pre.autosomes.vds", overwrite=args.overwrite)

    if not args.skip_post_process:
        # Post
        vds = hc.read(args.output + ".pre.autosomes.vds")
        release_dict = {
            'ge_': hc.read(final_exome_autosomes),
            'gg_': hc.read(final_genome_autosomes)
        }
        as_filter_attr = release_dict['ge_'].get_va_attributes('va.info.AS_FilterStatus')

        post_process_subset(vds, release_dict,
                            'va.info.ge_AS_FilterStatus',
                            as_filter_attr).write(args.output + ".autosomes.vds", overwrite=args.overwrite)

        vds = hc.read(args.output + ".autosomes.vds")
        sanity_check = run_sanity_checks(vds, pops, return_string=True)
        send_snippet('@konradjk', sanity_check, 'autosome_sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))

        vds = hc.read(args.output + ".autosomes.vds").filter_variants_intervals(autosomes_intervals).filter_variants_intervals(exome_calling_intervals)
        write_vcfs(vds, '', args.output + '.internal', args.output, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header)
        vds.export_samples(args.output + '.sample_meta.txt.bgz', 'sa.meta.*')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', help='Input VDS: exomes or genomes or file path', default='exomes')
    parser.add_argument('--release_only', help='Whether only releaseables should be included in subset (default: False)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--projects', help='File with projects to subset')
    parser.add_argument('--skip_pre_process', help='Skip pre-processing (assuming already done)', action='store_true')
    parser.add_argument('--skip_post_process', help='Skip pre-processing (assuming already done)', action='store_true')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    args = parser.parse_args()

    if args.input == 'exomes':
        from exomes_sites_vcf import *
    elif args.input == 'genomes':
        from genomes_sites_vcf import *

    main(args)

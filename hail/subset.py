#!/usr/bin/env python

__author__ = 'konrad'

import argparse
import gzip
from utils import *


def read_projects(project_file):
    projects = set()
    f = gzip.open(project_file) if project_file.endswith('gz') else open(project_file)
    for line in f:
        projects.add(line.strip())
    f.close()
    return projects


def main(args, pops):
    projects = read_projects(args.projects)

    hc = HailContext()

    create_sites_vds_annotations(
        preprocess_vds(args.input, release=args.release_only)
        .annotate_global_py('global.projects', projects, TSet(TString()))
        .filter_samples_expr('global.projects.contains(sa.meta.pid)', keep=True),
        pops,
        dbsnp_path=dbsnp_vcf,
        drop_samples=False
    ).write(args.output + ".pre.autosomes.vds")

    rf_vds = hc.read(rf_path)
    post_process_vds(hc, args.output + ".pre.autosomes.vds",
                     rf_vds,
                     RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                     'va.rf').write(args.output + ".autosomes.vds", overwrite=True)

    vds = hc.read(args.output + ".autosomes.vds")
    sanity_check = run_sanity_checks(vds, pops, return_string=send_to_slack)
    if send_to_slack: send_snippet('@konradjk', sanity_check, 'autosome_sanity_%s.txt' % date_time)

    vds = hc.read(args.output + ".autosomes.vds").filter_variants_intervals(autosomes_intervals).filter_variants_intervals(exome_calling_intervals)
    write_vcfs(vds, '', args.output + '.internal', args.output, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', help='Input VDS: exomes or genomes or file path', default='exomes')
    parser.add_argument('--release_only', help='Whether only releaseables should be included in subset (default: False)', action='store_true')
    parser.add_argument('--projects', help='File with projects to subset')
    parser.add_argument('--output', '-o', help='Output prefix')
    args = parser.parse_args()

    if args.input == 'exomes':
        from exomes_sites_vcf import *
        args.input = 'gs:///gnomad-exomes-raw/full/gnomad.exomes.all.vds'
    elif args.input == 'genomes':
        from genomes_sites_vcf import *
        args.input = ''

    main(args, pops)

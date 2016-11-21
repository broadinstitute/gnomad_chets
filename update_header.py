#!/usr/bin/env python

__author__ = 'konrad'

import sys
import argparse
import gzip
from collections import defaultdict
from vcf_parser import HeaderParser
import pipes

categories = {
    'AFR': 'African/African American',
    'AMR': 'Latino',
    'ASJ': 'Ashkenazi Jewish',
    'EAS': 'East Asian',
    'FIN': 'Finnish',
    'NFE': 'Non-Finnish European',
    'OTH': 'Other (population not assigned)',
    'SAS': 'South Asian',
    'Male': 'Male',
    'Female': 'Female',
    'Adj': 'Adjusted'
}

additional_descriptions = {
    'SOR': ('Strand Odds Ratio', '1'),
    'POPMAX': ('Population with max AF', 'A'),
    'AC_POPMAX': ('AC in the population with the max AF', 'A'),
    'AN_POPMAX': ('AN in the population with the max AF', 'A'),
    'AF_POPMAX': ('Maximum Allele Frequency across populations (excluding OTH)', 'A'),
    'AB_HIST': ('Histogram for Allele Balance; 100*AD[i_alt]/sum(AD); Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5', 'R'),
    'Hom': ('Count of homozygous individuals', 'A'),
    'Hemi': ('Count of hemizygous individuals', 'A'),
    'PROJECTMAX': ('Project with most individuals carrying alt allele', 'A'),
    'PROJECTMAX_NSamples': ('Number of samples in project with most individuals carrying alt allele', 'A'),
    'PROJECTMAX_NonRefSamples': ('Number of samples carrying alt allele in project with most individuals carrying alt allele', 'A'),
    'PROJECTMAX_PropNonRefSamples': ('Proportion of samples carrying alt allele in project with most individuals carrying alt allele', 'A')
}


def read_header(input):
    f = gzip.open(input) if input.endswith('gz') else open(input)
    header = HeaderParser()
    for line in f:
        if not line.startswith('##'): break
        header.parse_meta_data(line.strip())
    f.close()
    return header


def edit_line(line, label='Allele Count', data_type='Integer', number='A'):
    data = line['ID'].split('_')
    if len(data) == 2:
        ac, source = data
        if source not in categories:
            print >> sys.stderr, '%s not found' % source
            return line
        line['Description'] = '%s %s' % (categories[source], label)
    elif len(data) == 3:
        ac, source1, source2 = data
        if source1 not in categories or source2 not in categories:
            print >> sys.stderr, '%s or %s not found' % (source1, source2)
            return line
        line['Description'] = '%s among %s %ss' % (label, categories[source1], categories[source2])
    elif len(data) > 3:
        print >> sys.stderr, 'Got >3 values at %s' % line['ID']
    line['Type'] = data_type
    line['Number'] = number
    return line


def main(args):
    ref_header = read_header(args.reference)
    old_header = read_header(args.input)

    new_header = HeaderParser()

    new_header.add_fileformat('VCFv4.2')

    # for line in old_header.filter_lines:
    #     line['Description'] = [x for x in ref_header.filter_lines if x['ID'] == line['ID']][0]

    for line in old_header.info_lines:
        if line['ID'] in additional_descriptions:
            additional_data = additional_descriptions[line['ID']]
            if len(additional_data) == 2:
                line['Description'] = additional_data[0]
                line['Number'] = additional_data[1]
            else:
                line['Description'] = additional_data
        elif line['ID'].startswith('AC_'):
            line = edit_line(line)
        elif line['ID'].startswith('Hemi_'):
            line = edit_line(line, label='Hemizygote Count')
        elif line['ID'].startswith('Hom_'):
            line = edit_line(line, label='Homozygote Count')
        elif line['ID'].startswith('AF_'):
            line = edit_line(line, label='Allele Frequency', data_type='Float')
        elif line['ID'].startswith('AN_'):
            line = edit_line(line, label='Chromosome Count', number='1')
        elif line['ID'].startswith('GC_'):
            line = edit_line(line, label='Genotype Count', number='G')
        elif line['ID'] == 'GC':
            line['Number'] = 'G'
            line['Description'] = 'Counts of individuals carrying each genotype (Genotype Count)'
        elif line['ID'] in ref_header.extra_info and line['ID'] != 'CSQ':
            line['Description'] = ref_header.extra_info[line['ID']]['Description']
            line['Number'] = ref_header.extra_info[line['ID']]['Number']
        new_header.add_info(line['ID'], line['Number'], line['Type'], line['Description'])

    for line in old_header.contig_lines:
        new_header.add_contig(line['ID'], line['length'])

    if args.output.endswith('gz'):
        pipe = pipes.Template()
        pipe.append('bgzip -c /dev/stdin', '--')
        g = pipe.open(args.output, 'w')
    else:
        g = sys.stdout if args.output == 'stdout' else open(args.output, 'w')
    for line in new_header.print_header():
        if line.startswith('#CHROM'):
            g.write('##reference=file:///seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta\n')
        g.write(line + '\n')

    if args.full_vcf:
        f = gzip.open(args.input) if args.input.endswith('.gz') else open(args.input)
        for line in f:
            if line.startswith('#'): continue
            g.write(line)
        f.close()
    g.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--reference', '-r', help='Reference header', default='/humgen/atgu1/fs03/konradk/exac/v1/final/ExAC.r1.sites.vep.vcf.gz')
    parser.add_argument('--input', '-i', help='Header/VCF to annotate')
    parser.add_argument('--output', '-o', help='Output header')
    parser.add_argument('--full_vcf', '-f', help='Write full VCF (not just header)', action='store_true')
    args = parser.parse_args()
    main(args)

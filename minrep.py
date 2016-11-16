#!/usr/bin/env python

__author__ = 'konradjk'

import argparse
import gzip
import pipes
import sys
from minimal_representation import get_minimal_representation


def main(args):
    if args.vcf == 'stdin':
        f = sys.stdin
    else:
        f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)

    if args.output is None: args.output = args.vcf.replace('.vcf', '.minrep.vcf')
    if args.output == args.vcf:
        print >> sys.stderr, "VCF filename has no '.vcf' and no output file name was provided. Exiting."
        sys.exit(1)
    if not args.output.endswith('.gz'): args.output += '.gz'

    pipe = pipes.Template()
    pipe.append('bgzip -c /dev/stdin', '--')
    g = pipe.open(args.output, 'w')

    header = None
    for line in f:
        line = line.strip()

        # Reading header lines to get VEP and individual arrays
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                header = line.split('\t')
                header = dict(zip(header, range(len(header))))
            print >> g, line
            continue

        if header is None:
            print >> sys.stderr, "VCF file does not have a header line (CHROM POS etc.). Exiting."
            sys.exit(1)

        fields = line.split('\t')
        if fields[header['ALT']] == '*': continue

        alts = fields[header['ALT']].split(',')
        if len(alts) == 1:
            fields[header['POS']], fields[header['REF']], fields[header['ALT']] = map(str, get_minimal_representation(fields[header['POS']], fields[header['REF']], fields[header['ALT']]))
        elif '*' in alts and len(alts) == 2:
            star = alts.index('*')
            alt = alts[star - 1]
            fields[header['POS']], fields[header['REF']], new_alt = map(str, get_minimal_representation(fields[header['POS']], fields[header['REF']], alt))
            new_alts = [new_alt]
            new_alts.insert(star, '*')
            fields[header['ALT']] = ','.join(new_alts)
        print >> g, '\t'.join(fields)

    f.close()
    g.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file; may be gzipped', required=True)
    parser.add_argument('--output', '-o', help='Output VCF file; may be gzipped')
    args = parser.parse_args()
    main(args)
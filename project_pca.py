#!/usr/bin/env bash

from pyhail import *
from pyspark import SparkContext
import argparse
import os

vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
pca_vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2_gnomad_pca.cp5.vds'
meta_path = 'file:///mnt/lustre/konradk/exac/super_meta.txt.bgz'
output_path ='file:///mnt/lustre/konradk/exac/exac_full_with_pca'


def main(args):
    sc = SparkContext(appName='Hail')
    hc = HailContext(sc)

    vds = hc.read(args.vds)
    pca_vds = hc.read(args.pca)

    pc_number = 10
    pc_range = range(1, pc_number + 1)

    pca_vds.annotate_samples_expr('va.calldata = gs.callStats()')
    vds.annotate_samples_table(args.meta, 'sample', root='sa.meta')
    vds.annotate_variants_vds(pca_vds, root='va.pca')
    vds.filter_variants_expr('!isMissing(va.pca)')
    vds.annotate_variants_expr('va.PCs = [%s]' % ', '.join(['va.PC%s' % x for x in pc_range]))
    vds.annotate_samples_expr('sa.PCs = gs.map(g => let p = va.calldata.AF in if (p == 0 || p == 1) [%s] else (g.gt - 2 * p) / sqrt(2 * p * (1 - p)) * va.PCs).sum()' % ', '.join(['0.0']*pc_number))

    vds.write(os.path.join(args.output, '.vds'))
    vds.export_samples(os.path.join(args.output, '.samples.txt.bgz'), 'sample = s.id, sa.meta.sex, sa.PCs.*')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vds', '--input', '-i', help='VDS to annotate', default=vds_path)
    parser.add_argument('--pca', '-p', help='VDS with PCA data', default=pca_vds_path)
    parser.add_argument('--meta', '-m', help='Metadata file', default=meta_path)
    parser.add_argument('--output', '-o', help='Output file prefix', default=output_path)
    args = parser.parse_args()
    main(args)
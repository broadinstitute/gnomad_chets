#!/usr/bin/env bash

from pyhail import *
from pyspark import SparkContext
import argparse
import os
from utils import *

vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
hardcalls_vds_path = 'file:///mnt/lustre/konradk/exac/hardcalls/v2/exacv2.hardcalls.splitmulti.qc.vds'
pca_vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2_gnomad_pca.cp5.vds'
meta_path = 'file:///mnt/lustre/konradk/exac/super_meta.txt.bgz'
output_path ='file:///mnt/lustre/konradk/exac/exac_full_with_pca'



sc = SparkContext(appName='Hail')
hc = HailContext(sc)

vds = hc.read(vds_path)
pca_vds = hc.read(pca_vds_path)
hardcalls_vds = hc.read(hardcalls_vds_path)

pc_number = 10
pc_range = range(1, pc_number + 1)

new_pca_vds = pca_vds.annotate_variants_expr('va.calldata = gs.callStats(v)')
new_vds = hardcalls_vds.annotate_samples_table(meta_path, 'sample', root='sa.meta')\
    .annotate_variants_vds(new_pca_vds, root='va.pca')\
    .filter_variants_expr('!isMissing(va.pca)')\
    .annotate_variants_expr('va.PCs = [%s]' % ', '.join(['va.pca.pca_loadings.PC%s' % x for x in pc_range]))\
    .annotate_samples_expr('sa.PCs = gs.map(g => let p = va.pca.calldata.AF[1] in if (p == 0 || p == 1) [%s] else (g.gt - 2 * p) / sqrt(2 * p * (1 - p)) * va.PCs).sum()' % ', '.join(['0.0']*pc_number))

new_vds.write(os.path.join(output_path, '.vds'))
new_vds.export_samples(os.path.join(output_path, '.samples.txt.bgz'), 'sample = s.id, sa.meta.sex, sa.PCs.*')


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('--vds', '--input', '-i', help='VDS to annotate', default=vds_path)
#     parser.add_argument('--pca', '-p', help='VDS with PCA data', default=pca_vds_path)
#     parser.add_argument('--meta', '-m', help='Metadata file', default=meta_path)
#     parser.add_argument('--output', '-o', help='Output file prefix', default=output_path)
#     args = parser.parse_args()
#     main(args)
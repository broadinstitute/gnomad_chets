#!/usr/bin/env bash

__author__ = 'konrad'

from pyspark import SparkContext
from utils import *

root = 'file:///mnt/lustre/konradk/exac'
# vds_path = '%s/tiny.vds' % root
vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
meta_path = '%s/super_meta.txt.bgz' % root

sc = SparkContext(appName='Hail_exac_sites_vcf_y')
hc = HailContext(sc, log='sites_y.log')

pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']

criterion_pops = [('sa.meta.population', x) for x in pops]

a_based_annotations = ['va.info.AC', 'va.info.AC_Adj', 'va.info.AF', 'va.info.AF_Adj']
a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                            'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
a_based_annotations.extend(['va.info.AC_%s' % x for x in pops])
a_based_annotations.extend(['va.info.AF_%s' % x for x in pops])
r_based_annotations = ['va.info.GQ_HIST', 'va.info.DP_HIST']

# Dividing AC and AN by 2 on the y chromosome
correct_ac_an_command = ['va.info.AC = va.info.AC.map(x => (x/2).toInt), '
                         'va.info.AN = (va.info.AN/2).toInt']
for pop in pops + ['Adj']:
    correct_ac_an_command.append('va.info.AC_%(pop)s = va.info.AC_%(pop)s.map(x => (x/2).toInt), '
                                 'va.info.AN_%(pop)s = (va.info.AN_%(pop)s/2).toInt' % {'pop': pop})

correct_ac_an_command = ',\n'.join(correct_ac_an_command)

vds = hc.read(vds_path)

sites_vds = (vds.annotate_global_expr('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
             .annotate_samples_table('%s/super_meta.txt.bgz' % root, 'sample', impute=True, root='sa.meta')
             .filter_variants_intervals('%s/intervals/chrY.txt' % root)
             .filter_variants_expr('v.inYNonPar')
             .annotate_variants_loci('%s/sites/exac2.RF.allsites.txt.bgz' % root, 'Locus(chrom, pos)',
                                     code='va.info.RF_PROB = table.rfprob, va.info.RF_TYPE = table.type', impute=True)
             .filter_samples_expr('sa.meta.drop_status == "keep" && sa.meta.sex == "male"')
             .filter_genotypes('g.isHet', keep=False)
             .annotate_variants_expr('va.calldata.raw = gs.callStats(v)')
             .filter_alleles('va.calldata.raw.AC[aIndex] == 0', keep=False)  # change if default is no longer subset
             .konrad_special('va.info.GQ_HIST',
                             'gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.gq).hist(0, 100, 20)')
             .konrad_special('va.info.DP_HIST',
                             'gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.dp).hist(0, 100, 20)')
             .annotate_variants_expr('va.calldata.raw = gs.callStats(v)')
             .filter_to_adj()
             .projectmax()
             .annotate_variants_expr('va.calldata.Adj = gs.callStats(v)')
             .unfurl_callstats(criterion_pops, lower=True, gc=False)
             .filter_samples_all()
             .annotate_variants_expr('va.info.AC = va.calldata.raw.AC[1:], '
                                     'va.info.AN = va.calldata.raw.AN, '
                                     'va.info.AF = va.calldata.raw.AF[1:], '
                                     'va.info.AC_Adj = va.calldata.Adj.AC[1:], '
                                     'va.info.AN_Adj = va.calldata.Adj.AN, '
                                     'va.info.AF_Adj = va.calldata.Adj.AF[1:]')
             .annotate_variants_expr(correct_ac_an_command)
             .filter_star(a_based_annotations, r_based_annotations)
             .popmax(pops)
             .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
)

sites_vds.write('%s/sites/exac_v2_y_sites.vds' % root, overwrite=True)
sites_vds.export_vcf('%s/sites/exac_v2_y_internal.vcf.bgz' % root)

(sites_vds.annotate_variants_expr('va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
 .export_vcf('%s/sites/exac_v2_y_release.vcf.bgz' % root)
)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('--vds', '--input', '-i', help='VDS to annotate', default=vds_path)
#     parser.add_argument('--pca', '-p', help='VDS with PCA data', default=pca_vds_path)
#     parser.add_argument('--meta', '-m', help='Metadata file', default=meta_path)
#     parser.add_argument('--output', '-o', help='Output file prefix', default=output_path)
#     args = parser.parse_args()
#     main(args)
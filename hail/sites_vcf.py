#!/usr/bin/env bash

__author__ = 'konrad'

from pyspark import SparkContext
from utils import *

root = 'file:///mnt/lustre/konradk/exac'
# vds_path = '%s/tiny.vds' % root
vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
meta_path = '%s/super_meta.txt.bgz' % root

hc = HailContext(log='site_auto.log')

pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
sexes = ['Male', 'Female']

criterion_pops = [('sa.meta.population', x) for x in pops]
criterion_pops.extend([('sa.meta.sex', x) for x in sexes])

cuts = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS', 'Male', 'Female']

g_based_annotations = ['va.info.GC', 'va.info.GC_Adj']
g_based_annotations.extend(['va.info.GC_%s' % x for x in cuts])
a_based_annotations = ['va.info.AC', 'va.info.AC_Adj', 'va.info.AF', 'va.info.AF_Adj']
a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                            'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
a_based_annotations.extend(['va.info.AC_%s' % x for x in cuts])
a_based_annotations.extend(['va.info.AF_%s' % x for x in cuts])
r_based_annotations = ['va.info.GQ_HIST', 'va.info.DP_HIST', 'va.info.AB_HIST']

vds = hc.read(vds_path)

sites_vds = (vds.annotate_global_expr('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
             .annotate_samples_table('%s/super_meta.txt.bgz' % root, 'sample', impute=True, root='sa.meta')
             .filter_variants_intervals('%s/intervals/autosomes.txt' % root)
             .annotate_variants_loci('%s/sites/exac2.RF.allsites.txt.bgz' % root, 'Locus(chrom, pos)',
                                     code='va.info.RF_PROB = table.rfprob, va.info.RF_TYPE = table.type', impute=True)
             .filter_samples_expr('sa.meta.drop_status == "keep"')
             .annotate_variants_expr('va.calldata.raw = gs.callStats(v)')
             .filter_alleles('va.calldata.raw.AC[aIndex] == 0', keep=False)  # change if default is no longer subset
             .konrad_special('va.info.GQ_HIST',
                             'gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.gq).hist(0, 100, 20)')
             .konrad_special('va.info.DP_HIST',
                             'gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.dp).hist(0, 100, 20)')
             .annotate_variants_expr('va.calldata.raw = gs.callStats(v)')
             .filter_to_adj()
             .konrad_special('va.info.AB_HIST',
                             'gs.filter(g => (g.gtj == %s || g.gtk == %s) && g.isHet).map(g => 100*g.ad[%s]/g.dp).hist(0, 100, 20)')
             .projectmax()
             .annotate_variants_expr('va.calldata.Adj = gs.callStats(v)')
             .unfurl_callstats(criterion_pops, lower=True)
             .filter_samples_all()
             .annotate_variants_expr('va.info.AC = va.calldata.raw.AC[1:], '
                                     'va.info.AN = va.calldata.raw.AN, '
                                     'va.info.AF = va.calldata.raw.AF[1:], '
                                     'va.info.GC = va.calldata.raw.GC, '
                                     'va.info.AC_Adj = va.calldata.Adj.AC[1:], '
                                     'va.info.AN_Adj = va.calldata.Adj.AN, '
                                     'va.info.AF_Adj = va.calldata.Adj.AF[1:], '
                                     'va.info.GC_Adj = va.calldata.Adj.GC')
             .filter_star(a_based_annotations, r_based_annotations, g_based_annotations)
             .unfurl_hom(cuts, simple_hom=True)
             .popmax(pops)
             .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
             .repartition(1000)
)

print 'Writing VDS...'
sites_vds.write('%s/sites/exac_v2_sites.vds' % root)
print 'Done! Writing internal VCF...'
sites_vds.export_vcf('%s/sites/exac_v2_internal.vcf.bgz' % root)
print 'Done! Writing release VCF...'
(sites_vds.annotate_variants_expr('va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
 .export_vcf('%s/sites/exac_v2_release.vcf.bgz' % root)
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
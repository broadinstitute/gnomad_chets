#!/usr/bin/env bash

__author__ = 'konrad'

from pyspark import SparkContext
from utils import *

root = 'file:///mnt/lustre/konradk/exac'
# vds_path = '%s/tiny.vds' % root
vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
meta_path = '%s/super_meta.txt.bgz' % root

sc = SparkContext(appName='Hail_exac_sites_vcf_x')
hc = HailContext(sc, log='site_x.log')

pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
sexes = ['Male', 'Female']

cuts = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS', 'Male', 'Female']

g_based_annotations = ['va.info.GC', 'va.info.GC_Adj']
g_based_annotations.extend(['va.info.GC_%s' % x for x in sexes])
g_based_annotations.extend(['va.info.GC_%s_%s' % (y, x) for x in sexes for y in pops])

a_based_annotations = ['va.info.AC', 'va.info.AC_Adj', 'va.info.AF', 'va.info.AF_Adj']
a_based_annotations.extend(['va.info.AC_Male', 'va.info.AC_Female', 'va.info.AF_Male', 'va.info.AF_Female'])
a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                            'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
a_based_annotations.extend(['va.info.AC_%s_%s' % (y, x) for x in sexes for y in pops])
a_based_annotations.extend(['va.info.AF_%s_%s' % (y, x) for x in sexes for y in pops])

r_based_annotations = ['va.info.GQ_HIST', 'va.info.DP_HIST', 'va.info.AB_HIST']

generate_callstats_expression = []
for sex in sexes:
    for pop in pops:
        input_dict = {'sex': sex.lower(), 'sex_label': sex.capitalize(), 'pop': pop.lower(), 'pop_upper': pop.upper()}
        generate_callstats_expression.append('va.calldata.%(pop_upper)s_%(sex_label)s = gs.filter(g => sa.meta.population == "%(pop)s" && sa.meta.sex == "%(sex)s").callStats(v)' % input_dict)
    input_dict = {'sex': sex.lower(), 'sex_label': sex.capitalize()}
    generate_callstats_expression.append('va.calldata.%(sex_label)s = gs.filter(g => sa.meta.sex == "%(sex)s").callStats(v)' % input_dict)
generate_callstats_expression = ',\n'.join(generate_callstats_expression)

rearrange_callstats_expression = []
for metric in ['AC', 'AN', 'AF', 'GC']:
    for sex in sexes:
        start = ''
        end = '[1:]' if metric in ('AC', 'AF') else ''
        if sex == 'male' and metric in ('AC', 'AN'):
            start = 'let div_factor = if (v.inXNonPar) 2 else 1 in ('
            end += '/div_factor)'
            end += '.map(x => x.toInt)' if metric == 'AC' else '.toInt'
        for pop in pops:
            input_dict = {'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric, 'pop': pop, 'pop_upper': pop.upper(), 'start': start, 'end': end}
            rearrange_callstats_expression.append('va.info.%(metric)s_%(pop_upper)s_%(sex_label)s = %(start)sva.calldata.%(pop_upper)s_%(sex_label)s.%(metric)s%(end)s' % input_dict)
        input_dict = {'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric, 'start': start, 'end': end}
        rearrange_callstats_expression.append('va.info.%(metric)s_%(sex_label)s = %(start)sva.calldata.%(sex_label)s.%(metric)s%(end)s' % input_dict)

rearrange_callstats_expression.extend(['va.info.GC = va.calldata.raw.GC',
                                       'va.info.AC = (va.calldata.raw.AC - (va.calldata.hemi_raw.AC/2)).map(x => x.toInt)[1:]',
                                       'va.info.AN = (va.calldata.raw.AN - (va.calldata.hemi_raw.AN/2)).toInt',
                                       'va.info.GC_Adj = va.calldata.Adj.GC',
                                       'va.info.AC_Adj = (va.calldata.Adj.AC - (va.calldata.Hemi_Adj.AC/2)).map(x => x.toInt)[1:]',
                                       'va.info.AN_Adj = (va.calldata.Adj.AN - (va.calldata.Hemi_Adj.AN/2)).toInt'])
rearrange_callstats_expression = ',\n'.join(rearrange_callstats_expression)

ac_an_expression = []
for metric in ['AC', 'AN']:
    for pop in pops:
        input_dict = {'metric': metric, 'pop': pop, 'pop_upper': pop.upper()}
        ac_an_expression.append('va.info.%(metric)s_%(pop_upper)s = va.info.%(metric)s_%(pop_upper)s_Male + va.info.%(metric)s_%(pop_upper)s_Female' % input_dict)
ac_an_expression = ',\n'.join(ac_an_expression)

af_expression = []
for pop in pops:
    input_dict = {'pop': pop, 'pop_upper': pop.upper()}
    af_expression.append('va.info.AF_%(pop_upper)s = va.info.AC_%(pop_upper)s / va.info.AN_%(pop_upper)s' % input_dict)
af_expression = ',\n'.join(af_expression )

hom_hemi_expression = []
for sex in sexes:
    metric = 'Hom' if sex == 'female' else 'Hemi'
    for pop in pops:
        input_dict = {'pop': pop, 'pop_upper': pop.upper(), 'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric}
        hom_hemi_expression.append('va.info.%(metric)s_%(pop_upper)s = if (v.inXNonPar) range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(pop_upper)s_%(sex_label)s[(n * (n + 1) / 2).toInt - 1]) else NA: Array[Int]' % input_dict)
    input_dict = {'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric}
    hom_hemi_expression.append('va.info.%(metric)s = if (v.inXNonPar) range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(sex_label)s[(n * (n + 1) / 2).toInt - 1]) else NA: Array[Int]' % input_dict)
hom_hemi_expression = ',\n'.join(hom_hemi_expression)

vds = hc.read(vds_path)

sites_vds = (vds.annotate_global_expr('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
             .annotate_samples_table('%s/super_meta.txt.bgz' % root, 'sample', impute=True, root='sa.meta')
             .filter_variants_intervals('%s/intervals/chrX.txt' % root)
             .annotate_variants_loci('%s/sites/exac2.RF.allsites.txt.bgz' % root, 'Locus(chrom, pos)',
                                     code='va.info.RF_PROB = table.rfprob, va.info.RF_TYPE = table.type', impute=True)
             .filter_samples_expr('sa.meta.drop_status == "keep"')
             .filter_genotypes('sa.meta.sex == "male" && g.isHet && v.inXNonPar', keep=False)
             .annotate_variants_expr('va.calldata.raw = gs.callStats(v)')
             .filter_alleles('va.calldata.raw.AC[aIndex] == 0', keep=False)  # change if default is no longer subset
             .konrad_special('va.info.GQ_HIST',
                             'gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.gq).hist(0, 100, 20)')
             .konrad_special('va.info.DP_HIST',
                             'gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.dp).hist(0, 100, 20)')
             .annotate_variants_expr('va.calldata.raw = gs.callStats(v), '
                                     'va.calldata.hemi_raw = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(v)')
             .filter_to_adj()
             .konrad_special('va.info.AB_HIST',
                             'gs.filter(g => (g.gtj == %s || g.gtk == %s) && g.isHet).map(g => 100*g.ad[%s]/g.dp).hist(0, 100, 20)')
             .annotate_variants_expr('va.calldata.Adj = gs.callStats(v), '
                                     'va.calldata.Hemi_Adj = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(v)')
             .projectmax()
             .annotate_variants_expr(generate_callstats_expression)
             .filter_samples_all()
             .annotate_variants_expr(rearrange_callstats_expression)
             .annotate_variants_expr('va.info.AF = va.info.AC.map(x => if (va.info.AN > 0) x.toDouble/va.info.AN else NA: Double), '
                                     'va.info.AF_Adj = va.info.AC_Adj.map(x => if (va.info.AN_Adj > 0) x.toDouble/va.info.AN_Adj else NA: Double)')  # Got here
             .filter_star(a_based_annotations, r_based_annotations, g_based_annotations)
             .annotate_variants_expr(hom_hemi_expression)
             .annotate_variants_expr(ac_an_expression)
             .annotate_variants_expr(af_expression)
             .popmax(pops)
             .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
)

print 'Writing VDS...'
sites_vds.write('%s/sites/exac_v2_x_sites.vds' % root)
print 'Done! Writing internal VCF...'
sites_vds.export_vcf('%s/sites/exac_v2_x_internal.vcf.bgz' % root)
print 'Done! Writing release VCF...'
(sites_vds.annotate_variants_expr('va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
 .export_vcf('%s/sites/exac_v2_x_release.vcf.bgz' % root)
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
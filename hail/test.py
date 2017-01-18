
from utils import *

#Inputs

vds_path = "gs://gnomad/gnomad.10ksites.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats7.vds"
vep_config = "/vep/vep-gcloud.properties"

#Resources
lcr_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
decoy_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
autosomes_intervals = "gs://gnomad/autosomes.intervals"
dbsnp = "gs://gnomad-lfran/All_20160601.vcf.bgz"

#Outputs
out_root = "gs://gnomad-lfran/tmp"
sites_annotations_path = "%s/gnomad.sites.annotations.vds" % out_root
tmp_vds = "%s/test.rf-ann6.vds" % out_root
tmp_RF_ann_out = '%s/gnomad.rf.ann.txt.bgz' % out_root
tmp_RF_ann_exp = '%s/test.rf.ann.exp.txt.bgz' % out_root
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats9.vds"



hc = hail.HailContext(log='/site_auto.log')

pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']

rf_features = ['va.variantType',
                'va.info.QD',
                'va.info.MQ',
                'va.info.MQRankSum',
                'va.info.FS',
                'va.info.SOR',
                'va.info.InbreedingCoeff',
                'va.info.ReadPosRankSum',
                'va.stats.raw.nrq_median',
                'va.stats.raw.ab_median',
                'va.stats.raw.dp_median',
                'va.stats.raw.gq_median']

rf = hc.read(rf_path,sites_only=True)
rf.annotate_global_
rf.export_variants('gs://gnomad-lfran/')



global_expr = ['global.%s = variants.filter(x => isMissing(%s)).count()' % (a,a) for a in rf_features ]

rf_kt =rf.variants_keytable()

print(rf_kt.schema())

rf_kt = rf_kt.flatten().select(['va.info.MQRankSum'])

print(rf_kt.schema())

rf_df = rf_kt.to_dataframe()
rf_df.printSchema()
rf_df.show()

MQmedian = rf_df.approxQuantile("`va.info.MQRankSum`", [0.5], 0.05)

print(MQmedian)


# (
#     hc.read(rf_path,sites_only=True)
#     .annotate_global_expr_by_variant(global_expr)
#     .show_globals()
#
# )


# (hc.read("gs://gnomad/gnomad.raw_hardcalls.vds",sites_only=True)
#  .filter_variants_expr('v.altAlleles.exists(a => a.isSNP) && pcoin(0.92)',
#                        keep=False)
#  .filter_variants_expr('v.altAlleles.exists(a => !a.isSNP) && pcoin(0.6)',
#                        keep=False)
#  .export_variants(output= "gs://gnomad-lfran/tmp/variantTypeCheck.txt.bgz",
#                     condition='v = v, alts = v.altAlleles.map(a => a.alt).mkString(","), variantType = va.variantType')
#  )

#intervals = "gs://gnomad-lfran/tmp/maryam_variants.txt"

#vds = set_vcf_filters(hc, "gs://gnomad-lfran/tmp/maryam_variants.txt", rf_path, rf_ann_root, rf_snv_cutoff, rf_indel_cutoff, filters = {}, filters_to_keep = [], tmp_path = '/tmp')

# x=  (
#     hc.read(vds_path)
#         .annotate_global_expr_by_sample('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
#         .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=pyhail.TextTableConfig(impute=True))
#         .annotate_samples_expr(['sa.meta.population = sa.meta.predicted_pop',
#                                 'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
#         .filter_samples_expr('!isMissing(sa.meta.predicted_pop)')  # Could be cleaner
#         .filter_variants_intervals('gs://gnomad-lfran/tmp/test.interval')
#         .histograms('va.info')
#         .export_samples(output='file:///tmp/out_sampes.txt', condition='s.id,sa.meta.project_description')
#         .export_variants(output='file:///tmp/out.txt',
#                          condition='gq_hist_all = va.info.GQ_HIST_ALL,'
#                                    'gq_hist_alt = va.info.GQ_HIST_ALT,'
#                                    'dp_hist_all = va.info.DP_HIST_ALL,'
#                                    'dp_hist_alt = va.info.DP_HIST_ALT,'
#                                    'ab_hist_all = va.info.AB_HIST_ALL,'
#                                    'ab_hist_alt = va.info.AB_HIST_ALT'
#                          )
# )
# print(x.print_schema())

# res = annotate_non_split_from_split(hc, non_split_vds_path=vds_path,
#                               split_vds=hc.read(rf_path),
#                               annotations=['va.RF1'],
#                               annotation_exp_out_path=tmp_RF_ann_out)
#
#
#
# res.write(tmp_vds)
# rf = hc.read(rf_path)
# (
#     hc.read(tmp_vds)
#     .annotate_variants_expr('va.alts = v.altAlleles.map(a => a.alt)')
#     .split_multi()
#     .annotate_variants_vds(rf,code='va.rf_ann = vds.RF1')
#     .annotate_variants_expr('va.issame = va.rf_ann == va.RF1[va.aIndex-1]')
#     .export_variants(tmp_RF_ann_exp, 'v=v,va.alts=va.alts, ai1 = va.aIndex, va.rf_ann.prob = va.rf_ann.probability["TP"], va.RF1 =  va.RF1,va.issame=va.issame')
# )

#.annotate_variants_expr('va.rf_match = va.rf_ann.probability["TP"] == va.RF1[va.aIndex - 1].probability["TP"]')
#va.RF1.prob = va.RF1[va.aIndex - 1].probability["TP"],
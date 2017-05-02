from utils import *
from variantqc import *
from pprint import *

hc = hail.HailContext(log='/test.log')

(
    hc.read(full_genome_vds)
    .annotate_samples_table(genomes_meta, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
    .min_rep()
    .write('gs://gnomad-genomes-raw/full/gnomad.genomes.all.vds')

)


# rel = hc.read(final_genome_vds).filter_variants_intervals(Interval.parse('22'))
#
# nonpsych = hc.read("gs://gnomad-genomes/subsets/nonpsych/nonpsych.vds").filter_variants_intervals(Interval.parse('22'))
#
# both = nonpsych.annotate_variants_vds(rel, code = 'va.rel.RF = vds.info.AS_RF')
#
# print(both.query_variants('variants.filter(v => isDefined(va.rel.RF)).map(v => va.rel.RF == va.info.AS_RF).counter()'))







#print(vds.query_variants(['variants.filter(x => isMissing(%s)).count()' % f for f in rf_features]))


# allele_annotations = get_stats_expr("va.stats.qc_samples_raw", medians=True, samples_filter_expr='sa.meta.qc_sample')
# #allele_annotations.extend(
# #    get_stats_expr("va.stats.release_samples_raw", medians=True, samples_filter_expr='sa.meta.keep'))
#
# print(allele_annotations)
#
# (
#     hc.read("gs://gnomad/gnomad.raw_hardcalls.old2.vds")
#     .annotate_alleles_expr(allele_annotations)
#     .write(raw_hardcalls_path, overwrite=True)
# )
#
# hapmap = hc.read(hapmap_path)
# mills = hc.read(mills_path)
# omni = hc.read(omni_path)
#
# hardcalls_split = (
#     hc.read(raw_hardcalls_path)
#         .annotate_variants_expr('va.nonsplit_alleles = v.altAlleles.map(a => a.alt)')
#         .split_multi()
#         .annotate_variants_vds(hapmap, code='va.hapmap = isDefined(vds)')
#         .annotate_variants_vds(omni, code='va.omni = isDefined(vds)')
#         .annotate_variants_vds(mills, code='va.mills = isDefined(vds)')
#         .tdt(fam=fam_path)
#         .write(raw_hardcalls_split_path, overwrite=True)
# )

# vds = hc.read(rf_path).annotate_variants_table(mendel_path + ".lmendel",'SNP',code='va.mendel = table.N',config=hail.TextTableConfig(impute=True))
# fps = vds.query_variants(['variants.filter(v => va.tdt.pval < 0.001).count()',
#                           'variants.filter(v => va.mendel > 5).count()',
#                           'variants.filter(v => va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30).count()',
#                           'variants.filter(v => va.omni || va.mills || (va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1] == 2)).count()'])
#
# print(fps)

# (
#     hc.read(rf_path)
#     .annotate_variants_table(mendel_path + ".lmendel",'SNP',code='va.mendel = table.N',config=hail.TextTableConfig(impute=True))
#     .annotate_global_expr_by_variant('global.nBadTDT = variants.filter(v => va.tdt.pval < 0.001).count(), '
#                                      'global.nDoubleMendel = variants.filter(v => va.mendel > 5).count(),'
#                                      'global.nHF = variants.filter(v => va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30).count()')
#     .show_globals()
# )

## Show high confidence FP indels in mills / tx singletons
# (
#     hc.read("gs://gnomad/RF/gnomad.sites.RF.newStats12.vds")
#     .filter_variants_list("gs://gnomad/lfran/tmp/high_conf_fps.txt")
#     .annotate_global_expr_by_variant('global.fps_mills = variants.filter(v => va.mills).count(),'
#                                       'global.fps_tx_singletons = variants.filter(v => va.tdt.nTransmitted == 1 && '
#                                       'va.calldata.raw.AC[va.aIndex]==2).count()')
#     .show_globals()
# )

# Table type / wasSplit
# (
#     hc.read("gs://gnomad/gnomad.raw_hardcalls.vds",sites_only=True)
#     .annotate_variants_expr('va.nonsplit_alleles = v.altAlleles.map(a => a.alt)')
#     .split_multi()
#     .variants_keytable()
#     .aggregate_by_key(key_condition='type = va.variantType, '
#                                     'wasSplit = va.wasSplit, '
#                                     'hasStar = va.nonsplit_alleles.exists(a => a == "*")',
#                       agg_condition='n = va.count()')
#     .to_dataframe()
#     .show()
# )
#
# pprint(
#     hc.read("gs://gnomad/gnomad.raw_hardcalls.split.vds")
#     .filter_variants_expr('va.variantType == "mixed" && !va.wasSplit')
#     .head()
# )

send_message(channel='@laurent', message='Test2 is done processing!')

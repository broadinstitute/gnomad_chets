from utils import *
from variantqc import *
from pprint import *

hc = hail.HailContext(log='/test.log')

rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats12.vds"
mendel_path = "gs://gnomad/gnomad.raw_calls"
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"
raw_hardcalls_split_path = "gs://gnomad/gnomad.raw_hardcalls.split.vds"
fam_path = "gs://gnomad/gnomad.final.goodTrios.fam"
rf_ann_path = "gs://gnomad/RF/gnomad.sites.annotated_for_RF.vds"

for contig in ["X","autosomes"]:

    print("Computing %s" % contig)

    exomes = hc.read("gs://gnomad-exomes/sites/vds/gnomad.exomes.sites.release.%s.vds" % contig, sites_only=True).split_multi()
    exomes = exomes.annotate_samples_expr('sa = {}')
    genomes = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.%s.vds" % contig, sites_only=True).split_multi()
    genomes = genomes.annotate_samples_expr('sa = {}')

    both = genomes.join(exomes)
    both = both.annotate_variants_vds(exomes, 'va.exomes.pass = vds.filters.forall(f => f == "AC0") && vds.info.AS_FilterStatus[vds.aIndex - 1].forall(f => f == "AC0")')

    (
        both
            .variants_keytable().aggregate_by_key('genome_pass = va.filters.forall(f => f == "AC0") && va.info.AS_FilterStatus[va.aIndex - 1].forall(f => f == "AC0"),'
                                        'exome_pass = va.exomes.pass',
                          'n = va.count()')
        .to_dataframe()
        .show()
    )





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


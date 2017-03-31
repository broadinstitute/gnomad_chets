
from utils import *
import hail
from hail.representation import Interval

#Inputs
vds_path = "gs://gnomad/gnomad.10ksites.vds"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats24.vds"
vep_config = "/vep/vep-gcloud.properties"
mendel_path = "gs://gnomad/gnomad.raw_calls"
fam_path = "gs://gnomad/gnomad.final.goodTrios.fam"
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"

genome_pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']

hc = hail.HailContext(log='/test.log')

#print(hc.read('gs://gnomad-genomes/subsets/odonovan_quads/odonovan_quads.vds')
print(hc.read(final_genome_vds)
      .filter_variants_intervals(Interval.parse('5:11586-21586'))
      .query_variants(['variants.count()','variants.map(v => va.info.AC_ASJ_Male).collect()[1:10]']))


# exomes = hc.read(final_exome_vds)
# sanity_check_text = run_sanity_checks(exomes.filter_variants_intervals(IntervalTree.parse_all(['1-22'])),pops = POPS, return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'exomes merged file -- autosomes sanity check')
# sanity_check_text = run_sanity_checks(exomes.filter_variants_intervals(Interval.parse("X")),pops = POPS, contig = 'X', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'exomes merged file -- X sanity check')
# sanity_check_text = run_sanity_checks(exomes.filter_variants_intervals(Interval.parse("Y")),pops = POPS, contig='Y', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'exomes merged file -- Y sanity check')

# genomes = hc.read(final_genome_vds)
# sanity_check_text = run_sanity_checks(genomes.filter_variants_intervals(IntervalTree.parse_all(['1-22'])),pops = genome_pops, return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'genomes merged file -- autosomes sanity check')
# sanity_check_text = run_sanity_checks(genomes.filter_variants_intervals(Interval.parse("X")),pops = genome_pops, contig = 'X', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'genomes merged file -- X sanity check')
# sanity_check_text = run_sanity_checks(genomes.filter_variants_intervals(Interval.parse("Y")),pops = genome_pops, contig='Y', return_string=True, skip_star=True)
# send_snippet('#joint_calling', sanity_check_text, 'genomes merged file -- Y sanity check')

# exomes_vdses = [
# hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.autosomes.vds"),
# hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.X.vds"),
# hc.read("gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.Y.vds")
# ]
# exomes_vdses = merge_schemas(exomes_vdses)
# exomes_vdses[0].union(exomes_vdses[1:]).write(final_exome_vds)
#
# genomes_vdses = [
# hc.read("gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.autosomes.vds"),
# hc.read("gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.X.vds")
# ]
# genomes_vdses = merge_schemas(genomes_vdses)
# genomes_vdses[0].union(genomes_vdses[1:]).write(final_genome_vds)


#
# print(vds.variant_schema)
# print(vds.sample_schema)
# print(vds.query_samples(['samples.count()']))
# print(vds.query_variants(['variants.map(v => va.calldata.all_samples_raw.AN).stats().max',
#                           'variants.map(v => va.calldata.combined.AN).stats().max']))


def subset_tj(hc):
    genomes = False

    if(genomes):
        sample_lists = {'concordance': genomes_concordance_samples}
        projects = '"G68758","G87944","G77318","G29747","G77318","G77525","G29749","G87944","G87944","G31561","G77318","G89387","G87944","G84381","G29747","G89387","G31561","G29748","G29748","G26842"'
        projects_expr = '[%s].toSet.contains(sa.project_or_cohort)' % projects
        path = full_genome_vds
        out = "genomes"
    else:
        sample_lists = {'sa.concordance': exomes_concordance_samples,
                        'sa.hapmap': 'gs://gnomad-exomes-raw/gnomad_exomes_hapmap_samples_ids.txt'}
        projects = '"C1629","FINNMETSEQ","FINNMETSEQ","C1975","C1476","UK10K","C1975","UK10K","C1476","C1975","UK10K","UK10K","C1753","C1975","C1622","FINNMETSEQ","C1629","C1975","C1629","C1975"'
        projects_expr = '[%s].toSet.contains(sa.meta.pid) || sa.hapmap' % projects
        path = full_exome_vds
        out = "exomes"

    vds = hc.read(path)

    for root,sample_list in sample_lists.iteritems():
        vds = vds.annotate_samples_list(sample_list, root)

    filter_samples_expr = " || ".join(sample_lists.keys())

    vds = (
        vds
            .filter_samples_expr(filter_samples_expr)
            .filter_samples_expr(projects_expr)
            .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
            .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
            .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
    )
    print(vds.query_samples(['samples.count()']))
    vds.export_vcf("gs://gnomad-genomes/subsets/tj/duplicate_samples.%s.vcf.bgz" % out)



# paths = ["gs://gnomad/gnom.ad.vds","gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds"]
#
# for path in paths:
#     print("Running: %s\n\n" % path)
#     vds = (
#         hc.read(path,sites_only=True)
#             .annotate_variants_expr('va.v = v')
#             .split_multi()
#
#     )
#     print(vds.query_variants(['variants.map(v => v.start - va.v.start).counter()']))
#

#vds = hc.read(gnomad_path, sites_only=True)
# vds = hc.import_vcf("gs://gnomad-lfran/tmp/1var.vcf")
# #vds = vds.filter_variants_intervals(Interval.parse('1:1000000-1500000'))
# #vds = vds.filter_variants_intervals(Interval.parse('1:1087812-1087814'))
# #vds = vds.filter_variants_expr('v.altAlleles.forall(a => a.alt != "*")')
# vds = vds.vep(config=vep_config, csq=False, root='va.info.CSQ', force=True)
# #vds = vds.vep(config=vep_config, csq=False, root='va.info.CSQ', force=True)
# print(vds.variant_schema)
# print(vds.query_variants('variants.map(v => str(va.info.CSQ)).collect()'))
#
# vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)
# #vds = vds.vep(config=vep_config, csq=False, root='va.info.CSQ', force=True)
# print(vds.variant_schema)
# print(vds.query_variants('variants.map(v => str(va.info.CSQ)).collect()'))
# #vds.export_vcf("gs://gnomad-lfran/tmp/test-vep-attr.vcf")

# vds = hc.import_fasta("gs://gnomad-resources/Homo_sapiens_assembly19.fasta", filter_Ns=True, flanking_context=1, create_snv_alleles=True, create_deletion_size=4, create_insertion_size=2)
# print(vds.query_variants(['variants.count()',
#                           'variants.map(v => v.nAltAlleles).stats()',
#                           'variants.map(v => va.context).counter()']))


# vds = hc.read(gnomad_path)
# vds = filter_intervals(vds,['2:21229160-21229160','2:21280079-21280079'])
# vds = vds.annotate_samples_expr('sa.keepMe = gs.filter(g => v.start == 21280079 && g.isCalledNonRef).count() > 0')
# vds = vds.filter_samples_expr('sa.keepMe')
# vds.export_vcf('gs://gnomad-lfran/tmp/APOB_Lof.vcf')

# y_intervals = '/tmp/chrY.txt'
# with open(y_intervals, 'w') as f:
#     f.write('Y:1-1000000000')
#
#
#
#
# all = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.X.bad-multi.vds")
#
# extra_fields = {
#     "va.info.BaseQRankSum" : []
# }
#
# multi = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.X.vds")
# multi = copy_schema_attributes(all, multi)
# multi = multi.annotate_variants_expr('va.pass = va.filters.isEmpty')
#
# #print("\n----- ALL -----\n")
# #print_schema_attributes(all)
# #print("\n----- Multi -----\n")
# #print_schema_attributes(multi)
#
# res = all.annotate_variants_vds(multi,'va = if(isMissing(vds)) va else vds')
# res.write("gs://gnomad-genomes/sites/internal/gnomad.genomes.sites.X.vds")


#
#     print("All filters")
#     vds = hc.read(p)
#     (
#         vds
#             .split_multi()
#             .variants_keytable().aggregate_by_key(
#             key_condition='type = if(v.altAllele.isSNP) "snv" else if(v.altAllele.isIndel) "indel" else "other",'
#                           'filtered = !va.info.AS_FilterStatus[va.aIndex-1].isEmpty || !va.filters.isEmpty',
#             agg_condition='n = va.count()')
#             .to_dataframe()
#             .show()
#     )


#print(vds.query_variants("variants.map(v => v.altAllele.isSNP).counter()"))

#run_sanity_checks(hc.read("gs://gnomad-exomes/sites/gnomad.exomes.sites.autosomes.vds"),['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS'])
# vds = hc.read("gs://gnomad-genomes/sites/gnomad.genomes.sites.X.vds")
# print(vds.query_variants(['variants.count()',
#                           'variants.filter(v => v.inXNonPar).count()']))
# print(vds.query_variants(['variants.map(v => range(v.nAltAlleles)'
#                     '.map(i => (va.info.Hom_AFR[i] + va.info.Hom_AMR[i] + va.info.Hom_ASJ[i] + va.info.Hom_EAS[i] + '
#                     'va.info.Hom_FIN[i] + va.info.Hom_NFE[i] + va.info.Hom_OTH[i] - va.info.Hom[i]).abs).max).stats()',
#                     'variants.map(v => range(v.nAltAlleles)'
#                     '.map(i => (va.info.AC_AFR[i] + va.info.AC_AMR[i] + va.info.AC_ASJ[i] + va.info.AC_EAS[i] + '
#                     'va.info.AC_FIN[i] + va.info.AC_NFE[i] + va.info.AC_OTH[i] - va.info.AC[i]).abs).max).stats()',
#                     'variants.map(v => (va.info.AN_AFR + va.info.AN_AMR + va.info.AN_ASJ + va.info.AN_EAS '
#                     '+ va.info.AN_FIN + va.info.AN_NFE + va.info.AN_OTH - va.info.AN).abs).stats()'
#                    ]))

#print(hc.read(gnomad_path).query_variants('variants.map(v => v.contig).counter()'))
#print hc.read('gs://gnomad-lfran/tmp/gnomad.sites.tmp.2017-02-13_21-40.vds').query_variants(['variants.filter(v => isDefined(va.info.%s)).count()' % x for x in ('DS', 'END', 'MQ0', 'MQ', 'RAW_MQ')])

#vds = vds.filter_variants_intervals('file://' + y_intervals)
#vds = vds.filter_variants_intervals('gs://gnomad-lfran/tmp/1gene.intervals')
#print(vds.query_variants('variants.map(v => isDefined(va.info.DS)).counter()'))
# s =hc.read("gs://gnomad-lfran/tmp/gnomad.sites.tmp.2017-02-13_21-40.vds").variant_schema
# print(s)
# rf = [x for x in s.fields if x.name == "info"][0]
# print(rf)
#for f in rf.fields:




#pprint(hc.read(gnomad_path, sites_only=True).query_variants('variants.filter(v => pcoin(0.001)).map(v => v.start).stats()'))
# vds = hc.read(vds_path, sites_only=True)
# vds = vds.filter_variants_expr('pcoin(0.1)')
# vds.export_vcf('file:///test/test.before.vcf')
# vds = hc.import_vcf('/test.before.new.vcf')
# vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)
# print(vds.variant_schema)
# vds.export_vcf('file:///test/test.after.vcf')
#vds.write("gs://gnomad-lfran/tmp/vep.test.vds")

# vds = hc.import_vcf('/test.before.new.vcf')
# vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)
# vds.export_vcf('file:///test/test.after.vcf')

#vds = hc.read("gs://gnomad-lfran/tmp/vep.test.vds")
#pprint(vds.query_variants('variants.filter(v => pcoin(0.1)).map(v => va.info.CSQ).collect()'))

    #.export_vcf("gs://gnomad-lfran/tmp/vep.test.vcf")
#pprint(vds.query_variants('variants .filter(v => v.altAlleles.exists(a => a.alt.length >1032)).collect()'))
#pprint(vds.query_variants('variants.map(v => v.altAlleles.map(a => a.alt.length).max).stats().max'))

#mvs = vds.query_variants('variants.filter(v => va.mendel.length > 0).collect()')[0]

#print(mvs[:3])


#label_criteria = ['va.transmitted_singleton && v.altAllele.isIndel', 'va.mills']
#pprint(dict(zip(label_criteria, vds.query_variants(['variants.filter(x => %s).count()' % x for x in label_criteria]))))


#
# hc.read("gs://gnomad/gnomad.raw_hardcalls.tmp.vds").print_schema()
#
# pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
#
# rf_features = ['va.variantType',
#                 'va.info.QD',
#                 'va.info.MQ',
#                 'va.info.MQRankSum',
#                 'va.info.FS',
#                 'va.info.SOR',
#                 'va.info.InbreedingCoeff',
#                 'va.info.ReadPosRankSum',
#                 'va.stats.raw.nrq_median',
#                 'va.stats.raw.ab_median',
#                 'va.stats.raw.dp_median',
#                 'va.stats.raw.gq_median']
#
# rf = (hc.read(rf_path,sites_only=True)
#         .annotate_global_expr_by_variant('global.variantsByType = index(variants.map(v => va.variantType).counter(),key)')
#     )
# rf.show_globals()
# (
#     rf.filter_variants_expr('pcoin([1.0,1000000 / global.variantsByType[va.variantType].count].min)', keep=True)
#     .export_variants('gs://gnomad-lfran/gnomad.features_for_median.txt',
#                    condition='v=v,'
#                              'variantType=va.variantType,'
#                              'mqranksum=va.info.MQRankSum,'
#                              'readposranksum=va.info.ReadPosRankSum,'
#                              'ab_median=va.stats.raw.ab_median')
# )
#

# global_expr = ['global.%s = variants.filter(x => isMissing(%s)).count()' % (a,a) for a in rf_features ]
#
# rf_kt =rf.variants_keytable()
#
# print(rf_kt.schema())
#
# rf_kt = rf_kt.flatten().select(['va.info.MQRankSum'])
#
# print(rf_kt.schema())
#
# rf_df = rf_kt.to_dataframe()
# rf_df.printSchema()
# rf_df.show()
#
# MQmedian = rf_df.approxQuantile("`va.info.MQRankSum`", [0.5], 0.05)
#
# print(MQmedian)


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

send_message(channel='@laurent', message='Test is done processing!')
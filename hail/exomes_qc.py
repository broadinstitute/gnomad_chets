
from variantqc import *
import sys

# try:
#     hc
# except Exception, e:
#     hc = HailContext()

from hail import *
hc = HailContext()

root = 'gs://exac2'
meta_path = '%s/super_meta.txt.bgz' % root
autosome_intervals = '%s/intervals/autosomes.txt' % root
evaluation_intervals = '%s/intervals/exome_evaluation_regions.v1.intervals' % root
high_coverage_intervals = '%s/intervals/high_coverage.auto.interval_list' % root

raw_hardcall_vds_path = '%s/hardcalls/exacv2.raw.hardcalls.qc.vds' % root
raw_hardcalls_split_vds_path = '%s/hardcalls/exacv2.raw.hardcalls.splitmulti.qc.vds' % root

# raw_hardcall_vds = hc.read(raw_hardcall_vds_path)
# write_split(raw_hardcall_vds, raw_hardcalls_split_vds_path)
raw_hardcalls_split_vds = hc.read(raw_hardcalls_split_vds_path)

fam_path = '%s/variantqc/exac2.qctrios.fam' % root
tdt_vds_path = '%s/variantqc/v2_tdt.raw.vds' % root

# Use raw VDS for determing true positives
get_transmission_training_examples(raw_hardcalls_split_vds, tdt_vds_path, fam_path, autosome_intervals)

# omni_vds = hc.read('%s/1000G_omni2.5.b37.vds' % truth_dir)
# mills_vds = hc.read('%s/Mills_and_1000G_gold_standard.indels.b37.vds' % truth_dir)
# transmission_vds = hc.read(tdt_vds_path)

rf_variantqc_path = '%s/variantqc/exacv2_rf.vds' % root

# rf_vds = (annotate_for_random_forests(raw_hardcalls_split_vds, transmission_vds, omni_vds, mills_vds)
#           .random_forests(training='va.train', label='va.label', root='va.rf', features=rf_features, num_trees=500, max_depth=5))
# rf_vds.write(rf_variantqc_path)

# rf_vds = hc.read(rf_variantqc_path)
final_variantqc_path = '%s/variantqc/exacv2_variantqc.vds' % root

# TODO: Add 12878 and CHMI annotations to this command

# new_vds = (rf_vds.filter_variants_intervals(autosome_intervals)
#            .annotate_samples_fam(fam_path, root='sa.fam')
#            .annotate_variants_vds(hc.read('%s/variantqc/exac2_vqsr.vds' % root), code='va.info.VQSLOD = vds.info.VQSLOD')
#            .annotate_variants_intervals(evaluation_intervals, root='va.evaluation_interval')  # warning: this is not a boolean
#            .annotate_variants_intervals(high_coverage_intervals, root='va.high_coverage_interval')
#            .annotate_variants_table('%s/variantqc/v2.lmendel' % root, 'SNP', root='va.mendel', config=hail.TextTableConfig(impute=True))
#            .annotate_variants_table('%s/variantqc/validatedDN.cut.txt.bgz' % root, 'Variant(CHROM, POSITION.toInt, REF, ALT)', code='va.validated_denovo = table.DataSet', config=hail.TextTableConfig(impute=True))
#            .annotate_variants_expr('va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.oneHotAlleles(v)).sum()')
#            .filter_samples_expr('sa.meta.drop_status == "keep"')
#            .write(final_variantqc_path, overwrite=True)
# )
#
rf_vds = hc.read(final_variantqc_path)
#
# columns = '''chrom = v.contig,
# pos = v.start,
# ref = v.ref,
# alt = v.alt,
# evaluation_interval = !isMissing(va.evaluation_interval),
# high_coverage_interval = va.high_coverage_interval,
# vqslod = va.info.VQSLOD,
# type = va.variantType,
# rfprob = va.rf.probability["TP"],
# pass = va.pass,
# qd = va.info.QD,
# wassplit = va.wasSplit,
# mendel_errors = va.mendel.N,
# validated_denovo = va.validated_denovo,
# ac_unrelated = va.AC_unrelated[1].toInt,
# transmitted = va.tdt.nTransmitted,
# untransmitted = va.tdt.nUntransmitted,
# ac_orig = va.info.AC[va.aIndex - 1],
# ac_all_raw = va.calldata.all_samples_raw.AC[1],
# ac_qc_raw = va.calldata.qc_samples_raw.AC[1],
# an_qc_raw = va.calldata.release_samples_raw.AN,
# ac_release_raw = va.calldata.release_samples_raw.AC[1],
# nrq_median = va.stats.qc_samples_raw.nrq_median,
# gq_median = va.stats.qc_samples_raw.gq_median,
# dp_median = va.stats.qc_samples_raw.dp_median,
# ab_median = va.stats.qc_samples_raw.ab_median'''

# rf_vds.export_variants('%s/variantqc.txt.bgz' % root, columns)

truth_vds = (hc.read('%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vds' % truth_dir)
             .filter_variants_intervals(autosome_intervals)
)

(raw_split_hardcallvds.filter_variants_intervals(autosome_intervals)
 .filter_samples_expr('s.id == "C1975::NA12878"')
 .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')
 .concordance()
)

#     annotatevariants expr -c 'va.PASS = va.filters.contains("PASS")' \
#     concordance --right vds --samples s_concordance --variants v_concordance \
#     get -n s_concordance \
#     write -o $truthprefix.s_concordance.vds \
#     get -n v_concordance \
#     write -o $truthprefix.v_concordance.vds \
#     read -i $truthprefix.v_concordance.vds \
#     annotatevariants vds -r va.rf_new -i $root/variantqc/exacv2_rf.vds \
#     annotatevariants vds -r va.rf_all -i $root/sites/exac2.sites.RF2.vds \
#     annotatevariants vds -r va.rf_bal -i $root/sites/exac2.sites.RF3.vds \
#     annotatevariants intervals -r va.evaluation_interval -i $root/intervals/exome_evaluation_regions.v1.intervals \
#     annotatevariants intervals -r va.high_coverage_interval -i $root/intervals/high_coverage.auto.interval_list \
#     annotatevariants bed -r va.hc_12878 -i $truth_dir/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed \
#     annotateglobal expr -c 'global.gt_mappings = ["missing", "no_call", "homref", "het", "homvar"]' \
#     annotatevariants expr -c 'va.gt_arr = range(5).find(i => va.concordance[i].exists(x => x > 0))' \
#     annotatevariants expr -c '
#         va.truth_gt = global.gt_mappings[va.gt_arr],
#         va.called_gt = global.gt_mappings[range(5).find(i => va.concordance[va.gt_arr][i] > 0)]' \
#     exportvariants -c '
#         chrom = v.contig,
#         pos = v.start,
#         ref = v.ref,
#         alt = v.alt,
#         type = va.rf_new.variantType,
#         rfprob = va.rf_new.RF.probability["TP"],
#         rfprob_all = va.rf_all.RF.probability[2],
#         rfprob_bal = va.rf_bal.RF.probability[2],
#         mixed = va.right.isMixed,
#         evaluation_interval = va.evaluation_interval,
#         high_coverage_interval = va.high_coverage_interval,
#         hc_12878 = va.hc_12878,
#         wassplit_truth = va.left.wasSplit,
#         vqslod = va.right.info.VQSLOD,
#         pass = va.right.PASS,
#         wassplit = va.right.wasSplit,
#         truth_gt = va.truth_gt,
#         called_gt = va.called_gt,
#         concordance = va.concordance
#         ' \
#     -o ${truthprefix}_concordance.txt.bgz
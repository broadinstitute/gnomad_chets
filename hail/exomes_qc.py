
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

syndip_concordance_prefix = '%s/variantqc/syndip' % root
NA12878_concordance_prefix = '%s/variantqc/na12878' % root

concordance_annotations = ['chrom = v.contig',
                           'pos = v.start',
                           'ref = v.ref',
                           'alt = v.alt',
                           'concordance = va.concordance'
]

truth_concordance_annotations = list(concordance_annotations)
truth_concordance_annotations.extend(['type = va.rf.variantType',
                                      'wassplit = va.left.wasSplit',
                                      'vqslod = va.rf.info.VQSLOD',
                                      'truth_wassplit = va.right.wasSplit',
                                      'truth_gt = va.truth_gt',
                                      'called_gt = va.called_gt',
                                      'training = va.rf.train',
                                      'label = va.rf.label',
                                      'rfprob1 = va.rf.RF1.probability["TP"]'
                                      ])

compute_concordance(hc.read(raw_hardcalls_split_vds_path),
                    hc.read(syndip_path),
                    'CHMI_CHMI3_WGS1',
                    syndip_high_conf_regions_path,
                    syndip_concordance_prefix)

export_concordance(hc.read(syndip_concordance_prefix + ".v_concordance.vds"),
                   rf_vds,
                   syndip_concordance_prefix,
                   truth_concordance_annotations)

compute_concordance(hc.read(raw_hardcalls_split_vds_path),
                    hc.read(NA12878_path),
                    'G94982_NA12878',
                    NA12878_high_conf_exome_regions_path,
                    NA12878_concordance_prefix)

export_concordance(hc.read(NA12878_concordance_prefix + ".v_concordance.vds"),
                   rf_vds,
                   NA12878_concordance_prefix,
                   truth_concordance_annotations)

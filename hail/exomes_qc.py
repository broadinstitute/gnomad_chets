
from variantqc import *

root = 'file:///mnt/lustre/konradk/exac'
meta_path = '%s/super_meta.txt.bgz' % root
full_vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
full_v1_vds_path = 'file:///mnt/lustre/lfran/exac/exac_all.new.vds'

# All hardcalls
raw_hardcallvds_path = '%s/hardcalls/v2/exacv2.raw.hardcalls.qc.vds' % root
hardcallvds_path = '%s/hardcalls/v2/exacv2.hardcalls.qc.vds' % root
v1_hardcallvds_path = '%s/hardcalls/v1/exacv1.hardcalls.qc.vds' % root

# All splitmulti hardcalls
raw_split_hardcallvds_path = '%s/hardcalls/v2/exacv2.raw.hardcalls.splitmulti.qc.vds' % root
split_hardcallvds_path = '%s/hardcalls/v2/exacv2.hardcalls.splitmulti.qc.vds' % root
v1_split_hardcallvds_path = '%s/hardcalls/v1/exacv1.hardcalls.splitmulti.qc.vds' % root

full_vds = hc.read(full_vds_path)
full_v1_vds = hc.read(full_v1_vds_path)
# write_hardcalls(full_vds, raw_hardcallvds_path, meta_path, adj=False)
# write_hardcalls(full_vds, hardcallvds_path, meta_path)
# write_hardcalls(full_v1_vds, v1_hardcallvds_path, meta_path)
# write_split(raw_hardcallvds_path, raw_split_hardcallvds_path)
# write_split(hardcallvds_path, split_hardcallvds_path)
# write_split(v1_hardcallvds_path, v1_split_hardcallvds_path)

# TDT stuff
raw_fam_path = '%s/variantqc/exac2.qctrios.raw.fam' % root

# Run mendelerrors to remove bad trios
# (hc.read(split_hardcallvds_path)
#  .filter_variants_intervals('%s/intervals/autosomes.txt' % root)
#  .filter_variants_intervals('%s/intervals/exome_evaluation_regions.v1.intervals' % root)
#  .mendel_errors('%s/v2' % root, raw_fam_path))

fam_path = '%s/variantqc/exac2.qctrios.fam' % root

# Process exac2.qctrios.raw.fam into exac2.qctrios.fam by removing anyone in family with > 250 Mendelian errors
# Laurent integrated this into his fam script

tdt_vds_path = '%s/variantqc/v2_tdt.raw.vds' % root

# Use raw VDS for determing true positives
# get_transmitted_singletons(raw_split_hardcallvds_path, tdt_vds_path, fam_path, )

truth_dir = '%s/variantqc/truth_sets' % root
omni_vds = hc.read('%s/1000G_omni2.5.b37.splitmulti.vds' % truth_dir)
mills_vds = hc.read('%s/Mills_and_1000G_gold_standard.indels.b37.splitmulti.vds' % truth_dir)
transmission_vds = hc.read(tdt_vds_path)

rf_vds_path = '%s/variantqc/exacv2_rf.vds' % root
# run_random_forests(hardcallvds_path, rf_vds_path)


def temp():
    hardcallvds = hc.read(split_hardcallvds_path)
    rf_vds = hc.read(rf_vds_path)
    final_variantqc_path = '%s/variantqc/exacv2_variantqc.vds' % root

    new_vds = (hardcallvds.filter_variants_intervals('%s/intervals/autosomes.txt' % root)
               .annotate_samples_fam(fam_path, root='sa.fam')
               .filter_samples_expr('sa.meta.drop_status == "keep" || !isMissing(sa.fam.famID)')
               .annotate_variants_intervals('%s/intervals/exome_evaluation_regions.v1.intervals' % root, root='va.evaluation_interval')  # warning: this is not a boolean
               .annotate_variants_intervals('%s/intervals/high_coverage.auto.interval_list' % root, root='va.high_coverage_interval')
               .annotate_variants_table('%s/variantqc/v2.lmendel' % root, 'SNP', impute=True, root='va.mendel')
               .annotate_variants_table('%s/variantqc/validatedDN.cut.txt.bgz' % root, 'Variant(CHROM, POSITION.toInt, REF, ALT)', impute=True, code='va.validated_denovo = table.DataSet')
               .annotate_variants_vds(rf_vds, root='va.rf')
               .annotate_variants_expr('va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.oneHotAlleles(v)).sum(),'
                                       'va.pass = va.filters.contains("PASS")')
               .tdt(fam_path)
               .filter_samples_expr('sa.meta.drop_status == "keep"')
               .variant_qc()
               .write(final_variantqc_path)
    )

    columns = '''chrom = v.contig,
    pos = v.start,
    ref = v.ref,
    alt = v.alt,
    evaluation_interval = !isMissing(va.evaluation_interval),
    high_coverage_interval = va.high_coverage_interval,
    vqslod = va.info.VQSLOD,
    type = va.variantType,
    rfprob_all = va.rf.RF.probability["TP"],
    rfprob_bal = va.rf_old.RF.probability[2],
    pass = va.pass,
    is_mixed = va.isMixed,
    qd = va.info.QD,
    wassplit = va.wasSplit,
    mendel_errors = va.mendel.N,
    validated_denovo = va.validated_denovo,
    ac_unrelated = va.AC_unrelated[1].toInt,
    transmitted = va.tdt.nTransmitted,
    untransmitted = va.tdt.nUntransmitted,
    ac_orig = va.info.AC[va.aIndex - 1],
    ab_mean = va.info.ab_stats[va.aIndex - 1].mean,
    callrate = va.qc.callRate,
    af = va.qc.AF,
    ac = va.qc.AC'''

    rf_vds = hc.read('%s/sites/exac2.sites.RF2.vds' % root)

    hc.read('%s/v2_variantqc.vds' % root).annotate_variants_vds(rf_vds, root='va.rf_old').export_variants('%s/variantqc.txt.bgz' % root, columns)

    # TODO: remember when finalizing to .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0') and recalulate callStats

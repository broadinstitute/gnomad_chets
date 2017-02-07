
from variantqc import *

try:
    hc
except Exception, e:
    hc = HailContext()

# Actions
split = False
transmission = False
rf = True
finalize_rf = True
export_rf = True

syndip_compute = False
na12878_compute = False
syndip_export = True
na12878_export = True

root = 'gs://gnomad-exomes'
autosome_intervals = '%s/intervals/autosomes.txt' % root
evaluation_intervals = '%s/intervals/exome_evaluation_regions.v1.intervals' % root
high_coverage_intervals = '%s/intervals/high_coverage.auto.interval_list' % root
fam_path = '%s/variantqc/exac2.qctrios.fam' % root

raw_hardcall_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw.hardcalls.qc.vds'
raw_hardcalls_split_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw.hardcalls.splitmulti.qc.vds'

# print hc.read(raw_hardcall_vds_path).query_variants('variants.count()')[0]  # 15210030 sites
# print hc.read(raw_hardcalls_split_vds_path).query_variants('variants.count()')[0]  # 17215421 variants

if split:
    raw_hardcall_vds = hc.read(raw_hardcall_vds_path)
    write_split(raw_hardcall_vds, raw_hardcalls_split_vds_path)

sites_qc_vds_path = '%s/variantqc/exacv2.sites.qc.vds' % root

# Use raw VDS for determing true positives
if transmission:
    raw_hardcalls_split_vds = hc.read(raw_hardcalls_split_vds_path)
    transmission_mendel(raw_hardcalls_split_vds, sites_qc_vds_path, fam_path, autosome_intervals, mendel_path='%s/variantqc/exomes' % root)

extension = '.hardfilteronly'
rf_variantqc_path = '%s/variantqc/rf_explore/exacv2.rf%s.vds' % (root, extension)

# print hc.read(rf_variantqc_path).query_variants('variants.count()')[0]  # 16762994 autosomal variants

if rf:
    omni_vds = hc.read(omni_path)
    mills_vds = hc.read(mills_path)
    sites_qc_vds = hc.read(sites_qc_vds_path)

    rf_vds = sample_RF_training_examples(annotate_for_random_forests(sites_qc_vds, omni_vds, mills_vds))
    rf_vds = rf_vds.random_forests(training='va.train', label='va.label', root='va.rf', features=rf_features, num_trees=500, max_depth=5)#, perc_training=0.9))
    rf_vds.write(rf_variantqc_path, overwrite=True)

final_variantqc_path = '%s/variantqc/rf_explore/exacv2.variantqc%s.vds' % (root, extension)

if finalize_rf:
    rf_vds = hc.read(rf_variantqc_path)
    new_vds = (rf_vds.filter_variants_intervals(autosome_intervals)
               .annotate_variants_vds(hc.read('%s/variantqc/exac2_vqsr.vds' % root), code='va.info.VQSLOD = vds.info.VQSLOD')
               .annotate_variants_intervals(evaluation_intervals, root='va.evaluation_interval')  # warning: this is not a boolean
               .annotate_variants_intervals(high_coverage_intervals, root='va.high_coverage_interval')
               .annotate_variants_table('%s/variantqc/validatedDN.cut.txt.bgz' % root, 'Variant(CHROM, POSITION.toInt, REF, ALT)', code='va.validated_denovo = table.DataSet', config=hail.TextTableConfig(impute=True))
               .write(final_variantqc_path, overwrite=True)
    )

if export_rf:

    rf_vds = hc.read(final_variantqc_path)

    columns = '''chrom = v.contig,
    pos = v.start,
    ref = v.ref,
    alt = v.alt,
    evaluation_interval = !isMissing(va.evaluation_interval),
    high_coverage_interval = va.high_coverage_interval,
    vqslod = va.info.VQSLOD,
    type = va.variantType,
    train = if (va.train) va.label else NA: String,
    rfprob = va.rf.probability["TP"],
    pass = va.pass,
    qd = va.info.QD,
    wassplit = va.wasSplit,
    mendel_errors = va.mendel,
    validated_denovo = va.validated_denovo,
    ac_orig = va.info.AC[va.aIndex - 1],
    ac_all_raw = va.calldata.all_samples_raw.AC[1],
    ac_qc_raw = va.calldata.qc_samples_raw.AC[1],
    an_qc_raw = va.calldata.release_samples_raw.AN,
    ac_release_raw = va.calldata.release_samples_raw.AC[1],
    nrq_median = va.stats.qc_samples_raw.nrq_median,
    gq_median = va.stats.qc_samples_raw.gq_median,
    dp_median = va.stats.qc_samples_raw.dp_median,
    ab_median = va.stats.qc_samples_raw.ab_median'''

    rf_vds.export_variants('%s/variantqc/rf_explore/variantqc%s.txt.bgz' % (root, extension), columns)

# Truth sets
syndip_concordance_prefix = '%s/variantqc/syndip' % root
NA12878_concordance_prefix = '%s/variantqc/na12878' % root

concordance_annotations = ['chrom = v.contig',
                           'pos = v.start',
                           'ref = v.ref',
                           'alt = v.alt',
                           'indel = v.altAllele.isIndel',
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
                                      'rfprob = va.rf.rf.probability["TP"]'
                                      ])

if na12878_compute or syndip_compute:
    raw_hardcalls_split_vds = hc.read(raw_hardcalls_split_vds_path)

if na12878_export or syndip_export:
    rf_vds = hc.read(final_variantqc_path, sites_only=True)

if syndip_compute:
    compute_concordance(raw_hardcalls_split_vds,
                        hc.read(syndip_path).rename_samples('%s/variantqc/syndip.rename' % root),
                        'CHMI_CHMI3_Nex1',
                        syndip_high_conf_regions_path,
                        syndip_concordance_prefix)

if syndip_export:
    export_concordance(hc.read(syndip_concordance_prefix + ".v_concordance.vds")
                       .filter_variants_intervals(high_coverage_intervals)
                       .filter_variants_intervals(evaluation_intervals),
                       rf_vds,
                       truth_concordance_annotations,
                       syndip_concordance_prefix + extension)

if na12878_compute:
    compute_concordance(raw_hardcalls_split_vds,
                        hc.read(NA12878_path).rename_samples('%s/variantqc/na12878.rename' % root),
                        'C1975::NA12878',
                        NA12878_high_conf_exome_regions_path,
                        NA12878_concordance_prefix)

if na12878_export:
    export_concordance(hc.read(NA12878_concordance_prefix + ".v_concordance.vds")
                       .filter_variants_intervals(high_coverage_intervals)
                       .filter_variants_intervals(evaluation_intervals),
                       rf_vds,
                       truth_concordance_annotations,
                       NA12878_concordance_prefix + extension)

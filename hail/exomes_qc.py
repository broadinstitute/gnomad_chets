
from variantqc import *

try:
    hc
except Exception, e:
    hc = HailContext()

# Actions
hardcalls = False
split = False
transmission = False
rf = True
finalize_rf = True
export_rf = True

syndip_compute = False
na12878_compute = False
syndip_export = True
na12878_export = True

bucket = 'gs://gnomad-exomes'
autosome_intervals = '%s/intervals/autosomes.txt' % bucket
evaluation_intervals = '%s/intervals/exome_evaluation_regions.v1.intervals' % bucket
high_coverage_intervals = '%s/intervals/high_coverage.auto.interval_list' % bucket

root = '%s/variantqc' % bucket
fam_path = '%s/exac2.qctrios.fam' % root

full_vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
raw_hardcalls_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw.hardcalls.qc.vds'
raw_hardcalls_split_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw.hardcalls.splitmulti.qc.vds'

if hardcalls:
    # Beast mode
    (hc.read(full_vds_path)
     .annotate_samples_table('%s/super_meta.txt.bgz' % bucket, 'sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
     .annotate_variants_expr([get_variant_type_expr(),
                              'va.calldata.all_samples_raw = gs.callStats(g => v)',
                              'va.nAltAlleles = v.altAlleles.filter(a => a.alt != "*").length',
                              'va.nonsplit_alleles = v.altAlleles.map(a => a.alt)'])
     .annotate_alleles_expr(get_stats_expr("va.stats.all_samples_raw", medians=True))
     .histograms("va.hists.all_samples_raw")
     .annotate_samples_fam('%s/exac2.qctrios.fam' % root)
     .filter_samples_expr('sa.meta.drop_status == "keep" || (!isMissing(sa.fam.famID) && !("hard" ~ sa.meta.drop_condense)) || s.id == "C1975::NA12878" || s.id == "CHMI_CHMI3_Nex1"')
     .annotate_variants_expr('va.calldata.qc_samples_raw = gs.callStats(g => v)')
     .filter_alleles('va.calldata.qc_samples_raw.AC[aIndex] > 0', annotation='va.filtered_allele_indices = range(1, v.nAltAlleles + 1).map(i => !aIndices.toSet.contains(i))')
     .annotate_variants_expr('va.calldata.qc_samples_raw = gs.callStats(g => v), '
                             'va.calldata.release_samples_raw = gs.filter(g => sa.meta.drop_status == "keep").callStats(g => v)')
     .annotate_alleles_expr(get_stats_expr("va.stats.qc_samples_raw", medians=True))
     .histograms("va.hists.qc_samples_raw")
     .hardcalls()
     .min_rep() # Needs fixing before this can be run as is. Another strategy is to remove this, write a temp hardcalls here, then minrep
     .write(raw_hardcalls_vds_path))

if split:
    raw_hardcall_vds = hc.read(raw_hardcalls_vds_path)
    write_split(raw_hardcall_vds, raw_hardcalls_split_vds_path)

# print hc.read(raw_hardcalls_vds_path).query_variants('variants.count()')[0]  # 15210030 sites
# print hc.read(raw_hardcalls_split_vds_path).query_variants('variants.count()')[0]  # 17215421 variants

sites_qc_vds_path = '%s/gnomad.exomes.sites.qc.vds' % root

# Use raw VDS for determing true positives
if transmission:
    raw_hardcalls_split_vds = hc.read(raw_hardcalls_split_vds_path)
    transmission_mendel(raw_hardcalls_split_vds, sites_qc_vds_path, fam_path, autosome_intervals, mendel_path='%s/variantqc/exomes' % root)

rf_variantqc_path = '%s/gnomad.exomes.rf.vds' % root

# print hc.read(rf_variantqc_path).query_variants('variants.count()')[0]  # 16762994 autosomal variants

if rf:
    omni_vds = hc.read(omni_path)
    mills_vds = hc.read(mills_path)
    sites_qc_vds = hc.read(sites_qc_vds_path)

    rf_vds = sample_RF_training_examples(annotate_for_random_forests(sites_qc_vds, omni_vds, mills_vds))
    rf_vds = rf_vds.random_forests(training='va.train', label='va.label', root='va.rf', features=rf_features, num_trees=500, max_depth=5)#, perc_training=0.9))
    rf_vds.write(rf_variantqc_path, overwrite=True)

final_variantqc_path = '%s/gnomad.exomes.variantqc.vds' % root

if finalize_rf:
    rf_vds = hc.read(rf_variantqc_path)
    new_vds = (rf_vds.filter_variants_intervals(autosome_intervals)
               .annotate_variants_vds(hc.read('%s/gnomad.exomes.vqsr.vds' % root), code='va.info.VQSLOD = vds.info.VQSLOD')
               .annotate_variants_intervals(evaluation_intervals, root='va.evaluation_interval')  # warning: this is not a boolean
               .annotate_variants_intervals(high_coverage_intervals, root='va.high_coverage_interval')
               .annotate_variants_table('%s/validatedDN.cut.txt.bgz' % root, 'Variant(CHROM, POSITION.toInt, REF, ALT)', code='va.validated_denovo = table.DataSet', config=hail.TextTableConfig(impute=True))
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

    # rf_vds.export_variants('%s/gnomad.exomes.variantqc.txt.bgz' % root, columns)
    rf_vds.variants_keytable().to_dataframe().write.parquet('%s/gnomad.exomes.variantqc.parquet' % root)

# Truth sets
syndip_concordance_prefix = '%s/gnomad.exomes.syndip' % root
NA12878_concordance_prefix = '%s/gnomad.exomes.na12878' % root

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
                        hc.read(syndip_path).rename_samples('%s/syndip.rename' % root),
                        'CHMI_CHMI3_Nex1',
                        syndip_high_conf_regions_path,
                        syndip_concordance_prefix)

if syndip_export:
    export_concordance(hc.read(syndip_concordance_prefix + ".v_concordance.vds")
                       .filter_variants_intervals(high_coverage_intervals)
                       .filter_variants_intervals(evaluation_intervals),
                       rf_vds,
                       truth_concordance_annotations,
                       syndip_concordance_prefix)

if na12878_compute:
    compute_concordance(raw_hardcalls_split_vds,
                        hc.read(NA12878_path).rename_samples('%s/na12878.rename' % root),
                        'C1975::NA12878',
                        NA12878_high_conf_exome_regions_path,
                        NA12878_concordance_prefix)

if na12878_export:
    export_concordance(hc.read(NA12878_concordance_prefix + ".v_concordance.vds")
                       .filter_variants_intervals(high_coverage_intervals)
                       .filter_variants_intervals(evaluation_intervals),
                       rf_vds,
                       truth_concordance_annotations,
                       NA12878_concordance_prefix)

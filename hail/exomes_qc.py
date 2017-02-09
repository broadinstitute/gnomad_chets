
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
meta_path = '%s/super_meta.txt.bgz' % bucket
fam_path = '%s/exac2.qctrios.fam' % root

full_vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
raw_hardcalls_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw.hardcalls.qc.vds'
raw_hardcalls_split_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw.hardcalls.splitmulti.qc.vds'


def main():
    if hardcalls:
        full_vds = hc.read(full_vds_path)
        write_qc_hardcalls(full_vds, raw_hardcalls_vds_path, meta_path, fam_path)

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
        rf_vds.write(rf_variantqc_path)

    final_variantqc_path = '%s/gnomad.exomes.variantqc.vds' % root

    if finalize_rf:
        rf_vds = hc.read(rf_variantqc_path)
        (rf_vds.filter_variants_intervals(autosome_intervals)
         .annotate_variants_vds(hc.read('%s/gnomad.exomes.vqsr.vds' % root), code='va.info.VQSLOD = vds.info.VQSLOD')
         .annotate_variants_intervals(evaluation_intervals, root='va.evaluation_interval')  # warning: this is not a boolean
         .annotate_variants_intervals(high_coverage_intervals, root='va.high_coverage_interval')
         .annotate_variants_table('%s/validatedDN.cut.txt.bgz' % root, 'Variant(CHROM, POSITION.toInt, REF, ALT)', code='va.validated_denovo = table.DataSet', config=hail.TextTableConfig(impute=True))
         .write(final_variantqc_path)
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
                            hc.read(syndip_path).rename_samples('%s/gnomad.exomes.syndip.rename' % root),
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
                            hc.read(NA12878_path).rename_samples('%s/gnomad.exomes.na12878.rename' % root),
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


def write_qc_hardcalls(vds, output_path, meta_path, fam_path, adj=True):
    # Beast mode
    # Have not tested adj version yet
    ann_type = 'Adj' if adj else 'raw'
    vds = vds.annotate_samples_table(meta_path, 'sample', root='sa.meta', config=hail.TextTableConfig(impute=True))

    if adj: vds = vds.filter_genotypes(adj_criteria)

    # Get some variant information
    vds = (vds.annotate_variants_expr([get_variant_type_expr(),
                                       'va.calldata.all_samples_%s = gs.callStats(g => v)' % ann_type,
                                       'va.nAltAlleles = v.altAlleles.filter(a => a.alt != "*").length',
                                       'va.nonsplit_alleles = v.altAlleles.map(a => a.alt)']))

    # Annotate all smaples
    vds = (vds.annotate_alleles_expr(get_stats_expr("va.stats.all_samples_%s" % ann_type, medians=True))
           .histograms("va.hists.all_samples_%s" % ann_type))

    # Filter to QC samples
    vds = (vds.annotate_samples_fam(fam_path)
           .filter_samples_expr('sa.meta.drop_status == "keep" || (!isMissing(sa.fam.famID) && !("hard" ~ sa.meta.drop_condense)) || s.id == "C1975::NA12878" || s.id == "CHMI_CHMI3_Nex1"')
           .annotate_variants_expr('va.calldata.qc_samples_%s = gs.callStats(g => v)' % ann_type))

    # Filter out sites not in QC set
    vds = (vds.filter_alleles('va.calldata.qc_samples_%s.AC[aIndex] > 0' % ann_type,
                              annotation='va.filtered_allele_indices = range(1, v.nAltAlleles + 1).map(i => !aIndices.toSet.contains(i))'))

    # Calculate final QC stats
    vds = (vds.annotate_variants_expr(['va.calldata.qc_samples_%s = gs.callStats(g => v)' % ann_type,
                                    'va.calldata.release_samples_%s = gs.filter(g => sa.meta.drop_status == "keep").callStats(g => v)' % ann_type])
           .annotate_alleles_expr(get_stats_expr("va.stats.qc_samples_%s" % ann_type, medians=True))
           .histograms("va.hists.qc_samples_%s" % ann_type))

    # Current approach while smart shuffle not implemented
    temp_vds_path = output_path + '.temp.vds'
    vds.hardcalls().write(temp_vds_path)

    hc.read(temp_vds_path).min_rep().write(output_path)

    # Eventually this will prevail:
    # return (vds.hardcalls()
    #         .min_rep() # Needs fixing before this can be run as is. Another strategy is to remove this, write a temp hardcalls here, then minrep
    #         .write(output_path))


def write_split(input_vds, output_path):
    annotations = index_into_arrays([
        'va.stats.qc_samples_raw.gq',
        'va.stats.qc_samples_raw.dp',
        'va.stats.qc_samples_raw.nrq',
        'va.stats.qc_samples_raw.ab',
        'va.stats.qc_samples_raw.gq_median',
        'va.stats.qc_samples_raw.dp_median',
        'va.stats.qc_samples_raw.nrq_median',
        'va.stats.qc_samples_raw.ab_median'
    ])
    annotations.extend([
        'va.hasStar = va.nonsplit_alleles.exists(a => a == "*")',
        'va.wasMixed = va.variantType == "mixed"',
        'va.alleleType = if(v.altAllele.isSNP) "snv"'
        '   else if(v.altAllele.isInsertion) "ins"'
        '   else if(v.altAllele.isDeletion) "del"'
        '   else "complex"'])
    (input_vds
     .annotate_variants_expr(['va.nonsplit_alleles = v.altAlleles.map(a => a.alt)'])
     .split_multi()
     .annotate_variants_expr(annotations)
     .write(output_path))


def transmission_mendel(vds, output_vds_path, fam_path, autosomes_intervals, mendel_path=None):
    vds = vds.filter_variants_intervals(autosomes_intervals)

    if mendel_path is None: mendel_path = '/tmp/exomes'
    vds.mendel_errors(mendel_path, fam_path)

    (vds.tdt(fam_path)
     .annotate_variants_table(mendel_path + ".lmendel", 'SNP', code='va.mendel = table.N', config=hail.TextTableConfig(impute=True))
     .filter_samples_all()
     .write(output_vds_path))


if __name__ == '__main__':
    main()
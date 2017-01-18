__author__ = 'konrad'

try:
    hc
except NameError:
    from utils import *
    hc = HailContext()

# Input files
root = 'file:///mnt/lustre/konradk/exac'
meta_path = '%s/super_meta.txt.bgz' % root
full_vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
full_v1_vds_path = 'file:///mnt/lustre/lfran/exac/exac_all.new.vds'

adj_criteria = 'g.gq >= 20 && g.dp >= 10 & g.ad[1]/g.dp >= 0.25'


def prepare_files():
    CHROMS = map(str, range(1, 23))
    CHROMS.extend(['X', 'Y'])

    # Prepare some intervals files
    with open('%s/intervals/chrX.txt', 'w') as f:
        f.write('X:1-1000000000')
    with open('%s/intervals/chrY.txt', 'w') as f:
        f.write('Y:1-1000000000')
    with open('%s/intervals/autosomes.txt', 'w') as f:
        f.write('\n'.join(['%s:1-1000000000' % x for x in CHROMS]))


def write_hardcalls(input_path, output_path, meta_path, adj=True, metrics=True, partitions=10000, shuffle=True):
    get_variant_type_expression = '''va.variantType =
    let non_star = v.altAlleles.filter(a => a.alt != "*") in
        if (non_star.forall(a => a.isSNP))
            if (non_star.count() > 1)
                "multi-snv"
            else
                "snv"
        else if (non_star.forall(a => a.isIndel))
            if (non_star.count() > 1)
                "multi-indel"
            else
                "indel"
        else
            "mixed"'''

    hists = [
        'va.hists.%s.GQ_HIST = gs.filter(sa.meta.drop_status == "keep").map(g => g.gq).hist(0, 100, 20)',
        'va.hists.%s.DP_HIST = gs.filter(sa.meta.drop_status == "keep").map(g => g.dp).hist(0, 100, 20)',
        'va.hists.%s.AB_HIST = gs.filter(sa.meta.drop_status == "keep").map(g => 100*g.ad[%s]/g.dp).hist(0, 100, 20)'
    ]

    stats = ['va.stats.%s.gq = gs.filter(g => g.isCalledNonRef).map(g => g.gq).stats()',
             'va.stats.%s.dp = gs.filter(g => g.isCalledNonRef).map(g => g.dp).stats()',
             'va.stats.%s.nrq = gs.filter(g => g.isCalledNonRef).map(g => g.dosage[0]).stats()',
             'va.stats.%s.ab = gs.filter(g => g.isCalledNonRef).map(g => g.ad[1]/g.dp).stats()']

    template = ('%(destination)s = let sorted_vals = gs.filter(g => g.isCalledNonRef && !isMissing(%(metric)s)).map(g => %(metric)s).collect().sort() in '
                'if (sorted_vals.size == 0) NA: Double else '
                'if (sorted_vals.size %% 2 == 1) sorted_vals[(sorted_vals.size/2).toInt] else '
                '(sorted_vals[(sorted_vals.size/2).toInt] + sorted_vals[(sorted_vals.size/2).toInt - 1])/2.0')
    medians = [('g.gq', 'va.stats.raw.gq_median'), ('g.dp', 'va.stats.raw.dp_median'), ('g.dosage[0]', 'va.stats.raw.nrq_median'), ('g.ad[1]/g.dp', 'va.stats.raw.ab_median')]

    pre_adj_expression = [get_variant_type_expression]
    pre_adj_expression.extend([x % 'raw' for x in stats])
    pre_adj_expression.extend([x % 'raw' for x in hists])
    pre_adj_expression.extend([template % {'metric': metric, 'destination': destination} for (metric, destination) in medians])
    pre_adj_expression.append('va.calldata.allsamples_raw = gs.callStats(v)')
    pre_adj_expression = ','.join(pre_adj_expression)

    post_adj_expression = [x % 'Adj' for x in stats]
    post_adj_expression.extend([x % 'Adj' for x in hists])
    post_adj_expression.append('va.calldata.allsamples_Adj = gs.callStats(v)')

    if adj:
        if metrics:
            return (hc.read(input_path)
                    .annotate_samples_table(meta_path, 'sample', impute=True, root='sa.meta')
                    .annotate_alleles_expr(pre_adj_expression)
                    .filter_genotypes(adj_criteria)
                    .annotate_alleles_expr(post_adj_expression)
                    .hardcalls()
                    .repartition(numPartitions=partitions,shuffle=shuffle)
                    .write(output_path)
            )
        else:
            return (hc.read(input_path)
                    .annotate_samples_table(meta_path, 'sample', impute=True, root='sa.meta')
                    .filter_genotypes(adj_criteria)
                    .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')
                    .hardcalls()
                    .repartition(numPartitions=partitions,shuffle=shuffle)
                    .write(output_path)
            )
    else:
        return (hc.read(input_path)
                .annotate_samples_table(meta_path, 'sample', impute=True, root='sa.meta')
                .hardcalls()
                .repartition(numPartitions=partitions,shuffle=shuffle)
                .write(output_path)
        )


def write_split(input_path, output_path):
    return hc.read(input_path)\
        .split_multi()\
        .write(output_path)

# All hardcalls
raw_hardcallvds_path = '%s/hardcalls/v2/exacv2.raw.hardcalls.qc.vds' % root
hardcallvds_path = '%s/hardcalls/v2/exacv2.hardcalls.qc.vds' % root
v1_hardcallvds_path = '%s/hardcalls/v1/exacv1.hardcalls.qc.vds' % root

# All splitmulti hardcalls
raw_split_hardcallvds_path = '%s/hardcalls/v2/exacv2.raw.hardcalls.splitmulti.qc.vds' % root
split_hardcallvds_path = '%s/hardcalls/v2/exacv2.hardcalls.splitmulti.qc.vds' % root
v1_split_hardcallvds_path = '%s/hardcalls/v1/exacv1.hardcalls.splitmulti.qc.vds' % root

# write_hardcalls(full_vds_path, raw_hardcallvds_path, adj=False)
# write_hardcalls(full_vds_path, hardcallvds_path)
# write_hardcalls(full_v1_vds_path, v1_hardcallvds_path)
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


def get_transmitted_singletons(vds_path, output_vds_path):
    return (hc.read(vds_path)
            .filter_variants_intervals('%s/intervals/autosomes.txt' % root)
            .tdt(fam_path)
            .filter_variants_expr('va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1]==2')
            .filter_samples_all()
            .write(output_vds_path))

tdt_vds_path = '%s/variantqc/v2_tdt.raw.vds' % root

# Use raw VDS for determing true positives
# get_transmitted_singletons(raw_split_hardcallvds_path, tdt_vds_path)


def run_random_forests(vds_path, output_vds_path):
    vds = hc.read(vds_path)

    features = ['va.variantType',
                'va.info.QD',
                'va.info.MQ',
                'va.info.MQRankSum',
                'va.info.FS',
                'va.info.SOR',
                'va.info.InbreedingCoeff',
                'va.info.ReadPosRankSum',
                'va.stats.nrq_median',
                'va.stats.ab_median',
                'va.stats.dp_median',
                'va.stats.gq_median']

    truth_dir = '%s/variantqc/truth_sets' % root
    omni_vds = hc.read('%s/1000G_omni2.5.b37.splitmulti.vds' % truth_dir)
    mills_vds = hc.read('%s/Mills_and_1000G_gold_standard.indels.b37.splitmulti.vds' % truth_dir)
    transmission_vds = hc.read(tdt_vds_path)


    new_vds = (vds.annotate_variants_vds(omni_vds, code='va.omni = isDefined(vds)')
               .annotate_variants_vds(mills_vds, code='va.mills = isDefined(vds)')
               .annotate_variants_vds(transmission_vds, code='va.transmitted_singleton = isDefined(vds)')
               .annotate_variants_expr('va.TP = va.omni || va.mills || va.transmitted_singleton, '
                                          'va.FP = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30')
               .annotate_variants_expr('va.label = if(!isMissing(va.FP) && va.FP) "FP" else if(va.TP) "TP" else NA: String, '
                                       'va.train = v.contig != "20" && (va.TP || va.FP)')
               .annotate_global_expr('global.ac_hist_indels_mills = variants.filter(v => va.mills).map(v => log10(va.info.AF)).hist(-6, 0, 20),'
                                     'global.ac_hist_indels_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isIndel).map(v => log10(va.info.AF)).hist(-6, 0, 20),'
                                     'global.ac_hist_snps_omni = variants.filter(v => va.omni).map(v => log10(va.info.AF)).hist(-6, 0, 20),'
                                     'global.ac_hist_snps_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isSNP).map(v => log10(va.info.AF)).hist(-6, 0, 20),'
                                     'global.indels_mills = variants.filter(v => va.mills).count(),'
                                     'global.indels_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isIndel).count(),'
                                     'global.snps_omni = variants.filter(v => va.omni).count(),'
                                     'global.snps_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isSNP).count(),'
                                     'global.nTP = variants.filter(x => va.label == "TP").count(), '
                                     'global.nFP = variants.filter(x => va.label == "FP").count()')
               .show_globals()
               .random_forests(training='va.train', label='va.label', root='va.RF', features=features)
               .write(output_vds_path, overwrite=True)
    )

    new_vds.export_variants('exac2_rf.txt.bgz', 'chrom = v.contig, pos = v.start, ref = v.ref, alt = v.altAllele, type = va.variantType, label = va.label, rfprob = va.RF.probability["TP"]')


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
__author__ = 'konrad'

from utils import *

ab_cutoff = 0.2
adj_criteria = 'g.gq >= 20 && g.dp >= 10 && (' \
               '!g.isHet || ' \
               '(g.gtj == 0 && g.ad[1]/g.dp >= %(ab)s) || ' \
               '(g.gtj > 0 && g.ad[0]/g.dp >= %(ab)s && g.ad[1]/g.dp >= %(ab)s)' \
               ')' % {'ab': ab_cutoff}

rf_features = ['va.variantType',
            'va.info.QD',
            'va.info.MQ',
            'va.info.MQRankSum',
            'va.info.FS',
            'va.info.SOR',
            'va.info.InbreedingCoeff',
            'va.info.ReadPosRankSum',
            'va.stats.raw.nrq_median',
            'va.stats.raw.ab_median',
            'va.stats.raw.dp_median',
            'va.stats.raw.gq_median']


def write_hardcalls(vds, output_path, meta_path, adj=True, metrics=True, partitions=10000, shuffle=True):

    out = vds.annotate_samples_table(meta_path, 'sample', root='sa.meta', config=hail.TextTableConfig(impute=True))

    if metrics:
        pre_adj_expression = get_variant_type_expr("va.variantType")
        pre_adj_expression.append('va.calldata.all_samples_raw = gs.callStats(g => v)')
        out = (out
               .annotate_variants_expr(pre_adj_expression)
               .annotate_alleles_expr(get_stats_expr("va.stats.raw", medians=True))
               .histograms("va.hists.raw"))

    if adj:
        out = (
            out.filter_genotypes(adj_criteria)
            .annotate_variants_expr('va.calldata.all_samples_Adj = gs.callStats(g => v)')
            .filter_alleles('va.calldata.allsamples_Adj.AC[aIndex] == 0', subset=True, keep=False)
        )

        if metrics:
            out = (out
                   .annotate_variants_expr('va.calldata.all_samples_Adj = gs.callStats(g => v)')
                   .annotate_alleles_expr(get_stats_expr("va.stats.Adj", medians=True))
                   .histograms("va.hists.Adj"))

    return (out.hardcalls()
            .repartition(partitions, shuffle=shuffle)
            .write(output_path))


def write_split(input_vds, output_path):
    a_indexed = ['va.calldata.allsamples_raw',
                 'va.stats.raw.gq',
                 'va.stats.raw.dp',
                 'va.stats.raw.nrq',
                 'va.stats.raw.ab',
                 'va.stats.raw.gq_median',
                 'va.stats.raw.dp_median',
                 'va.stats.raw.nrq_median',
                 'va.stats.raw.ab_median']
    return (input_vds
            .split_multi()
            .annotate_variants_expr(index_into_arrays(a_indexed))
            .write(output_path))


def get_transmitted_singletons(vds, output_vds_path, fam_path, autosomes_intervals):
    return (vds
            .filter_variants_intervals(autosomes_intervals)
            .tdt(fam_path)
            .filter_variants_expr('va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1] == 2')
            .filter_samples_all()
            .write(output_vds_path))


def annotate_for_random_forests(vds, transmission_vds=None, omni_vds=None, mills_vds=None, sample=True):

    if transmission_vds is not None:
        vds = vds.annotate_variants_vds(transmission_vds, code='va.transmitted_singleton = isDefined(vds)')

    if omni_vds is not None:
        vds = vds.annotate_variants_vds(omni_vds, code='va.omni = isDefined(vds)')

    if mills_vds is not None:
        vds = vds.annotate_variants_vds(mills_vds, code='va.mills = isDefined(vds)')

    # Calculating median for each feature to impute as needed
    sample_text = '&& pcoin([1.0, 1000000 / global.variantsByType[va.variantType].count].min)' if sample else ''
    global_features_expr_median = []
    for feature in features_for_median:
        for variant_type in variant_types:
            global_features_expr_median.append(
                'global.median.`%s`.`%s` = variants'
                '.filter(v => va.variantType == "%s" %s)'
                '.map(v => %s).collect().median()' % (feature, variant_type, variant_type, sample_text, feature))

    # Reformat into a dict
    median_to_dicts = []
    for variant_type in variant_types:
        median_to_dicts.append('{variantType : "%s", %s }' % (variant_type,
                                                ",\n".join(['`%s` : global.median.`%s`.`%s`' % (feature, feature, variant_type) for feature in features_for_median])))

    dict_median_expression = 'global.median = index([ ' + ",\n".join(median_to_dicts) + '], variantType)'

    variants_features_imputation = ['%(f)s = if(isDefined(%(f)s)) %(f)s else global.median[va.variantType].`%(f)s`'
                                    % {'f': feature} for feature in features_for_median]

    global_missing_features_expr_before = ['global.missing.before.%s = variants.filter(x => isMissing(%s)).count()' % (a, a) for a in rf_features]
    global_missing_features_expr_after = ['global.missing.after.%s = variants.filter(x => isMissing(%s)).count()' % (a, a) for a in rf_features]

    return (vds
            # .annotate_global_expr_by_variant('global.variantsByType = index(variants.map(v => va.variantType).counter(),key)')
            # .annotate_global_expr_by_variant(global_missing_features_expr_before)
            .annotate_global_expr_by_variant(global_features_expr_median)
            .annotate_global_expr_by_variant(dict_median_expression)
            .annotate_variants_expr(variants_features_imputation)
            # .annotate_global_expr_by_variant(global_missing_features_expr_after)
            .annotate_variants_expr('va.TP = va.omni || va.mills || va.transmitted_singleton, '
                                    'va.FP = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30')
            .annotate_variants_expr('va.label = if(!isMissing(va.FP) && va.FP) "FP" else if(va.TP) "TP" else NA: String, '
                                    'va.train = v.contig != "20" && (va.TP || va.FP)')
            .annotate_global_expr_by_variant('global.nTP = variants.filter(x => va.label == "TP").count(), '
                                             'global.nFP = variants.filter(x => va.label == "FP").count()')
            .show_globals()
    )


__author__ = 'konrad'

from utils import *
import sys

ab_cutoff = 0.2
adj_criteria = 'g.gq >= 20 && g.dp >= 10 && (' \
               '!g.isHet || ' \
               '(g.gtj == 0 && g.ad[1]/g.dp >= %(ab)s) || ' \
               '(g.gtj > 0 && g.ad[0]/g.dp >= %(ab)s && g.ad[1]/g.dp >= %(ab)s)' \
               ')' % {'ab': ab_cutoff}

rf_features = ['va.alleleType',
              'va.nAltAlleles',
               'va.wasMixed',
               'va.hasStar',
            'va.info.MQRankSum',
            'va.info.SOR',
            'va.info.InbreedingCoeff',
            'va.info.ReadPosRankSum',
            'va.stats.qc_samples_raw.nrq_median',
            'va.stats.qc_samples_raw.ab_median',
            'va.stats.qc_samples_raw.dp_median',
            'va.stats.qc_samples_raw.gq_median']

features_for_median = [
    'va.info.MQRankSum',
    'va.info.ReadPosRankSum',
    'va.stats.qc_samples_raw.ab_median'
]

variant_types = ['snv', 'multi-snv', 'indel', 'multi-indel', 'mixed']


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
    a_indexed = [
        'va.stats.qc_samples_raw.gq',
        'va.stats.qc_samples_raw.dp',
        'va.stats.qc_samples_raw.nrq',
        'va.stats.qc_samples_raw.ab',
        'va.stats.qc_samples_raw.gq_median',
        'va.stats.qc_samples_raw.dp_median',
        'va.stats.qc_samples_raw.nrq_median',
        'va.stats.qc_samples_raw.ab_median'
    ]
    return (input_vds
            .split_multi()
            .annotate_variants_expr(index_into_arrays(a_indexed))
            .write(output_path))


def transmission_mendel(vds, output_vds_path, fam_path, autosomes_intervals):
    return (vds
            .filter_variants_intervals(autosomes_intervals)
            .tdt(fam_path)
            .mendel_errors('va.mendel', fam_path)
            .filter_samples_all()
            .write(output_vds_path))


def annotate_for_random_forests(vds, omni_vds, mills_vds, sample=True):

    vds_schema = [f.name for f in vds.variant_schema.fields]

    if "tdt" not in vds_schema or "mendel" not in vds_schema:
        print >> sys.stderr, "va.tdt or va.mendel missing"
        sys.exit(2)

    vds = vds.annotate_variants_vds(omni_vds, code='va.omni = isDefined(vds)')
    vds = vds.annotate_variants_vds(mills_vds, code='va.mills = isDefined(vds)')


    vds = (vds.annotate_variants_expr('va.transmitted_singleton = va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1] == 2,'
                                      'va.transmission_disequilibrated = va.tdt.pval < 0.001,'
                                      'va.mendel_excess = va.mendel >= 10,'
                                      'va.failing_hard_filters = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30,'
                                      'va.TP = va.omni || va.mills || va.transmitted_singleton'))

    # Variants per type (for downsampling)
    variant_counts = dict([(x['key'], x['count']) for x in vds.query_variants('variants.map(v => va.variantType).counter()')[0]])

    # Missing features before imputation
    # vds.query_variants(['variants.filter(x => isMissing(%s)).count()' % (a, a) for a in rf_features])

    # Prepare query for "median per feature per variant type"
    sample_text = '&& pcoin([1.0, 1000000 / global.variantsByType[va.variantType].count].min)' if sample else ''
    feature_medians_expr = []
    for feature in features_for_median:
        for variant_type in variant_types:
            feature_medians_expr.append(
                'variants.filter(v => va.variantType == "%s" %s)'
                '.map(v => %s).collect().median()' % (variant_type, sample_text, feature))

    # Process query into dict
    feature_medians_query = vds.query_variants(feature_medians_expr)
    i = 0
    from collections import defaultdict
    feature_medians = defaultdict(dict)
    for feature in features_for_median:
        for variant_type in variant_types:
            feature_medians[feature][variant_type] = feature_medians_query[i]
            i += 1

    variants_features_imputation = ['%(f)s = if(isDefined(%(f)s)) %(f)s else global.median["%(f)s"][va.variantType]'
                                    % {'f': feature} for feature in features_for_median]

    #Get number of training examples
    training_counts = vds.query_variants(['variants.filter(v => va.transmission_disequilibrated).count()',
                              'variants.filter(v => va.mendel_excess).count()',
                              'variants.filter(v => va.failing_hard_filters).count()',
                              'variants.filter(v => va.omni || va.mills || va.transmitted_singleton).count()'])

    print(training_counts)

    ntraining = min(sum(training_counts[0:2]),training_counts[3])
    training_probs = {
        'tp': training_counts[4] / ntraining,
        'tdt': training_counts[0] / (1 / 3 * ntraining),
        'mendel': training_counts[1] / (1 / 3 * ntraining),
        'hard': training_counts[2] / (1 / 3 * ntraining)
    }


    vds = (vds
            .annotate_global_py('global.variantsByType', variant_counts, TDict(TLong()))
            .annotate_global_py('global.median', feature_medians, TDict(TDict(TDouble())))
            .annotate_variants_expr(variants_features_imputation)
            .annotate_variants_expr('va.FP = va.transmission_disequilibrated || va.mendel_excess || va.failing_hard_filters')
            .annotate_variants_expr('va.label = if(!isMissing(va.FP) && va.FP) "FP" else if(va.TP) "TP" else NA: String, '
                                    'va.train = (va.TP && pcoin(%(tp).3f)) || '
                                    '(va.transmission_disequilibrated && pcoin(%(tdt).3f)) ||'
                                    '(va.mendel_excess && pcoin(%(mendel).3f)) ||'
                                    '(va.failing_hard_filters && pcoin(%(hard).3f))' % training_probs)
    )
    print(vds.query_variants(['variants.filter(x => va.label == "TP").count()',
                             'variants.filter(x => va.label == "FP").count()',
                              'variants.filter(x => va.label == "TP" && va.train).count()',
                              'variants.filter(x => va.label == "FP" && va.train).count()']))
    return vds


def filter_for_concordance(vds, high_conf_regions):
    return(
        vds.filter_variants_intervals(lcr_path, keep=False)
        .filter_variants_intervals(decoy_path, keep=False)
        .filter_variants_intervals(high_conf_regions, keep=True)
    )


def compute_concordance(vds, truth_vds, sample, high_conf_regions, out_prefix):
    truth = filter_for_concordance(truth_vds, high_conf_regions=high_conf_regions)

    (s_concordance, v_concordance) = (filter_for_concordance(vds, high_conf_regions=high_conf_regions)
                                      .filter_variants_intervals(lcr_path, keep=False)
                                      .filter_variants_intervals(decoy_path, keep=False)
                                      .filter_samples_expr('s.id == "%s"' % sample, keep=True)
                                      .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count()>0', keep=True)
                                      .concordance(right=truth)
                                      )
    s_concordance.write(out_prefix + ".s_concordance.vds")
    v_concordance.write(out_prefix + ".v_concordance.vds")


# def compute_concordance(vds, rf_vds, sample, truth_vds, high_conf_regions, out_prefix, out_annotations, recompute=True):
#     if recompute:
#         recompute_concordance(vds, sample, truth_vds, high_conf_regions, out_prefix)
#
#     else:
#         v_concordance = hc.read(out_prefix + ".v_concordance.vds")
#     export_concordance(v_concordance, rf_vds, out_prefix, out_annotations)


def export_concordance(conc_vds, rf_vds, out_annotations, out_prefix):
    (
        conc_vds.annotate_variants_vds(rf_vds, root='va.rf')
        .annotate_global_py('global.gt_mappings', ["missing", "no_call" ,"homref" ,"het" ,"homvar"], TArray(TString()))
        .annotate_variants_expr('va.gt_arr = range(5).find(i => va.concordance[i].exists(x => x > 0))')
        .annotate_variants_expr('va.truth_gt =  global.gt_mappings[va.gt_arr],'
                                'va.called_gt = global.gt_mappings[range(5).find(i => va.concordance[va.gt_arr][i] >0)]')
        .annotate_variants_expr('va.variantType = if(isDefined(va.rf.variantType)) va.rf.variantType '
                                'else if(v.altAlleles.forall(x => x.isSNP)) "snv" '
                                'else if(v.altAlleles.forall(x => x.isIndel)) "indel"'
                                'else "mixed"')
        .export_variants(out_prefix + ".stats.txt.bgz", ",".join(out_annotations))
     )
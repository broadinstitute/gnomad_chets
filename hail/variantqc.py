__author__ = 'konrad'

from utils import *
import sys
from collections import defaultdict

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


def sample_RF_training_examples(vds,
                                tp_criteria="va.omni || va.mills || va.transmitted_singleton",
                                fp_criteria="va.failing_hard_filters",
                                fp_to_tp=1.0):

    vds = vds.annotate_variants_expr(["va.TP = " + tp_criteria,
                                      "va.FP = " + fp_criteria])

    # Get number of training examples
    training_criteria = ['va.TP', 'va.FP']
    training_counts = dict(zip(training_criteria, vds.query_variants(
        ['variants.filter(v => %s).count()' % criterion for criterion in training_criteria])))

    # Get titvs of each training criterion
    pprint(
        dict(zip(training_criteria, vds.query_variants(['variants.filter(v => %s && v.altAllele.isTransition).count()/'
                                                        'variants.filter(v => %s && v.altAllele.isTransversion).count()' %
                                                        (criterion, criterion) for criterion in training_criteria]))))

    # Balancing FPs to match TP rate
    print("\nTraining examples:")
    pprint(training_counts)

    training_probs = {'va.TP': 1,
                      'va.FP': float(fp_to_tp) * training_counts['va.TP'] / training_counts['va.FP']}

    print("Probability of using training example:")
    pprint(training_probs)

    training_selection = ' || '.join(['%s && pcoin(%.3f)' % (crit, prob) for crit, prob in training_probs.items()])

    vds = (vds
           .annotate_variants_expr(
        'va.label = if(!isMissing(va.FP) && va.FP) "FP" else if(va.TP) "TP" else NA: String, '
        'va.train = %s' % training_selection)
           )

    label_criteria = ['va.label == "TP"', 'va.label == "FP"', 'va.label == "TP" && va.train',
                      'va.label == "FP" && va.train']
    print("\nNumber of training examples used:")
    pprint(
        dict(zip(label_criteria, vds.query_variants(['variants.filter(x => %s).count()' % x for x in label_criteria]))))

    return vds


def annotate_for_random_forests(vds, omni_vds, mills_vds, sample=True):

    vds_schema = [f.name for f in vds.variant_schema.fields]

    if "tdt" not in vds_schema:
        print >> sys.stderr, "va.tdt missing"
        sys.exit(2)

    vds = vds.annotate_variants_vds(omni_vds, code='va.omni = isDefined(vds)')
    vds = vds.annotate_variants_vds(mills_vds, code='va.mills = isDefined(vds)')

    vds = (vds.annotate_variants_expr('va.transmitted_singleton = va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1] == 2,'
                                      'va.transmission_disequilibrated = va.tdt.pval < 0.001,'
                                      # 'va.mendel_excess = va.mendel >= 10,'
                                      'va.failing_hard_filters = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30'
                                      ))

    # Variants per type (for downsampling)
    variant_counts = dict([(x['key'], x['count']) for x in vds.query_variants('variants.map(v => va.variantType).counter()')[0]])
    print("\nCount by variant type:")
    pprint(variant_counts)

    vds = vds.annotate_global_py('global.variantsByType', variant_counts, TDict(TLong()))

    # Missing features before imputation
    # vds.query_variants(['variants.filter(x => isMissing(%s)).count()' % (a, a) for a in rf_features])

    # Prepare query for "median per feature per variant type"
    sample_text = '&& pcoin([1.0, 1000000 / global.variantsByType[va.variantType]].min)' if sample else ''
    feature_medians_expr = []
    for feature in features_for_median:
        for variant_type in variant_types:
            feature_medians_expr.append(
                'variants.filter(v => va.variantType == "%s" %s)'
                '.map(v => %s).collect().median()' % (variant_type, sample_text, feature))

    # Process query into dict
    feature_medians_query = vds.query_variants(feature_medians_expr)
    i = 0
    feature_medians = defaultdict(dict)
    for feature in features_for_median:
        for variant_type in variant_types:
            feature_medians[feature][variant_type] = feature_medians_query[i]
            i += 1
    print("\nMedians per feature per variant type")
    pprint(dict(feature_medians))

    variants_features_imputation = ['%(f)s = if(isDefined(%(f)s)) %(f)s else global.median["%(f)s"][va.variantType]'
                                    % {'f': feature} for feature in features_for_median]

    vds = (vds
           .annotate_global_py('global.median', feature_medians, TDict(TDict(TDouble())))
           .annotate_variants_expr(variants_features_imputation)
    )

    return vds


def filter_for_concordance(vds, high_conf_regions=None):
    vds = (vds.filter_variants_intervals(lcr_path, keep=False)
           .filter_variants_intervals(decoy_path, keep=False))
    if high_conf_regions is None:
        return vds
    else:
        return vds.filter_variants_intervals(high_conf_regions, keep=True)


def compute_concordance(vds, truth_vds, sample, high_conf_regions, out_prefix):
    truth = filter_for_concordance(truth_vds, high_conf_regions=high_conf_regions)

    vds = filter_for_concordance(vds, high_conf_regions=high_conf_regions)

    if sample.startswith('gs://'):
        vds = vds.filter_samples_list(sample).filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0', keep=True)
    elif sample != '':
        vds = vds.filter_samples_expr('s == "%s"' % sample, keep=True).filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0', keep=True)

    (s_concordance, v_concordance) = vds.concordance(right=truth)  # TODO: make sure truth is already minrepped
    s_concordance.write(out_prefix + ".s_concordance.vds")
    v_concordance.write(out_prefix + ".v_concordance.vds")


def export_concordance(conc_vds, rf_vds, out_annotations, out_prefix, single_sample=True):
    vds = conc_vds.annotate_variants_vds(rf_vds, root='va.rf')

    if single_sample:
        vds = (vds.annotate_global_py('global.gt_mappings', ["missing", "no_call" ,"homref" ,"het" ,"homvar"], TArray(TString()))
               .annotate_variants_expr('va.gt_arr = range(5).find(i => va.concordance[i].exists(x => x > 0))')
               .annotate_variants_expr('va.called_gt =  global.gt_mappings[va.gt_arr],'
                                       'va.truth_gt = global.gt_mappings[range(5).find(i => va.concordance[va.gt_arr][i] > 0)]'))
    else:
        vds = (vds.annotate_variants_expr('va.correct = va.concordance[3][3] + va.concordance[4][4], '
                                          'va.both_ref = va.concordance[2][2], '
                                          'va.wrong_call = va.concordance[2][3] + va.concordance[3][2] + va.concordance[2][4] + '
                                          'va.concordance[4][2] + va.concordance[3][4] + va.concordance[4][3], '
                                          'va.missing_called = va.concordance[0][2] + va.concordance[0][3] + va.concordance[0][4], '
                                          'va.missing_gt_called = va.concordance[1][2] + va.concordance[1][3] + va.concordance[1][4], '
                                          'va.missing_truth = va.concordance[2][0] + va.concordance[3][0] + va.concordance[4][0], '
                                          'va.missing_gt_truth = va.concordance[2][1] + va.concordance[3][1] + va.concordance[4][1], '
                                          'va.total = va.concordance[2:].map(x => x[2:].sum).sum - va.concordance[2][2]'))

    return (vds.annotate_variants_expr('va.variantType = if(isDefined(va.rf.variantType)) va.rf.variantType '
                                       'else if(v.altAlleles.forall(x => x.isSNP)) "snv" '
                                       'else if(v.altAlleles.forall(x => x.isIndel)) "indel"'
                                       'else "mixed"')
            .export_variants(out_prefix + ".stats.txt.bgz", ",".join(out_annotations))
    )
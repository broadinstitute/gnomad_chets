from utils import *
import sys
from collections import defaultdict

vqsr_features = ['va.info.MQRankSum',
                 'va.info.SOR',
                 'va.info.InbreedingCoeff',
                 'va.info.ReadPosRankSum',
                 'va.info.FS',
                 'va.info.QD',
                 'va.info.MQ',
                 'va.info.DP'
                 ]

rf_features = ['va.alleleType',
               'va.nAltAlleles',
               'va.wasMixed',
               'va.hasStar',
               'va.info.MQRankSum',
               'va.info.SOR',
               'va.info.InbreedingCoeff',
               'va.info.ReadPosRankSum',
               'va.stats.qc_samples_raw.qd',
               'va.stats.qc_samples_raw.pab.max'
               ]

features_for_median = [
    'va.info.MQRankSum',
    'va.info.ReadPosRankSum',
    'va.stats.qc_samples_raw.ab_median',
    'va.stats.qc_samples_raw.best_ab'
]

variant_types = ['snv', 'multi-snv', 'indel', 'multi-indel', 'mixed']


def sample_RF_training_examples(vds,
                                tp_criteria="va.omni || va.mills || va.transmitted_singleton",
                                fp_criteria="va.failing_hard_filters",
                                fp_to_tp=1.0):

    training_classes_expr = ["va.TP = (%s) && !((%s).orElse(false))" % (tp_criteria, fp_criteria),
                                      "va.FP = (%s) && !((%s).orElse(false))" % (fp_criteria, tp_criteria)]

    logger.info("Training classes defined as: " + ",".join(training_classes_expr))

    vds = vds.annotate_variants_expr(training_classes_expr)

    # Get number of training examples
    training_criteria = ['va.TP', 'va.FP']
    training_counts = dict(zip(training_criteria, vds.query_variants(
        ['variants.filter(v => %s).count()' % criterion for criterion in training_criteria])))

    # Get titvs of each training criterion
    logger.info(pformat(
        dict(zip(training_criteria, vds.query_variants(['variants.filter(v => %s && v.altAllele.isTransition).count()/'
                                                        'variants.filter(v => %s && v.altAllele.isTransversion).count()' %
                                                        (criterion, criterion) for criterion in training_criteria]))))
    )

    # Balancing FPs to match TP rate
    logger.info("\nTraining examples:\n%s" % pformat(training_counts))

    if fp_to_tp > 0:
        if(training_counts['va.TP'] < training_counts['va.FP']):
            training_probs = {'va.TP': 1,
                              'va.FP': float(fp_to_tp) * training_counts['va.TP'] / training_counts['va.FP']}
        else:
            training_probs = {'va.TP': float(fp_to_tp) * training_counts['va.FP'] / training_counts['va.TP'],
                              'va.FP': 1}

        logger.info("Probability of using training example:\n%s" % pformat(training_probs))

        training_selection = ' || '.join(['%s && pcoin(%.3f)' % (crit, prob) for crit, prob in training_probs.items()])
    else:
        logger.info("Using all training examples.")
        training_selection = 'va.TP || va.FP'

    vds = (vds
           .annotate_variants_expr(
        'va.label = if(!isMissing(va.FP) && va.FP) "FP" else orMissing(va.TP, "TP"), '
        'va.train = %s' % training_selection)
           )

    label_criteria = ['va.label == "TP"', 'va.label == "FP"', 'va.label == "TP" && va.train',
                      'va.label == "FP" && va.train']
    logger.info("\nNumber of training examples used:\n%s" % pformat(
        dict(zip(label_criteria, vds.query_variants(['variants.filter(x => %s).count()' % x for x in label_criteria]))))
                )

    return vds


def annotate_for_random_forests(vds, omni_vds, mills_vds, sample=True):

    vds_schema = [f.name for f in vds.variant_schema.fields]

    if "tdt" not in vds_schema:
        logger.fatal("va.tdt missing")
        sys.exit(2)

    vds = vds.annotate_variants_vds(omni_vds, expr='va.omni = isDefined(vds)')
    vds = vds.annotate_variants_vds(mills_vds, expr='va.mills = isDefined(vds)')

    vds = (vds.annotate_variants_expr('va.transmitted_singleton = va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1] == 2,'
                                      'va.transmission_disequilibrated = va.tdt.pval < 0.001,'
                                      # 'va.mendel_excess = va.mendel >= 10,'
                                      'va.failing_hard_filters = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30'
                                      ))

    # Variants per type (for downsampling)
    variant_counts =vds.query_variants('variants.map(v => va.variantType).counter()')
    logger.info("\nCount by variant type:\n%s" % pformat(variant_counts))

    vds = vds.annotate_global('global.variantsByType', variant_counts, TDict(TString(),TLong()))

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
    logger.info("\nMedians per feature per variant type\n%s" % pformat(dict(feature_medians)))

    variants_features_imputation = ['%(f)s = if(isDefined(%(f)s)) %(f)s else global.median["%(f)s"][va.variantType]'
                                    % {'f': feature} for feature in features_for_median]

    vds = (vds
           .annotate_global('global.median', feature_medians, TDict(TString(),TDict(TString(),TDouble())))
           .annotate_variants_expr(variants_features_imputation)
    )

    return vds


def filter_for_concordance(vds, samples, high_conf_regions=None):
    vds = filter_low_conf_regions(vds, high_conf_regions=high_conf_regions)

    vds = vds.filter_samples_list(samples)
    if not vds.was_split():
        vds = vds.annotate_variants_expr('va.altAlleles = v.altAlleles')
        vds = vds.split_multi()

    vds = vds.filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0', keep=True)

    return vds


def compute_concordance(left_vds, right_vds, out_prefix, high_conf_regions = None, samples = None, overwrite = False,
                        left_name = 'left', right_name = 'right'):

    if samples is None:
        samples = list(set(left_vds.sample_ids).intersection(set(right_vds.sample_ids)))
        logger.debug("Found %d samples:\n%s\n..." % (len(samples), "\n".join(samples[:20])))
    elif not isinstance(samples, list):
        samples = [samples]

    left_vds = filter_for_concordance(left_vds, samples, high_conf_regions=high_conf_regions)
    right_vds = filter_for_concordance(right_vds, samples, high_conf_regions=high_conf_regions)

    global_concordance, s_concordance, v_concordance = left_vds.concordance(right=right_vds)  # TODO: make sure truth is already minrepped

    s_concordance.write(out_prefix + ".s_concordance.kt", overwrite=overwrite)

    vds = hail.VariantDataset.from_table(v_concordance)
    vds = vds.annotate_variants_vds(left_vds, root = 'va.%s' % left_name)
    vds = vds.annotate_variants_vds(right_vds, root='va.%s' % right_name)

    if len(samples) == 1:
        vds = (
            vds.annotate_global('global.gt_mappings', ["missing", "no_call", "homref", "het", "homvar"],
                                             TArray(TString()))
                .annotate_variants_expr('va.gt_arr = range(5).find(i => va.concordance[i].exists(x => x > 0))')
                .annotate_variants_expr(['va.%s_gt =  global.gt_mappings[va.gt_arr]' % left_name,
                                        'va.%s_gt = global.gt_mappings[range(5).find(i => va.concordance[va.gt_arr][i] > 0)]' % right_name,
                                         'va = drop(va, gt_arr)'])
        )

    vds = vds.annotate_variants_expr('va.variantType = if(v.altAllele.isSNP) "snv" '
                                                     'else if(v.altAllele.isIndel) "indel"'
                                                     'else "mixed"')

    vds.write(out_prefix + ".v_concordance.vds", overwrite=overwrite)


def export_concordance(vds, out_annotations, out_prefix, single_sample=True):

    if single_sample:
        vds = (vds.annotate_global('global.gt_mappings', ["missing", "no_call" ,"homref" ,"het" ,"homvar"], TArray(TString()))
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


def annotate_with_additional_rf_files(vds, rf_ann_files):
    """
    A helper function that takes a list of files in the format and add the RF annotation to the VDS

    :param VariantDataset vds: Input VDS
    :param rf_ann_files list of str: List of files in the format:  name|file_path|rf_root
    :return: Annotated VDS
    :rtype: VariantDataset
    """
    logger.info("Annotating with the following files: %s" % ",".join(rf_ann_files))

    out_metrics = []

    rf_ann_files = [re.compile("\\s*\\|\\s*").split(x) for x in rf_ann_files]

    for f in rf_ann_files:
        expr_dict = {'name': f[0], 'path': re.sub('^va\.', 'vds.', f[2])}
        out_metrics.extend([x % expr_dict for x in ['rfpred_%(name)s = va.rf_%(name)s.prediction',
                                                    'rfprob_%(name)s = va.rf_%(name)s.probability["TP"]',
                                                    'train_%(name)s = va.rf_%(name)s.train',
                                                    'label_%(name)s = va.rf_%(name)s.label']])
        vds = vds.annotate_variants_vds(vds.hc.read(f[1], drop_samples=True), expr=
        'va.rf_%(name)s.prediction = %(path)s.prediction,'
        'va.rf_%(name)s.probability = %(path)s.probability,'
        'va.rf_%(name)s.train = vds.train,'
        'va.rf_%(name)s.label = vds.label' % expr_dict)

    return vds, out_metrics

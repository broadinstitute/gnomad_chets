from utils import *
from collections import defaultdict
import argparse
from rf import *


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


def annotate_for_random_forests(vds, omni_vds, mills_vds, sample=True, fam_file=None):

    vds_schema = [f.name for f in vds.variant_schema.fields]

    if "tdt" not in vds_schema:
        if fam_file is not None:
            vds = vds.tdt(Pedigree.read(fam_file))

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


def main(args):

    hc = HailContext(log='/variantqc.log')

    if args.debug:
        logger.setLevel(logging.DEBUG)

    mendel_path = args.mendel_path if args.mendel_path else args.output
    rf_ann_path = args.rf_ann_path if args.rf_ann_path else args.output + ".annotated_for_rf.vds"
    rf_model_path = (args.rf_path if args.rf_path else args.output) + '.rf.model'
    rf_train_path = (args.rf_path if args.rf_path else args.output) + '.rf.training_sites.bgz'
    rf_path = args.rf_path if args.rf_path else args.output + '.rf.vds'
    rf_model = None
    vds = None

    if args.genomes:
        fam_path = genomes_fam_path
        hardcalls_path = full_genome_hardcalls_split_vds_path
        rf = hc.read(hardcalls_path, drop_samples=True)
    else:
        fam_path = exomes_fam_path
        hardcalls_path = full_exome_hardcalls_split_vds_path
        rf = (hc.read(hardcalls_path, drop_samples=True)
              .annotate_variants_vds(hc.read(vqsr_vds_path).split_multi(),
                                     expr='va.info.VQSLOD = vds.info.VQSLOD,'
                                          'va.info.POSITIVE_TRAIN_SITE = vds.info.POSITIVE_TRAIN_SITE,'
                                          'va.info.NEGATIVE_TRAIN_SITE = vds.info.NEGATIVE_TRAIN_SITE'))

    if args.annotate_for_rf:
        rf = rf.filter_variants_expr('va.calldata.qc_samples_raw.AC > 0')
        rf = rf.annotate_variants_table(hc.import_table(mendel_path + ".lmendel", impute=True).key_by('SNP'),
                                        expr='va.mendel = table.N')

        rf = annotate_for_random_forests(rf, hc.read(omni_vds_path), hc.read(mills_vds_path), fam_path)
        print(rf.variant_schema)

        rf.write(rf_ann_path, overwrite=args.overwrite)

    if args.train_rf:

        features = vqsr_features if args.vqsr_features else rf_features

        if args.vqsr_training:
            tp_criteria = "va.info.POSITIVE_TRAIN_SITE"
            fp_criteria = "va.info.NEGATIVE_TRAIN_SITE"
        else:
            tp_criteria = "va.omni || va.mills || va.info.POSITIVE_TRAIN_SITE"
            fp_criteria = "va.failing_hard_filters"
            if not args.no_transmitted_singletons:
                tp_criteria += " || va.transmitted_singleton"

        vds = sample_RF_training_examples(hc.read(rf_ann_path, drop_samples=True),
                                          fp_to_tp=args.fp_to_tp, tp_criteria=tp_criteria, fp_criteria=fp_criteria)

        if args.train_on_vqsr_sites:
            vds = vds.annotate_variants_expr(['va.train = isDefined(va.info.VQSR_NEGATIVE_TRAIN_SITE) || isDefined(va.info.VQSR_POSITIVE_TRAIN_SITE)',
                                              'va.label = ...'])  # TODO: wat.
        else:
            vds = vds.annotate_variants_expr(['va.train = va.train && '
                                              'va.calldata.qc_samples_raw.AC > 0 && '
                                              'v.contig != "22"'])

        rf_model = train_rf(vds, rf_features=features, num_trees=args.num_trees, max_depth=args.max_depth)
        save_model(rf_model, rf_model_path, overwrite=args.overwrite)
        vds.export_variants(rf_train_path, 'v, va.train')

    if args.apply_rf:

        if not rf_model:
            rf_model = load_model(rf_model_path)

        if not vds:
            vds = hc.read(rf_ann_path, drop_samples=True)
            vds = vds.annotate_variants_table(
                hc.import_table(rf_train_path,
                                no_header=True,
                                types={'f0': TVariant(), 'f1': TBoolean()})
                .key_by('f0'),
                expr='va.train = table')
            vds = vds.annotate_variants_expr('va.label = NA: String')

        features = vqsr_features if args.vqsr_features else rf_features

        vds = apply_rf_model(vds, rf_model, features)
        vds.write(rf_path, overwrite=args.overwrite)

    if args.write:
        # Output metrics for RF evaluation
        out_metrics = [
            'chrom = v.contig',
            'pos = v.start',
            'ref = v.ref',
            'alt = v.alt',
            'multi = va.wasSplit',
            'vqslod = va.info.VQSLOD',
            'pass = va.filters.contains("PASS")',
            'ti = v.altAllele.isTransition.toInt',
            'tv = v.altAllele.isTransversion.toInt',
            'ins = v.altAllele.isInsertion.toInt',
            'trans = va.tdt.nTransmitted',
            'untrans = va.tdt.nUntransmitted',
            'type = va.variantType',
            'qd = va.info.QD',
            'train_vqsr = va.info.NEGATIVE_TRAIN_SITE || va.info.POSITIVE_TRAIN_SITE',
            'label_vqsr = if(va.info.POSITIVE_TRAIN_SITE) "TP" else orMissing(isDefined(va.info.NEGATIVE_TRAIN_SITE), "FP")',
            'ac_origin = va.info.AC[va.aIndex]',
            'an_origin = va.info.AN',
            'ac_unrelated = va.AC_unrelated',
            'mendel_err = va.mendel',
            'qd2 = va.stats.qc_samples_raw.qd'
        ]

        rf_out = (
            hc.read(rf_ann_path, drop_samples=True)
            .filter_variants_table(KeyTable.import_interval_list(lcr_intervals_path), keep=False)
            .filter_variants_table(KeyTable.import_interval_list(decoy_intervals_path), keep=False)
            .filter_variants_expr('va.calldata.qc_samples_raw.AC > 0')
        )

        if args.add_default_rf:
            rf_out, additional_metrics = annotate_with_additional_rf_files(rf_out, ['current|{}|va.rf'.format(rf_path)])  # hack AF
            out_metrics.extend(additional_metrics)

        if args.rf_ann_files:
            rf_out, additional_metrics = annotate_with_additional_rf_files(rf_out, args.rf_ann_files)
            out_metrics.extend(additional_metrics)

        # Get number of SNVs and indels to sample
        nvariants = rf_out.query_variants(["variants.filter(x => x.altAllele.isSNP).count()",
                                           "variants.filter(x => x.altAllele.isIndel).count()"])

        logger.info("Number of SNVs: {}, number of Indels: {}".format(*nvariants))

        (
            rf_out
            .annotate_variants_table(hc.import_table(mendel_path + ".lmendel", impute=True).key_by('SNP'), expr='va.mendel = table.N')
            .filter_variants_expr('(v.altAllele.isSNP && pcoin(2500000.0 / {})) || '
                                  '(v.altAllele.isIndel && pcoin(2500000.0 / {})) ||'
                                  '(va.mendel > 0 && va.calldata.all_samples_raw.AC == 1)'.format(*nvariants))
            .export_variants(rf_path + ".va.txt.bgz", ",".join(out_metrics))
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--add_default_rf', help='A Konrad Konvenience', action='store_true')
    parser.add_argument('--rf_ann_files', help='RF files to annotate results with in pipe-delimited format: name|location|rf_root', nargs='+')

    actions = parser.add_argument_group('Actions')
    actions.add_argument('--compute_mendel', help='Computes Mendel errors', action='store_true')
    actions.add_argument('--annotate_for_rf', help='Creates an annotated VDS with features for RF', action='store_true')
    actions.add_argument('--train_rf', help='Trains RF model', action='store_true')
    actions.add_argument('--apply_rf', help='Applies RF model to the data', action='store_true')
    actions.add_argument('--write', help='Writes the results of the current RF model with the given name.', action='store_true')

    paths = parser.add_argument_group('Paths')
    paths.add_argument('--mendel_path', help='Overrides the default mendel path prefix ($output)')
    paths.add_argument('--rf_ann_path', help='Overrides the default rf annotation path ($output + .annotated_for_rf.vds)')
    paths.add_argument('--rf_model_path', help='Overrides the default rf model prefix ($output)')
    paths.add_argument('--rf_path', help='Overrides the default rf path ($output + .rf.vds)')

    rf_params = parser.add_argument_group('Random Forest parameters')
    rf_params.add_argument('--fp_to_tp', help='Sets to ratio of TPs to FPs for creating the RF model.', default=1.0, type=float)
    rf_params.add_argument('--num_trees', help='Number of trees in the RF model.', default=500)
    rf_params.add_argument('--max_depth', help='Maxmimum tree depth in the RF model.', default=5)

    training_params = parser.add_argument_group('Training data parameters')
    training_params.add_argument('--vqsr_features', help='Use VQSR features only (+ snv/indel variant type)', action='store_true')
    training_params.add_argument('--vqsr_training', help='Use VQSR training examples', action='store_true')
    training_params.add_argument('--no_transmitted_singletons', help='Do not use transmitted singletons for training.', action='store_true')
    training_params.add_argument('--train_on_vqsr_sites', help='Use VQSR training sites', action='store_true')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

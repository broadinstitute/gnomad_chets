from utils import *
import pyspark.sql
from pyspark.ml.feature import *
from pyspark.ml.classification import *
from pyspark.ml import *
from pyspark.sql import Row
from pyspark.sql import SparkSession
import pyspark


def run_rf_test(vds, output='/Users/laurent/tmp'):
    """
    Runs a dummy test RF on a given VDS:
    1. Creates variant annotations and labels to run model on
    2. Trains a RF pipeline model (including median imputation of missing values in created annotations)
    3. Saves the RF pipeline model
    4. Applies the model to the VDS and prints features importance

    :param VariantDataset vds: Input VDS
    :param str output: Output files prefix to save the RF model
    :return: VDS after applying RF model
    :rtype: VariantDataset
    """
    vds = vds.annotate_variants_expr(['va.train = pcoin(0.9), '
                                      'va.feature1 = pcoin(0.1)',
                                      'va.feature2 = rnorm(0.0, 1.0)',
                                      'va.feature3 = orMissing(pcoin(0.5), rnorm(0.0, 1.0))'])
    vds = vds.annotate_variants_expr('va.label = if(va.feature1 && va.feature2>0) "TP" else "FP"')
    rf_features = ['va.feature1', 'va.feature2', 'va.feature3']


    logger.info('Feature3 defined values before imputation: {}'.format(vds.query_variants(
        'variants.map(v => isDefined(va.feature3)).counter()')))
    logger.info('Feature3 median: {}'.format(vds.query_variants(
        'variants.map(v => va.feature3).collect().median()')))
    vds = impute_features_median(vds, rf_features)
    logger.info('Feature3 defined values after imputation: {}'.format(vds.query_variants(
        'variants.map(v => isDefined(va.feature3)).counter()')))

    rf_model = train_rf(vds, rf_features)
    save_model(rf_model, out=output + '/rf.model', overwrite=True)
    rf_model = load_model(output + '/rf.model')
    return apply_rf_model(vds, rf_model, rf_features)


def impute_features_median(vds, features, relative_error=0.01):
    """
    Imputes numeric fields with missing values using the approximate median of non-missing values.
    Non-numeric fields are ignored.
    If a column only has NAs, nothing is done.

    :param VariantDataset vds: input VDS
    :param list of str features: list of features to impute. If none given, all numerical features with missing data are imputed
    :param float relative_error: The relative error on the median approximation
    :return: VDS with missing values imputed with median
    :rtype: VariantDataset
    """

    schema = vds.variant_schema
    # Select numerical fields and create valid SQL name for all fields to avoid problem with dot-delimited paths
    SSQL_features = {f: toSSQL(f) for f in features if annotation_type_is_numeric(get_ann_type(f,schema))}

    logger.info("Imputing features {} with median.".format(",".join(SSQL_features.keys())))

    kt = vds.split_multi().variants_table()
    kt = kt.annotate(['{1} = {0}'.format(f,c) for f,c in SSQL_features.iteritems()])
    kt = kt.select([c for c in SSQL_features.values()])
    df = kt.to_dataframe()

    quantiles = {}
    for f,c in SSQL_features.iteritems():
        col_no_na = df.select(c).dropna()
        if col_no_na.first() is not None:
            quantiles[f] = col_no_na.approxQuantile(c, [0.5], relative_error)[0]

    vds = vds.annotate_variants_expr(['{0} = orElse({0},{1})'.format(f,quantiles[f]) for f in SSQL_features.keys()])
    return vds


def vds_to_rf_df(vds, rf_features, label='va.label'):
    """

    Creates a Dataframe ready to compute Random forest from a VDS.

    :param VariantDataset vds: Input VDS
    :param list of str rf_features: Features (variant annotations) to use to train the RF model
    :param str label: Variant annotation containing the label to predict
    :return: Dataframe with columns containing the features listed, the label and the variant
    :rtype: DataFrame
    """

    cols = rf_features + [label]
    kt = vds.split_multi().variants_table()

    # Rename everything to avoid problem with dot-delimited paths
    kt = kt.annotate(['%s = %s' % (toSSQL(x), x) for x in cols] +
                     ['variant = str(v)'])
    kt = kt.select([toSSQL(x) for x in cols] + ['variant'])

    # Create dataframe
    # 1) drop rows with missing values (not supported for RF)
    # 2) replace missing labels with standard value since StringIndexer doesn't handle missing values
    df = kt.to_dataframe().dropna(subset=[toSSQL(x) for x in rf_features]).fillna('NA', subset=[toSSQL(label)])

    return df


def get_features_importance(rf_pipeline, rf_index=-2, assembler_index=-3):
    """
    Extract the features importance from a Pipeline model containing a RandomForestClassifier stage.

    :param PipelineModel rf_pipeline: Input pipeline
    :param int rf_index: index of the RandomForestClassifier stage
    :param int assembler_index: index of the VectorAssembler stage
    :return: feature importance for each feature in the RF model
    :rtype: dict of str: float
    """

    feature_names = [x[:-len("_indexed")] if x.endswith("_indexed") else x for x in
                     rf_pipeline.stages[assembler_index].getInputCols()]
    feature_importance = {fromSSQL(new_name): importance for
                          (new_name, importance) in zip(feature_names, rf_pipeline.stages[rf_index].featureImportances)}
    return feature_importance


def get_labels(rf_pipeline):
    """
    Returns the labels from the StringIndexer stage at index 0 from an RF pipeline model

    :param PipelineModel rf_pipeline: Input pipeline
    :return: labels
    :rtype: list of str
    """
    return rf_pipeline.stages[0].labels


def apply_rf_model(vds, rf_pipeline, rf_features, label='va.label', va_root='va.rf', globals_root='globals.rf'):
    """
    Applies a Random Forest (RF) pipeline model to a VDS and annotate the RF probabilities and predictions as a variant annotation.

    :param VariantDataset vds: Input VDS
    :param PipelineModel rf_pipeline: Random Forest pipeline model
    :param list of str rf_features: List of features in the pipeline. !Should match the model list of features!
    :param str label: Variant annotation containing the labels. !Should match the model labels!
    :param str va_root: Root of output for va RF annotation (predictions and probabilities)
    :param str globals_root: Root of output for global RF annotation (features importance)
    :return: VDS with RF annotations
    :rtype: VariantDataset
    """

    logger.info("Applying RF model to VDS")

    df = vds_to_rf_df(vds, rf_features, label=label)

    feature_importance = get_features_importance(rf_pipeline)

    transformed = rf_pipeline.transform(df)

    logger.info("Annotating dataset with results")

    # Required for RDD.toDF() !
    spark = SparkSession(vds.hc.sc)

    kt = KeyTable.from_dataframe(
        transformed.rdd.map(
            lambda row:
            Row(variant=row['variant'],
                probability=row["probability"].toArray().tolist(),
                prediction=row["predictedLabel"])
        ).toDF()
    ).persist()

    probability_to_dict_expr = 'probability = index([{%s}], label).mapValues(x => x.prob)' % "},{".join(
        ['label: "%s", prob: probability[%d]' % (l, i) for (i, l) in enumerate(get_labels(rf_pipeline))])

    kt = kt.annotate(['variant = Variant(variant)',
                      probability_to_dict_expr]).key_by('variant')

    vds = vds.annotate_variants_table(kt,
                                      expr="%s.prediction = table.prediction, %s.probability = table.probability" % (
                                          va_root, va_root))
    vds = vds.annotate_global(globals_root, feature_importance, TDict(TString(), TDouble()))

    return vds


def save_model(rf_pipeline, out, overwrite=False):
    """
    Saves a Random Forest pipeline model.

    :param PipelineModel rf_pipeline: Pipeline to save
    :param str out: Output path
    :param bool overwrite: If set, will overwrite existing file(s) at output location
    :return: Nothing
    :rtype: NoneType
    """
    logger.info("Saving model to %s" % out)
    if overwrite:
        rf_pipeline.write().overwrite().save(out)
    else:
        rf_pipeline.save(out)


def load_model(input_path):
    """
    Loads a Random Forest pipeline model.

    :param str input_path: Location of model to load
    :return: Random Forest pipeline model
    :rtype: PipelineModel
    """
    logger.info("Loading model from %s" % input_path)
    return pyspark.ml.PipelineModel.load(input_path)


def train_rf(vds, rf_features, training='va.train', label='va.label', num_trees=500, max_depth=5):
    """
    Trains a Random Forest (RF) pipeline model.

    :param VariantDataset vds: Input VDS
    :param list of str rf_features: List of VDS annotations to be used as features
    :param str training: Boolean VDS annotation indicating which sites should be used for training
    :param str label: Annotation containing the label to predict
    :param int num_trees: Number of trees to use
    :param int max_depth: Maximum tree depth
    :return: Random Forest pipeline model
    :rtype: PipelineModel
    """
    logger.info("Training RF model using:\n"
                "features: %s\n"
                "training: %s\n"
                "labels: %s\n"
                "num_trees: %d\n"
                "max_depth: %d" % (",".join(rf_features), training, label, num_trees, max_depth))

    df = vds_to_rf_df(vds, rf_features + [training], label=label)
    df = df.drop('variant')

    SSQL_training = toSSQL(training)
    SSQL_label = toSSQL(label)

    label_indexer = StringIndexer(inputCol=SSQL_label, outputCol=SSQL_label + "_indexed").fit(df)
    labels = label_indexer.labels
    logger.info("Found labels: %s" % labels)

    string_features = [x[0] for x in df.dtypes if x[0] != SSQL_label and x[0] != SSQL_training and x[1] == 'string']
    if string_features:
        logger.info("Indexing string features: %s", ",".join(string_features))
    string_features_indexers = [StringIndexer(inputCol=x, outputCol=x + "_indexed").fit(df)
                                for x in string_features]

    assembler = VectorAssembler(inputCols=[x[0] + "_indexed" if x[1] == 'string' else x[0]
                                           for x in df.dtypes if x[0] != SSQL_label and x[0] != SSQL_training],
                                outputCol="features")

    rf = RandomForestClassifier(labelCol=SSQL_label + "_indexed", featuresCol="features",
                                maxDepth=max_depth, numTrees=num_trees)

    label_converter = IndexToString(inputCol='prediction', outputCol='predictedLabel', labels=labels)

    pipeline = Pipeline(stages=[label_indexer] + string_features_indexers +
                               [assembler, rf, label_converter])

    #Train model
    logger.info("Training RF model")
    training_df = df.filter(SSQL_training).drop(SSQL_training)
    rf_model = pipeline.fit(training_df)

    feature_importance = get_features_importance(rf_model)

    logger.info(
        "RF features importance:\n%s" % "\n".join(["%s: %s" % (f, i) for (f, i) in feature_importance.iteritems()]))

    return rf_model

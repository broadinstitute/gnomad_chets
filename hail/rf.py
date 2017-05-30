from utils import *
from hail.expr import *
from pyspark.ml.feature import *
from pyspark.ml.classification import *
from pyspark.ml import *
from pyspark.sql import Row
from pyspark.sql import SparkSession
import pyspark
from pyspark import SparkContext

def run_rf_test(vds, output = '/Users/laurent/tmp'):
    vds = vds.annotate_variants_expr('va.train = pcoin(0.9), va.feature1 = pcoin(0.1), va.feature2 = rnorm(0.0, 1.0)')
    vds = vds.annotate_variants_expr('va.label = if(va.feature1 && va.feature2>0) "TP" else "FP"')
    rf_features = ['va.feature1','va.feature2']

    rf_model = train_rf(vds, rf_features)
    save_model(rf_model, out = output +'/rf.model', overwrite=True)
    rf_model = load_model(output + '/rf.model')
    return apply_rf_model(vds, rf_model, rf_features)


def df_type_is_numeric(t):
    return (isinstance(t, pyspark.sql.types.DoubleType) or
            isinstance(t, pyspark.sql.types.DecimalType) or
            isinstance(t, pyspark.sql.types.FloatType) or
            isinstance(t, pyspark.sql.types.DoubleType) or
            isinstance(t, pyspark.sql.types.ByteType) or
            isinstance(t, pyspark.sql.types.IntegerType) or
            isinstance(t, pyspark.sql.types.LongType) or
            isinstance(t, pyspark.sql.types.ShortType)
            )


def impute_features_median(df):
    0 #TODO


#Replaces `.` with `_`, since Spark ML doesn't support column names with `.`
def toSSQL(str):
    return str.replace('.','_')


def vds_to_rf_df(vds, rf_features, label='va.label'):

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


def get_features_importance(rf_model, rf_index = -2, assembler_index = -3):
    feature_names = [x[:-len("_indexed")] if x.endswith("_indexed") else x for x in rf_model.stages[assembler_index].getInputCols()]
    feature_importance = {toSSQL(new_name): importance for
                          (new_name, importance) in zip(feature_names, rf_model.stages[rf_index].featureImportances)}
    return feature_importance


def get_labels(rf_model):
    return rf_model.stages[0].labels


def apply_rf_model(vds, rf_model, rf_features, root='va.rf', label='va.label'):
    logger.info("Applying RF model to VDS")

    df = vds_to_rf_df(vds, rf_features, label=label)

    feature_importance = get_features_importance(rf_model)

    transformed = rf_model.transform(df)

    logger.info("Annotating dataset with results")

    # Required for RDD.toDF() !
    spark = SparkSession(vds.hc.sc)

    kt = hail.KeyTable.from_dataframe(
        transformed.rdd.map(
            lambda row:
            Row(variant=row['variant'],
                probability=row["probability"].toArray().tolist(),
                prediction=row["predictedLabel"])
        ).toDF()
    ).persist()

    probability_to_dict_expr = 'probability = index([{%s}], label).mapValues(x => x.prob)' % "},{".join(
        ['label: "%s", prob: probability[%d]' % (l, i) for (i, l) in enumerate(get_labels(rf_model))])

    kt = kt.annotate(['variant = Variant(variant)',
                      probability_to_dict_expr]).key_by('variant')

    vds = vds.annotate_variants_table(kt, expr = "%s.prediction = table.prediction, %s.probability = table.probability" % (
    root, root))
    vds = vds.annotate_global('global.%s' % (root[3:]), feature_importance, TDict(TString(), TDouble()))

    return vds


def save_model(rf_model, out, overwrite = False):
    logger.info("Saving model to %s" % out)
    if overwrite:
        rf_model.write().overwrite().save(out)
    else:
        rf_model.save(out)


def load_model(input):
    logger.info("Loading model from %s" % input)
    return pyspark.ml.PipelineModel.load(input)


def train_rf(vds, rf_features, training='va.train', label='va.label', num_trees=500, max_depth=5):
    logger.info("Training RF model using:\n"
                "features: %s\n"
                "training: %s\n"
                "labels: %s\n"
                "num_trees: %d\n"
                "max_depth: %d" %( ",".join(rf_features), training, label, num_trees, max_depth ))

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
    string_features_indexers = [StringIndexer(inputCol= x, outputCol= x + "_indexed").fit(df)
                                 for x in string_features]

    assembler = VectorAssembler(inputCols= [ x[0] + "_indexed" if x[1] == 'string' else x[0]
                                             for x in df.dtypes if x[0] != SSQL_label and x[0] != SSQL_training],
                                outputCol="features")

    rf = RandomForestClassifier(labelCol=SSQL_label + "_indexed", featuresCol="features",
                                maxDepth=max_depth, numTrees=num_trees)

    label_converter = IndexToString(inputCol='prediction', outputCol='predictedLabel', labels=labels)

    pipeline = Pipeline(stages = [label_indexer] + string_features_indexers +
                                 [assembler, rf, label_converter])

    #rTain model on training sites
    logger.info("Training RF model")
    training_df = df.filter(SSQL_training).drop(SSQL_training)
    rf_model = pipeline.fit(training_df)

    feature_importance = get_features_importance(rf_model)

    logger.info("RF features importance:\n%s" % "\n".join(["%s: %s" % (f,i) for (f,i) in feature_importance.iteritems()]))

    return rf_model

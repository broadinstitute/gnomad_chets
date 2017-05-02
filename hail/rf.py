from utils import *
from hail.expr import *
from pyspark.ml.feature import *
from pyspark.ml.classification import *
from pyspark.ml import *
from pyspark.sql import Row
from pyspark.sql import SparkSession
import pyspark
from pyspark import SparkContext

TRAIN = 'train'
LABEL = 'label'

def run_rf_test(hc, vds):
    vds = vds.annotate_variants_expr('va.train = pcoin(0.9), va.feature1 = pcoin(0.1), va.feature2 = rnorm(0.0, 1.0)')
    vds = vds.annotate_variants_expr('va.label = if(va.feature1 && va.feature2>0) "TP" else "FP"')

    return run_rf(hc, vds, ['va.feature1','va.feature2'], out = '/Users/laurent/tmp/rf.model')

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
    0


def vds_to_rf_df(vds, rf_features, training='va.train', label='va.label'):

    kt = vds.split_multi().variants_keytable()

    # Rename everything to avoid problem with dot-delimited paths
    new_to_old_features = {'f%d' % i: old for (i, old) in enumerate(rf_features)}
    kt = kt.annotate(['%s = %s' % (new, old) for (new, old) in new_to_old_features.iteritems()] +
                     ['%s = %s' % (TRAIN, training), '%s = %s' % (LABEL, label), 'variant = str(v)'])
    kt = kt.select(new_to_old_features.keys() + ['variant', TRAIN, LABEL])

    # Create dataframe
    # 1) drop rows with missing values (not supported for RF)
    # 2) replace missing labels with standard value since StringIndexer doesn't handle missing values
    df = kt.to_dataframe().dropna(subset=new_to_old_features.keys()).fillna('NA', subset=[LABEL])

    return df, new_to_old_features

def run_rf(vds, rf_features, training='va.train', label='va.label', root='va.rf', num_trees=500, max_depth=5, out = None, overwrite = False):

    df, new_to_old_features = vds_to_rf_df(vds, rf_features, training=training, label=label)

    #Select training sites
    training_df = df.filter(TRAIN).drop(TRAIN).drop('variant')

    label_indexer = StringIndexer(inputCol=LABEL, outputCol=LABEL + "_indexed").fit(df)
    labels = label_indexer.labels

    logger.info("Found labels: %s" % labels)

    string_features_indexers = [ StringIndexer(inputCol= x[0], outputCol= x[0] + "_indexed")
                                 for x in training_df.dtypes if x[0] != LABEL and x[1] == 'string' ]

    assembler = VectorAssembler(inputCols= [ x[0] + "_indexed" if x[1] == 'string' else x[0]
                                             for x in training_df.dtypes if x[0] != LABEL ],
                                outputCol="features")

    features_indexer = VectorIndexer(inputCol="features", outputCol="features_indexed", maxCategories=10)

    rf = RandomForestClassifier(labelCol=LABEL + "_indexed", featuresCol="features_indexed",
                                maxDepth=max_depth, numTrees=num_trees)

    label_converter = IndexToString(inputCol='prediction', outputCol='predictedLabel', labels=labels)

    pipeline = Pipeline(stages = [label_indexer] + string_features_indexers +
                                 [assembler, features_indexer, rf, label_converter])

    logger.info("Training RF model")

    rf_model = pipeline.fit(training_df)

    feature_names = [ x[:-len("_indexed")] if x.endswith("_indexed") else x for x in assembler.getInputCols()]
    feature_importance = { new_to_old_features[new_name]: importance for
                           (new_name, importance) in zip(feature_names, rf_model.stages[-2].featureImportances) }

    logger.info("RF features importance:\n%s" % "\n".join(["%s: %s" % (f,i) for (f,i) in feature_importance.iteritems()]))

    if out:
        logger.info("Saving model to %s" % out)
        if overwrite:
            rf_model.write().overwrite().save(out)
        else:
            rf_model.save(out)

    logger.info("Applying RF model to entire data")

    transformed = rf_model.transform(df)

    logger.info("Annotating dataset with results")

    #Required for RDD.toDF() !
    spark = SparkSession(vds.hc.sc)

    kt = vds.hc.dataframe_to_keytable(
        transformed.rdd.map(
            lambda row:
            Row(variant=row['variant'],
                probability=row["probability"].toArray().tolist(),
                prediction=row["predictedLabel"])
        ).toDF()
    )

    probability_to_dict_expr = 'probability = index([{%s}], label).mapValues(x => x.prob)' % "},{".join(
        ['label: "%s", prob: probability[%d]' %(l,i) for (i,l) in enumerate(labels)])

    kt = kt.annotate(['variant = Variant(variant)',
                      probability_to_dict_expr]).key_by('variant')

    vds = vds.annotate_variants_keytable(kt, "%s.prediction = table.prediction, %s.probability = table.probability" % (root, root))
    vds = vds.annotate_global_py('global.%s' % (root[3:]), feature_importance, TDict(TString(), TDouble()))

    return vds














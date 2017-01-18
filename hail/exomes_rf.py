from pyhail import *

hc = HailContext('/exome.rf.log')
rf_vds = hc.read('gs://exac2/exacv2_annotated.vds')
features = ['va.info.QD',
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
(
    rf_vds.annotate_variants_expr(
        'va.stats.raw.nrq_median = va.stats.raw.nrq_median[va.aIndex - 1],'
        'va.stats.raw.ab_median = va.stats.raw.ab_median[va.aIndex - 1],'
        'va.stats.raw.dp_median = va.stats.raw.dp_median[va.aIndex - 1],'
        'va.stats.raw.gq_median = va.stats.raw.gq_median[va.aIndex - 1]')
        .random_forests(training='va.train', label='va.label', root='va.rf', features=features, num_trees=500,
                          max_depth=5)
        .write('gs://exac2/exacv2_annotated.rf.vds')
)

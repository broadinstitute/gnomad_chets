
from variantqc import *
import time
from resources import *

try:
    hc
except NameError:
    hc = HailContext(log='/variantqc.log')

#Inputs
gnomad_path = "gs://gnomad/gnom.ad.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
fam_path = "gs://gnomad/gnomad.final.goodTrios.fam"
qcsamples_path = "gs://gnomad/gnomad.qcsamples.txt"

#Outputs
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"
raw_hardcalls_split_path = "gs://gnomad/gnomad.raw_hardcalls.split.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats12.vds"
mendel_path = "gs://gnomad/gnomad.raw_calls"
date_time = time.strftime("%Y-%m-%d_%H-%M")
tmp_vds = "gs://gnomad-lfran/temp." + date_time + ".hardcalls.vds"
tmp_vds2 = "gs://gnomad-lfran/temp." + date_time + ".rf.vds"

#Actions
create_hardcalls_vds=False
compute_mendel = False
run_rf=True
write_results=True

#Create hardcalls file with raw annotations
if(create_hardcalls_vds):

    variant_annotations = get_variant_type_expr()

    allele_annotations = get_stats_expr("va.stats.raw", medians=True)
    allele_annotations.append("va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.nNonRefAlleles).sum()")

    hardcalls_vds = (
        hc.read(gnomad_path)
        .annotate_samples_fam(fam_path)
        .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
        .annotate_variants_expr(variant_annotations)
        .annotate_variants_expr("va.calldata.raw = gs.callStats(g => v) ")
        .annotate_variants_expr("va.calldata.qc_samples_raw = gs.filter(g => sa.meta.qc_sample).callStats(g => v) ")
        .annotate_variants_expr("va.calldata.release_samples_raw = gs.filter(g => sa.meta.keep).callStats(g => v) ")
        .annotate_alleles_expr(allele_annotations)
        .hardcalls()
        .write(tmp_vds, overwrite=True)
    )

    hardcalls_vds = (
         hc.read(tmp_vds)
         .repartition(npartition=2500, shuffle=False)
         .write(raw_hardcalls_path)
     )

    hapmap = hc.read(hapmap_path)
    mills = hc.read(mills_path)
    omni = hc.read(omni_path)

    hardcalls_split = (
        hc.read(raw_hardcalls_path)
        .split_multi()
        .annotate_variants_vds(hapmap, code='va.hapmap = isDefined(vds)')
        .annotate_variants_vds(omni, code='va.omni = isDefined(vds)')
        .annotate_variants_vds(mills, code='va.mills = isDefined(vds)')
        .tdt(fam=fam_path)
        .write(raw_hardcalls_split_path,overwrite=True)
    )

if(compute_mendel):
    (
        hc.read(raw_hardcalls_split_path)
        .mendel_errors(mendel_path,fam=fam_path)
    )

#Run random forests
if(run_rf):

    features1 = ['va.variantType',
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
    features2 = ['va.variantType',
                'va.info.MQ',
                'va.info.MQRankSum',
                'va.info.SOR',
                'va.info.InbreedingCoeff',
                'va.info.ReadPosRankSum',
                'va.stats.raw.nrq_median',
                'va.stats.raw.ab_median',
                'va.stats.raw.dp_median',
                'va.stats.raw.gq_median']

    features_for_median = [
        'va.info.MQRankSum',
        'va.info.ReadPosRankSum',
        'va.stats.raw.ab_median'
    ]

    variant_types = ['snv', 'multi-snv', 'indel', 'multi-indel', 'mixed']
    global_features_expr_median = []
    for feature in features_for_median:
        for variant_type in variant_types:
            global_features_expr_median.append(
                'global.median.`%s`.`%s` = variants'
                '.filter(v => va.variantType == "%s" && pcoin([1.0,1000000 / global.variantsByType[va.variantType].count].min))'
                '.map(v => %s).collect().median()\n' % (
                    feature,variant_type, variant_type,feature))

    x = []
    for variant_type in variant_types:
        x.append('{variantType : "%s", ' % variant_type
                 + ",\n".join(
            ['`%s` : global.median.`%s`.`%s`' % (feature, feature, variant_type) for feature in features_for_median]) +
                 "}")

    global_features_expr_median2 = 'global.median = index([ ' + ",\n".join(x) + '], variantType)'

    variants_features_imputation = ['%s = if(isDefined(%s)) %s else global.median[va.variantType].`%s`'
                                    % (feature, feature, feature, feature) for feature in features_for_median]

    global_missing_features_expr = ['global.missing_%s = variants.filter(x => isMissing(%s)).count()' % (a, a) for a in
                                    features1]
    rf = (
        hc.read(raw_hardcalls_split_path, sites_only=True)
            .filter_variants_expr('va.calldata.qc_samples_raw.AC[va.aIndex] > 0')
            .annotate_variants_expr(
            'va.transmitted_singleton = va.tdt.nTransmitted == 1 && va.calldata.raw.AC[va.aIndex]==2,'
            'va.stats.raw.nrq_median = va.stats.raw.nrq_median[va.aIndex - 1],'
            'va.stats.raw.ab_median = va.stats.raw.ab_median[va.aIndex - 1],'
            'va.stats.raw.dp_median = va.stats.raw.dp_median[va.aIndex - 1],'
            'va.stats.raw.gq_median = va.stats.raw.gq_median[va.aIndex - 1]')
            .annotate_global_expr_by_variant('global.variantsByType = index(variants.map(v => va.variantType).counter(),key)')
            .annotate_global_expr_by_variant(global_missing_features_expr)
            .annotate_global_expr_by_variant(global_features_expr_median)
            .annotate_global_expr_by_variant(global_features_expr_median2)
    )
    rf.show_globals()

    rf = (
        annotate_for_random_forests(
            rf.annotate_variants_expr(variants_features_imputation)

        )
            .random_forests(training='va.train', label='va.label', root='va.RF1', features=features1, num_trees=500,
                            max_depth=5)
            .write(rf_path)
    )

    # rf = (
    #     hc.read(tmp_vds2)
    #     .random_forests(training='va.train', label='va.label', root='va.RF2', features=features2, num_trees=500, max_depth=5)
    #     .write(rf_path, overwrite=True)
    # )

if(write_results):
    #Output metrics for RF evaluation
    out_metrics = [
        'chrom = v.contig',
        'pos = v.start',
        'ref = v.ref',
        'alt = v.alt',
        'multi = va.wasSplit',
        'vqslod = va.info.VQSLOD',
        'ti = v.altAllele.isTransition.toInt',
        'tv = v.altAllele.isTransversion.toInt',
        'ins = v.altAllele.isInsertion.toInt',
        'trans = va.tdt.nTransmitted',
        'untrans = va.tdt.nUntransmitted',
        'type = va.variantType',
        'training = va.train',
        'label = va.label',
        'rfpred1 = va.RF1.prediction',
        'rfprob1 = va.RF1.probability["TP"]',
        'qd = va.info.QD',
        'negtrain = va.info.NEGATIVE_TRAIN_SITE',
        'postrain = va.info.POSITIVE_TRAIN_SITE',
        'ac_origin = va.info.AC[va.aIndex-1]',
        'an_origin = va.info.AN',
        'ac_unrelated = va.AC_unrelated[va.aIndex-1]',
        'mendel_err = va.mendel'
    ]

    (
        hc.read(rf_path)
        .filter_variants_intervals(lcr_path, keep=False)
        .filter_variants_intervals(decoy_path, keep=False)
        .annotate_variants_table(mendel_path + ".lmendel",'SNP',code='va.mendel = table.N',config=hail.TextTableConfig(impute=True))
        .filter_variants_expr('v.altAllele.isSNP && pcoin(0.92) && !(va.mendel>0 && va.calldata.raw.AC[va.aIndex] == 1 )',keep=False)
        .filter_variants_expr('!v.altAllele.isSNP && pcoin(0.6) && !(va.mendel>0 && va.calldata.raw.AC[va.aIndex] == 1 )',keep=False)
        .export_variants(rf_path + ".va.txt.bgz", ",".join(out_metrics))
    )




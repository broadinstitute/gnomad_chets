
from variantqc import *
import time
from resources import *
import re

try:
    hc
except NameError:
    hc = HailContext(log='/variantqc.log')

#Inputs
gnomad_path = "gs://gnomad/gnom.ad.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
fam_path = "gs://gnomad/gnomad.final.goodTrios.fam"
qcsamples_path = "gs://gnomad/gnomad.qcsamples.txt"
iteration_re = r"newStats(\d+)"
other_rf_ann_files = ["gs://gnomad/RF/gnomad.sites.RF.newStats17.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats18.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats19.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats20.vds"]

#Outputs
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"
raw_hardcalls_split_path = "gs://gnomad/gnomad.raw_hardcalls.split.vds"
rf_ann_path = "gs://gnomad/RF/gnomad.sites.annotated_for_RF_unbalanced.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats23.vds"
mendel_path = "gs://gnomad/gnomad.raw_calls"
date_time = time.strftime("%Y-%m-%d_%H-%M")
#date_time = "2017-01-31_16-47"
tmp_vds = "gs://gnomad-lfran/temp." + date_time + ".hardcalls.vds"
tmp_vds2 = "gs://gnomad-lfran/temp." + date_time + ".rf.vds"

#Actions
create_hardcalls_vds=False
compute_mendel = False
annotate_for_rf=False
run_rf=True
write_results=True

#Create hardcalls file with raw annotations
if(create_hardcalls_vds):

    variant_annotations = get_variant_type_expr()

    allele_annotations = get_stats_expr("va.stats.raw", medians=True)
    allele_annotations.extend(get_stats_expr("va.stats.qc_samples_raw", medians=True, samples_filter_expr='sa.meta.qc_sample'))
    allele_annotations.extend(get_stats_expr("va.stats.release_samples_raw", medians=True, samples_filter_expr='sa.meta.keep'))
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
        .annotate_variants_expr(['va.nAltAlleles = v.altAlleles.filter(a => a.alt != "*").count()'])
        .hardcalls()
        .write(tmp_vds, overwrite=True)
    )

    hardcalls_vds = (
         hc.read(tmp_vds)
         .repartition(num_partitions=4000, shuffle=False)
         .write(raw_hardcalls_path)
     )

    hapmap = hc.read(hapmap_path)
    mills = hc.read(mills_path)
    omni = hc.read(omni_path)

    hardcalls_split = (
        hc.read(raw_hardcalls_path)
            .annotate_variants_expr(['va.nonsplit_alleles = v.altAlleles.map(a => a.alt)'])
            .split_multi()
            .annotate_variants_expr([
            'va.hasStar = va.nonsplit_alleles.exists(a => a == "*")',
            'va.wasMixed = va.variantType == "mixed"',
            'va.alleleType = if(v.altAllele.isSNP) "snv"'
            '   else if(v.altAllele.isInsertion) "ins"'
            '   else if(v.altAllele.isDeletion) "del"'
            '   else "complex"'])
            .annotate_variants_vds(hapmap, code='va.hapmap = isDefined(vds)')
            .annotate_variants_vds(omni, code='va.omni = isDefined(vds)')
            .annotate_variants_vds(mills, code='va.mills = isDefined(vds)')
            .tdt(fam=fam_path)
            .write(raw_hardcalls_split_path, overwrite=True)
    )

if(compute_mendel):
    (
        hc.read(raw_hardcalls_split_path)
        .mendel_errors(mendel_path,fam=fam_path)
    )

#Random forests
if(annotate_for_rf):
    rf = (
        hc.read(raw_hardcalls_split_path, sites_only=True)
            .filter_variants_expr('va.calldata.qc_samples_raw.AC[va.aIndex] > 0')
            .annotate_variants_expr([
            'va.stats.qc_samples_raw.nrq_median = va.stats.qc_samples_raw.nrq_median[va.aIndex - 1]',
            'va.stats.qc_samples_raw.ab_median = va.stats.qc_samples_raw.ab_median[va.aIndex - 1]',
            'va.stats.qc_samples_raw.dp_median = va.stats.qc_samples_raw.dp_median[va.aIndex - 1]',
            'va.stats.qc_samples_raw.gq_median = va.stats.qc_samples_raw.gq_median[va.aIndex - 1]'])
    )
    rf = rf.annotate_variants_table(mendel_path + ".lmendel", 'SNP', code='va.mendel = table.N',
                                    config=hail.TextTableConfig(impute=True))
    # rf.show_globals()

    rf = annotate_for_random_forests(rf, hc.read(omni_path), hc.read(mills_path), balance=False)
    print(rf.variant_schema)

    rf.write(rf_ann_path)

if(run_rf):

    print("Starting RF with features: \n")
    pprint(rf_features)

    vds = sample_RF_training_examples(hc.read(rf_ann_path, sites_only=True),
                                      fp_to_tp = 0.15)

    (
        vds
        .annotate_variants_expr(['va.train = va.train && va.calldata.qc_samples_raw.AC[va.aIndex] > 0'])
        .random_forests(training='va.train', label='va.label', root='va.RF1', features=rf_features, num_trees=500,
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
    out_num = re.search(iteration_re,rf_path).groups()[0]
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
        'training = va.train',
        'label = va.label',
        'rfpred%s = va.RF1.prediction' % out_num,
        'rfprob%s = va.RF1.probability["TP"]' % out_num,
        'qd = va.info.QD',
        'negtrain = va.info.NEGATIVE_TRAIN_SITE',
        'postrain = va.info.POSITIVE_TRAIN_SITE',
        'ac_origin = va.info.AC[va.aIndex-1]',
        'an_origin = va.info.AN',
        'ac_unrelated = va.AC_unrelated[va.aIndex-1]',
        'mendel_err = va.mendel'
    ]

    rf_out = (
        hc.read(rf_path, sites_only=True)
        .filter_variants_intervals(lcr_path, keep=False)
        .filter_variants_intervals(decoy_path, keep=False)
        .filter_variants_expr('va.calldata.qc_samples_raw.AC[va.aIndex] > 0')
    )

    for f in other_rf_ann_files:
        out_num = re.search(iteration_re,f).groups()[0]
        out_metrics.append('rfpred%s = va.RF%s.prediction,rfprob%s = va.RF%s.probability["TP"]' % (out_num,out_num,out_num,out_num))
        rf_out = rf_out.annotate_variants_vds(hc.read(f, sites_only=True), code=
                                              'va.RF%s.prediction = vds.RF1.prediction,'
                                              'va.RF%s.probability = vds.RF1.probability' % (out_num,out_num))

    #Get number of SNVs and indels to sample
    nvariants = rf_out.query_variants(["variants.filter(x => x.altAllele.isSNP).count()",
                           "variants.filter(x => x.altAllele.isIndel).count()"])
    print(nvariants)

    (
        rf_out
        .annotate_variants_table(mendel_path + ".lmendel",'SNP',code='va.mendel = table.N',config=hail.TextTableConfig(impute=True))
        .filter_variants_expr('(v.altAllele.isSNP && pcoin(2500000.0 / %d) )|| '
                              '(v.altAllele.isIndel && pcoin(2500000.0 / %d) )||'
                              '(va.mendel>0 && va.calldata.raw.AC[va.aIndex] == 1 )' % (nvariants[0],nvariants[1]))
        .export_variants(rf_path + ".va.txt.bgz", ",".join(out_metrics))
    )

#.filter_variants_expr('pcoin([1.0,1000000 / global.variantsByType[va.variantType]].min) || (va.mendel>0 && va.calldata.qc_samples_raw.AC[va.aIndex] == 1 )')


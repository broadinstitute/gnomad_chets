
from variantqc import *
from resources import *
import re
import argparse
import rf

#Inputs
iteration_re = r"newStats(\d+)"
other_rf_ann_files = ["gs://gnomad/RF/gnomad.sites.RF.newStats17.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats18.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats19.vds",
                      "gs://gnomad/RF/gnomad.sites.RF.newStats20.vds"]


def main(args):

    hc = hail.HailContext(log='/variantqc.log')

    if args.debug:
        logger.setLevel(logging.DEBUG)

    input = args.input if args.input else args.output
    raw_hardcalls_path = input + ".raw_hardcalls.vds"
    raw_hardcalls_split_path = input + ".raw_hardcalls.split.vds"
    mendel_path = input
    rf_ann_path = input + ".annotated_for_rf.vds"
    rf_path = input + '.rf.vds'


    #Create hardcalls file with raw annotations
    if args.write_hardcalls:

        variant_annotations = get_variant_type_expr()

        allele_annotations = get_stats_expr("va.stats.all_samples_raw", medians=True)
        allele_annotations.extend(get_stats_expr("va.stats.qc_samples_raw", medians=True, samples_filter_expr='sa.meta.qc_sample'))
        allele_annotations.extend(get_stats_expr("va.stats.release_samples_raw", medians=True, samples_filter_expr='sa.meta.keep'))
        allele_annotations.append("va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.nNonRefAlleles).sum()")

        (
            hc.read(full_genome_vds)
            .annotate_samples_fam(genomes_fam)
            .annotate_samples_table(genomes_meta, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
            .annotate_variants_expr(variant_annotations)
            .annotate_variants_expr("va.calldata.raw = gs.callStats(g => v) ")
            .annotate_variants_expr("va.calldata.qc_samples_raw = gs.filter(g => sa.meta.qc_sample).callStats(g => v) ")
            .annotate_variants_expr("va.calldata.release_samples_raw = gs.filter(g => sa.meta.keep).callStats(g => v) ")
            .annotate_alleles_expr(allele_annotations)
            .annotate_variants_expr(['va.nAltAlleles = v.altAlleles.filter(a => !a.isStar).length'])
            .hardcalls()
            .write(args.output + '.tmp.raw_hardcalls.vds', overwrite=True)
        )

        (
             hc.read(args.output + '.tmp.raw_hardcalls.vds')
             .min_rep()
             .repartition(num_partitions=4000)
             .write(raw_hardcalls_path, overwrite=args.overwrite)
         )

        hapmap = hc.read(hapmap_path)
        mills = hc.read(mills_path)
        omni = hc.read(omni_path)

        (
            hc.read(raw_hardcalls_path)
                .annotate_variants_expr(['va.nonsplit_alleles = v.altAlleles.map(a => a.alt)',
                                         'va.hasStar = v.altAlleles.exists(a => a.isStar)'])
                .split_multi()
                .annotate_variants_expr([
                'va.wasMixed = va.variantType == "mixed"',
                'va.alleleType = if(v.altAllele.isSNP) "snv"'
                '   else if(v.altAllele.isInsertion) "ins"'
                '   else if(v.altAllele.isDeletion) "del"'
                '   else "complex"'])
                .annotate_variants_vds(hapmap, code='va.hapmap = isDefined(vds)')
                .annotate_variants_vds(omni, code='va.omni = isDefined(vds)')
                .annotate_variants_vds(mills, code='va.mills = isDefined(vds)')
                .tdt(fam=genomes_fam)
                .write(raw_hardcalls_split_path, overwrite=args.overwrite)
        )

    if args.compute_mendel:
        (
            hc.read(raw_hardcalls_split_path)
            .mendel_errors(args.output,fam=genomes_fam)
        )
        mendel_path = args.output

    #Random forests
    if args.annotate_for_rf:
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

        rf = annotate_for_random_forests(rf, hc.read(omni_path), hc.read(mills_path))
        print(rf.variant_schema)

        rf.write(rf_ann_path)

        rf_ann_path = args.output + ".annotated_for_rf.vds"

    if args.rf:

        features = vqsr_features if args.vqsr_features else rf_features

        print("Starting RF with features: \n")
        pprint(features)

        vds = sample_RF_training_examples(hc.read(rf_ann_path, sites_only=True),
                                          fp_to_tp = args.fp_to_tp)

        if args.train_on_vsqsr_sites:
            vds = vds.annotate_variants_expr(['va.train = isDefined(va.info.VQSR_NEGATIVE_TRAIN_SITE) || isDefined(va.info.VQSR_POSITIVE_TRAIN_SITE)',
                                              'va.label = ...'])
        else:
            vds = vds.annotate_variants_expr(['va.train = va.train && va.calldata.qc_samples_raw.AC[va.aIndex] > 0'])
        vds = rf.run_rf(hc, vds, rf_features = features, num_trees =args.num_trees, max_depth = args.max_depth)
        vds.write(rf_path)

        rf_path = args.output + '.rf.vds'

    if args.write_results:
        #Output metrics for RF evaluation
        out_num = re.search(iteration_re, rf_path).groups()[0]
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
                                  '(va.mendel>0 && va.calldata.all_samples_raw.AC[va.aIndex] == 1 )' % (nvariants[0],nvariants[1]))
            .export_variants(rf_path + ".va.txt.bgz", ",".join(out_metrics))
        )

#.filter_variants_expr('pcoin([1.0,1000000 / global.variantsByType[va.variantType]].min) || (va.mendel>0 && va.calldata.qc_samples_raw.AC[va.aIndex] == 1 )')
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write_hardcalls', help='Creates a hardcalls vds', action='store_true')
    parser.add_argument('--compute_mendel', help='Computes Mendel errors', action='store_true')
    parser.add_argument('--annotate_for_rf', help='Creates an annotated VDS with features for RF', action='store_true')
    parser.add_argument('--rf', help='Run RF model', action='store_true')
    parser.add_argument('--write_results', help='Skip writing results', action='store_true')
    parser.add_argument('--fp_to_tp', help='Sets to ratio of TPs to FPs for creating the RF model.', default=1.0)
    parser.add_argument('--num_trees', help='Number of tress in the RF model.', default=500)
    parser.add_argument('--max_depth', help='Maxmimum tree depth in the RF model.', default=5)
    parser.add_argument('--vqsr_features', help='Use VQSR features only (+ snv/indel variant type)', action='store_true')
    parser.add_argument('--train_on_vsqsr_sites', help='Use VQSR training sites',
                        action='store_true')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--input', '-i', help='Input prefix (if different from output)', required=False)
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    args = parser.parse_args()
    main(args)

from variantqc import annotate_with_additional_rf_files
from resources import *
import argparse
from utils import *


# Annotations to output
concordance_annotations = ['chrom = v.contig',
'pos = v.start',
'ref = v.ref',
'alt = v.alt',
'snp = v.altAllele.isSNP',
'ins = v.altAllele.isInsertion',
'del = v.altAllele.isDeletion',
'concordance = va.concordance'
                           ]

truth_concordance_annotations = concordance_annotations + ['type = va.variantType',
                                      'wassplit = va.left.wasSplit',
                                      'vqslod = va.left.info.VQSLOD',
                                      'truth.wassplit = va.right.wasSplit',
                                      'truth_gt = va.truth_gt',
                                      'called_gt = va.called_gt'
                                      ]

omes_concordance_annotations = concordance_annotations + [
    'type = va.variantType',
    'gg_nalleles = va.gg.altAlleles.length',
    'gg_ac_raw = va.gg.calldata.raw.AC[va.gg.aIndex]',
    'gg_ac_release_samples_raw = va.gg.calldata.release_samples_raw.AC[va.gg.aIndex]',
    'gg_is_minor_allele = va.gg.calldata.raw.AC.exists(x => x > va.gg.calldata.raw.AC[va.gg.aIndex])',
    'gg_vqslod = va.gg.info.VQSLOD',
    'gg_vqsr_label = if(va.gg.info.POSITIVE_TRAIN_SITE) "TP" else '
    '    if(va.gg.info.NEGATIVE_TRAIN_SITE) "FP"'
    '    else NA:String',
    'ge_nalleles = va.ge.altAlleles.length',
    'ge_ac_raw = va.ge.calldata.raw.AC[va.ge.aIndex]',
    'ge_is_minor_allele = va.ge.calldata.raw.AC.exists(x => x > va.ge.calldata.raw.AC[va.ge.aIndex])',
    'ge_vqslod = va.ge.info.VQSLOD',
    'ge_vqsr_label = if(va.ge.info.POSITIVE_TRAIN_SITE) "TP" else '
    '    if(va.ge.info.NEGATIVE_TRAIN_SITE) "FP"'
    '    else NA:String',
    'n_discordant = va.nDiscordant'
]


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

    global_concordance, s_concordance, v_concordance = left_vds.concordance(right=right_vds)

    s_concordance.write(out_prefix + ".s_concordance.kt", overwrite=overwrite)

    vds = VariantDataset.from_table(v_concordance)
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


def export_truth_concordance(vds, rf_ann_files, output):
    """

    Exports concordance between a VDS and a truth sample

    :param VariantDataset vds: Concordance VDS
    :param list of str rf_ann_files: RF annotation files, in format name|path|va.rf root
    :param str output: Output file path
    """
    out_annotations = truth_concordance_annotations
    if rf_ann_files:
        vds, additional_out_metrics = annotate_with_additional_rf_files(vds, rf_ann_files)
        out_annotations.extend(additional_out_metrics)

    vds.export_variants(output, ",".join(out_annotations))


def main(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hc = HailContext(log='/concordance.log')

    raw_hardcalls_split_path = full_genome_hardcalls_split_vds_path if (args.genomes) else full_exome_hardcalls_split_vds_path
    if args.compute_syndip_concordance:

        sample_map = {'CHMI_CHMI3_WGS1' : 'CHMI_CHMI3_Nex1'} if (args.genomes) else {'CHMI_CHMI3_WGS1' : 'CHMI_CHMI3_WGS1'}

        compute_concordance(left_vds=hc.read(raw_hardcalls_split_path),
                            right_vds= hc.read(syndip_vds_path).rename_samples(sample_map),
                            out_prefix = args.output + ".syndip",
                            high_conf_regions = syndip_high_conf_regions_bed_path,
                            samples = sample_map.values(),
                            overwrite = args.overwrite,
                            left_name='called',
                            right_name='truth')

    if args.compute_na12878_concordance:
        high_conf_regions = NA12878_high_conf_regions_bed_path if (args.genomes) else NA12878_high_conf_exome_regions_bed_path
        sample_map = {'INTEGRATION': 'G94982_NA12878'} if (args.genomes) else {'INTEGRATION': 'C1975::NA12878'}
        compute_concordance(left_vds= hc.read(raw_hardcalls_split_path),
                            right_vds= hc.read(NA12878_vds_path).rename_samples(sample_map),
                            samples = sample_map.values,
                            high_conf_regions = high_conf_regions,
                            out_prefix = args.output + ".NA12878",
                            overwrite = args.overwrite,
                            left_name='called',
                            right_name='truth'
                            )

    if args.write_syndip_concordance:
        export_truth_concordance(hc.read(args.output + ".syndip.v_concordance.vds"), args.rf_ann_files, args.output + ".syndip.stats.txt.bgz")

    if args.write_na12878_concordance:
        export_truth_concordance(hc.read(args.output + ".NA12878.v_concordance.vds"), args.rf_ann_files, args.output + ".NA12878.stats.txt.bgz")

    if args.compute_omes_concordance:
        exomes = hc.read(full_exome_hardcalls_split_vds_path).filter_samples_list(read_list_data(exomes_qc_pass_samples_list_path))
        exomes = rename_samples(exomes, exomes_to_combined_IDs_tsv_path)

        genomes = hc.read(full_genome_hardcalls_split_vds_path).filter_samples_list(read_list_data(genomes_qc_pass_samples_list_path))
        genomes = rename_samples(genomes, genomes_to_combined_IDs_tsv_path)

        compute_concordance(left_vds= genomes,
                            right_vds= exomes,
                            high_conf_regions = exomes_high_conf_regions_intervals_path,
                            out_prefix = args.output + ".omes",
                            overwrite = args.overwrite,
                            left_name="gg",
                            right_name="ge")

    if args.write_omes_concordance:
        vds = hc.read(args.output + ".omes.v_concordance.vds")
        exomes_vqsr = hc.read(vqsr_vds_path)
        vds = vds.annotate_variants_vds(exomes_vqsr,
                                        'va.ge.info.VQSLOD = vds.info.VQSLOD, '
                                        'va.ge.info.POSITIVE_TRAIN_SITE = vds.info.POSITIVE_TRAIN_SITE,'
                                        'va.ge.info.NEGATIVE_TRAIN_SITE = vds.info.NEGATIVE_TRAIN_SITE')

        out_metrics = omes_concordance_annotations
        if args.rf_ann_files:
            vds, additional_out_metrics = annotate_with_additional_rf_files(vds, args.rf_ann_files)
            if args.exomes:
                additional_out_metrics = ["ge_" + a for a in additional_out_metrics]
            else:
                additional_out_metrics = ["gg_" + a for a in additional_out_metrics]
            out_metrics.extend(additional_out_metrics)

        vds = vds.filter_variants_expr('va.concordance[3].exists(x => x>0) || va.concordance[4].exists(x => x>0) || va.concordance[0:2].map(x => x[3:4].sum).sum >0')
        vds.export_variants(args.output + ".omes.concordance.txt.bgz", ",".join(out_metrics))


    if args.write_omes_interval_concordance:
        vds = hc.read(args.output + ".omes.v_concordance.vds")
        vds = vds.filter_variants_expr(
            'va.concordance[3].exists(x => x>0) || va.concordance[4].exists(x => x>0) || va.concordance[0:2].map(x => x[3:4].sum).sum >0')
        intervals_kt = KeyTable.import_interval_list(exomes_high_conf_regions_intervals_path)
        intervals_kt = intervals_kt.annotate('interval_name = str(interval)')
        vds = vds.annotate_variants_table(intervals_kt, root = 'va.interval')
        concordance_kt = vds.variants_table()
        concordance_kt = concordance_kt.aggregate_by_key(['interval = va.interval'], ['n = va.count()',
                                                                           'n_missing_gg = va.map(x => x.concordance).filter(c => c[0].exists(x => x > 0)).count()',
                                                                           'n_missing_ge = va.map(x => x.concordance).filter(c => c.exists(x => x[0] > 0)).count()'])
        concordance_kt.export(args.output + '.omes_concordance_intervals.txt.bgz')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Computes concordance from exomes', action='store_true')
    parser.add_argument('--genomes', help='Computes concordance from genomes', action='store_true')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')

    compute = parser.add_argument_group('Computing concordance')
    compute.add_argument('--compute_syndip_concordance', help='Computes syndip concordance', action='store_true')
    compute.add_argument('--compute_na12878_concordance', help='Computes NA12878 concordance', action='store_true')
    compute.add_argument('--compute_omes_concordance', help='Computes exomes/genomes concordance', action='store_true')

    export = parser.add_argument_group('Exporting concordance stats')
    export.add_argument('--write_syndip_concordance', help='Write syndip concordance', action='store_true')
    export.add_argument('--write_na12878_concordance', help='Write NA12878 concordance', action='store_true')
    export.add_argument('--write_omes_concordance', help='Write exomes/genomes concordance', action='store_true')
    export.add_argument('--write_omes_interval_concordance', help='Write exomes/genomes concordance by interval',
                        action='store_true')
    export.add_argument('--rf_ann_files',
                        help='Annotation files to annotate RF results in pipe-delimited format: name|location|rf_root',
                        nargs='+')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

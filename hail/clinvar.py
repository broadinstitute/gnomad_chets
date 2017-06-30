from utils import *
import argparse
from variantqc import annotate_with_additional_rf_files


def get_pop_expr(strats, metrics=['AN', 'AF', 'AC'], index_into_array=False):
    res = {}
    postfix = '[vds.aIndex -1]' if index_into_array else ''
    for metric in metrics:
        res.update({'{}_{}'.format(metric, x): 'vds.info.{}_{}{}'.format(metric, x, postfix) for x in
                    strats})
    return res


def main(args):
    hc = HailContext(log='/clinvar.log')

    if args.write_clinvar:
        logger.info("Creating clinvar VDS")
        clinvar_kt = hc.import_table(clinvar_tsv_path, types='chrom: String', impute=True, missing='NA')
        clinvar_kt = clinvar_kt.annotate('v = Variant(str(chrom),pos,ref,alt)').key_by('v')
        vds = VariantDataset.from_keytable(clinvar_kt)
        vds = vds.annotate_variants_expr('va = drop(va, chrom, pos, ref, alt)')
        vds = vds.repartition(1000)
        vds.write(clinvar_vds_path, overwrite=True)

    if args.write:

        out_metrics = ['chrom = v.contig',
                       'pos = v.start',
                       'ref = v.ref',
                       'alt = v.alt',
                       'snv = v.altAllele.isSNP',
                       'segdup = va.segdup',
                       'lcr = va.lcr',
                       'clinical_significance = va.clinical_significance',
                       'molecular_consequence = va.molecular_consequence',
                       'pathogenic = va.pathogenic',
                       'benign = va.benign',
                       'conflicted = va.conflicted',
                       'review_status = va.review_status',
                       'symbol = va.symbol',
                       'hgvs_p = va.hgvs_p',
                       'hgvs_c = va.hgvs_c',
                       'gold_stars = va.gold_stars',
                       'all_submitters = va.all_submitters',
                       'all_traits = va.all_traits',
                       'inheritance_modes = va.inheritance_modes',
                       'age_of_onset = va.age_of_onset',
                       'prevalence = va.prevalence',
                       'disease_mechanism = va.disease_mechanism',
                       'origin = va.origin']

        exac_expr = {'AC': 'vds.info.AC_Adj[vds.aIndex -1]',
                     'filtered': '!(vds.filters.isEmpty || vds.filters.contains("PASS"))',
                     'wasSplit': 'vds.wasSplit'}

        gnomad_expr = {'AC': 'vds.info.AC',
                       'wasSplit': 'vds.wasSplit',
                       'filtered': 'vds.info.AS_FilterStatus.contains("RF")',
                       }

        if args.export_INFO:
            ge_expr =  get_pop_expr(EXOME_POPS + ['raw', 'Male', 'Female', 'POPMAX'])
            ge_expr.update(gnomad_expr)
            gg_expr = get_pop_expr(GENOME_POPS + ['raw', 'Male', 'Female', 'POPMAX'])
            gg_expr.update(gnomad_expr)
            exac_expr.update(get_pop_expr(EXAC_POPS + ['Adj', 'MALE', 'FEMALE', 'POPMAX'], metrics=['AC'], index_into_array=True))
            exac_expr.update(get_pop_expr(EXAC_POPS + ['Adj', 'MALE', 'FEMALE', 'POPMAX'], metrics=['AN']))

        else:
            exac_expr.update({'VQSLOD': 'vds.info.VQSLOD',
                              'vqsr_label': 'if(vds.info.POSITIVE_TRAIN_SITE) "TP" else '
                                            '    if(vds.info.NEGATIVE_TRAIN_SITE) "FP"'
                                            '    else NA:String',
                              'AC0': 'vds.info.AC_Adj[vds.aIndex -1] == 0'})

            gnomad_expr.update({
                x: 'vds.info.{}'.format(x) for x in ['VQSLOD',
                                                     'MQRankSum',
                                                     'SOR',
                                                     'InbreedingCoeff',
                                                     'ReadPosRankSum',
                                                     'GQ_MEDIAN',
                                                     'DP_MEDIAN',
                                                     'AB_MEDIAN',
                                                     'DREF_MEDIAN']
            })

            gnomad_expr.update({
                'vqsr_label': 'if(vds.info.VQSR_POSITIVE_TRAIN_SITE) "TP" else '
                              '    if(vds.info.VQSR_NEGATIVE_TRAIN_SITE) "FP"'
                              '    else NA:String',
                'rfprob': 'vds.info.AS_RF',
                'rf_label': 'if (vds.info.AS_RF_POSITIVE_TRAIN.toSet.contains(vds.aIndex -1)) "TP" '
                            '    else if (vds.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(vds.aIndex -1)) "FP"'
                            '    else NA:String',
                'AC0': 'vds.info.AC == 0',
            })
            ge_expr = gg_expr = gnomad_expr

        clinvar = hc.read(clinvar_vds_path)
        clinvar = clinvar.annotate_variants_table(KeyTable.import_interval_list(lcr_intervals_path), root='va.lcr')
        clinvar = clinvar.annotate_variants_table(KeyTable.import_interval_list(decoy_intervals_path), root='va.segdup')

        if args.exomes_new_rf_metrics_file:
            clinvar = clinvar.annotate_variants_vds(hc.read(args.exomes_new_rf_metrics_file),
                                                    expr = ",".join(['va.ge.max_pab = vds.stats.qc_samples_raw.pab.max',
                                                                     'va.ge.as_qd = vds.stats.qc_samples_raw.qd']))
            out_metrics.extend(['ge_max_pab = va.ge.max_pab','ge_as_qd = va.ge.as_qd'])

        if args.genomes_new_rf_metrics_file:
            clinvar = clinvar.annotate_variants_vds(hc.read(args.genomes_new_rf_metrics_file),
                                                    expr = ",".join(['va.gg.max_pab = vds.stats.qc_samples_raw.pab.max',
                                                                     'va.gg.as_qd = vds.stats.qc_samples_raw.qd']))
            out_metrics.extend(['gg_max_pab = va.gg.max_pab', 'gg_as_qd = va.gg.as_qd'])

        clinvar = clinvar.annotate_variants_vds(hc.read(final_exome_split_vds_path),
                                                expr=",".join(['va.ge.{} = {}'.format(k, v) for k, v in ge_expr.iteritems()]))

        clinvar = clinvar.annotate_variants_vds(hc.read(final_genome_split_vds_path),
                                                expr=",".join(['va.gg.{} = {}'.format(k, v) for k, v in gg_expr.iteritems()]))

        clinvar = clinvar.annotate_variants_vds(hc.read(final_exac_sites_vds_path).split_multi(),
                                                expr=",".join(['va.exac.{} = {}'.format(k, v) for k, v in exac_expr.iteritems()]))


        out_metrics.extend(['ge_{0} = va.ge.{0}'.format(k) for k in ge_expr.keys()])
        out_metrics.extend(['gg_{0} = va.gg.{0}'.format(k) for k in gg_expr.keys()])
        out_metrics.extend(['exac_{0} = va.exac.{0}'.format(k) for k in exac_expr.keys()])

        if args.rf_ann_files:
            clinvar, additional_out_metrics = annotate_with_additional_rf_files(clinvar, args.rf_ann_files)
            out_metrics.extend(additional_out_metrics)

        clinvar.export_variants(args.write, ",".join(out_metrics))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write_clinvar', help='Creates the clinvar vds', action='store_true')
    parser.add_argument('--export_INFO', help='Creates the clinvar vds', action='store_true')
    parser.add_argument('--exomes_new_rf_metrics_file', help='Optional file containing new exomes RF metrics: pab.max, qd')
    parser.add_argument('--genomes_new_rf_metrics_file',
                        help='Optional file containing new genomes RF metrics: va.stats.qc_samples_raw.pab.max, va.stats.qc_samples_raw.as_qd')
    parser.add_argument('--write', help='Writes the clinvar variants along with their annotations to the path provided')
    parser.add_argument('--rf_ann_files',
                        help='RF files to annotate results with in pipe-delimited format: name|location|rf_root',
                        nargs='+')
    args = parser.parse_args()
    main(args)

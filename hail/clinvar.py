from utils import *
import re
import argparse
from hail.representation import *
from variantqc import *

def main(args):

    hc = hail.HailContext(log='/clinvar.log')

    if args.write_clinvar:
        logger.info("Creating clinvar VDS")
        clinvar_kt = hc.import_keytable(clinvar_variants, config=hail.TextTableConfig(types='chrom: String', impute=True, missing='NA'))
        clinvar_kt = clinvar_kt.annotate('v = Variant(str(chrom),pos,ref,alt)').key_by('v')
        vds = hail.VariantDataset.from_keytable(clinvar_kt)
        vds = vds.annotate_variants_expr('va = drop(va, chrom, pos, ref, alt)')
        vds = vds.repartition(1000)
        vds.write(clinvar_vds, overwrite = True)


    if args.write:
        clinvar = hc.read(clinvar_vds)
        clinvar = clinvar.annotate_variants_table(hail.KeyTable.import_interval_list(lcr_path), root = 'va.lcr')
        clinvar = clinvar.annotate_variants_table(hail.KeyTable.import_interval_list(decoy_path), root= 'va.segdup')

        clinvar = clinvar.annotate_variants_vds(hc.read("gs://gnomad-genomes/variantqc/rf_new_stats.annotated_for_rf.vds"),
                                                expr = ",".join(['va.gg.max_pab = vds.stats.qc_samples_raw.pab.max',
                                                                 'va.gg.as_qd = vds.stats.qc_samples_raw.qd']))

        clinvar = clinvar.annotate_variants_vds(hc.read(final_exome_vds).split_multi(),
                                                expr=",".join(['va.ge.ac = vds.info.AC[vds.aIndex -1]',
                                                               'va.ge.vqslod = vds.info.VQSLOD',
                                                               'va.ge.vqsr_label = if(vds.info.VQSR_POSITIVE_TRAIN_SITE) "TP" else '
                                                               '    if(vds.info.VQSR_NEGATIVE_TRAIN_SITE) "FP"'
                                                               '    else NA:String',
                                                               'va.ge.rfprob = vds.info.AS_RF[vds.aIndex - 1]',
                                                               'va.ge.filtered = vds.info.AS_FilterStatus[vds.aIndex - 1].contains("RF")',
                                                               'va.ge.rf_label = if (vds.info.AS_RF_POSITIVE_TRAIN.toSet.contains(vds.aIndex -1)) "TP" '
                                                               '    else if (vds.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(vds.aIndex -1)) "FP"'
                                                               '    else NA:String',
                                                               'va.ge.AC0 = vds.info.AC[vds.aIndex -1] == 0',
                                                               'va.ge.wasSplit = vds.wasSplit',
                                                               'va.ge.MQRankSum = vds.info.MQRankSum',
                                                               'va.ge.SOR = vds.info.SOR',
                                                               'va.ge.InbreedingCoeff = vds.info.InbreedingCoeff',
                                                               'va.ge.ReadPosRankSum = vds.info.ReadPosRankSum',
                                                               'va.ge.GQ_MEDIAN = vds.info.GQ_MEDIAN[vds.aIndex -1]',
                                                               'va.ge.DP_MEDIAN = vds.info.DP_MEDIAN[vds.aIndex -1]',
                                                               'va.ge.AB_MEDIAN = vds.info.AB_MEDIAN[vds.aIndex -1]',
                                                               'va.ge.DREF_MEDIAN = vds.info.DREF_MEDIAN[vds.aIndex -1]'
                                                               ])
                                                )
        clinvar = clinvar.annotate_variants_vds(hc.read(final_genome_vds).split_multi(),
                                                expr=",".join(['va.gg.ac = vds.info.AC[vds.aIndex -1]',
                                                               'va.gg.vqslod = vds.info.VQSLOD',
                                                               'va.gg.vqsr_label = if(vds.info.VQSR_POSITIVE_TRAIN_SITE) "TP" else '
                                                               '    if(vds.info.VQSR_NEGATIVE_TRAIN_SITE) "FP"'
                                                               '    else NA:String',
                                                               'va.gg.rfprob = vds.info.AS_RF[vds.aIndex - 1]',
                                                               'va.gg.filtered = vds.info.AS_FilterStatus[vds.aIndex - 1].contains("RF")',
                                                               'va.gg.rf_label = if (vds.info.AS_RF_POSITIVE_TRAIN.toSet.contains(vds.aIndex -1)) "TP" '
                                                               '    else if (vds.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(vds.aIndex -1)) "FP"'
                                                               '    else NA:String',
                                                               'va.gg.AC0 = vds.info.AC[vds.aIndex -1] == 0',
                                                               'va.gg.wasSplit = vds.wasSplit',
                                                               'va.gg.MQRankSum = vds.info.MQRankSum',
                                                               'va.gg.SOR = vds.info.SOR',
                                                               'va.gg.InbreedingCoeff = vds.info.InbreedingCoeff',
                                                               'va.gg.ReadPosRankSum = vds.info.ReadPosRankSum',
                                                               'va.gg.GQ_MEDIAN = vds.info.GQ_MEDIAN[vds.aIndex -1]',
                                                               'va.gg.DP_MEDIAN = vds.info.DP_MEDIAN[vds.aIndex -1]',
                                                               'va.gg.AB_MEDIAN = vds.info.AB_MEDIAN[vds.aIndex -1]',
                                                               'va.gg.DREF_MEDIAN = vds.info.DREF_MEDIAN[vds.aIndex -1]'])
                                                )
        clinvar = clinvar.annotate_variants_vds(hc.read(final_exac_sites_vds).split_multi(),
                                                expr=",".join(['va.exac.ac = vds.info.AC_Adj[vds.aIndex -1]',
                                                               'va.exac.vqslod = vds.info.VQSLOD',
                                                               'va.exac.vqsr_label = if(vds.info.POSITIVE_TRAIN_SITE) "TP" else '
                                                               '    if(vds.info.NEGATIVE_TRAIN_SITE) "FP"'
                                                               '    else NA:String',
                                                               'va.exac.filtered = !(vds.filters.isEmpty || vds.filters.contains("PASS"))',
                                                               'va.exac.AC0 = vds.info.AC_Adj[vds.aIndex -1] == 0',
                                                               'va.exac.wasSplit = vds.wasSplit'])
                                                )

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
                       'origin = va.origin',
                       'ac_ge = va.ge.ac',
                       'rfprob_ge = va.ge.rfprob',
                       'rf_label_ge = va.ge.rf_label',
                       'vqslod_ge = va.ge.vqslod',
                       'vqsr_label_ge = va.ge.vqsr_label',
                       'filtered_ge = va.ge.filtered',
                       'ac0_ge = va.ge.AC0',
                       'ac_gg = va.gg.ac',
                       'rfprob_gg = va.gg.rfprob',
                       'rf_label_gg = va.gg.rf_label',
                       'vqslod_gg = va.gg.vqslod',
                       'vqsr_label_gg = va.gg.vqsr_label',
                       'filtered_gg = va.gg.filtered',
                       'ac0_gg = va.gg.AC0',
                       'ac_exac = va.exac.ac',
                       'vqslod_exac = va.exac.vqslod',
                       'vqsr_label_exac = va.exac.vqsr_label',
                       'filtered_exac = va.exac.filtered',
                       'ac0_exac = va.gg.AC0',
                       'was_split_gg = va.gg.wasSplit',
                       'was_split_ge = va.ge.wasSplit',
                       'was_split_exac = va.exac.wasSplit',
                       'ge_MQRankSum = va.ge.MQRankSum',
                       'ge_SOR = va.ge.SOR',
                       'ge_InbreedingCoeff = va.ge.InbreedingCoeff',
                       'ge_ReadPosRankSum = va.ge.ReadPosRankSum',
                       'ge_GQ_MEDIAN = va.ge.GQ_MEDIAN',
                       'ge_DP_MEDIAN = va.ge.DP_MEDIAN',
                       'ge_AB_MEDIAN = va.ge.AB_MEDIAN',
                       'ge_DREF_MEDIAN = va.ge.DREF_MEDIAN',
                       'gg_MQRankSum = va.gg.MQRankSum',
                       'gg_SOR = va.gg.SOR',
                       'gg_InbreedingCoeff = va.gg.InbreedingCoeff',
                       'gg_ReadPosRankSum = va.gg.ReadPosRankSum',
                       'gg_GQ_MEDIAN = va.gg.GQ_MEDIAN',
                       'gg_DP_MEDIAN = va.gg.DP_MEDIAN',
                       'gg_AB_MEDIAN = va.gg.AB_MEDIAN',
                       'gg_DREF_MEDIAN = va.gg.DREF_MEDIAN',
                       'gg_max_pab = va.gg.max_pab',
                       'gg_as_qd = va.gg.max_pab'
                       ]

        if args.rf_ann_files:
            clinvar, additional_out_metrics = annotate_with_additional_rf_files(clinvar, args.rf_ann_files)
            out_metrics.extend(additional_out_metrics)

        clinvar.export_variants(args.write, ",".join(out_metrics))




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write_clinvar', help='Creates the clinvar vds', action='store_true')
    parser.add_argument('--write', help='Writes the clinvar variants along with their annotations to the path provided')
    parser.add_argument('--rf_ann_files', help='RF files to annotate results with in pipe-delimited format: name|location|rf_root', nargs='+')
    args = parser.parse_args()
    main(args)
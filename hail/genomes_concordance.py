from pyhail import *

try:
    hc
except NameError:
    hc = HailContext(log='/variantqc.log')

#magic(hc)

#Inputs
raw_hardcalls_split_path = "gs://gnomad/gnomad.raw_hardcalls.split.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats7.vds"
syndip_path = "gs://gnomad/truth-sets/hybrid.m37m.vds"
NA12878_path = " gs://gnomad/truth-sets/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vds"
exomes_hardcalls_path = "gs://exac2/exacv2.raw.hardcalls.splitmulti.qc.concordance_samples.vds"
exomes_to_combined_IDs = "gs://gnomad/exac_to_combined.IDs.txt"
exomes_concordance_samples = "gs://exac2/exac_samples_in_gnomad.txt"
genomes_to_combined_IDs = "gs://gnomad/gnomad_to_combined.IDs.txt"
genomes_concordance_samples = "gs://gnomad/gnomad_samples_in_exac.txt"

#Resources
lcr_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
decoy_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
syndip_high_conf_regions_path = "gs://gnomad/truth-sets/hybrid.m37m.bed"
NA12878_high_conf_regions_path = "gs://gnomad/truth-sets/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed"
exomes_high_conf_regions_path = "gs://exac2/high_coverage.auto.interval_list"

#Outputs
syndip_concordance_prefix = "gs://gnomad/truth-sets/gnomad_hybrid"
NA12878_concordance_prefix = "gs://gnomad/truth-sets/gnomad_NA12878"
exomes_concordance_prefix = "gs://gnomad/gnomad_exac_raw"

#Annotations to output
concordance_annotations = ['chrom = v.contig',
'pos = v.start',
'ref = v.ref',
'alt = v.alt',
'snp = v.altAllele.isSNP',
'ins = v.altAllele.isInsertion',
'del = v.altAllele.isDeletion',
'concordance = va.concordance'
                           ]

truth_concordance_annotations = list(concordance_annotations)
truth_concordance_annotations.extend(['type = va.rf.variantType',
                                      'wassplit = va.left.wasSplit',
                                      'vqslod = va.rf.info.VQSLOD',
                                      'qd = va.rf.info.QD',
                                      'truth.wassplit = va.right.wasSplit',
                                      'truth_gt = va.truth_gt',
                                      'called_gt = va.called_gt',
                                      'training = va.rf.train',
                                      'label = va.rf.label',
                                      'rfpred1 = va.rf.RF1.prediction',
                                      'rfprob1 = va.rf.RF1.probability["TP"]',
                                      'rfpred2 = va.rf.RF2.prediction',
                                      'rfprob2 = va.rf.RF2.probability["TP"]'
                                      ])

exomes_concordance_annotations = list(concordance_annotations)
exomes_concordance_annotations.extend(['gnomad.multi = va.left.wasSplit',
                                       'gnomad.vqslod = va.left.info.VQSLOD',
                                       'gnomad.qd = va.left.info.QD',
                                       'gnomad.ac = va.left.info.AC[va.right.aIndex-1]',
                                       'gnomad.an = va.left.info.AN',
                                       'exac.multi = va.right.wasSplit',
                                       'exac.vqslod = va.right.info.VQSLOD',
                                       'exac.ac = va.right.info.AC[va.left.aIndex-1]',
                                       'exac.an = va.right.info.AN',
                                       'exac.qd = va.right.info.QD',
                                       'gnomad.training = va.rf.train',
                                       'gnomad.label = va.rf.label',
                                       'gnomad.rfpred1 = va.rf.RF1.prediction',
                                       'gnomad.rfprob1 = va.rf.RF1.probability["TP"]',
                                       'gnomad.rfpred2 = va.rf.RF2.prediction',
                                       'gnomad.rfprob2 = va.rf.RF2.probability["TP"]'
                                       ])

#Actions
compute_syndip_concordance=False
export_syndip_concordance=True
compute_NA12878_concordance=False
export_NA12878_concordance=True
compute_Exomes_concordance=True
export_Exomes_concordance=True

def filter_for_concordance(vds,high_conf_regions):
    return(
        vds.filter_variants_intervals(lcr_path, keep=False)
            .filter_variants_intervals(decoy_path, keep=False)
            .filter_variants_intervals(high_conf_regions, keep=True)
    )

def compute_concordance(vds, rf_vds, sample, truth_path, high_conf_regions, out_prefix, out_annotations, recompute=True):
    if(recompute):
        truth = filter_for_concordance( hc.read(truth_path), high_conf_regions=high_conf_regions)

        (s_concordance, v_concordance) = (filter_for_concordance(vds, high_conf_regions=high_conf_regions)
                                          .filter_variants_intervals(lcr_path, keep=False)
                                          .filter_variants_intervals(decoy_path, keep=False)
                                          .filter_samples_expr('s.id == "%s"' % sample, keep=True)
                                          .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count()>0', keep=True)
                                          .concordance(right=truth)
                                          )
        s_concordance.write(out_prefix + ".s_concordance.vds")
        v_concordance.write(out_prefix + ".v_concordance.vds")

    else:
        v_concordance = hc.read(out_prefix + ".v_concordance.vds")

    (
        v_concordance.annotate_variants_vds(rf_vds, root='va.rf')
        .annotate_global_expr_by_variant('global.gt_mappings = ["missing","no_call","homref","het","homvar"]')
        .annotate_variants_expr('va.gt_arr = range(5).find(i => va.concordance[i].exists(x => x > 0))')
        .annotate_variants_expr('va.truth_gt =  global.gt_mappings[va.gt_arr],'
                                'va.called_gt = global.gt_mappings[range(5).find(i => va.concordance[va.gt_arr][i] >0)]')
        .export_variants(out_prefix + ".stats.txt.bgz", ",".join(out_annotations))
     )

if(compute_syndip_concordance or export_syndip_concordance):
    compute_concordance(hc.read(raw_hardcalls_split_path),
                        hc.read(rf_path),
                        'CHMI_CHMI3_WGS1',
                        syndip_path,
                        syndip_high_conf_regions_path,
                        syndip_concordance_prefix,
                        truth_concordance_annotations,
                        compute_syndip_concordance)

if(compute_NA12878_concordance or export_NA12878_concordance):
    compute_concordance(hc.read(raw_hardcalls_split_path),
                        hc.read(rf_path),
                        'G94982_NA12878',
                        NA12878_path,
                        NA12878_high_conf_regions_path,
                        NA12878_concordance_prefix,
                        truth_concordance_annotations,
                        compute_NA12878_concordance)

if (compute_Exomes_concordance or export_Exomes_concordance):
    if(compute_Exomes_concordance):
        exomes = (
            filter_for_concordance(hc.read(exomes_hardcalls_path), high_conf_regions=exomes_high_conf_regions_path)
            .filter_samples_list(exomes_concordance_samples)
            .rename_samples(exomes_to_combined_IDs)
            .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count()>0', keep=True)
        )

        (s_concordance, v_concordance) = (
            filter_for_concordance(hc.read(raw_hardcalls_split_path), high_conf_regions=exomes_high_conf_regions_path)
            .filter_samples_list(genomes_concordance_samples)
            .rename_samples(genomes_to_combined_IDs)
            .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count()>0', keep=True)
            .concordance(right=exomes)
                                          )
        s_concordance.write(exomes_concordance_prefix + ".s_concordance.vds")
        v_concordance.write(exomes_concordance_prefix + ".v_concordance.vds")

    else:
        v_concordance = hc.read(exomes_concordance_prefix + ".v_concordance.vds")

    if(export_Exomes_concordance):
        (
            v_concordance.annotate_variants_vds(hc.read(rf_path), root='va.rf')
                .filter_variants_expr('va.concordance[3].exists(x => x>0) || va.concordance[4].exists(x => x>0) || va.concordance[0:2].map(x => x[3:4].sum).sum >0')
            .export_variants(exomes_concordance_prefix + ".stats.txt.bgz", ",".join(exomes_concordance_annotations))
        )
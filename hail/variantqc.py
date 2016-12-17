__author__ = 'konrad'

import pyhail

from utils import *

adj_criteria = 'g.gq >= 20 && g.dp >= 10 & g.ad[1]/g.dp >= 0.25'


def write_hardcalls(vds, output_path, meta_path, adj=True, metrics=True, partitions=10000, shuffle=True):

    out = vds.annotate_samples_table(meta_path, 'sample', impute=True, root='sa.meta')

    if(metrics):
        pre_adj_expression = (
            get_variant_type_expr("va.variantType")
                .extend(get_stats_expr("va.stats.raw", medians=True))
                .extend(get_hists_expr("va.hists.raw"))
                .append('va.calldata.allsamples_raw = gs.callStats(v)')
        )
        out = out.annotate_alleles_expr(pre_adj_expression)

    if(adj):
        out = (
            out.filter_genotypes(adj_criteria)
                .filter_alleles('gs.filter(g => g.gtj + g.gtk).sum() > 0', annotation='')  # TODO: Add annotation
               )

        if(metrics):
            post_adj_expression = (
                get_stats_expr("va.stats.Adj", medians=True)
                    .extend(get_hists_expr("va.hists.Adj"))
                    .append('va.calldata.allsamples_Adj = gs.callStats(v)')
            )
            out = out.annotate_alleles_expr(post_adj_expression)

    return(out.hardcalls()
           .repartition(numPartitions=partitions,shuffle=shuffle)
           .write(output_path))


def get_transmitted_singletons(vds, output_vds_path, fam_path, autosomes_intervals):
    return (vds
            .filter_variants_intervals(autosomes_intervals)
            .tdt(fam_path)
            .filter_variants_expr('va.tdt.nTransmitted == 1 && va.info.AC[va.aIndex - 1]==2')
            .filter_samples_all()
            .write(output_vds_path))


# Use raw VDS for determing true positives
# get_transmitted_singletons(raw_split_hardcallvds_path, tdt_vds_path)


def annotate_for_random_forests(vds, transmission_vds = None , omni_vds = None, mills_vds = None ):

    new_vds = vds

    if(transmission_vds is not None):
        new_vds = new_vds.annotate_variants_vds(transmission_vds, code='va.transmitted_singleton = isDefined(vds)')

    if(omni_vds is not None):
        new_vds = new_vds.annotate_variants_vds(omni_vds, code='va.omni = isDefined(vds)')

    if (mills_vds is not None):
        new_vds = new_vds.annotate_variants_vds(mills_vds, code='va.mills = isDefined(vds)')

    return(new_vds
               .annotate_variants_expr('va.TP = va.omni || va.mills || va.transmitted_singleton, '
                                          'va.FP = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30')
               .annotate_variants_expr('va.label = if(!isMissing(va.FP) && va.FP) "FP" else if(va.TP) "TP" else NA: String, '
                                       'va.train = v.contig != "20" && (va.TP || va.FP)')
               .annotate_global_expr_by_variant('global.ac_hist_indels_mills = variants.filter(v => va.mills).map(v => log10(va.calldata.raw.AF[va.aIndex])).hist(-6, 0, 20),'
                                     'global.ac_hist_indels_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isIndel).map(v => log10(va.calldata.raw.AF[va.aIndex])).hist(-6, 0, 20),'
                                     'global.ac_hist_snps_omni = variants.filter(v => va.omni).map(v => log10(va.calldata.raw.AF[va.aIndex])).hist(-6, 0, 20),'
                                     'global.ac_hist_snps_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isSNP).map(v => log10(va.calldata.raw.AF[va.aIndex])).hist(-6, 0, 20),'
                                     'global.indels_mills = variants.filter(v => va.mills).count(),'
                                     'global.indels_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isIndel).count(),'
                                     'global.snps_omni = variants.filter(v => va.omni).count(),'
                                     'global.snps_tx_singleton = variants.filter(v => va.transmitted_singleton && v.altAllele.isSNP).count(),'
                                     'global.nTP = variants.filter(x => va.label == "TP").count(), '
                                     'global.nFP = variants.filter(x => va.label == "FP").count()')
               .show_globals()
    )

    #new_vds.export_variants('exac2_rf.txt.bgz', 'chrom = v.contig, pos = v.start, ref = v.ref, alt = v.altAllele, type = va.variantType, label = va.label, rfprob = va.RF.probability["TP"]')



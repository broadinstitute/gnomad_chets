from utils import *
from hail import *
from pprint import pprint

hc = HailContext(log='hail.log')

vds = hc.read(final_genome_autosomes)
vds = vds.filter_variants_intervals(Interval.parse("22"))

a_anns = [u'va.info.AC',
 u'va.info.AF',
 u'va.info.GQ_HIST_ALT',
 u'va.info.DP_HIST_ALT',
 u'va.info.AB_HIST_ALT',
 u'va.info.AC_AFR',
 u'va.info.AC_AMR',
 u'va.info.AC_ASJ',
 u'va.info.AC_EAS',
 u'va.info.AC_FIN',
 u'va.info.AC_NFE',
 u'va.info.AC_OTH',
# u'va.info.AC_SAS',
 u'va.info.AC_Male',
 u'va.info.AC_Female',
 u'va.info.AF_AFR',
 u'va.info.AF_AMR',
 u'va.info.AF_ASJ',
 u'va.info.AF_EAS',
 u'va.info.AF_FIN',
 u'va.info.AF_NFE',
 u'va.info.AF_OTH',
# u'va.info.AF_SAS',
 u'va.info.AF_Male',
 u'va.info.AF_Female',
 u'va.info.AC_raw',
 u'va.info.AF_raw',
 u'va.info.Hom_AFR',
 u'va.info.Hom_AMR',
 u'va.info.Hom_ASJ',
 u'va.info.Hom_EAS',
 u'va.info.Hom_FIN',
 u'va.info.Hom_NFE',
 u'va.info.Hom_OTH',
# u'va.info.Hom_SAS',
 u'va.info.Hom_Male',
 u'va.info.Hom_Female',
 u'va.info.Hom_raw',
 u'va.info.Hom',
 u'va.info.POPMAX',
 u'va.info.AC_POPMAX',
 u'va.info.AN_POPMAX',
 u'va.info.AF_POPMAX',
 u'va.info.DP_MEDIAN',
 u'va.info.DREF_MEDIAN',
 u'va.info.GQ_MEDIAN',
 u'va.info.AB_MEDIAN',
 u'va.info.AS_RF',
 u'va.info.AS_FilterStatus']

vds = (vds.split_multi()
       .annotate_global_py('global.csqs', CSQ_ORDER, TArray(TString()))
       .annotate_variants_expr(index_into_arrays(a_anns, vep_root='va.vep'))
       .annotate_variants_expr([
    'va.alt = v.altAlleles[0]',
    'va.filters = va.filters.mkString("|")',
    'va.info.AS_FilterStatus = va.info.AS_FilterStatus.mkString("|")',
    'va.vep = va.vep.transcript_consequences.map(x => drop(x, domains))'])
       .annotate_variants_bed('gs://exac2/gm12878_hmm.bed', 'va.chrom_state')
       .annotate_variants_bed('gs://gnomad-resources/annotations/finucane_et_al/TSS_Hoffman.bed', 'va.tss')
       .annotate_variants_bed('gs://gnomad-resources/annotations/finucane_et_al/PromoterFlanking_Hoffman.bed',
                              'va.promoter')
       .annotate_variants_bed('gs://gnomad-resources/annotations/finucane_et_al/Enhancer_Hoffman.bed', 'va.enhancer')
       .annotate_variants_bed('gs://gnomad-resources/annotations/finucane_et_al/Coding_UCSC.bed', 'va.coding')
       )

vds = (vds.annotate_variants_expr(
    'va.worst_csq = '
    'let csq = global.csqs.find(c => va.vep.flatMap(x => x.consequence_terms).toSet().contains(c)) in '
    'if (va.vep.filter(x => x.lof == "HC").length > 0)'
    '   csq + "-HC" '
    'else '
    '   if (va.vep.filter(x => x.lof == "LC").length > 0)'
    '       csq + "-LC" '
    '   else '
    '       if (va.vep.filter(x => x.polyphen_prediction == "probably_damaging").length > 0)'
    '           csq + "-probably_damaging"'
    '       else'
    '           if (va.vep.filter(x => x.polyphen_prediction == "possibly_damaging").length > 0)'
    '               csq + "-possibly_damaging"'
    '           else'
    '               csq'
).annotate_variants_expr('va.info = drop(va.info, AS_RF_POSITIVE_TRAIN, AS_RF_NEGATIVE_TRAIN, CSQ)'))

kt = vds.variants_keytable().flatten()

non_coding =  kt.query(['tss = va.tss.filter(x => x).map(x => va.info.AF).stats()'])

print(non_coding)
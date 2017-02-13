
from hail import *

hc = HailContext()

root = 'gs://gnomad-public/truth-sets'

truth_sets = [
    '1000G_omni2.5.b37.vcf.bgz',
    'hapmap_3.3.b37.vcf.bgz',
    'Mills_and_1000G_gold_standard.indels.b37.vcf.bgz',
    'NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.bgz'
]

for truth_vcf in truth_sets:
    vds = truth_vcf.replace('.vcf.bgz', '.vds')
    store_gq = truth_vcf.startswith('NA12878')  # 12878 VCF is missing PLs
    hc.import_vcf('%s/vcf/%s' % (root, truth_vcf), store_gq=store_gq).split_multi().write('%s/%s' % (root, vds))

hc.import_vcf('gs://gnomad/truth-sets/hybrid.m37m.vcf.gz', force_bgz=True).rename_samples('gs://gnomad/truth-sets/hybrid.m37m.gnomad_rename').split_multi().write('gs://gnomad/truth-sets/hybrid.m37m.vds')

hc.import_vcf('gs://gnomad-exomes/variantqc/ExAC.merged.sites_only.vcf.ICfiltered.recalibrated.vcf.bgz').split_multi().write('gs://gnomad-exomes/variantqc/gnomad.exomes.vqsr.vds')

hc.import_vcf('gs://gnomad-exomes/variantqc/ExAC.merged.sites_only.vcf.ICfiltered.recalibrated.vcf.bgz').write('gs://gnomad-exomes/variantqc/gnomad.exomes.vqsr.unsplit.vds')

hc.import_vcf('gs://gnomad-public/cpg.vcf.bgz').write('gs://gnomad-public/cpg.vds')

# Raw counts:
# print hc.import_vcf('gs://exac2/variantqc/ExAC.merged.sites_only.vcf.ICfiltered.recalibrated.vcf.bgz').query_variants('variants.count()')[0]  # 18444471 sites
# print hc.read('gs://exac2/variantqc/exac2_vqsr.vds').query_variants('variants.count()')[0]  # 21352671 variants
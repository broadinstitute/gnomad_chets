
from utils import *
from pyhail import *

#Inputs

#vds_path = "gs://gnomad/gnom.ad.vds"
vds_path =  "gs://gnomad/gnomad.10ksites.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats7.vds"
vep_config = "/vep/vep-gcloud.properties"

#Resources
lcr_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
decoy_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
autosomes_intervals = "gs://gnomad/autosomes.txt"
dbsnp = "gs://gnomad-lfran/All_20160601.vcf.bgz"

#Outputs
out_root = "gs://gnomad-lfran/tmp"
out_vds = "%s/gnomad.sites.annotations.vds" % out_root
out_internal_vcf_prefix = "%s/gnomad.sites.internal" % out_root
out_external_vcf_prefix = "%s/gnomad.sites" % out_root
tmp_vds = "gs://gnomad-lfran/tmp/gnomad.sites.tmp.vds"
tmp_RF_ann_out = 'gs://gnomad-lfran/tmp/gnomad.rf.ann.txt.bgz'


hc = HailContext(log='/site_auto.log')

pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']

vep = hc.read(vep_path)

vds = ( hc.read(vds_path)
        .annotate_global_expr_by_sample('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
        .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=pyhail.TextTableConfig(impute=True))
        .filter_variants_intervals(autosomes_intervals)
        .filter_samples_expr('!isMissing(sa.meta.predicted_pop)') #Could be cleaner
        .annotate_variants_intervals(decoy_path, 'va.decoy')
        .annotate_variants_expr(get_add_filter_annotation('segdup','va.decoy'))
        .annotate_variants_intervals(lcr_path, 'va.lcr')
        .annotate_variants_expr(get_add_filter_annotation('LCR','va.lcr'))
        )

sites_vds = create_sites_vds_annotations(vds, pops, 1000, shuffle=False)
sites_vds.write(tmp_vds)
sites_vds = (
    # annotate_non_split_from_split(hc,non_split_vds_path=tmp_vds,
    #                               split_vds=hc.read(rf_path),
    #                               annotations=['va.RF'],
    #                               annotation_exp_out_path=tmp_RF_ann_out)
    hc.read(tmp_vds)
    .annotate_variants_vds(vep, code='va.info.CSQ = vds.csq')
    .vep(config=vep_config,csq=True,root='va.info.CSQ')
)

print 'Writing VDS...'
sites_vds.write(out_vds)
print 'Done! Writing VCF...'
for i in range(1,22):
    (hc.read(out_vds)
        .filter_variants_expr('v.contig == %d' % i)
        .export_vcf(out_internal_vcf_prefix + "%d.vcf.bgz" % i)
        .annotate_variants_expr('va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
        .export_vcf(out_external_vcf_prefix + "%d.vcf.bgz" % i)
     )
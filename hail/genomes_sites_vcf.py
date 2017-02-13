from utils import *
import time
from resources import *

# Inputs

# vds_path = "gs://gnomad/gnom.ad.vds"
vds_path = "gs://gnomad/gnomad.10ksites.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats12.vds"
vep_config = "/vep/vep-gcloud.properties"

# Resources
autosomes_intervals = "gs://gnomad/autosomes.txt"

# Outputs
date_time = time.strftime("%Y-%m-%d_%H-%M")
#date_time = '2017-01-30_22-18'
#date_time = '2016-12-20_21-20'
out_root = "gs://gnomad-lfran/tmp"
out_vds_prefix = "%s/gnomad.sites.annotations.%s" % (out_root, date_time)
out_internal_vcf_prefix = "%s/gnomad.sites.internal" % out_root
out_external_vcf_prefix = "%s/gnomad.sites" % out_root
tmp_vds_prefix = "gs://gnomad-lfran/tmp/gnomad.sites.tmp." + date_time
tmp_RF_ann_out = 'gs://gnomad-lfran/tmp/gnomad.rf.ann.tmp.' + date_time + 'txt.bgz'

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
rf_snv_cutoff = 0.5
rf_indel_cutoff = 0.75

#Actions
preprocess_autosomes = True
postprocess_autosomes = True
write_autosomes = True
preprocess_X = True
postprocess_X = True
write_X = True
preprocess_Y = True
postprocess_Y = True
write_Y = True


hc = HailContext(log='/site_auto.log')


def preprocess_vds(vds_path):
    print("Preprocessing %s\n" % vds_path)
    return (
        hc.read(vds_path)
            .annotate_global_py('global.pops',map(lambda x: x.lower(), pops), TArray(TString()))
            .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
            .annotate_samples_expr(['sa.meta.population = sa.meta.final_pop',
                                    'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
            .filter_samples_expr('!isMissing(sa.meta.predicted_pop)')  # Could be cleaner
            .filter_variants_intervals('gs://gnomad-lfran/tmp/test.interval')
            .annotate_variants_intervals(decoy_path, 'va.decoy')
            .annotate_variants_intervals(lcr_path, 'va.lcr')
    )

if preprocess_autosomes:
    (
        create_sites_vds_annotations(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf,
            npartitions=1000,
            shuffle=False)
        .write(tmp_vds_prefix + ".vds")
    )

if postprocess_autosomes:
    post_process_vds(hc, out_vds_prefix + ".vds", rf_path, rf_snv_cutoff, rf_indel_cutoff, vep_config).write(out_vds_prefix + ".vds")

if write_autosomes:
    for i in range(1, 23):
        write_vcfs(out_vds_prefix + ".vds", i, out_internal_vcf_prefix, out_external_vcf_prefix)

if preprocess_X:
    (
        create_sites_vds_annotations_X(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf,
            npartitions=100,
            shuffle=False)
        .write(tmp_vds_prefix + ".X.vds")
    )

if postprocess_X:
    post_process_vds(hc, out_vds_prefix + ".X.vds", rf_path, rf_snv_cutoff, rf_indel_cutoff, vep_config).write(out_vds_prefix + ".X.vds")

if write_X:
    write_vcfs(out_vds_prefix + ".X.vds", "X", out_internal_vcf_prefix, out_external_vcf_prefix)

if preprocess_Y:
    (
        create_sites_vds_annotations_Y(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf,
            npartitions=10,
            shuffle=False)
        .write(tmp_vds_prefix + ".Y.vds")
    )

if postprocess_Y:
    post_process_vds(hc, out_vds_prefix + ".Y.vds", rf_path, rf_snv_cutoff, rf_indel_cutoff, vep_config).write(out_vds_prefix + ".Y.vds")

if write_Y:
    write_vcfs(out_vds_prefix + ".Y.vds", "Y", out_internal_vcf_prefix, out_external_vcf_prefix)
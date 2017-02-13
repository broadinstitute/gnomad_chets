from variantqc import *
import time

# Inputs

bucket = 'gs://gnomad-exomes'
autosomes_intervals = '%s/intervals/autosomes.txt' % bucket
# evaluation_intervals = '%s/intervals/exome_evaluation_regions.v1.intervals' % bucket
# high_coverage_intervals = '%s/intervals/high_coverage.auto.interval_list' % bucket
meta_path = 'gs://gnomad-exomes-raw/super_meta.txt.bgz'

root = '%s/sites' % bucket

vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'

rf_path = '%s/variantqc/gnomad.exomes.variantqc.vds' % bucket
vep_config = "/vep/vep-gcloud.properties"

# Outputs
out_vds_prefix = "%s/gnomad.sites" % root
out_internal_vcf_prefix = "%s/gnomad.sites.internal" % root
out_external_vcf_prefix = "%s/gnomad.sites" % root

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
rf_snv_cutoff = 0.1
rf_indel_cutoff = 0.2

#Actions
preprocess_autosomes = False
postprocess_autosomes = True
write_autosomes = True
preprocess_X = False
postprocess_X = False
write_X = False
preprocess_Y = False
postprocess_Y = False
write_Y = False

hc = HailContext()


def preprocess_vds(vds_path):
    print("Preprocessing %s\n" % vds_path)
    return (hc.read(vds_path)
            .annotate_global_py('global.pops', map(lambda x: x.lower(), pops), TArray(TString()))
            .annotate_samples_table(meta_path, 'sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
            .annotate_samples_expr(['sa.meta.project_description = sa.meta.description'])  # Could be cleaner
            .filter_variants_intervals('gs://gnomad-lfran/tmp/test.interval')
            .annotate_variants_intervals(decoy_path, 'va.decoy')
            .annotate_variants_intervals(lcr_path, 'va.lcr')
    )


if preprocess_autosomes:
    (
        create_sites_vds_annotations(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf)
        .write(out_vds_prefix + ".pre.vds")
    )

# TODO: fix the prefixes here
if postprocess_autosomes:
    post_process_vds(hc, out_vds_prefix + ".test.vds", rf_path, 'va.rf', 'va.train', 'va.label', rf_snv_cutoff, rf_indel_cutoff, vep_config).write(out_vds_prefix + ".vds")
    # .repartition(1000, shuffle=False)

if write_autosomes:
    for i in range(1, 23):
        write_vcfs(hc.read(out_vds_prefix + ".vds"), i, out_internal_vcf_prefix, out_external_vcf_prefix)

if preprocess_X:
    (
        create_sites_vds_annotations_X(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf)
        .write(out_vds_prefix + ".pre.X.vds")
    )

if postprocess_X:
    post_process_vds(hc, out_vds_prefix + ".pre.X.vds", rf_path, 'va.rf', 'va.train', 'va.label', rf_snv_cutoff, rf_indel_cutoff, vep_config).repartition(100, shuffle=False).write(out_vds_prefix + ".X.vds")

if write_X:
    write_vcfs(hc.read(out_vds_prefix + ".X.vds"), "X", out_internal_vcf_prefix, out_external_vcf_prefix)

if preprocess_Y:
    (
        create_sites_vds_annotations_Y(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf)
        .write(out_vds_prefix + ".pre.Y.vds")
    )

if postprocess_Y:
    post_process_vds(hc, out_vds_prefix + ".pre.Y.vds", rf_path, 'va.rf', 'va.train', 'va.label', rf_snv_cutoff, rf_indel_cutoff, vep_config).repartition(10, shuffle=False).write(out_vds_prefix + ".Y.vds")

if write_Y:
    write_vcfs(hc.read(out_vds_prefix + ".Y.vds"), "Y", out_internal_vcf_prefix, out_external_vcf_prefix)

send_message()
from utils import *
import time
from resources import *

# Inputs

vds_path = "gs://gnomad/gnom.ad.vds"
#vds_path = "gs://gnomad/gnomad.10ksites.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats24.vds"
vep_config = "/vep/vep-gcloud.properties"
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"

# Resources
autosomes_intervals = "gs://gnomad/autosomes.txt"

# Outputs
#date_time = time.strftime("%Y-%m-%d_%H-%M")
#date_time = '2017-02-13_21-40'
date_time = '2017-02-15_17-40'
out_root = "gs://gnomad-genomes/sites"
out_vds_prefix = "%s/gnomad.genomes.sites" % out_root
out_internal_vcf_prefix = "%s/vcf/gnomad.genomes.sites.internal" % out_root
out_external_vcf_prefix = "%s/vcf/gnomad.genomes.sites" % out_root
tmp_vds_prefix = "gs://gnomad-lfran/tmp/gnomad.sites.tmp." + date_time

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
RF_SNV_CUTOFF=0.4
RF_INDEL_CUTOFF=0.4

#Actions
preprocess_autosomes = False
postprocess_autosomes = False
write_autosomes = False
write_autosomes_matt = False
run_autosomes_checks = False
preprocess_X = False
postprocess_X = True
write_X = True
run_X_checks = True
preprocess_Y = False
postprocess_Y = False
write_Y = False
run_Y_checks = False


hc = HailContext(log='/site_auto.log')

def preprocess_vds(vds_path):
    print("Preprocessing %s\n" % vds_path)
    return (
        hc.read(vds_path)
            .annotate_global_py('global.pops',map(lambda x: x.lower(), pops), TArray(TString()))
            .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
            .annotate_samples_expr(['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop',
                                    'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
            .filter_samples_expr('sa.meta.keep')
            .annotate_variants_intervals(decoy_path, 'va.decoy')
            .annotate_variants_intervals(lcr_path, 'va.lcr')
            .annotate_variants_expr('va.info = drop(va.info, MQ0, RAW_MQ)')
    )

if preprocess_autosomes:
    (
        create_sites_vds_annotations(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf)
        .write(tmp_vds_prefix + ".vds")
    )

if postprocess_autosomes:
    post_process_vds(hc, tmp_vds_prefix + ".fixed.vds",
                     hc.read(rf_path, sites_only=True),
                     'va.RF1', 'va.train', 'va.label' ,
                     vep_config).write(out_vds_prefix + ".vds")

if write_autosomes:
    for i in range(1, 23):
        write_vcfs(hc.read(out_vds_prefix + ".vds"), i, out_internal_vcf_prefix, out_external_vcf_prefix, append_to_header=additional_vcf_header)

if write_autosomes_matt:
    vds = hc.read(out_vds_prefix + ".vds")
    vds = vds.annotate_variants_expr(
            'va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
    vds = vds.repartition(32,shuffle=False)
    vds.export_vcf(out_internal_vcf_prefix + ".autosomes.vcf.bgz", append_to_header=additional_vcf_header,parallel=True)

if run_autosomes_checks:
    run_sanity_checks(hc.read(out_vds_prefix + ".vds"),pops)

if preprocess_X:
    x_intervals = '%s/chrX.txt' % '/tmp'
    with open(x_intervals, 'w') as f:
        f.write('X:1-1000000000')

    vds = preprocess_vds(vds_path)
    vds = vds.filter_variants_intervals('file://' + x_intervals)

    vds.repartition(3000).write(tmp_vds_prefix + ".X.repartition.vds")

    (
        create_sites_vds_annotations_X(
            hc.read(tmp_vds_prefix + ".X.repartition.vds"),
            pops,
            dbsnp_path=dbsnp_vcf
        )
        .write(tmp_vds_prefix + ".X.vds")
    )

if postprocess_X:
    post_process_vds(hc, tmp_vds_prefix + ".X.vds",
                     hc.read(rf_path, sites_only=True),
                     'va.RF1','va.train', 'va.label',
                     vep_config).write(out_vds_prefix + ".X.vds")

if write_X:
    write_vcfs(hc.read(out_vds_prefix + ".X.vds"), "X",
               out_internal_vcf_prefix, out_external_vcf_prefix,
               append_to_header=additional_vcf_header,
               drop_fields=['Hemi_raw'])

if run_X_checks:
    run_sanity_checks(hc.read(out_vds_prefix + ".X.vds"),pops,X=True)

if preprocess_Y:
    (
        create_sites_vds_annotations_Y(
            preprocess_vds(vds_path),
            pops,
            dbsnp_path=dbsnp_vcf
        )
        .write(tmp_vds_prefix + ".Y.vds")
    )

if postprocess_Y:
    post_process_vds(hc, tmp_vds_prefix + ".Y.vds",
                     hc.read(rf_path, sites_only=True),
                     'va.RF1','va.train', 'va.label',
                     vep_config).write(out_vds_prefix + ".Y.vds")

if write_Y:
    write_vcfs(hc.read(out_vds_prefix + ".Y.vds"), "Y", out_internal_vcf_prefix, out_external_vcf_prefix, append_to_header=additional_vcf_header)

if run_Y_checks:
    run_sanity_checks(hc.read(out_vds_prefix + ".X.vds"), pops)
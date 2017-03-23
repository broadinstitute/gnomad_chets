from utils import *
import time
from resources import *
import traceback

# Inputs

vds_path = full_genome_vds
#vds_path = "gs://gnomad/gnomad.10ksites.vds"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats24.vds"
raw_hardcalls_path = "gs://gnomad/gnomad.raw_hardcalls.vds"

# Resources
autosomes_intervals = "gs://gnomad/autosomes.txt"

# Outputs
date_time = time.strftime("%Y-%m-%d_%H-%M")

print("Date time for temp files: %s\n" % date_time)

#Last X: date_time = '2017-02-21_12-20'
date_time = "2017-02-27_19-05"
out_root = "gs://gnomad-genomes/sites"
out_vds_prefix = "%s/gnomad.genomes.sites" % out_root
out_internal_vds_prefix = "%s/internal/gnomad.genomes.sites" % out_root
out_vds_prefix = out_internal_vds_prefix
out_external_vds_prefix = "gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites"
out_external_vds_prefix = "gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites"
out_internal_vcf_prefix = "%s/vcf/gnomad.genomes.sites.internal" % out_root
out_external_vcf_prefix = "%s/vcf/gnomad.genomes.sites" % out_root
tmp_vds_prefix = "gs://gnomad-lfran/tmp/gnomad.sites.multi.tmp." + date_time

out_matt_prefix = "gs://gnomad-genomes/sites/vcf/gnomad.genomes.sites"

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
RF_SNV_CUTOFF=0.4
RF_INDEL_CUTOFF=0.4

#Actions
preprocess_autosomes = False
postprocess_autosomes = False
write_autosomes = False
write_autosomes_matt = False
write_autosomes_coding = False
write_autosomes_public_vds = False
run_autosomes_checks = False
run_pre_calculate_metrics = False


preprocess_X = False
postprocess_X = False
write_X = True
write_X_coding = False
write_X_matt = False
write_X_public_vds = False
run_X_checks = False

vcf_filters = {
        'RF': 'isMissing(va.info.AS_FilterStatus) || '
              '(va.info.AS_FilterStatus.forall(x => !x.isEmpty) && va.info.AS_FilterStatus.exists(x => x.contains("RF")))',
        'AC0': '(va.info.AS_FilterStatus.forall(x => !x.isEmpty) && va.info.AS_FilterStatus.exists(x => x.contains("AC0")))',
        'SEGDUP': 'va.decoy',
        'LCR': 'va.lcr'
    }

def preprocess_vds(vds, vqsr_vds=None, release=True):
    print("Preprocessing %s\n" % vds_path)
    pre_vds = (vds
               .annotate_global_py('global.pops',map(lambda x: x.lower(), pops), TArray(TString()))
               .annotate_samples_table(genomes_meta, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
               .annotate_samples_expr(['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop',
                                       'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
               .filter_variants_expr('v.nAltAlleles > 1')
               .annotate_variants_intervals(decoy_path, 'va.decoy')
               .annotate_variants_intervals(lcr_path, 'va.lcr')
               .annotate_variants_expr('va.info = drop(va.info, MQ0, RAW_MQ)')
               )
    if release:
        return vds.filter_samples_expr('sa.meta.keep')
    else:
        return vds

if __name__ == '__main__':

    hc = hail.HailContext(log='/site_auto.log')

    #try:
    if preprocess_autosomes:
        (
            create_sites_vds_annotations(
                preprocess_vds(hc.read(vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
                .write(tmp_vds_prefix + ".vds")
        )

    if postprocess_autosomes:
        vds = post_process_vds(hc, tmp_vds_prefix + ".vds",
                         hc.read(rf_path, sites_only=True),
                         RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.RF1')
        vds = vds.annotate_variants_vds(hc.read('gs://gnomad-genomes/sites/internal/gnomad.genomes.sites.autosomes.vds'),'va.vep = vds.vep, va.info.CSQ = vds.info.CSQ')
        vds.write(out_vds_prefix + ".autosomes.vds", overwrite=True)

    if run_autosomes_checks:
        print(run_sanity_checks(hc.read(out_vds_prefix + ".autosomes.vds"), pops))

    if write_autosomes_matt:
        write_vcfs(hc.read(out_internal_vds_prefix + ".autosomes.vds"), '', out_matt_prefix + ".internal.autosomes.vds", out_matt_prefix + ".autosomes.vds",
                   RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   append_to_header=additional_vcf_header,
                   nchunks=500, export_internal=False)

    if write_autosomes:
        for i in range(1, 23):
            write_vcfs(hc.read(out_vds_prefix + ".autosomes.vds"), i, out_internal_vcf_prefix, out_external_vcf_prefix,
                       RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                       append_to_header=additional_vcf_header,
                       export_internal=False)

    if write_autosomes_coding:
        vds = hc.read(out_vds_prefix + ".autosomes.vds")
        vds = vds.filter_variants_intervals(exome_calling_intervals)
        write_vcfs(vds, '', out_internal_vcf_prefix + ".coding", out_external_vcf_prefix+ ".coding",
                   RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   append_to_header=additional_vcf_header,
                   export_internal=False)

    if write_autosomes_public_vds:
        vds = hc.read(out_vds_prefix + ".autosomes.vds")
        write_public_vds(hc,vds, out_internal_vds_prefix + ".autosomes.vds", out_external_vds_prefix + "release.autosomes.vds")

    if run_pre_calculate_metrics:
        vds = hc.read(out_vds_prefix + ".autosomes.vds")
        pre_calculate_metrics(vds, out_root + '/genome_precalculated_metrics.txt')
        send_snippet('#exac_browser', open(out_root + '/genome_precalculated_metrics.txt').read())

    if preprocess_X:
    #    vds.repartition(3000).write(tmp_vds_prefix + ".X.repartition.vds")
        vds = hc.read("gs://gnomad-lfran/tmp/gnomad.sites.tmp.2017-02-15_17-40.X.repartition.vds")
        vds = vds.filter_variants_expr('v.nAltAlleles > 1')
        (
            create_sites_vds_annotations_X(
                vds,
                pops,
                dbsnp_path=dbsnp_vcf
            )
                .write(tmp_vds_prefix + ".X.vds")
        )

    if postprocess_X:
        vds = post_process_vds(hc, tmp_vds_prefix + ".X.vds",
                         hc.read(rf_path, sites_only=True),
                         RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.RF1')

        vds = vds.annotate_variants_vds(hc.read('gs://gnomad-genomes/sites/internal/gnomad.genomes.sites.X.vds'),
                                        'va.vep = vds.vep, va.info.CSQ = vds.info.CSQ')
        vds.write(out_vds_prefix + ".X.vds", overwrite=True)

    if run_X_checks:
        print(run_sanity_checks(hc.read(out_vds_prefix + ".X.vds"), pops, contig="X"))

    if write_X_matt:
        write_vcfs(hc.read(out_internal_vds_prefix + ".X.vds"), 'X', out_matt_prefix + ".internal", out_matt_prefix,
                   RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   append_to_header=additional_vcf_header,
                   nchunks=50, export_internal=False)

    if write_X:
        write_vcfs(hc.read(out_vds_prefix + ".X.vds"), "X",
                   out_internal_vcf_prefix, out_external_vcf_prefix,
                   RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   append_to_header=additional_vcf_header,
                   export_internal=False)

    if write_X_coding:
        vds = hc.read(out_vds_prefix + ".X.vds")
        vds = vds.filter_variants_intervals(exome_calling_intervals)
        write_vcfs(vds, 'X', out_internal_vcf_prefix + ".coding", out_external_vcf_prefix+ ".coding",
                   RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   append_to_header=additional_vcf_header,
                   export_internal=False)

    if write_X_public_vds:
        vds = hc.read(out_vds_prefix + ".X.vds")
        write_public_vds(hc,vds, out_internal_vds_prefix + ".X.vds", out_external_vds_prefix + ".X.vds")

    send_message(channel='@laurent', message='Genomes are done processing!')

    #except Exception, e:
    #    send_message(channel='@laurent', message='Genomes failed processing :facepalm:\n```%s```' % traceback.print_stack())
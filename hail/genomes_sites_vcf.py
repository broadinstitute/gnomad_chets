from utils import *
import time
from resources import *

# Inputs

# vds_path = "gs://gnomad/gnom.ad.vds"
vds_path = "gs://gnomad/gnomad.10ksites.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats7.vds"
vep_config = "/vep/vep-gcloud.properties"

# Resources
autosomes_intervals = "gs://gnomad/autosomes.txt"

# Outputs
date_time = time.strftime("%Y-%m-%d_%H-%M")
#date_time = '2016-12-20_21-20'
out_root = "gs://gnomad-lfran/tmp"
out_vds_prefix = "%s/gnomad.sites.annotations.%s" % (out_root, date_time)
out_internal_vcf_prefix = "%s/gnomad.sites.internal" % out_root
out_external_vcf_prefix = "%s/gnomad.sites" % out_root
tmp_vds_prefix = "gs://gnomad-lfran/tmp/gnomad.sites.tmp." + date_time
tmp_RF_ann_out = 'gs://gnomad-lfran/tmp/gnomad.rf.ann.tmp.' + date_time + 'txt.bgz'
intervals_tmp = '/tmp'

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']

#Actions
preprocess_autosomes = True
postprocess_autosomes = True
write_autosomes = True
preprocess_X = False
postprocess_X = False
write_X = False
preprocess_Y = False
postprocess_Y = False
write_Y = False


hc = HailContext(log='/site_auto.log')


def preprocess_vds(vds_path):
    print("Preprocessing %s\n" % vds_path)
    return (
        hc.read(vds_path)
            .annotate_global_expr_by_sample('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
            .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=pyhail.TextTableConfig(impute=True))
            .annotate_samples_expr(['sa.meta.population = sa.meta.predicted_pop',
                                    'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
            .filter_samples_expr('!isMissing(sa.meta.predicted_pop)')  # Could be cleaner
            .filter_variants_intervals('gs://gnomad-lfran/tmp/test.interval')
            .annotate_variants_intervals(decoy_path, 'va.decoy')
            .annotate_variants_intervals(lcr_path, 'va.lcr')
    )


def post_process_vds(vds_path):
    print("Postprocessing %s\n" % vds_path)
    vep = hc.read(vep_path)

    filters = {
        'RF' : 'isMissing(va.info.AS_FilterStatus) || va.info.AS_FilterStatus.exists(x => x != "PASS")',
        'SEGDUP' : 'va.decoy',
        'LCR': 'va.lcr'
    }

    return (
            set_vcf_filters(hc, vds_path, rf_path, 'va.RF1', 0.5, 0.75, filters=filters,
                        filters_to_keep=['InbreedingCoefficient'], tmp_path='/tmp')
            .annotate_variants_vds(vep, code='va.info.CSQ = vds.csq')
            .vep(config=vep_config, csq=True, root='va.info.CSQ')
    )


def write_vcfs(vds_path, contig):
    print 'Writing VCFs for %s' % vds_path
    interval_path = '%s/%s.txt' % (intervals_tmp, str(contig))
    with open(interval_path, 'w') as f:
        f.write('%s:1-1000000000' % str(contig))

    (hc.read(vds_path)
     .filter_variants_intervals('file://' + interval_path)
     .export_vcf(out_internal_vcf_prefix + ".%s.vcf.bgz" % str(contig))
     .annotate_variants_expr(
        'va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
     .export_vcf(out_external_vcf_prefix + ".%s.vcf.bgz" % str(contig))
     )

if preprocess_autosomes:
    (
        create_sites_vds_annotations(
            preprocess_vds(vds_path),
            pops,
            tmp_path=intervals_tmp,
            dbsnp_path=dbsnp_vcf,
            npartitions=1000,
            shuffle=False)
        .write(tmp_vds_prefix + ".vds")
    )

if postprocess_autosomes:
    post_process_vds(tmp_vds_prefix + ".vds").write(out_vds_prefix + ".vds")

if write_autosomes:
    for i in range(1,22):
        write_vcfs(out_vds_prefix + ".vds", i)

if preprocess_X:
    (
        create_sites_vds_annotations_X(
            preprocess_vds(vds_path),
            pops,
            tmp_path=intervals_tmp,
            dbsnp_path=dbsnp_vcf,
            npartitions=100,
            shuffle=False)
        .write(tmp_vds_prefix + ".X.vds")
    )

if postprocess_X:
    post_process_vds(tmp_vds_prefix + ".X.vds").write(out_vds_prefix + ".X.vds")

if write_X:
    write_vcfs(out_vds_prefix + ".X.vds", "X")

if preprocess_Y:
    (
        create_sites_vds_annotations_Y(
            preprocess_vds(vds_path),
            pops,
            tmp_path=intervals_tmp,
            dbsnp_path=dbsnp_vcf,
            npartitions=10,
            shuffle=False)
        .write(tmp_vds_prefix + ".Y.vds")
    )

if postprocess_Y:
    post_process_vds(tmp_vds_prefix + ".Y.vds").write(out_vds_prefix + ".Y.vds")

if write_Y:
    write_vcfs(out_vds_prefix + ".Y.vds", "Y")
from pyhail import *
from utils import *
import time

# Inputs

# vds_path = "gs://gnomad/gnom.ad.vds"
vds_path = "gs://gnomad/gnomad.10ksites.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"
vep_path = "gs://gnomad/gnomad.splitmulti.vep.vds"
rf_path = "gs://gnomad/RF/gnomad.sites.RF.newStats7.vds"
vep_config = "/vep/vep-gcloud.properties"

# Resources
lcr_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
decoy_path = "gs://gnomad-lfran/annotations/LCR.interval_list"
autosomes_intervals = "gs://gnomad/autosomes.txt"
dbsnp = "gs://gnomad-lfran/All_20160601.vcf.bgz"

# Outputs
out_root = "gs://gnomad-lfran/tmp"
out_vds_prefix = "%s/gnomad.sites.annotations" % out_root
out_internal_vcf_prefix = "%s/gnomad.sites.internal" % out_root
out_external_vcf_prefix = "%s/gnomad.sites" % out_root
date_time = time.strftime("%Y-%m-%d_%H-%M")
tmp_vds_prefix = "gs://gnomad-lfran/tmp/gnomad.sites.tmp." + date_time
tmp_RF_ann_out = 'gs://gnomad-lfran/tmp/gnomad.rf.ann.tmp.' + date_time + 'txt.bgz'
intervals_tmp = '/tmp'

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']

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
            .annotate_global_expr_by_sample('global.pops=["%s"]' % '", "'.join(map(lambda x: x.lower(), pops)))
            .annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=pyhail.TextTableConfig(impute=True))
            .annotate_samples_expr(['sa.meta.population = sa.meta.predicted_pop',
                                    'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
            .filter_samples_expr('!isMissing(sa.meta.predicted_pop)')  # Could be cleaner
    )


def post_process_vds(vds_path):
    print("Postprocessing %s\n" % vds_path)
    vep = hc.read(vep_path)
    return (
        annotate_non_split_from_split(hc, non_split_vds_path=tmp_vds,
                                      split_vds=hc.read(rf_path),
                                      annotations=['va.RF1'],
                                      annotation_exp_out_path=tmp_RF_ann_out)
            .annotate_variants_intervals(decoy_path, 'va.decoy')
            .annotate_variants_expr(get_add_filter_annotation('segdup', 'va.decoy'))
            .annotate_variants_intervals(lcr_path, 'va.lcr')
            .annotate_variants_expr(get_add_filter_annotation('LCR', 'va.lcr'))
            .annotate_variants_vds(vep, code='va.info.CSQ = vds.csq')
            .vep(config=vep_config, csq=True, root='va.info.CSQ')
    )


def write_vcfs(vds_path, contig):
    print 'Writing VCFs for %s' % vds_path
    interval_path = '%s/%s.txt' % (intervals_tmp, str(contig))
    with open(interval_path, 'w') as f:
        f.write('%d:1-1000000000' % str(contig))

    (hc.read(vds_path)
     .filter_variants_intervals('file://' + interval_path)
     .export_vcf(out_internal_vcf_prefix + "%d.vcf.bgz" % str(contig))
     .annotate_variants_expr(
        'va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
     .export_vcf(out_external_vcf_prefix + "%d.vcf.bgz" % str(contig))
     )

if preprocess_autosomes:
    (
        create_sites_vds_annotations(
            preprocess_vds(vds_path),
            pops,
            tmp_path=intervals_tmp,
            dbsnp_path=dbsnp,
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
        create_sites_vds_annotationsX(
            preprocess_vds(vds_path),
            pops,
            tmp_path=intervals_tmp,
            dbsnp_path=dbsnp,
            npartitions=1000,
            shuffle=False)
        .write(tmp_vds_prefix + ".X.vds")
    )

if postprocess_X:
    post_process_vds(tmp_vds_prefix + ".X.vds").write(out_vds_prefix + ".X.vds")

if write_X:
    write_vcfs(out_vds_prefix + ".X.vds", "X")

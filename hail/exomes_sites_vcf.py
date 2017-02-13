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

vep_path = "%s/gnomad.splitmulti.vep.vds" % root
rf_path = '%s/variantqc/gnomad.exomes.variantqc.vds' % bucket
vep_config = "/vep/vep-gcloud.properties"

# Outputs
out_vds_prefix = "%s/gnomad.sites" % root
out_internal_vcf_prefix = "%s/gnomad.sites.internal" % root
out_external_vcf_prefix = "%s/gnomad.sites" % root
intervals_tmp = '/tmp'

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']

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


def post_process_vds(vds_path):
    print("Postprocessing %s\n" % vds_path)
    vep = hc.read(vep_path)

    filters = {
        'RF': 'isMissing(va.info.AS_FilterStatus) || va.info.AS_FilterStatus.forall(x => x != "PASS")',
        'SEGDUP': 'va.decoy',
        'LCR': 'va.lcr'
    }

    return (
        set_vcf_filters(hc, vds_path, rf_path, 'va.RF',
                        rf_snv_cutoff=0.1, rf_indel_cutoff=0.2, filters=filters,
                        filters_to_keep=['InbreedingCoefficient'], tmp_path='/tmp')
        .annotate_variants_vds(vep, code='va.info.CSQ = vds.csq')
        .vep(config=vep_config, csq=True, root='va.info.CSQ')
    )


def write_vcfs(vds_path, contig):
    print 'Writing VCFs for %s' % vds_path
    interval_path = '%s/%s.txt' % (intervals_tmp, str(contig))
    with open(interval_path, 'w') as f:
        f.write('%s:1-1000000000' % str(contig))

    vds = hc.read(vds_path).filter_variants_intervals('file://' + interval_path)
    vds.export_vcf(out_internal_vcf_prefix + ".%s.vcf.bgz" % str(contig))

    (
        vds.annotate_variants_expr(
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
        .write(out_vds_prefix + ".vds")
    )

if postprocess_autosomes:
    post_process_vds(out_vds_prefix + ".vds").write(out_vds_prefix + ".vds")

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
        .write(out_vds_prefix + ".X.vds")
    )

if postprocess_X:
    post_process_vds(out_vds_prefix + ".X.vds").write(out_vds_prefix + ".X.vds")

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
        .write(out_vds_prefix + ".Y.vds")
    )

if postprocess_Y:
    post_process_vds(out_vds_prefix + ".Y.vds").write(out_vds_prefix + ".Y.vds")

if write_Y:
    write_vcfs(out_vds_prefix + ".Y.vds", "Y")
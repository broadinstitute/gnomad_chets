from variantqc import *
from hail import *
import time

# Inputs

bucket = 'gs://gnomad-exomes'
autosomes_intervals = '%s/intervals/autosomes.txt' % bucket
# evaluation_intervals = '%s/intervals/exome_evaluation_regions.v1.intervals' % bucket
# high_coverage_intervals = '%s/intervals/high_coverage.auto.interval_list' % bucket
date_time = time.strftime("%Y-%m-%d_%H:%M")

root = '%s/sites' % bucket

vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
vqsr_vds_path = 'gs://gnomad-exomes/variantqc/gnomad.exomes.vqsr.unsplit.vds'
rf_path = '%s/variantqc/gnomad.exomes.rf.vds' % bucket

# Outputs
out_vds_prefix = "%s/internal/gnomad.exomes.sites" % root
out_internal_vcf_prefix = "%s/internal/gnomad.exomes.sites.internal" % root
out_external_vcf_prefix = "%s/vcf/gnomad.exomes.sites" % root

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
RF_SNV_CUTOFF = 0.1
RF_INDEL_CUTOFF = 0.2
send_to_slack = True

#Actions
run_all = False
run_pre = False
preprocess_autosomes = run_all or run_pre or False
preprocess_X = run_all or run_pre or False
preprocess_Y = run_all or run_pre or False
postprocess = False
write = False
run_pre_calculate_metrics = False


def preprocess_vds(vds, vqsr_vds, vds_pops=pops, release=True):
    print("Preprocessing %s\n" % vds_path)
    annotations = ['culprit', 'POSITIVE_TRAIN_SITE', 'NEGATIVE_TRAIN_SITE', 'VQSLOD']
    pre_vds = (vds
               .annotate_global_py('global.pops', map(lambda x: x.lower(), vds_pops), TArray(TString()))
               .annotate_samples_table(exomes_meta, 'sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
               .annotate_samples_expr(['sa.meta.project_description = sa.meta.description'])  # Could be cleaner
               .annotate_variants_intervals(decoy_path, 'va.decoy')
               .annotate_variants_intervals(lcr_path, 'va.lcr')
               .annotate_variants_vds(vqsr_vds, code=', '.join(['va.info.%s = vds.info.%s' % (a, a) for a in annotations]))
    )
    return pre_vds.filter_samples_expr('sa.meta.drop_status == "keep"') if release else pre_vds

if __name__ == '__main__':
    hc = HailContext()

    if preprocess_autosomes:
        (
            create_sites_vds_annotations(
                preprocess_vds(hc.read(vds_path), hc.read(vqsr_vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
            .write(out_vds_prefix + ".pre.autosomes.vds")
        )

    if preprocess_X:
        (
            create_sites_vds_annotations_X(
                preprocess_vds(hc.read(vds_path), hc.read(vqsr_vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
            .write(out_vds_prefix + ".pre.X.vds")
        )

    if preprocess_Y:
        (
            create_sites_vds_annotations_Y(
                preprocess_vds(hc.read(vds_path), hc.read(vqsr_vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
            .write(out_vds_prefix + ".pre.Y.vds")
        )

    if postprocess:
        auto_vds = hc.read(out_vds_prefix + ".pre.autosomes.vds")
        x_vds = hc.read(out_vds_prefix + ".pre.X.vds")
        y_vds = hc.read(out_vds_prefix + ".pre.Y.vds")
        vds = auto_vds.union([x_vds, y_vds])  # TODO: union schema

        rf_vds = hc.read(rf_path)
        post_process_vds(vds, rf_vds, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.rf').write(out_vds_prefix + ".post.vds", overwrite=True)

        vds = hc.read(out_vds_prefix + ".post.vds")
        sanity_check = run_sanity_checks(vds, pops, return_string=send_to_slack)
        if send_to_slack: send_snippet('#joint_calling', sanity_check, 'sanity_%s.txt' % date_time)

    if write:
        vds = hc.read(out_vds_prefix + ".post.vds").filter_variants_intervals(IntervalTree.read(exome_calling_intervals))
        write_vcfs(vds, '', out_internal_vcf_prefix, out_external_vcf_prefix, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header)
        write_public_vds(hc, vds, out_vds_prefix + ".internal.vds", out_external_vcf_prefix.replace('vcf', 'vds') + ".vds")

    if run_pre_calculate_metrics:
        vds = hc.read(out_external_vcf_prefix.replace('vcf', 'vds') + ".autosomes.vds")
        pre_calculate_metrics(vds, "exome_precalculated_metrics.txt")
        send_snippet('#exac_browser', open('exome_precalculated_metrics.txt').read())

    send_message(channel='@konradjk', message='Exomes are done processing!')

# zcat gnomad.exomes.sites.autosomes.vcf.bgz | head -250 | grep "^##" > header
# zcat gnomad.exomes.sites.X.vcf.bgz | head -250 | grep "^##" | while read i; do grep -F "$i" header; if [[ $? != 0 ]]; then echo $i >> header; fi; done
# Optional: nano header to move CSQ, contigs, and reference below X specific annotations
# cat header <(zcat gnomad.exomes.sites.autosomes.vcf.bgz | grep -v "^##") <(zcat gnomad.exomes.sites.X.vcf.bgz | grep -v "^#") <(zcat gnomad.exomes.sites.Y.vcf.bgz | grep -v "^#") | bgzip -c > gnomad.exomes.sites.vcf.gz
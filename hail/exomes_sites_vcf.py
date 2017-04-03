from variantqc import *
from hail import *
import time
import argparse

# Inputs

bucket = 'gs://gnomad-exomes'
autosomes_intervals = '%s/intervals/autosomes.txt' % bucket
date_time = time.strftime("%Y-%m-%d_%H:%M")

root = '%s/sites' % bucket

vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
vqsr_vds_path = 'gs://gnomad-exomes/variantqc/gnomad.exomes.vqsr.unsplit.vds'
rf_path = '%s/variantqc/gnomad.exomes.rf.vds' % bucket

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
RF_SNV_CUTOFF = 0.1
RF_INDEL_CUTOFF = 0.2


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


def main(args):
    out_vds_prefix = args.output
    out_internal_vcf_prefix = args.output + '.internal'
    out_external_vcf_prefix = args.output.replace('internal', 'vcf')
    if args.debug: logger.setLevel(logging.DEBUG)
    hc = HailContext()

    if not args.skip_preprocess_autosomes:
        (
            create_sites_vds_annotations(
                preprocess_vds(hc.read(vds_path), hc.read(vqsr_vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
            .write(out_vds_prefix + ".pre.autosomes.vds")
        )

    if not args.skip_preprocess_X:
        (
            create_sites_vds_annotations_X(
                preprocess_vds(hc.read(vds_path), hc.read(vqsr_vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
            .write(out_vds_prefix + ".pre.X.vds")
        )

    if not args.skip_preprocess_Y:
        (
            create_sites_vds_annotations_Y(
                preprocess_vds(hc.read(vds_path), hc.read(vqsr_vds_path)),
                pops,
                dbsnp_path=dbsnp_vcf)
            .write(out_vds_prefix + ".pre.Y.vds")
        )

    if not args.skip_postprocess:
        auto_vds = hc.read(out_vds_prefix + ".pre.autosomes.vds")
        x_vds = hc.read(out_vds_prefix + ".pre.X.vds")
        y_vds = hc.read(out_vds_prefix + ".pre.Y.vds")
        vds = auto_vds.union([x_vds, y_vds])  # TODO: union schema

        rf_vds = hc.read(rf_path)
        post_process_vds(vds, rf_vds, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.rf').write(out_vds_prefix + ".post.vds", overwrite=True)

        vds = hc.read(out_vds_prefix + ".post.vds")
        sanity_check = run_sanity_checks(vds, pops, return_string=True)
        if args.slack_channel: send_snippet(args.slack_channel, sanity_check, 'sanity_%s.txt' % date_time)

    if not args.skip_write:
        vds = hc.read(out_vds_prefix + ".post.vds").filter_variants_intervals(IntervalTree.read(exome_calling_intervals))
        write_vcfs(vds, '', out_internal_vcf_prefix, out_external_vcf_prefix, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header)
        write_public_vds(hc, vds, out_vds_prefix + ".internal.vds", out_external_vcf_prefix.replace('vcf', 'vds') + ".vds")

    if not args.skip_pre_calculate_metrics:
        vds = hc.read(out_external_vcf_prefix.replace('vcf', 'vds') + ".autosomes.vds")
        pre_calculate_metrics(vds, "exome_precalculated_metrics.txt")
        send_snippet('#exac_browser', open('exome_precalculated_metrics.txt').read())

    if args.slack_channel: send_message(channel=args.slack_channel, message='Exomes are done processing!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--skip_preprocess_autosomes', help='Skip pre-processing autosomes (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_X', help='Skip pre-processing X (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_Y', help='Skip pre-processing Y (assuming already done)', action='store_true')
    parser.add_argument('--skip_postprocess', help='Skip merge and post-process (assuming already done)', action='store_true')
    parser.add_argument('--skip_write', help='Skip writing data (assuming already done)', action='store_true')
    parser.add_argument('--skip_pre_calculate_metrics', help='Skip pre-calculating metrics (assuming already done)', action='store_true')
    parser.add_argument('--expr', help='''Additional expression (e.g. "!sa.meta.remove_for_non_tcga)"''')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    args = parser.parse_args()
    main(args)
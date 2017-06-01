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
rf_path = '%s/variantqc/gnomad.exomes.rf.vds' % bucket

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
RF_SNV_CUTOFF = 0.1
RF_INDEL_CUTOFF = 0.2


def preprocess_vds(vds, meta_kt, vqsr_vds, vds_pops=pops, release=True):
    print("Preprocessing %s\n" % vds_path)
    annotations = ['culprit', 'POSITIVE_TRAIN_SITE', 'NEGATIVE_TRAIN_SITE', 'VQSLOD']
    pre_vds = (vds
               .annotate_global('global.pops', map(lambda x: x.lower(), vds_pops), TArray(TString()))
               .annotate_samples_table(meta_kt, root='sa.meta')
               .annotate_samples_expr(['sa.meta.project_description = sa.meta.description'])  # Could be cleaner
               .annotate_variants_intervals(decoy_path, 'va.decoy')
               .annotate_variants_intervals(lcr_path, 'va.lcr')
               .annotate_variants_vds(vqsr_vds, code=', '.join(['va.info.%s = vds.info.%s' % (a, a) for a in annotations]))
    )
    return pre_vds.filter_samples_expr('sa.meta.drop_status == "keep"') if release else pre_vds


def main(args):
    if args.debug: logger.setLevel(logging.DEBUG)
    hc = HailContext()

    if not args.skip_preprocess_autosomes or args.skip_preprocess_X or args.skip_preprocess_Y:
        vds = hc.read(vds_path)
        vqsr_vds = hc.read(vqsr_vds_path)
        meta_kt = hc.import_table(exomes_meta, impute=True).key_by('sample')
        vds = preprocess_vds(vds, meta_kt, vqsr_vds)
        if args.expr:
            vds = vds.filter_samples_expr(args.expr)
        logger.info('Found %s samples', vds.query_samples('samples.count()'))

        if not args.skip_preprocess_autosomes:
            (
                create_sites_vds_annotations(vds, pops, dbsnp_path=dbsnp_vcf)
                .write(args.output + ".pre.autosomes.vds")
            )

        if not args.skip_preprocess_X:
            (
                create_sites_vds_annotations_X(vds, pops, dbsnp_path=dbsnp_vcf)
                .write(args.output + ".pre.X.vds")
            )

        if not args.skip_preprocess_Y:
            (
                create_sites_vds_annotations_Y(vds, pops, dbsnp_path=dbsnp_vcf)
                .write(args.output + ".pre.Y.vds")
            )

    if not args.skip_merge:
        # Combine VDSes
        vdses = [hc.read(args.output + ".pre.autosomes.vds"), hc.read(args.output + ".pre.X.vds"), hc.read(args.output + ".pre.Y.vds")]
        vdses = merge_schemas(vdses)
        vds = vdses[0].union(vdses[1:])
        vds.write(args.output + '.pre.vds', overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.vep.vds", overwrite=args.overwrite)
         )

    if not args.skip_postprocess:
        vds = hc.read(args.output + ".pre.vep.vds")
        rf_vds = hc.read(rf_path)
        post_process_vds(vds, rf_vds, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.rf').write(args.output + ".post.vds", overwrite=args.overwrite)

        vds = hc.read(args.output + ".post.vds")
        sanity_check = run_sites_sanity_checks(vds, pops)
        if args.slack_channel: send_snippet(args.slack_channel, sanity_check, 'sanity_%s.txt' % date_time)

    if not args.skip_write:
        exome_intervals = KeyTable.import_interval_list(exome_calling_intervals)
        vds = hc.read(args.output + ".post.vds").filter_variants_table(exome_intervals)
        write_vcfs(vds, '', args.output, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header)
        write_public_vds(vds, args.output + ".vds", overwrite=args.overwrite)

    if not args.skip_pre_calculate_metrics:
        vds = hc.read(args.output + ".vds")
        pre_calculate_metrics(vds, "exome_precalculated_metrics.txt")
        send_snippet('#exac_browser', open('exome_precalculated_metrics.txt').read())

    if args.slack_channel: send_message(channel=args.slack_channel, message='Exomes are done processing!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--skip_preprocess_autosomes', help='Skip pre-processing autosomes (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_X', help='Skip pre-processing X (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_Y', help='Skip pre-processing Y (assuming already done)', action='store_true')
    parser.add_argument('--skip_postprocess', help='Skip merge and post-process (assuming already done)', action='store_true')
    parser.add_argument('--skip_merge', help='Skip merging data (assuming already done)', action='store_true')
    parser.add_argument('--skip_vep', help='Skip VEPping data (assuming already done)', action='store_true')
    parser.add_argument('--skip_write', help='Skip writing data (assuming already done)', action='store_true')
    parser.add_argument('--skip_pre_calculate_metrics', help='Skip pre-calculating metrics (assuming already done)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    parser.add_argument('--expr', help='''Additional expression (e.g. "!sa.meta.remove_for_non_tcga)"''')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    args = parser.parse_args()
    main(args)
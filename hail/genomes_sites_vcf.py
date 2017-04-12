from utils import *
import time
from resources import *
import argparse
import os

# Inputs
vds_path = full_genome_vds

# Outputs
date_time = time.strftime("%Y-%m-%d_%H-%M")

#Config
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
RF_SNV_CUTOFF=0.4
RF_INDEL_CUTOFF=0.4

vcf_filters = {
        'RF': 'isMissing(va.info.AS_FilterStatus) || '
              '(va.info.AS_FilterStatus.forall(x => !x.isEmpty) && va.info.AS_FilterStatus.exists(x => x.contains("RF")))',
        'AC0': '(va.info.AS_FilterStatus.forall(x => !x.isEmpty) && va.info.AS_FilterStatus.exists(x => x.contains("AC0")))',
        'SEGDUP': 'va.decoy',
        'LCR': 'va.lcr'
    }

def preprocess_vds(vds, vqsr_vds=None, vds_pops=pops, release=True):
    print("Preprocessing %s\n" % vds_path)
    vds = (vds
               .annotate_global_py('global.pops',map(lambda x: x.lower(), vds_pops), TArray(TString()))
               .annotate_samples_table(genomes_meta, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
               .annotate_samples_expr(['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop',
                                       'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
               .annotate_variants_intervals(decoy_path, 'va.decoy')
               .annotate_variants_intervals(lcr_path, 'va.lcr')
               .annotate_variants_expr('va.info = drop(va.info, MQ0, RAW_MQ)')
               )
    if release:
        vds = vds.filter_samples_expr('sa.meta.keep')

    return vds

def main(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hc = hail.HailContext(log='/site_auto.log')

    vds = None

    if not args.skip_preprocess_autosomes:
        vds = preprocess_vds(hc.read(vds_path))
        if args.expr:
            vds = vds.filter_samples_expr(args.expr)

        samples_info = vds.query_samples(['samples.count()',
                                          'samples.map(s => sa.meta.population).counter()'])

        logger.info("Pre-proocessing done, found: %d samples.\n Populations: %s" %
                    (samples_info[0],
                     ", ".join(["%s: %s" % (k, v) for k, v in samples_info[1].iteritems()])))

        (
            create_sites_vds_annotations(
                vds,
                pops,
                dbsnp_path=dbsnp_vcf)
                .write(args.output + ".pre.autosomes.vds", overwrite=args.overwrite)
        )

    if not args.skip_preprocess_X:
        vds = preprocess_vds(hc.read(vds_path))
        if args.expr:
            vds = vds.filter_samples_expr(args.expr)
        (
            create_sites_vds_annotations_X(
                vds ,
                pops,
                dbsnp_path=dbsnp_vcf
            )
                .write(args.output + ".pre.X.vds", overwrite=args.overwrite)
        )

    if not args.skip_merge:
        # Combine VDSes
        vdses = [hc.read(args.output + ".pre.autosomes.vds"), hc.read(args.output + ".pre.X.vds")]
        vdses = merge_schemas(vdses)
        sites_vds = vdses[0].union(vdses[1:])
        sites_vds.write(args.output + '.pre.vds', overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.vep.vds", overwrite=args.overwrite)
         )

    if not args.skip_postprocess:
        vds = post_process_vds(hc.read(args.output + '.pre.vep.vds'),
                         hc.read(genomes_rf_path, sites_only=True),
                         RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.RF1')
        #vds = vds.annotate_variants_vds(hc.read('gs://gnomad-genomes/sites/internal/gnomad.genomes.sites.autosomes.vds'),'va.vep = vds.vep, va.info.CSQ = vds.info.CSQ')
        vds.write(args.output + ".vds", overwrite=True)

    if not args.skip_sanity_checks:
        sanity_check_text = run_sites_sanity_checks(hc.read(args.output + ".vds"), pops, return_string=True, skip_star=False)
        if args.slack_channel:
            send_snippet(args.slack_channel, sanity_check_text,
                         'sanity_%s_%s.txt' % (os.path.basename(args.output), date_time))
        else:
            logger.info(sanity_check_text)

    if not (args.skip_write_public_vcf or args.skip_write_internal_vcf) :
        if vds is None:
            vds = hc.read(args.output + ".vds")

        public_out = args.output
        if args.skip_write_public_vcf:
            public_out = None
        internal_out = args.output + ".internal"
        if args.skip_write_internal_vcf:
            internal_out = None

        if args.write_vcf_per_chrom:
            for contig in range(1, 23):
                write_vcfs(vds, contig, public_out, internal_out, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                           append_to_header=additional_vcf_header)
            write_vcfs(vds, 'X', public_out, internal_out, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                       append_to_header=additional_vcf_header)
        else:
            write_vcfs(vds, '', public_out, internal_out, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                       append_to_header=additional_vcf_header)

    if not args.skip_write_public_vds:
        if vds is None:
            vds = hc.read(args.output + ".vds")
        write_public_vds(vds, args.output + ".public.vds")


    if not args.skip_pre_calculate_metrics:
        vds = hc.read(args.output + ".vds")
        pre_calculate_metrics(vds, args.output + '.genome_precalculated_metrics.txt')
        if args.slack_channel:
            send_snippet('#exac_browser', open(args.output + '.genome_precalculated_metrics.txt').read())

    if args.slack_channel:
        send_message(channel=args.slack_channel, message='Genomes are done processing!')
    else:
        send_message(channel='@laurent', message='Genomes are done processing!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--skip_preprocess_autosomes', help='Skip pre-processing autosomes (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_X', help='Skip pre-processing X (assuming already done)', action='store_true')
    parser.add_argument('--skip_merge', help='Skip merge step (assuming already done)', action='store_true')
    parser.add_argument('--skip_vep', help='Skip VEP (assuming already done)', action='store_true')
    parser.add_argument('--skip_postprocess', help='Skip merge and post-process (assuming already done)', action='store_true')
    parser.add_argument('--skip_sanity_checks', help='Skip sanity checks', action='store_true')
    parser.add_argument('--skip_write_public_vds', help='Skip writing public VDS (assuming already done)', action='store_true')
    parser.add_argument('--skip_write_public_vcf', help='Skip writing final VCF(s) (assuming already done)', action='store_true')
    parser.add_argument('--skip_write_internal_vcf', help='Skip writing final VCF(s) (assuming already done)',
                        action='store_true')
    parser.add_argument('--skip_pre_calculate_metrics', help='Skip pre-calculating metrics (assuming already done)', action='store_true')
    parser.add_argument('--write_vcf_per_chrom',
                        help='If set, generates a VCF for each chromosome. Otherwise, creates a single VCF.',
                        action='store_true')
    parser.add_argument('--expr', help='''Additional expression (e.g. "!sa.meta.remove_for_non_tcga)"''')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    args = parser.parse_args()
    main(args)
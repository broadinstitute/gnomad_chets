__author__ = 'gtiao'


from hail import *
from tabulate import tabulate
import argparse
import sys
import copy
from resources import additional_vcf_header_path, exome_calling_intervals_path
from utils import *
from sites_vcf import *
import logging

# NOTES: input is hard-coded to be public release sites files for v2.0.1; output prefixes are hard-coded also


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("make_release_2.0.2_sites")
logger.setLevel(logging.INFO)


def remove_filters_set_flags(vds, remove_filters_dict):
    flags = ['va.info.%s = va.filters.contains("%s")' % (x, y) for x, y in remove_filters_dict.items()]
    flags.extend(['va.old_filter = va.filters', 'va.filters = va.filters.filter(x => !["{}"].toSet.contains(x))'.format('","'.join(remove_filters_dict.values()))])
    vds = vds.annotate_variants_expr(flags)
    return vds

remove_filters_dict = {'lcr': 'LCR', 'segdup': 'SEGDUP'}


# Remove lcr and segdup annotations from header definitions:
for anno in remove_filters_dict.values():
    del FILTERS_DESC[anno]

va_attr = {
    'lcr': {"Number": "0", "Description": "In a low complexity region"},
    'segdup': {"Number": "0", "Description": "In a segmental duplication region"}
}

rf_ann_expr = ['va.info.AS_RF_NEGATIVE_TRAIN = isDefined(va.info.AS_RF_NEGATIVE_TRAIN) && va.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(va.aIndex)',
    'va.info.AS_RF_POSITIVE_TRAIN = isDefined(va.info.AS_RF_POSITIVE_TRAIN) && va.info.AS_RF_POSITIVE_TRAIN.toSet.contains(va.aIndex)']

# Function to check for annotations to remove
def discover_drop_annotations(vds):
    # Consider only INFO fields to drop
    schema = flatten_struct(vds.variant_schema, root="va")
    info_annotations = [k[8:] for k in schema.keys() if k.startswith('va.info')]

    anno_counts = vds.query_variants(['variants.filter(v => isDefined(va.info.{}).count()'.format(ann) for ann in info_annotations])
    return [ann for ann, count in zip(info_annotations, anno_counts) if count == 0]



def main(args):
    hc = HailContext(log='/hail.sites_vcf.log')

    # TODO: update paths in resources file
    # vds_path = final_genome_vds_path if args.genomes else final_exome_vds_path
    vds_path = 'gs://gnomad-public/release/2.0.1/vds/genomes/gnomad.genomes.r2.0.1.sites.vds' if args.genomes else 'gs://gnomad-public/release/2.0.1/vds/exomes/gnomad.exomes.r2.0.1.sites.vds'
    out_external_vcf_prefix = 'gs://gnomad/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites' if args.genomes else 'gs://gnomad/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites'
    RF_SNV_CUTOFF = 0.4 if args.genomes else 0.1
    RF_INDEL_CUTOFF = 0.4 if args.genomes else 0.2
    overwrite_vds = True if args.overwrite_vds else False

    # Import exome/genome sites release VDS; re-write filter columns and introduce new flags
    vds = hc.read(vds_path)
    if args.coding_only:
        exome_intervals = KeyTable.import_interval_list(exome_calling_intervals_path)
        vds = vds.filter_variants_table(exome_intervals)
        out_external_vcf_prefix = 'gs://gnomad/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.coding_only' if args.genomes else 'gs://gnomad/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.coding_only'
    vds = remove_filters_set_flags(vds, remove_filters_dict)

    # Set header annotation with new LCR and segdup filters
    logger.info("Converting requested FILTER values into flags...")
    for ann in remove_filters_dict.keys():
        vds = vds.set_va_attributes("va.info.%s" % ann, va_attr[ann])

    # Count new columns for specific contig
    logger.info("Counting FILTER values converted into flags...")
    kt = vds.variants_table()
    check_cols = ['n{0} =  va.filter(x => x.info.{0}).count()'.format(x) for x in remove_filters_dict.keys()]
    kt = kt.aggregate_by_key('Old_Filter = va.old_filter.mkString("|")',
                             ['Filter = va.flatMap(x => x.filters).collect().toSet().mkString("|")',
                              'Count = va.count()'] + check_cols)
    kt = kt.persist()
    out_external_report = out_external_vcf_prefix + ".filter_counts_post_removal.txt"
    kt.export(out_external_report)

    # Print report to screen
    report_pd = kt.to_pandas()
    print "Counts of filter annotations, before and after removal of flags:"
    print tabulate(report_pd, headers='keys', tablefmt='psql')

    # Write new release sites vds (if desired)
    if args.write_new_sites_vds:
        logger.info("Dropping unused annotations and writing new sites-only vds...")
        drop_anno = discover_drop_annotations(vds)
        new_vds_path = out_external_vcf_prefix + '.vds'
        if drop_anno:
            vds.annotate_variants_expr('va.info = drop(va.info, %s)' % ", ".join(drop_anno)).write(new_vds_path,
                                                                                                   overwrite_vds)
        else:
            vds.write(new_vds_path, overwrite_vds)

    # Export VCFs for release
    # NOTE: out_external_vcf_prefix is supplied where out_internal_vcf_prefix is normally supplied, to avoid the PROJECTMAX
    # operations in the code (the release VDS being used here has no PROJECTMAX annotations)

    # Write out VCFs
    if args.write_vcf_per_chrom:
        logger.info("Writing new VCFs by chromosome...")
        contigs = [str(c) for c in range(1, 23)] + ["X"]
        if args.exomes:
            contigs.append("Y")
    elif args.coding_only:
        logger.info("Writing new coding-only VCFs by autosomes and X...")
        contigs = ["1-22", "X"]
    else:
        logger.info("Writing new VCF with all contigs...")
        contigs = [""]

    for contig in contigs:
        vds_contig = vds.filter_intervals(Interval.parse(str(contig)))
        drop_anno = discover_drop_annotations(vds_contig)
        if args.write_vcf_per_chrom or args.coding_only:
            logger.info("Dropping %s annotations from unsplit vds for contig(s) %s", (drop_anno, contig))
        else:
            logger.info("Dropping %s annotations from unsplit vds", drop_anno)
        write_vcfs(vds, contig, out_external_vcf_prefix, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                   drop_fields=drop_anno,
                   as_filter_status_fields=['va.info.AS_FilterStatus'], append_to_header=additional_vcf_header_path)

    # Write new release sites vds with multiallelics split (if desired)
    if args.write_new_sites_vds_split:
        logger.info("Splitting multi-allelics and writing new sites-only vds...")
        new_vds_path = out_external_vcf_prefix + '.split.vds'
        vds = split_vds_and_annotations(vds, rf_ann_expr)
        drop_anno = discover_drop_annotations(vds)
        if drop_anno:
            logger.info("Dropping %s annotations from multiallelic vds", drop_anno)
            vds.annotate_variants_expr('va.info = drop(va.info, %s)' % ", ".join(drop_anno)).write(new_vds_path,
                                                                                                   overwrite_vds)
        else:
            vds.write(new_vds_path, overwrite_vds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--write_new_sites_vds', help='Write out a sites-only release VDS', action='store_true')
    parser.add_argument('--write_new_sites_vds_split',
                        help='Write out a sites-only release VDS with multi-allelics split', action='store_true')
    parser.add_argument('--write_vcf_per_chrom', help='Write out sites VCF by chromosome', action='store_true')
    parser.add_argument('--coding_only', help='Output results in coding regions only', action='store_true')
    parser.add_argument('--overwrite_vds', help='Overwrite pre-existing VDS', action='store_true')

    args = parser.parse_args()

    if int(args.coding_only) + int(args.exomes) != 1:
        sys.exit('Error: One and only one of --coding_only or --exomes must be specified')
    if int(args.coding_only) + int(args.write_vcf_per_chrom) > 1:
        sys.exit('Error: Cannot use both  --coding_only and --write_vcf_per_chrom')
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    main(args)
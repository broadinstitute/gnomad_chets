__author__ = 'gtiao'


from hail import *
from tabulate import tabulate
import argparse
import sys
import copy
from resources import additional_vcf_header_path
from utils import flatten_struct, set_filters_attributes, ADJ_GQ, ADJ_DP, ADJ_AB, FILTERS_DESC
from sites_py import write_vcfs
import logging


#NOTES: input is hard-coded to be public release sites files for v2.0.1; output prefixes are hard-coded also


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("make_release_2.0.2_sites")
logger.setLevel(logging.INFO)


def remove_filters_set_flags(vds, remove_filters_dict):
    flags = ['va.info.%s = va.filters.contains("%s")' % (x, y) for x, y in remove_filters_dict.items()]
    vds = vds.annotate_variants_expr(flags + ['va.old_filter = va.filters', 'va.filters = va.filters.filter(x => !["{}"].toSet.contains(x))'.format('","'.join(remove_filters_dict.values()))])
    return vds

remove_filters_dict = {'lcr': 'LCR', 'segdup': 'SEGDUP'}


#Remove lcr and segdup annotations from header definitions:
for anno in remove_filters_dict.values():
    del FILTERS_DESC[anno]

va_attr = {
    'lcr': {"Number": "0", "Description": "In a low complexity region"},
    'segdup': {"Number": "0", "Description": "In a segmental duplication region"}
}


def main(args):
    hc = HailContext(log='/hail.sites_vcf.log')

    #Set parameters
    vds_path = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.vds' if args.genomes else 'gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.vds'
    out_external_vcf_prefix = 'gs://gnomad/release_2.0.2/gnomad.genomes.r2.0.2.sites' if args.genomes else 'gs://gnomad/release_2.0.2/gnomad.exomes.r2.0.2.sites'
    RF_SNV_CUTOFF = .4 if args.genomes else 0.1
    RF_INDEL_CUTOFF = .4 if args.genomes else 0.2


    #Import exome/genome sites release VDS; re-write filter columns and introduce new flags
    vds = hc.read(vds_path)
    vds = remove_filters_set_flags(vds, remove_filters_dict)

    #Set header annotation with new LCR and segdup filters
    for ann in remove_filters_dict.keys():
        vds = vds.set_va_attributes("va.info.%s" % ann, va_attr[ann])

    #Count new columns for specific contig
    table = vds.variants_table()
    check_cols = ['n%s =  va.filter(x => x.info.%s).count()' % (x, x) for x in remove_filters_dict.keys()]
    table = table.aggregate_by_key('Old_Filter = va.old_filter.mkString("|")', ['Filter = va.flatMap(x => x.filters).collect().toSet().mkString("|")', 'Count = va.count()'] + check_cols)
    out_external_report = out_external_vcf_prefix + ".filter_counts_post_removal.txt"
    table.export(out_external_report)

    #Print report to screen
    df = table.to_pandas()
    print tabulate(df, headers='keys', tablefmt='psql')

    #Write new release sites vds (if desired)
    if args.write_new_sites_vds:
        new_vds_path = out_external_vcf_prefix + '.vds'
        vds.annotate_variants_expr('va.info = drop(va.info, old_filter)').write(new_vds_path, overwrite=True)

    #Export VCFs for release
    #NOTE: out_external_vcf_prefix is supplied where out_internal_vcf_prefix is normally supplied, to avoid the PROJECTMAX
    #operations in the code (the release VDS being used here has no PROJECTMAX annotations)


    #TODO: write function to discover empty info annotations; amend drop_fields to contain all fields to drop;
    # and remove old_filter
    if args.write_vcf_per_chrom:
        for contig in range(1, 23):
            write_vcfs(vds, contig, out_external_vcf_prefix, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, drop_fields = ['old_filter'],
                as_filter_status_fields=['va.info.AS_FilterStatus'], append_to_header = additional_vcf_header_path)
        write_vcfs(vds, 'X', out_external_vcf_prefix, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, drop_fields = ['old_filter'],
                as_filter_status_fields=['va.info.AS_FilterStatus'], append_to_header = additional_vcf_header_path)
        if args.exomes:
            write_vcfs(vds, 'Y', out_external_vcf_prefix, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, drop_fields = ['old_filter'],
                as_filter_status_fields=['va.info.AS_FilterStatus'], append_to_header = additional_vcf_header_path)
    else:
        write_vcfs(vds, '', out_external_vcf_prefix, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, drop_fields = ['old_filter'],
            as_filter_status_fields=['va.info.AS_FilterStatus'], append_to_header = additional_vcf_header_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--write_new_sites_vds', help='Skip pre-processing autosomes (assuming already done)', action='store_false')
    parser.add_argument('--write_vcf_per_chrom', help='Write out sites VCF by chromosome', action='store_true', default=False)

    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    main(args)
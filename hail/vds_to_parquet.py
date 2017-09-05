__author__ = 'konradk'

from utils import *
import argparse

try:
    hc = HailContext(log="/hail.log")
except Exception:
    pass


def main(args):
    vds_path = final_exome_split_vds_path if args.exomes else final_genome_split_vds_path

    vds = hc.read(vds_path)
    vds = process_consequences(vds)
    vds = vds.annotate_variants_expr('va.worst_csq_suffix = va.vep.worst_csq_suffix, va.worst_csq = va.vep.worst_csq, va.info = drop(va.info, CSQ), va.pass = va.filters.isEmpty')
    vds = vds.annotate_variants_expr('va = drop(va, vep)')
    vds = vds.annotate_variants_expr('va.filters = va.filters.toArray().sort().mkString("|"), va.info.AS_FilterStatus = va.info.AS_FilterStatus.toArray().sort().mkString("|")')

    kt = vds.variants_table()
    kt = kt.annotate('v = {chrom: v.contig, pos: v.start, ref: v.ref, alt: v.altAlleles[0].alt}')
    kt = kt.flatten()
    kt = kt.rename({x.name: x.name.replace('.', '_').strip('`') for x in kt.schema.fields})
    df = kt.to_dataframe().write

    if args.overwrite:
        df = df.mode('overwrite')
    df.parquet(args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--exomes', help='Input gnomAD exomes', action='store_true')
    parser.add_argument('--genomes', help='Input gnomAD genomes', action='store_true')
    parser.add_argument('--output', help='Output file')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --input, --exomes, or --genomes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
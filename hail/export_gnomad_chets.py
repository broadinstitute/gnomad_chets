import argparse
from compound_hets_utils import *
from hail import *
from collections import OrderedDict
import re


def write_or_export_kt(kt, output, overwrite, part = None, export=False, non_overlapping_gene_intervals=None): #TODO -- make sure that the export of intervals is clean

    filename = '{}{}{}'.format(
        output,
        '.part{}'.format(part) if part is not None else '',
        '' if export else '.kt'
    )

    if export:
        export_variant_pairs(kt, filename)
        if non_overlapping_gene_intervals:
            with hadoop_write("{}.genes.part{}.list".format(output, part)) as gene_list_file:
                for row in non_overlapping_gene_intervals:
                    gene_list_file.write('{}\t{}\t{}\t{}\t{}\n'.format(row.gene, row.interval.start.contig, row.interval.start.position, row.interval.end.position, row.n))
    else:
        kt.write(filename, overwrite=overwrite)


def get_sorted_gene_intervals(vds, tmp_file=None):

    def cmp_interval(x, y):
        if (x.interval.start.contig == y.interval.start.contig):
            return cmp(x.interval.start.position, y.interval.start.position)
        else:
            xnum = re.search(r"^\d+", x.interval.start.contig)
            ynum = re.search(r"^\d+", y.interval.start.contig)

            if xnum:
                if ynum:
                    return cmp(int(x.interval.start.contig), int(y.interval.start.contig))
                else:
                    return -1
            elif ynum:
                return 1
        return cmp(x, y)

    kt = None

    if tmp_file:
        try:
            kt = vds.hc.read_table(tmp_file)
        except:
            pass

    if kt is None:
        kt = vds.variants_table()
        kt = kt.annotate('gene = va.all_genes').explode('gene')
        kt = kt.select(['gene', 'v'])
        kt = kt.repartition(100)
        kt = kt.aggregate_by_key(["chrom = v.contig","gene = gene"], [
            'start = v.takeBy(v => v.start, 1).map(v => Locus(v.contig, v.start))[0]',
            'end = v.takeBy(v => -v.start, 1).map(v => Locus(v.contig, v.start))[0]',
            'n = v.count()'
        ])
        kt = kt.filter('end.position - start.position > 0')
        kt = kt.annotate('interval = Interval(start,end)')
        kt = kt.select(['gene','interval','n'])

        if tmp_file:
            kt.write(tmp_file)
            kt = vds.hc.read_table(tmp_file)

    gene_intervals = kt.collect()
    longest_interval = (0,'')
    max_variants = (0, '')
    for x in gene_intervals:
        interval_length = x.interval.end.position - x.interval.start.position
        if interval_length > longest_interval[0]:
            longest_interval = (interval_length, x.gene)
        if x.n > max_variants[0]:
            max_variants = (x.n, x.gene)

    logger.info('Found {} genes. Longest gene: {} ({}bp). Gene with most variants: {} ({} variants).'.format(
        len(gene_intervals),
        longest_interval[1],
        longest_interval[0],
        max_variants[1],
        max_variants[0]))

    gene_intervals.sort(cmp=cmp_interval)

    logger.info("Found {} gene intervals.".format(len(gene_intervals)))

    return gene_intervals


def get_non_overlapping_gene_intervals(gene_intervals):

    non_overlapping_gene_intervals = [gene_intervals.pop(0)]
    overlapping_intervals = []

    for gene_interval in gene_intervals:
        if gene_interval.interval.start.contig > non_overlapping_gene_intervals[-1].interval.end.contig or gene_interval.interval.start.position > non_overlapping_gene_intervals[-1].interval.end.position:
            non_overlapping_gene_intervals.append(gene_interval)
        else:
            overlapping_intervals.append(gene_interval)

    return non_overlapping_gene_intervals, overlapping_intervals


def get_gene_variant_pairs(vds, gene_intervals):
    gene_intervals_kt = KeyTable.from_py(hc=vds.hc,
                                         rows_py=gene_intervals,
                                         schema=TStruct(['gene','interval'],[TString(), TInterval()]),
                                         key_names = ['interval'],
                                         num_partitions = 1)

    vds = vds.filter_variants_table(gene_intervals_kt)
    vds = vds.annotate_variants_table(gene_intervals_kt, expr='va = {gene: table}')

    return vds.phase_em(va_keys=['va.gene'],
                        sa_keys=['sa.pop'],
                        num_partitions=len(gene_intervals)*2,
                        per_sample=False)


def export_variant_pairs(kt, output):

    variant_cols = OrderedDict([('chrom', 'v{}.contig'),
                   ('pos', 'v{}.start'),
                   ('ref', 'v{}.ref'),
                   ('alt', 'v{}.alt')
                    ])

    # kt = conditional_column_swap(kt,
    #                              swap_expr='v1.start > v2.start',
    #                              columns=[('v1','v2')],
    #                              gt_counts_col="genotype_counts",
    #                              hc_counts_col="haplotype_counts")


    kt = kt.annotate(['{0}{1} = {2}'.format(name, n, expr.format(n)) for n in ["1", "2"] for name, expr in variant_cols.iteritems()])
    kt, gc_flatten_cols = flatten_haplotype_counts(kt, hc_col="", out_prefix='g', numeric_naming=True)
    kt, hc_flatten_cols = flatten_haplotype_counts(kt, gc_col="", out_prefix='h', numeric_naming=True)
    kt = kt.drop(['v1','v2','va1', 'va2', 'genotype_counts', 'haplotype_counts'])
    kt.export(output + '.tsv.bgz', parallel=True)


def main(args):

    hc = HailContext(log='/gnomad_compound_hets.log')

    if args.write_pairs_kt or args.export_pairs_no_kt:

        data_type = "exomes" if args.exomes else "genomes"
        vds = get_gnomad_data(hc,
                              data_type,
                              hardcalls="adj",
                              split=True,
                              release_samples=True)

        vds = vds.annotate_samples_expr('sa = {{pop: {} }}'.format('sa.meta.population' if args.exomes else 'if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop'))

        if args.chrom20:
            vds = vds.filter_intervals(Interval.parse("20:1-10000000"))

        vds = vds.annotate_variants_vds(get_gnomad_public_data(hc, data_type=data_type, split=True),
                                        root='va')

        vds = vds.annotate_global('global.coding_csq',
                                  set(CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT),
                                  TSet(TString()))
        vds = vds.filter_variants_expr('isDefined(va.vep.transcript_consequences) && va.vep.transcript_consequences.exists(t => t.biotype == "protein_coding" && t.consequence_terms.exists(c => global.coding_csq.contains(c)))')
        vds = vds.annotate_variants_expr('va.all_genes = va.vep.transcript_consequences.filter(t => t.biotype == "protein_coding" && t.consequence_terms.exists(c => global.coding_csq.contains(c))).map(t => t.gene_symbol).toSet')
        if args.pops:
            vds = vds.filter_samples_expr('["{}"].toSet().contains(sa.pop)'.format('","'.join(args.pops.split(","))))
        vds = vds.persist()

        gene_intervals = get_sorted_gene_intervals(vds, "{}.gene_intervals.kt".format(args.output))
        if args.input_genes:
            input_genes = hc.import_table(args.input_genes).query('gene.collect()')
            gene_intervals = [g for g in gene_intervals if g.gene in input_genes]
            logger.info("{} gene intervals lefts after filter for input gene list".format(len(gene_intervals)))

        non_overlapping_gene_intervals, gene_intervals = get_non_overlapping_gene_intervals(gene_intervals)

        logger.info("Processing {} genes. {} genes overlapping those intervals".format(len(non_overlapping_gene_intervals),
                                                                                       len(gene_intervals)))

        variant_pairs = get_gene_variant_pairs(vds, non_overlapping_gene_intervals)
        part = None
        if args.write_by_part:
            write_or_export_kt(variant_pairs, args.output, args.overwrite, 0, args.export_pairs_no_kt, non_overlapping_gene_intervals)
            part = 1

        while(gene_intervals):
            non_overlapping_gene_intervals, gene_intervals = get_non_overlapping_gene_intervals(gene_intervals)
            logger.info("Processing {} genes. {} genes overlapping those intervals".format(len(non_overlapping_gene_intervals),
                                                                                               len(gene_intervals)))
            kt = get_gene_variant_pairs(vds, non_overlapping_gene_intervals)

            if args.write_by_part:
                write_or_export_kt(kt, args.output, args.overwrite, part, args.export_pairs_no_kt, non_overlapping_gene_intervals)
                part += 1
            else:
                variant_pairs = variant_pairs.union(kt)

        if not args.write_by_part:
            write_or_export_kt(variant_pairs, args.output, args.overwrite, part, args.export_pairs_no_kt)

    if args.export_pairs:
        export_variant_pairs(hc.read_table(args.output + ".kt"), args.output + ".tsv.bgz")


    if args.slack_channel:
        send_message(args.slack_channel, 'gnomAD compound hets {} is done processing!'.format(args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Use exomes', required=False, action='store_true')
    parser.add_argument('--genomes', help='Use genomes', required=False, action='store_true')
    parser.add_argument('--write_pairs_kt', help='Writes variant pairs KeyTable.', required=False, action='store_true')
    parser.add_argument('--write_by_part', help='Writes variant pairs KeyTable or export (if using --export_pairs_no_kt) by part.', required=False, action='store_true')
    parser.add_argument('--export_pairs_no_kt', help='Writes variant pairs as a tsv without saving the KT.', required=False, action='store_true')
    parser.add_argument('--export_pairs', help='Writes variant pairs tsv.', required=False, action='store_true')
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrites existing results.', required=False,
                        action='store_true')
    parser.add_argument('--chrom20', help='Process chrom 20 only', required=False, action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--input_genes', help='File containing a list of genes to extract.', required=False)
    parser.add_argument('--pops', help='Comma-separated list of populations (not blank spaces allowed) to extract (if not set => all populations).', required=False)
    args = parser.parse_args()

    if args.exomes and args.genomes:
        sys.exit("Only one of --exomes, --genomes can be specified.")

    if args.export_pairs_no_kt and args.export_pairs:
        sys.exit("Only one of --export_pairs and --export_pairs_no_kt can be specified.")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
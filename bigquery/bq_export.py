from gnomad_hail import *
from bq_utils import *


def export_genotypes(data_type: str, export_missing_genotypes: bool, output_dir: str, max_freq: Optional[float] = None, least_consequence: str = None, variant_index: bool = True) -> None:
    mt = get_gnomad_data(data_type, non_refs_only=True)
    mt = mt.select_cols().select_rows().add_row_index()

    vep = hl.read_table(annotations_ht_path(data_type, 'vep'))
    if least_consequence is not None:
        vep_consequences = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1]))
        vep = vep.filter(
            vep.vep.transcript_consequences.any(
                lambda x: (x.biotype == 'protein_coding') & hl.any(lambda csq: vep_consequences.contains(csq), x.consequence_terms)
            )
        )
    else:
        vep = vep.filter(
            vep.vep.transcript_consequences.any(lambda x: x.biotype == 'protein_coding')
        )
    vep = vep.persist()
    logger.info(f"Found {vep.count()} variants with a VEP coding transcript and a consequence worst or equal to {least_consequence}.")

    select_expr = hl.is_defined(vep[mt.row_key])

    if max_freq is not None:
        freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
        freq = freq.filter(freq.freq[0].AF < max_freq)
        freq = freq.persist()
        logger.info(f"Found {freq.count()} variants with AF < {max_freq}")
        select_expr = select_expr & hl.is_defined(vep[mt.row_key])

    mt = mt.filter_rows(select_expr)
    ht = mt.entries()
    if export_missing_genotypes:
        ht = ht.filter(ht.is_missing)
    else:
        ht = ht.filter(hl.is_defined(ht.GT))

    ht = ht.key_by()
    if variant_index:
        select_expr = {'v': ht.row_idx}
    else:
        select_expr = {
            'chrom': ht.locus.contig,
            'pos': ht.locus.position,
            'ref': ht.alleles[0],
            'alt': ht.alleles[1]
        }
    select_expr.update({
        's': ht.s,
        'is_het': ht.GT.is_het(),
        'is_adj': ht.adj,
        'dp': ht.DP,
        'gq': ht.GQ,
        'ad0': ht.AD[0],
        'ad1': ht.AD[1],
        'pl0': ht.PL[0],
        'pl1': ht.PL[1],
        'pl2': ht.PL[2],
        'pid': ht.PID,
        'pgt': ht.PGT
    })
    ht = ht.select(**select_expr)

    ht.to_spark().write.parquet('{}/gnomad_{}{}_genotypes.parquet'.format(output_dir, data_type,'missing_' if export_missing_genotypes else '' ), mode='overwrite')


def export_variants(data_type: str, export_subsets_freq: bool, export_all_sex_freq: bool, output_dir: str, add_vep: bool) -> None:

    def format_freq(ht: hl.Table, subset: str):
        def filter_freqs(x: hl.expr.ArrayExpression):
            if export_all_sex_freq:
                return ~x[0].contains('platform') & ~x[0].contains('downsampling')
            else:
                return ~x[0].contains('platform') & ~x[0].contains('downsampling') & ~(x[0].contains('sex') & x[0].contains('pop'))

        def get_freq_struct(x: hl.expr.ArrayExpression):
            freq_struct = hl.struct(
                adj=x[0]['group'] == 'adj',
                pop=x[0].get('subpop', x[0].get('pop', 'all')),
                sex=x[0].get('sex', 'all'),
                ac=x[1].AC,
                an=x[1].AN,
                af=x[1].AF,
                hom=x[1].homozygote_count
            )

            if export_subsets_freq:
                freq_struct = freq_struct.annotate(
                    subset=subset
                )
            return freq_struct

        return  ht.annotate(freq=hl.zip(freq_meta, freq[ht.key].freq)
                     .filter(filter_freqs)
                     .map(get_freq_struct))

    ht = get_gnomad_data(data_type, non_refs_only=True).rows()
    ht = ht.select().add_index()  # TODO: Add interesting annotations
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    freq_meta = hl.literal(freq.freq_meta.collect()[0])
    ht = format_freq(ht, 'all')

    if export_subsets_freq:
        subsets = {'control': 'controls',
                   'neuro': 'non_neuro',
                   'topmed': 'non_topmed'}
        if data_type == 'exomes':
            subsets.update({'tcga': 'non_cancer'})
        for file_name, subset_name in subsets.items():
            subset_ht = format_freq(hl.read_table(annotations_ht_path(data_type, f'frequencies_{file_name}')), subset_name)
            ht = ht.annotate(freq=ht.freq.extend(subset_ht[ht.key].freq))

    if add_vep:
        vep = hl.read_table(annotations_ht_path(data_type, 'vep'))
        vep = vep.select(
            transcript_consequences=vep.vep.transcript_consequences,
            most_severe_consequence=vep.vep.most_severe_consequence
        )
        ht = ht.annotate(**vep[ht.key])

    rf = hl.read_table(annotations_ht_path(data_type, 'rf'))
    rf_expr = {f: rf[ht.key][f] for f in rf.row_value if not f.endswith('rank')}
    rf_expr['rf_probability'] = rf[ht.key]['rf_probability']['TP']
    ht = ht.annotate(**rf_expr)

    lcr_intervals = hl.import_locus_intervals(lcr_intervals_path)
    decoy_intervals = hl.import_locus_intervals(decoy_intervals_path)
    segdup_intervals = hl.import_locus_intervals(segdup_intervals_path)
    ht = ht.annotate(lcr=hl.is_defined(lcr_intervals[ht.locus]))  # TODO: combine annotations if Hail bug is fixed
    ht = ht.annotate(decoy=hl.is_defined(decoy_intervals[ht.locus]))
    ht = ht.annotate(segdup=hl.is_defined(segdup_intervals[ht.locus]))
    ht = ht.annotate(nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar()))

    export_ht_for_bq(ht, f'{output_dir}/gnomad_{data_type}_variants.parquet')


def main(args):

    hl.init(log='/bq.log')
    data_types = []
    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:

        if args.export_metadata:
            meta = get_gnomad_meta(data_type=data_type, full_meta=True)
            export_ht_for_bq(meta, f'{args.output_dir}/gnomad_{data_type}_meta.parquet')

        if args.export_genotypes:
            export_genotypes(data_type,
                             False,
                             args.output_dir,
                             args.max_freq,
                             args.least_consequence,
                             True)

        if args.export_missing_genotypes:
            export_genotypes(data_type,
                             True,
                             args.output_dir,
                             args.max_freq,
                             args.least_consequence,
                             True)

        if args.export_variants:
            export_variants(data_type, args.export_subset_freq, args.export_all_sex_freq, args.output_dir, True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--exomes', help='Run on exomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--export_metadata', help='Export samples metadata', action='store_true')
    parser.add_argument('--export_genotypes', help='Export non-ref genotypes', action='store_true')
    parser.add_argument('--export_variants', help='Export variants', action='store_true')
    parser.add_argument('--output_dir', help='Output root. default: gs://gnomad-tmp/bq', default='gs://gnomad-tmp/bq')

    var_exp = parser.add_argument_group('Export variants', description='Options related to exporting gnomAD variants')
    var_exp.add_argument('--export_subset_freq', help='If set, exports subset frequencies', action='store_true')
    var_exp.add_argument('--export_all_sex_freq', help='If set, exports sex-specific frequencies for each pop', action='store_true')

    gt_exp = parser.add_argument_group('Export genotypes', description='Options related to exporting gnomAD genotypes')
    gt_exp.add_argument('--max_freq', help='If specified, maximum global adj AF for genotypes table to emit. (default: 0.02)', default=0.02, type=float)
    gt_exp.add_argument('--least_consequence', help='When exporting genotypes, includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: 3_prime_UTR_variant)',
                        default='3_prime_UTR_variant')
    gt_exp.add_argument('--export_missing_genotypes', help='Export missing genotypes (missing genotypes export to a different (_missing) file.', action='store_true')

    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if not args.export_genotypes and not args.export_variants and not args.export_metadata and not args.export_transcripts:
        sys.exit('Error: At least one of --export_metadata, --export_variants or --export_genotypes or --export_transcripts must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)





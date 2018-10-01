from gnomad_hail import *
import hail as hl


def export_for_bq(ht: hl.Table, destination: str, lower: bool = True, mode: str = 'overwrite'):

    # Special case locus and alleles
    ht = ht.key_by()
    if 'locus' in ht.row:
        ht = ht.transmute(
            chrom=ht.locus.contig,
            pos=ht.locus.position
        )

    if 'alleles' in ht.row:
        ht = ht.transmute(
            ref=ht.alleles[0],
            alt=hl.delimit(ht.alleles[1:], ",")
        )

    # Flatten all structs
    while [c for c in ht.row.items() if isinstance(c[1].dtype, hl.tstruct)]:
        ht = ht.flatten()

    # Rename any field with '.' (and lower case if set)
    names = {name: name.replace('.', '_') for name in ht.row}
    if lower:
        names = {old_name: new_name.lower() for old_name, new_name in names.items()}

    ht = ht.rename(names)

    ht.to_spark().write.parquet(destination, mode=mode)


def export_genotypes(data_type: str, output_dir: str, max_freq: Optional[float] = None, least_consequence: str = None, variant_index: bool = True) -> None:
    mt = get_gnomad_data(data_type, non_refs_only=True)
    mt = mt.select_cols().select_rows().add_row_index().add_col_index()

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
        freq = freq.filter(freq.freq[0].AF[1] < max_freq)
        freq = freq.persist()
        logger.info(f"Found {freq.count()} variants with AF < {max_freq}")
        select_expr = select_expr & hl.is_defined(vep[mt.row_key])

    mt = mt.filter_rows(select_expr)
    ht = mt.entries()
    ht = ht.filter(ht.is_missing | hl.is_defined(ht.GT))
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
        's_id': ht.col_idx,
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

    ht.to_spark().write.parquet(f'{output_dir}/gnomad_{data_type}_genotypes.parquet', mode='overwrite')


def export_variants(data_type: str, output_dir: str) -> None:
    ht = get_gnomad_data(data_type, non_refs_only=True).rows()
    ht = ht.select().add_index()  # TODO: Add interesting annotations
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    freq_meta = hl.literal(freq.freq_meta.collect()[0])
    ht = ht.annotate(freq=hl.zip(freq_meta, freq[ht.key].freq)[:10].map(lambda x: hl.struct(
        adj=x[0]['group'] == 'adj',
        pop=x[0].get('pop', 'all'),
        ac=x[1].AC[1],
        an=x[1].AN,
        af=x[1].AF[1],
        hom=x[1].homozygote_count[1]
    )))
    export_for_bq(ht, f'{output_dir}/gnomad_{data_type}_variants.parquet')


def export_transcripts(data_type: str, output_dir: str) -> None:
    ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    ht = ht.add_index()
    ht = ht.key_by()
    ht = ht.select(
        v=ht.idx,
        transcript_consequences=ht.vep.transcript_consequences,
        most_severe_consequence=ht.vep.most_severe_consequence
    )
    ht = ht.explode('transcript_consequences')
    ht = ht.flatten()
    ht = ht.transmute(
        transcript_consequences_domains=hl.fold(lambda a, b: a + "," + b, "", ht['transcript_consequences.domains'].map(lambda x: x.db + ":" + x.name))[1:]
    )
    ht = ht.explode('transcript_consequences.consequence_terms')
    ht = ht.repartition(1000, shuffle=False)

    export_for_bq(ht, f'{args.output_dir}/gnomad_{data_type}_transcripts.parquet')


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
            gnomad = get_gnomad_data(data_type, non_refs_only=True).cols().add_index()
            meta = meta.annotate(s_id=gnomad[meta.key].idx)
            export_for_bq(meta, f'{args.output_dir}/gnomad_{data_type}_meta.parquet')

        if args.export_genotypes:
            export_genotypes(data_type,
                             args.output_dir,
                             args.max_freq,
                             args.least_consequence,
                             True)

        if args.export_variants:
            export_variants(data_type, args.output_dir)

        if args.export_transcripts:
            export_transcripts(data_type, args.output_dir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--exomes', help='Run on exomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--max_freq', help='If specified, maximum global adj AF for genotypes table to emit. (default: 0.05)', default=0.05, type=float)
    parser.add_argument('--least_consequence', help='Includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: 3_prime_UTR_variant)',
                        default='3_prime_UTR_variant')
    parser.add_argument('--export_metadata', help='Export samples metadata', action='store_true')
    parser.add_argument('--export_genotypes', help='Export genotypes', action='store_true')
    parser.add_argument('--export_variants', help='Export variants', action='store_true')
    parser.add_argument('--export_transcripts', help='Export transcripts', action='store_true')

    parser.add_argument('--output_dir', help='Output root. default: gs://gnomad-tmp/bq', default='gs://gnomad-tmp/bq')
    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if not args.export_genotypes and not args.export_variants and not args.export_metadata and not args.export_transcripts:
        sys.exit('Error: At least one of --export_metadata, --export_variants or --export_genotypes or --export_transcripts must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)





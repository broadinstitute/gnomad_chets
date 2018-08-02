from gnomad_hail import *
import hail as hl


def export_for_bq(ht: hl.Table, destination: str, lower: bool):

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
            alt=hl.delimit(ht.alleles[:1], ",")
        )

    # Flatten all structs
    while [c for c in ht.row.items() if isinstance(c[1].dtype, hl.tstruct)]:
        ht = ht.flatten()

    # Rename any field with '.' (and lower case if set)
    names = {name: name.replace('.', '_') for name in ht.row}
    if lower:
        names = {old_name: new_name.lower() for old_name, new_name in names.items}

    ht.rename(names)

    ht.to_spark().write.parquet(destination)


def export_genotypes(data_type: str, output_dir: str, max_freq: Optional[float] = None) -> None:
    mt = get_gnomad_data(data_type, non_refs_only=True)
    mt = mt.select_cols().select_rows()

    vep = hl.read_table(annotations_ht_path(data_type, 'vep'))
    vep = vep.filter(
        vep.vep.transcript_consequences.any(lambda x: x.biotype == 'protein_coding')
    )
    vep = vep.persist()
    logger.info(f"Found {vep.count()} variants with a VEP coding transcript.")

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
    ht = ht.select(
        chrom=ht.locus.contig,
        pos=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
        is_het=ht.GT.is_het(),
        is_adj=ht.adj,
        dp=ht.DP,
        gq=ht.GQ,
        ad=ht.AD,
        ad1=ht.AD[1],
        pl0=ht.PL[0],
        pl1=ht.PL[1],
        pl2=ht.PL[2],
        pid=ht.PID,
        pgt=ht.PGT
    )

    ht.to_spark().write.parquet(f'{output_dir}/gnomad_{data_type}_genotypes.parquet')


def main(args):

    hl.init(log='/create_rank.log')
    data_types = []
    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:

        if args.export_metadata:
            export_for_bq(get_gnomad_meta(data_type=data_type, full_meta=True), f'{args.output_dir}/gnomad_{data_type}_meta.parquet')

        if args.export_genotypes:
            export_genotypes(data_type, args.max_freq, args.output_dir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--exomes', help='Run on exomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--max_freq', help='If specified, maximum global adj AF for genotypes table to emit.', type=float)
    parser.add_argument('--export_metadata', help='Export samples metadata', action='store_true')
    parser.add_argument('--export_genotypes', help='Export genotypes', action='store_true')
    parser.add_argument('--export_variants', help='Export variants', action='store_true')

    parser.add_argument('--output_dir', help='Output root. default: gs://gnomad-tmp/bq', default='gs://gnomad-tmp/bq')
    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if not args.export_genotypes and not args.export_variants and not args.export_metadata:
        sys.exit('Error: At least onr of --export_metadata, --export_variants or --export_genotypes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)





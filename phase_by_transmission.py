from gnomad_hail import *
import hail as hl


def explode_trio_matrix(tm: hl.MatrixTable, col_keys: List[str] = ['s'], keep_trio_cols: bool = True, keep_trio_entries: bool = False) -> hl.MatrixTable:
    """Splits a trio MatrixTable back into a sample MatrixTable.

    Example
    -------
    >>> # Create a trio matrix from a sample matrix
    >>> pedigree = hl.Pedigree.read('data/case_control_study.fam')
    >>> trio_dataset = hl.trio_matrix(dataset, pedigree, complete_trios=True)

    >>> # Explode trio matrix back into a sample matrix
    >>> exploded_trio_dataset = explode_trio_matrix(trio_dataset)

    Notes
    -----
    The resulting MatrixTable column schema is the same as the proband/father/mother schema,
    and the resulting entry schema is the same as the proband_entry/father_entry/mother_entry schema.
    If the `keep_trio_cols` option is set, then an additional `source_trio` column is added with the trio column data.
    If the `keep_trio_entries` option is set, then an additional `source_trio_entry` column is added with the trio entry data.

    Note
    ----
    This assumes that the input MatrixTable is a trio MatrixTable (similar to the result of :meth:`.methods.trio_matrix`)
    Its entry schema has to contain 'proband_entry`, `father_entry` and `mother_entry` all with the same type.
    Its column schema has to contain 'proband`, `father` and `mother` all with the same type.

    Parameters
    ----------
    tm : :class:`.MatrixTable`
        Trio MatrixTable (entries have to be a Struct with `proband_entry`, `mother_entry` and `father_entry` present)
    col_keys : :obj:`list` of str
        Column key(s) for the resulting sample MatrixTable
    keep_trio_cols: bool
        Whether to add a `source_trio` column with the trio column data (default `True`)
    keep_trio_entries: bool
        Whether to add a `source_trio_entries` column with the trio entry data (default `False`)

    Returns
    -------
    :class:`.MatrixTable`
        Sample MatrixTable"""

    select_entries_expr = {'__trio_entries': hl.array([tm.proband_entry, tm.father_entry, tm.mother_entry])}
    if keep_trio_entries:
        select_entries_expr['source_trio_entry'] = hl.struct(**tm.entry)
    tm = tm.select_entries(**select_entries_expr)

    tm = tm.key_cols_by()
    select_cols_expr = {'__trio_members': hl.zip_with_index(hl.array([tm.proband, tm.father, tm.mother]))}
    if keep_trio_cols:
        select_cols_expr['source_trio'] = hl.struct(**tm.col)
    tm = tm.select_cols(**select_cols_expr)

    mt = tm.explode_cols(tm.__trio_members)

    mt = mt.transmute_entries(
        **mt.__trio_entries[mt.__trio_members[0]]
    )

    mt = mt.key_cols_by()
    mt = mt.transmute_cols(**mt.__trio_members[1])

    if col_keys:
        mt = mt.key_cols_by(*col_keys)

    return mt


def main(args):
    data_type = 'exomes' if args.exomes else 'genomes'

    if args.pbt_tm:
        mt = get_gnomad_data(data_type, split=False)
        meta = mt.cols()
        hq_samples = meta.aggregate(hl.agg.filter(meta.meta.high_quality, hl.agg.collect(meta.s)))
        ped = hl.Pedigree.read(fam_path(data_type), delimiter='\\t').filter_to(hq_samples)
        ped_samples = hl.literal(set([s for trio in ped.complete_trios() for s in [trio.s, trio.pat_id, trio.mat_id]]))

        mt = mt.filter_cols(ped_samples.contains(mt.s))
        mt = mt.select_cols().select_rows()
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

        tm = hl.trio_matrix(mt, ped, complete_trios=True)
        tm = hl.experimental.phase_trio_matrix_by_transmission(tm)
        tm.write(pbt_phased_trios_mt_path(data_type, split=False, trio_matrix=True), overwrite=args.overwrite)

    if args.pbt_explode:
        tm = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type, split=False, trio_matrix=True))

        tm = tm.annotate_entries(trio_adj=tm.proband_entry.adj & tm.father_entry.adj & tm.mother_entry.adj)
        pmt = explode_trio_matrix(tm, keep_trio_entries=True)
        pmt = pmt.transmute_entries(trio_adj=pmt.source_trio_entry.trio_adj)
        pmt.write(pbt_phased_trios_mt_path(data_type, split=False), overwrite=args.overwrite)

        pmt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type, split=False))
        pmt = pmt.rename({'PBT_GT': 'PGT'}) # ugly but supported by hl.split_multi_hts
        pmt = hl.split_multi_hts(pmt)
        pmt = pmt.rename({'PGT': 'PBT_GT'})
        pmt.write(pbt_phased_trios_mt_path(data_type), overwrite=args.overwrite)

    if args.phase_multi_families:
        pbt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
        # Keep samples that:
        # 1. There are more than one entry in the Matrix (i.e. they are part of multiple trios)
        # 2. In all their entries, the parents are the same (there are only two exceptions to this, so best to ignore these and focus on parents/multi-offspring families)
        nt_samples = pbt.cols()
        nt_samples = nt_samples.group_by('s').aggregate(trios=hl.agg.collect(nt_samples.source_trio))
        nt_samples = nt_samples.filter((hl.len(nt_samples.trios) > 1) &
                                       nt_samples.trios[1:].any(
                                           lambda x: (x.mother.s != nt_samples.trios[0].mother.s) | (x.father.s != nt_samples.trios[0].father.s)),
                                       keep=False)
        pbt = pbt.filter_cols(hl.is_defined(nt_samples[pbt.col_key]))

        # Group cols for these samples, keeping all GTs in an array
        # Compute the consensus GT (incl. phase) + QC metrics based on (a) phased genotypes have priority, (b) genotypes with most votes
        pbt = pbt.group_cols_by('s').aggregate(
            PBT_GTs=hl.agg.filter(hl.is_defined(pbt.GT), hl.agg.collect(pbt.GT))
        )
        gt_counter = hl.sorted(hl.array(pbt.PBT_GTs.group_by(lambda x: x).map_values(lambda x: hl.len(x))), key=lambda x: x[0].phased * 100 + x[1], reverse=True)
        phased_gt_counts = gt_counter.filter(lambda x: x[0].phased).map(lambda x: x[1])
        pbt = pbt.annotate_entries(
            consensus_gt=gt_counter.map(lambda x: x[0]).find(lambda x: True),
            phase_concordance=phased_gt_counts.find(lambda x: True) / hl.sum(phased_gt_counts),
            discordant_gts=hl.len(hl.set(pbt.PBT_GTs.map(lambda x: hl.cond(x.phased, hl.call(x[0],x[1]), x)))) > 1
        )
        pbt.write('gs://gnomad/projects/compound_hets/pbt_multi_families.mt')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--exomes', help='Run on exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--pbt_tm', help='Creates a PBT-phased trio-matrix.',
                        action='store_true')
    parser.add_argument('--pbt_explode', help='Creates a PBT-phased MT by exploding the pbt_mt.',
                        action='store_true')
    parser.add_argument('--phase_multi_families', help='Computes consensus phase from PBT in families with multiple offspring.',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)




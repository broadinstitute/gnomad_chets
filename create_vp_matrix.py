from gnomad_hail import *
import hail as hl
from resources import *

def create_variant_pair_mt(mt: hl.MatrixTable, row_groups: List[str]):
    """
    Create a variant-pair MatrixTable containing all variant-pairs that appear both in the same individual within the row group.
    E.g., if the row group is a gene, this creates a variant-pair MT with all variant-pairs in a given gene.

    Note that this produces an empty shell of an MT with no data in the rows, columns or entries.

    :param MatrixTable mt: Input MatrixTable
    :param list of str row_groups: Row annotations for delimiting variant-pairs territory
    :return: Variant-pair MT
    :rtype: MatrixTable
    """

    mt = mt.filter_entries(mt.GT.is_non_ref())
    et = mt.select_cols().select_rows(*row_groups).entries()
    et = et.group_by(*row_groups, *mt.col_key).aggregate(
        vgt=(hl.agg.collect(hl.struct(locus=et.locus, alleles=et.alleles)))
    )

    et = et.annotate(
        vgt=hl.range(0, hl.len(et.vgt))
                      .flatmap(lambda i1: hl.range(i1+1, hl.len(et.vgt))
                               .map(lambda i2: hl.struct(v1=et.vgt[i1], v2=et.vgt[i2])))
    )

    et = et.explode(et.vgt)
    et = et.key_by(*mt.col_key, locus1=et.vgt.v1.locus, alleles1=et.vgt.v1.alleles, locus2=et.vgt.v2.locus, alleles2=et.vgt.v2.alleles)
    et = et.select()
    et = et.distinct()

    vp_mt = et.to_matrix_table(row_key=['locus2','alleles2','locus1','alleles1'],
                               col_key=[*mt.col_key])
    return vp_mt


def filter_freq_and_csq(mt: hl.MatrixTable, data_type: str, max_freq: float, least_consequence: str):
    """
    Filters MatrixTable to include variants that:
    1. Have a global AF <= `max_freq`
    2. Have a consequence at least as severe as `least_consequence` (based on ordering from CSQ_ORDER)

    :param MatrixTable mt: Input MT
    :param str data_type: One of 'exomes' or 'genomes'
    :param float max_freq: Max. AF to keep
    :param str least_consequence: Least consequence to keep.
    :return: Filtered MT
    :rtype: MatrixTable
    """

    # mt = mt.select_cols(pop=mt.meta.pop)
    mt = mt.select_rows('a_index')

    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    vep_consequences = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1]))

    mt = mt.select_rows(
        vep=(
            hl.set(
                vep_ht[mt.row_key].vep.transcript_consequences
                    .filter(
                    lambda tc: (tc.biotype == 'protein_coding') &
                               (tc.allele_num == mt.a_index) &
                               (tc.consequence_terms.any(lambda c: vep_consequences.contains(c)))
                )
                    .map(lambda x: x.gene_id)
            )
        ),
        af=hl.float32(freq[mt.row_key].freq[0].AF)
    )

    mt = mt.filter_rows(hl.is_defined(mt.vep) & (hl.len(mt.vep) > 0) & (mt.af > 0) & (mt.af <= max_freq))
    mt = mt.explode_rows(mt.vep)
    mt = mt.rename({'vep': 'gene_id'})
    return mt


def get_counts_agg_expr(mt: hl.MatrixTable):
    return (hl.case(missing_false=True)
    #0x
    .when(hl.is_missing(mt.gt1) & ~mt.missing1,
          hl.case(missing_false=True)
          .when(hl.is_missing(mt.gt2) & ~mt.missing2,[1,0,0,0,0,0,0,0,0])
          .when(mt.gt2.is_het(), [0,1,0,0,0,0,0,0,0])
          .when(mt.gt2.is_hom_var(), [0,0,1,0,0,0,0,0,0])
          .default([0,0,0,0,0,0,0,0,0]))
    #1x
    .when(mt.gt1.is_het(),
          hl.case(missing_false=True)
          .when(hl.is_missing(mt.gt2) & ~mt.missing2,[0,0,0,1,0,0,0,0,0])
          .when(mt.gt2.is_het(), [0,0,0,0,1,0,0,0,0])
          .when(mt.gt2.is_hom_var(), [0,0,0,0,0,1,0,0,0])
          .default([0,0,0,0,0,0,0,0,0]))
    #2x
    .when(mt.gt1.is_hom_var(),
          hl.case(missing_false=True)
          .when(hl.is_missing(mt.gt2) & ~mt.missing2,[0,0,0,0,0,0,1,0,0])
          .when(mt.gt2.is_het(), [0,0,0,0,0,0,0,1,0])
          .when(mt.gt2.is_hom_var(), [0,0,0,0,0,0,0,0,1])
          .default([0,0,0,0,0,0,0,0,0]))
    .default([0,0,0,0,0,0,0,0,0]))


def annotate_vp_ht(ht, data_type):

    def get_freq_expr(ht,  freq_ht, freq_dict, by_pop):
        pop_index = freq_dict[ht.pop] if by_pop else 0
        freq_expr = freq_ht[ht.key].freq[pop_index]
        return hl.struct(ac=freq_expr.AC, an=freq_expr.AN, af=freq_expr.AF)

    ht = ht.persist()
    # Annotate freq and CpG information
    methyation_ht = hl.read_table(methylation_sites_mt_path())
    freq_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    rf_ht = hl.read_table(annotations_ht_path(data_type, 'rf'))
    freq_meta = hl.eval(freq_ht.globals.freq_meta)
    freq_dict = {f['pop']: i for i, f in enumerate(freq_meta[:10]) if 'pop' in f}
    freq_dict['all'] = 0
    freq_dict = hl.literal(freq_dict)
    ht_ann = ht.key_by('locus1', 'alleles1').select('locus2', 'alleles2', 'pop')
    ht_ann = ht_ann.annotate(
        filters1=rf_ht[ht_ann.key].filters,
        global_freq1=get_freq_expr(ht_ann, freq_ht, freq_dict, False),
        pop_freq1=hl.or_missing(hl.is_defined(ht_ann.pop), get_freq_expr(ht_ann, freq_ht, freq_dict, True)),
        cpg1=methyation_ht[ht_ann.locus1].MEAN > 0.6,
        snv1=hl.is_snp(ht_ann.alleles1[0], ht_ann.alleles1[1])
    )
    ht_ann.write(f'gs://gnomad-tmp/compound_hets/{data_type}_ann1.ht', overwrite=True)
    ht_ann = hl.read_table(f'gs://gnomad-tmp/compound_hets/{data_type}_ann1.ht')
    ht_ann = ht_ann.key_by('locus2', 'alleles2')
    ht_ann = ht_ann.annotate(
        filters2=rf_ht[ht_ann.key].filters,
        global_freq2=get_freq_expr(ht_ann, freq_ht, freq_dict, False),
        pop_freq2=hl.or_missing(hl.is_defined(ht_ann.pop), get_freq_expr(ht_ann, freq_ht, freq_dict, True)),
        cpg2=methyation_ht[ht_ann.locus2].MEAN > 0.6,
        snv2=hl.is_snp(ht_ann.alleles2[0], ht_ann.alleles2[1]))
    ht_ann.write(f'gs://gnomad-tmp/compound_hets/{data_type}_ann2.ht', overwrite=True)
    ht_ann = hl.read_table(f'gs://gnomad-tmp/compound_hets/{data_type}_ann2.ht').key_by('locus1', 'alleles1', 'locus2', 'alleles2', 'pop')
    ht = ht.key_by('locus1', 'alleles1', 'locus2', 'alleles2', 'pop')
    return ht.annotate(**ht_ann[ht.key])


def remove_pbt_dups(pbt_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    PBT MT contains a column per trio-individual.
    For samples involved in multiple trios (in this case, a parent with multiple offspring),
    this filters samples to only keep each sample ones, picking randomly amongst possible duplicates
    """
    pbt_cols = pbt_mt.cols()
    pbt_dups = pbt_cols.group_by('s').aggregate(trio_ids=hl.agg.collect(pbt_cols.trio_id))
    pbt_dups = pbt_dups.filter(hl.len(pbt_dups.trio_ids) > 1)
    pbt_dups = pbt_dups.persist()
    pbt_samples_to_remove = pbt_dups.transmute(trio_id=pbt_dups.trio_ids[1:]).explode('trio_id').key_by('s', 'trio_id')
    return pbt_mt.filter_cols(hl.is_missing(pbt_samples_to_remove[pbt_mt.col_key]))


def create_full_vp(data_type, path_args, args):
    vp_mt = hl.read_matrix_table(vp_list_mt_path(*path_args))
    if args.pbt:
        gnomad = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
        gnomad = gnomad.key_cols_by(s=gnomad.s, trio_id=gnomad.source_trio.id)
        gnomad = gnomad.select_entries(gt=hl.case()
                                       .when(gnomad.PBT_GT.is_non_ref(), gnomad.PBT_GT)
                                       .when(gnomad.GT.is_non_ref(), hl.call(gnomad.GT[0], gnomad.GT[1]))
                                       .or_missing(),
                                       missing=hl.is_missing(gnomad.GT),
                                       adj=gnomad.adj,
                                       trio_adj=gnomad.trio_adj).select_cols().select_rows()
    else:
        gnomad = hl.read_matrix_table(gnomad_adj_missing_path(data_type))

    gnomad.describe()
    vp_mt = vp_mt.key_rows_by('locus2', 'alleles2')
    gnomad_joined = gnomad[vp_mt.row_key, vp_mt.col_key]
    if args.pbt:
        vp_mt = vp_mt.annotate_entries(gt2=gnomad_joined.gt, missing2=gnomad_joined.missing, adj2=gnomad_joined.adj, trio_adj2=gnomad_joined.trio_adj)
    else:
        vp_mt = vp_mt.annotate_entries(gt2=gnomad_joined.gt, missing2=gnomad_joined.missing, adj2=gnomad_joined.adj)
    vp_mt.write(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp.mt', overwrite=True)
    vp_mt = hl.read_matrix_table(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp.mt')
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
    gnomad_joined = gnomad[vp_mt.row_key, vp_mt.col_key]
    if args.pbt:
        vp_mt = vp_mt.annotate_entries(gt1=gnomad_joined.gt, missing1=gnomad_joined.missing, adj1=gnomad_joined.adj, trio_adj1=gnomad_joined.trio_adj)
    else:
        vp_mt = vp_mt.annotate_entries(gt1=gnomad_joined.gt, missing1=gnomad_joined.missing, adj1=gnomad_joined.adj)
    vp_mt.write(full_mt_path(*path_args), overwrite=args.overwrite)


def create_vp_summary(data_type, path_args, args):
    mt = hl.read_matrix_table(full_mt_path(data_type, False, args.least_consequence, args.max_freq, args.chrom))
    meta = get_gnomad_meta(data_type).select('pop', 'release')
    mt = mt.annotate_cols(**meta[mt.col_key])
    mt = mt.filter_cols(mt.release)

    if args.pbt:
        pbt_samples = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom)).cols().key_by('s')
        mt = mt.filter_cols(hl.is_missing(pbt_samples[mt.col_key]))

    mt = mt.annotate_cols(pop=[mt.pop, 'all']).explode_cols('pop')
    ht = mt.group_cols_by(mt.pop).aggregate(
        gt_counts=hl.agg.group_by(mt.adj1 & mt.adj2,
                                  hl.agg.counter(get_counts_agg_expr(mt))).map_values(
            lambda x: hl.fold(lambda i, j: i + j[0] * j[1], [0, 0, 0, 0, 0, 0, 0, 0, 0], hl.zip(x.keys(), x.values().map(lambda v: hl.int32(v))))
        )
    ).entries()
    ht = ht.annotate(gt_counts=hl.fold(lambda i, j: hl.struct(raw=i[0] + j[1], adj=hl.cond(j[0], i[1] + j[1], i[1])),
                                       hl.struct(raw=[0, 0, 0, 0, 0, 0, 0, 0, 0], adj=[0, 0, 0, 0, 0, 0, 0, 0, 0]),
                                       hl.zip(ht.gt_counts.keys(), ht.gt_counts.values())
                                       )
                     )
    ht.write(f'gs://gnomad-tmp/compound_hets/ht_sites{"_pbt" if args.pbt else ""}_by_pop.ht', overwrite=True)
    ht = hl.read_table(f'gs://gnomad-tmp/compound_hets/ht_sites{"_pbt" if args.pbt else ""}_by_pop.ht')
    ht = ht.key_by('locus1', 'alleles1', 'locus2', 'alleles2', 'pop')
    ht = ht.repartition(1000, shuffle=False)
    ht.write(vp_count_ht_path(*path_args), overwrite=args.overwrite)


def create_pbt_summary(data_type, path_args, args):
    pbt = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))
    pbt = remove_pbt_dups(pbt)

    # Remove probands
    pbt_meta = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type)).cols()
    pbt_probands = pbt_meta.filter(pbt_meta.s == pbt_meta.source_trio.proband.s)
    pbt = pbt.key_cols_by('s')
    pbt = pbt.filter_cols(hl.is_missing(pbt_probands[pbt.col_key])).persist()
    logger.info(f"Found {pbt.count_cols()} unique parents in PBT file.")

    # Compute counts, grouped by pop
    meta = get_gnomad_meta('exomes').select('pop')
    pbt = pbt.annotate_cols(**meta[pbt.col_key])

    # Assumes only non-ref genotypes are present
    def get_hap_counter(adj, same):
        if same:
            cond = hl.range(0, pbt.gt1.ploidy).any(lambda a: pbt.gt1[a] == pbt.gt2[a])
        else:
            cond = hl.range(0, pbt.gt1.ploidy).any(lambda a: pbt.gt1[a] != pbt.gt2[a])

        if adj:
            cond = cond & pbt.trio_adj1 & pbt.trio_adj2

        return hl.agg.count_where(cond)

    pbt = pbt.filter_entries(pbt.gt1.phased & pbt.gt2.phased & (pbt.gt1.ploidy == pbt.gt2.ploidy))
    pbt = pbt.annotate_cols(pop=[pbt.pop, 'all']).explode_cols('pop')
    pbt = pbt.group_cols_by(pbt.pop).aggregate(
        adj=hl.struct(
            same_hap=get_hap_counter(True, True),
            diff_hap=get_hap_counter(True, False)
        ),
        raw=hl.struct(
            same_hap=get_hap_counter(False, True),
            diff_hap=get_hap_counter(False, False)
        )
    )

    trio_same_hap = hl.struct(
        raw=hl.case()
            .when((pbt.raw.same_hap > 0) & (pbt.raw.diff_hap == 0), True)
            .when((pbt.raw.same_hap == 0) & (pbt.raw.diff_hap > 0), True)
            .or_missing(),
        adj=hl.case()
            .when((pbt.adj.same_hap > 0) & (pbt.adj.diff_hap == 0), True)
            .when((pbt.adj.same_hap == 0) & (pbt.adj.diff_hap > 0), False)
            .or_missing()
    )

    pbt = pbt.filter_entries(pbt.raw.same_hap + pbt.raw.diff_hap > 0).entries()
    pbt.write(f'gs://gnomad-tmp/compound_hets/{data_type}_pbt_counts.ht', overwrite=args.overwrite)

    # Add freq and CpG annotations
    ht = hl.read_table(f'gs://gnomad-tmp/compound_hets/{data_type}_pbt_counts.ht')
    ht = annotate_vp_ht(ht, data_type)
    ht = ht.repartition(1000, shuffle=False)
    ht.write(pbt_phase_count_ht_path(*path_args), overwrite=args.overwrite)


def create_pbt_trio_ht(data_type, args):
    pbt = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))
    meta = get_gnomad_meta(data_type)[pbt.s]
    pbt = pbt.annotate_cols(
        sex=meta.sex,
        pop=meta.pop)
    pbt = pbt.annotate_cols(
        role=hl.case()
            .when(pbt.s == pbt.trio_id, 0)
            .when(pbt.sex == 'male', 1)
            .when(pbt.sex == 'female', 2)
            .or_missing()
    )

    # create trio matrix
    # ped = hl.Pedigree.read(fam_path(data_type), delimiter="\t")
    entries = hl.sorted(
        hl.agg.collect(
            hl.struct(
                role=pbt.role,
                entry=pbt.entry.annotate(
                    sex=pbt.sex,
                    pop=pbt.pop,
                    same_hap=hl.or_missing(
                        pbt.gt1.ploidy == pbt.gt2.ploidy,
                        hl.range(0, pbt.gt1.ploidy).any(lambda a: pbt.gt1[a] == pbt.gt2[a])
                    )
                )
            )
        ),
        key=lambda x: x.role).map(lambda x: x.entry)
    tm = pbt.group_cols_by(pbt.trio_id).aggregate(
        child=entries[0],
        father=entries[1],
        mother=entries[2]
    )
    tm.write("gs://gnomad-tmp/pbt_vp_tm.tmp.mt", overwrite=True)

    # annotate variant-phase phase per trio / overall
    tm = hl.read_matrix_table("gs://gnomad-tmp/pbt_vp_tm.tmp.mt")

    # Assumes no hom ref genotypes
    tm = tm.annotate_entries(
        same_hap=(
            hl.case(missing_false=True)
                .when(tm.father.same_hap == tm.mother.same_hap, tm.father.same_hap)
                .when(hl.is_defined(tm.father.same_hap) & hl.is_missing(tm.mother.same_hap), tm.father.same_hap)
                .when(hl.is_missing(tm.father.same_hap) & hl.is_defined(tm.mother.same_hap), tm.mother.same_hap)
                .or_missing()
        ),
        adj1=tm.child.trio_adj1,
        adj2=tm.child.trio_adj2,
        pop=(
            hl.case(missing_false=True)
                .when(tm.father.pop == tm.mother.pop, tm.father.pop)
                .when((tm.father.gt1.is_het() | tm.father.gt2.is_het()) & ~(tm.mother.gt1.is_het() | tm.mother.gt2.is_het()), tm.father.pop)
                .when(~(tm.father.gt1.is_het() | tm.father.gt2.is_het()) & (tm.mother.gt1.is_het() | tm.mother.gt2.is_het()), tm.mother.pop)
                .or_missing()
        ),
        child=tm.child.drop('trio_adj1', 'trio_adj2'),
        father=tm.father.drop('trio_adj1', 'trio_adj2'),
        mother=tm.mother.drop('trio_adj1', 'trio_adj2')
    )

    tm.write(pbt_trio_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom), overwrite=args.overwrite)

    # Create entries table
    tm = hl.read_matrix_table(pbt_trio_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))
    et = tm.entries()
    et = et.filter(hl.is_defined(et.same_hap))
    et = et.flatten()

    # Add annotations
    et = annotate_vp_ht(et, data_type)

    et.write(pbt_trio_et_path(data_type, True, args.least_consequence, args.max_freq, args.chrom), overwrite=args.overwrite)


def main(args):

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]

    if args.create_mini_mt:
        if args.pbt:
            mt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
            mt = mt.filter_entries(mt.GT.is_non_ref())
            mt = mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)
            # mt = mt.annotate_entries(GT=mt.PBT_GT)
        else:
            mt = get_gnomad_data(data_type, non_refs_only=True)
            mt = mt.filter_cols(mt.meta.high_quality)
            mt = mt.select_entries('GT', 'adj', 'is_missing')

        if args.chrom:
            print(f"Selecting chrom {args.chrom}")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.chrom)])

        mt = filter_freq_and_csq(mt, data_type, args.max_freq, args.least_consequence)
        mt.write(mini_mt_path(*path_args), overwrite=args.overwrite)

    if args.create_vp_list:
        mt = hl.read_matrix_table(mini_mt_path(*path_args))
        mt = create_variant_pair_mt(mt, ['gene_id'])
        mt.write(vp_list_mt_path(*path_args), overwrite=args.overwrite)

    if args.create_gnomad_adj_missing:
        gnomad = get_gnomad_data(data_type).select_cols().select_rows() if not args.pbt else hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
        gnomad = gnomad.select_entries(gt=hl.or_missing(gnomad.GT.is_non_ref(), gnomad.GT), missing=hl.is_missing(gnomad.GT), adj=gnomad.adj).select_cols().select_rows()
        gnomad.write(gnomad_adj_missing_path(data_type), overwrite=args.overwrite)

    if args.create_full_vp:
        create_full_vp(data_type, path_args, args)

    if args.create_vp_summary:
        create_vp_summary(data_type, path_args, args)

    if args.create_pbt_summary:
        create_pbt_summary(data_type, path_args, args)

    if args.create_pbt_trio_ht:
        create_pbt_trio_ht(data_type, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--pbt', help='Runs on PBT-phased data instead of the entire gnomAD. Note that the PBT_GT will be renamed as GT',
                        action='store_true')
    parser.add_argument('--create_mini_mt', help='Creates a filtered, minimal MT that is then used to create the VP MT.',
                        action='store_true')
    parser.add_argument('--create_vp_list', help='Creates an empty VP MT containing all variant pairs but no other data.', action='store_true')
    parser.add_argument('--create_gnomad_adj_missing', help='Creates a gnomAD MT with only missing and adj fields.', action='store_true')
    parser.add_argument('--create_full_vp', help='Creates the VP MT.', action='store_true')
    parser.add_argument('--create_vp_summary', help='Creates a summarised VP table, with counts in release samples only. If --pbt is specified, then only sites present in PBT samples are used and counts exclude PBT samples.',
                        action='store_true')
    parser.add_argument('--create_pbt_summary', help='Creates a summarised PBT table, with counts of same/diff hap in unique parents. Note that --pbt flag has no effect on this.',
                        action='store_true')
    parser.add_argument('--create_pbt_trio_ht', help='Creates a HT with one line per trio/variant-pair (where trio is non-ref). Note that --pbt flag has no effect on this.',
                        action='store_true')
    parser.add_argument('--least_consequence', help=f'Includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--max_freq', help=f'If specified, maximum global adj AF for genotypes table to emit. (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chrom', help='Only run on given chromosome')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)




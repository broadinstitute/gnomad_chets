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
    et = et.select(dummy_entry=hl.null(hl.tbool))
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

def create_full_vp(data_type, path_args, args):
    vp_mt = hl.read_matrix_table(vp_list_mt_path(*path_args))
    if args.pbt:
        mt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
        mt = mt.key_cols_by('s', trio_id=mt.source_trio.id)
        mt = extract_pbt_probands(mt, data_type)
        mt = mt.key_cols_by('s')
        mt = mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)
        mt = mt.select_entries(gt=hl.case()
                                       .when(mt.PBT_GT.is_non_ref(), mt.PBT_GT)
                                       .when(mt.GT.is_non_ref(), hl.call(mt.GT[0], mt.GT[1]))
                                       .or_missing(),
                                       missing=hl.is_missing(mt.GT),
                                       adj=mt.adj,
                                       trio_adj=mt.trio_adj).select_cols().select_rows()
    else:
        mt = hl.read_matrix_table(gnomad_adj_missing_path(data_type))

    vp_mt = vp_mt.key_rows_by('locus2', 'alleles2')
    vp_mt = vp_mt.unfilter_entries()
    mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
    if args.pbt:
        vp_mt = vp_mt.select_entries(gt2=mt_joined.gt, missing2=mt_joined.missing, adj2=mt_joined.adj, trio_adj2=mt_joined.trio_adj)
    else:
        vp_mt = vp_mt.select_entries(gt2=mt_joined.gt, missing2=mt_joined.missing, adj2=mt_joined.adj)
    vp_mt.write(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp.mt', overwrite=True)
    vp_mt = hl.read_matrix_table(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp.mt')
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
    mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
    if args.pbt:
        vp_mt = vp_mt.annotate_entries(gt1=mt_joined.gt, missing1=mt_joined.missing, adj1=mt_joined.adj, trio_adj1=mt_joined.trio_adj)
    else:
        vp_mt = vp_mt.annotate_entries(gt1=mt_joined.gt, missing1=mt_joined.missing, adj1=mt_joined.adj)
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


def extract_pbt_probands(pbt_mt: hl.MatrixTable, data_type: str):

    # Keep a single proband from each family with > 1  proband.
    meta = get_gnomad_meta(data_type)
    hq_samples = hl.literal(meta.aggregate(hl.agg.filter(meta.high_quality & (meta.project_id != 'C978'), hl.agg.collect(meta.s))))
    fam_ht = hl.import_fam(fam_path(data_type), delimiter='\\t')
    fam_ht = fam_ht.filter(
        hq_samples.contains(fam_ht.id) &
        hq_samples.contains(fam_ht.pat_id) &
        hq_samples.contains(fam_ht.mat_id)
    )
    fam_ht = fam_ht.key_by('pat_id').distinct()
    fam_ht = fam_ht.key_by('mat_id').distinct()
    fam_ht = fam_ht.annotate(s=[fam_ht.id, fam_ht.pat_id, fam_ht.mat_id]).explode('s')
    fam_ht = fam_ht.key_by('s', 'id')

    pbt_mt = pbt_mt.filter_cols(hl.is_defined(fam_ht[pbt_mt.col_key]) & (pbt_mt.s == pbt_mt.trio_id)).key_cols_by('s').persist()
    logger.info(f"Found {pbt_mt.count_cols()} probands.")
    return pbt_mt


def create_pbt_summary(data_type, path_args, args):

    pbt = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))

    # Compute counts, grouped by pop
    meta = get_gnomad_meta('exomes').select('pop', 'project_id')
    pbt = pbt.key_cols_by('s')
    pbt = pbt.annotate_cols(**meta[pbt.col_key])
    pbt = pbt.filter_cols(pbt.project_id != 'C978') # TODO remove if re-generated entirely

    def get_trio_phase_expr(pbt) -> hl.expr.StructExpression:
        return hl.bind(
            lambda haps: hl.struct(
                adj=hl.struct(
                    same_hap_samples=haps.filter(lambda x: x[0] & x[1]).map(lambda x: x[2]),
                    chet_samples=haps.filter(lambda x: ~x[0] & x[1]).map(lambda x: x[2]),
                    n_same_hap=hl.len(haps.filter(lambda x: x[0] & x[1])),
                    n_chet=hl.len(haps.filter(lambda x: ~x[0] & x[1]))
                ),
                raw=hl.struct(
                    same_hap_samples=haps.filter(lambda x: x[0]).map(lambda x: x[2]),
                    chet_samples=haps.filter(lambda x: ~x[0]).map(lambda x: x[2]),
                    n_same_hap=hl.len(haps.filter(lambda x: x[0])),
                    n_chet=hl.len(haps.filter(lambda x: ~x[0]))
                )
            ),
            hl.agg.collect(hl.tuple([hl.range(0, pbt.gt1.ploidy).any(lambda a: pbt.gt1[a] == pbt.gt2[a]), pbt.trio_adj1 & pbt.trio_adj2, pbt.s]))
        )

    pbt = pbt.filter_entries(pbt.gt1.phased & pbt.gt2.phased & (pbt.gt1.ploidy == pbt.gt2.ploidy) & pbt.gt1.is_het() & pbt.gt2.is_het())
    pbt = pbt.annotate_cols(pop=[pbt.pop, 'all']).explode_cols('pop')
    pbt = pbt.group_cols_by(pbt.pop).aggregate(**get_trio_phase_expr(pbt))

    pbt = pbt.filter_entries(pbt.raw.n_same_hap + pbt.raw.n_chet > 0).entries()
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
                entry=hl.struct(
                    s=pbt.s,
                    gt1=pbt.entry.gt1,
                    missing1=pbt.entry.missing1,
                    adj1=pbt.entry.adj1,
                    gt2=pbt.entry.gt1,
                    missing2=pbt.entry.missing1,
                    adj2=pbt.entry.adj1,
                    sex=pbt.sex,
                    pop=pbt.pop,
                    chet=hl.or_missing(
                        pbt.gt1.ploidy == pbt.gt2.ploidy,
                        hl.range(0, pbt.gt1.ploidy).any(lambda a: pbt.gt1[a] != pbt.gt2[a])
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
        chet=(
            hl.case(missing_false=True)
                .when(hl.is_defined(tm.father.chet) & tm.father.chet == tm.mother.chet, tm.father.chet)
                .when(hl.is_defined(tm.father.chet) & hl.is_missing(tm.mother.chet), tm.father.chet)
                .when(hl.is_missing(tm.father.chet) & hl.is_defined(tm.mother.chet), tm.mother.chet)
                .or_missing()
        ),
        adj1=tm.child.adj1 & tm.father.adj1 & tm.mother.adj1,
        adj2=tm.child.adj2 & tm.father.adj2 & tm.mother.adj2,
        pop=(
            hl.case(missing_false=True)
                .when(tm.father.pop == tm.mother.pop, tm.father.pop)
                .when(hl.is_defined(tm.father.chet) & hl.is_missing(tm.mother.chet), tm.father.pop)
                .when(hl.is_defined(tm.mother.chet) & hl.is_missing(tm.father.chet), tm.mother.pop)
                .or_missing()
        )
    )

    tm.write(pbt_trio_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom), overwrite=args.overwrite)

    # Create entries table
    tm = hl.read_matrix_table(pbt_trio_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))
    et = tm.entries()
    et = et.filter(hl.is_defined(et.chet))
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
            mt = mt.key_cols_by('s', trio_id=mt.source_trio.id)
            mt = extract_pbt_probands(mt, data_type)
            mt = mt.filter_entries(mt.GT.is_non_ref())
            mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
            mt = mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)
            mt = filter_freq_and_csq(mt, data_type, max_freq=1.0, least_consequence=args.least_consequence)
            # mt = mt.annotate_entries(GT=mt.PBT_GT)
        else:
            mt = get_gnomad_data(data_type, non_refs_only=True)
            mt = mt.filter_cols(mt.meta.high_quality)
            mt = mt.select_entries('GT', 'adj', 'is_missing')
            mt = filter_freq_and_csq(mt, data_type, args.max_freq, args.least_consequence)

        if args.chrom:
            print(f"Selecting chrom {args.chrom}")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.chrom)])

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




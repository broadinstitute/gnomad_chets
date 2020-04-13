from gnomad.utils.vep import CSQ_ORDER
from gnomad.utils.slack import try_slack
import hail as hl
from resources import *
from gnomad_qc.v2.resources import annotations_ht_path, pbt_phased_trios_mt_path, get_gnomad_meta, fam_path, methylation_sites_ht_path, get_gnomad_data
import argparse
import logging
from typing import List

logger = logging.getLogger("create_vp_matrix")
BAD_THAI_TRIOS_PROJECT_ID = 'C978'


def create_variant_pair_ht(mt: hl.MatrixTable, row_groups: List[str]):
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
    et = et.key_by(locus2=et.vgt.v2.locus, alleles2=et.vgt.v2.alleles, locus1=et.vgt.v1.locus, alleles1=et.vgt.v1.alleles)
    et = et.distinct()

    return et


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

    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    vep_consequences = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1]))

    mt = mt.select_rows(
        vep=(
            hl.set(
                vep_ht[mt.row_key].vep.transcript_consequences
                    .filter(
                    lambda tc: (tc.biotype == 'protein_coding') &
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


def create_full_vp(data_type, path_args, args):
    vp_ht = hl.read_table(vp_list_ht_path(*path_args))
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
        mt = hl.read_matrix_table(adj_missing_mt_path(data_type))

    vp_ht = vp_ht.key_by('locus2', 'alleles2')
    vp_ht = vp_ht.select(locus1=vp_ht.locus1, alleles1=vp_ht.alleles1)
    vp_mt = mt.annotate_rows(v1=vp_ht.index(mt.row_key, all_matches=True))
    if args.pbt:
        vp_mt = vp_mt.select_entries(gt2=vp_mt.gt, missing2=vp_mt.missing, adj2=vp_mt.adj, trio_adj2=vp_mt.trio_adj)
    else:
        vp_mt = vp_mt.select_entries(gt2=vp_mt.gt, missing2=vp_mt.missing, adj2=vp_mt.adj)

    vp_mt = vp_mt.explode_rows(vp_mt.v1)
    vp_mt = vp_mt.transmute_rows(**vp_mt.v1)
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
    mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
    if args.pbt:
        vp_mt = vp_mt.annotate_entries(gt1=mt_joined.gt, missing1=mt_joined.missing, adj1=mt_joined.adj, trio_adj1=mt_joined.trio_adj)
    else:
        vp_mt = vp_mt.annotate_entries(gt1=mt_joined.gt, missing1=mt_joined.missing, adj1=mt_joined.adj)

    vp_mt.write(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp.mt', overwrite=True)
    vp_mt = hl.read_matrix_table(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp.mt')
    vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')

    vp_mt.write(full_mt_path(*path_args), overwrite=args.overwrite)


def create_vp_summary(data_type, path_args, args):
    mt = hl.read_matrix_table(full_mt_path(data_type, False, args.least_consequence, args.max_freq, args.chrom))
    meta = get_gnomad_meta(data_type).select('pop', 'release')
    mt = mt.annotate_cols(**meta[mt.col_key])
    mt = mt.filter_cols(mt.release)

    if args.pbt:
        pbt_samples = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom)).cols().key_by('s')
        mt = mt.filter_cols(hl.is_missing(pbt_samples[mt.col_key]))

    # mt = mt.annotate_cols(pop=[mt.pop, 'all']).explode_cols('pop')
    mt = mt.select_entries('adj1', 'adj2', gt_array=get_counts_agg_expr(mt))
    ht = mt.annotate_rows(
        gt_counts=hl.agg.group_by(
            mt.pop,
            hl.struct(
                raw=hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.gt_array),
                adj=hl.agg.filter(mt.adj1 & mt.adj2, hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.gt_array))
            )
        )
    ).rows()

    ht = ht.select(
        'locus2',
        'alleles2',
        gt_counts=hl.bind(
            lambda x: hl.dict(
                hl.zip(ht.gt_counts.keys(), ht.gt_counts.values()).append(
                    (
                        'all',
                        hl.fold(
                            lambda i, j: hl.struct(raw=i.raw + j.raw, adj=i.adj+j.adj),
                            x[0],
                            x[1:]
                        )
                    )
                )
            ),
            ht.gt_counts.values()
        )
    )
    ht.describe()
    ht.write(f'gs://gnomad-tmp/compound_hets/ht_sites{"_pbt" if args.pbt else ""}_by_pop.ht', overwrite=True)
    ht = hl.read_table(f'gs://gnomad-tmp/compound_hets/ht_sites{"_pbt" if args.pbt else ""}_by_pop.ht')
    ht = ht.key_by('locus1', 'alleles1', 'locus2', 'alleles2')
    ht = ht.repartition(1000, shuffle=False)
    ht.write(vp_count_ht_path(*path_args), overwrite=args.overwrite)


def extract_pbt_probands(pbt_mt: hl.MatrixTable, data_type: str):

    # Keep a single proband from each family with > 1  proband.
    meta = get_gnomad_meta(data_type)
    hq_samples = hl.literal(meta.aggregate(hl.agg.filter(meta.high_quality & (meta.project_id != BAD_THAI_TRIOS_PROJECT_ID), hl.agg.collect(meta.s))))
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

def create_vp_ann(path_args, args):
    data_type = path_args[0]
    vp_ht = hl.read_matrix_table(full_mt_path(*path_args)).rows()
    # vp_ht = hl.read_table(vp_list_ht_path(*path_args))

    def get_freq_expr(ht,  freq_ht, freq_dict, by_pop):
        pop_index = freq_dict[ht.pop] if by_pop else 0
        freq_expr = freq_ht[ht.key].freq[pop_index]
        return hl.struct(ac=freq_expr.AC, an=freq_expr.AN, af=freq_expr.AF)

    # Annotate freq, VEP and CpG information
    methyation_ht = hl.read_table(methylation_sites_ht_path())
    freq_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    rf_ht = hl.read_table(annotations_ht_path(data_type, 'rf'))
    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    vep_ht = vep_ht.annotate(
        vep=vep_ht.vep.select(
            'ancestral',
            'most_severe_consequence',
            'transcript_consequences'
        )
    )
    freq_meta = hl.eval(freq_ht.globals.freq_meta)
    freq_dict = {f['pop']: i for i, f in enumerate(freq_meta[:10]) if 'pop' in f}
    freq_dict['all'] = 0
    freq_dict = hl.literal(freq_dict)
    ht_ann = vp_ht.key_by('locus2', 'alleles2').select('locus1', 'alleles1')
    ht_ann = ht_ann.annotate(
        filters2=rf_ht[ht_ann.key].filters,
        freq2=freq_dict.map_values(lambda i: freq_ht[ht_ann.key].freq[i]),
        cpg2=methyation_ht[ht_ann.locus2].MEAN > 0.6,
        snv2=hl.is_snp(ht_ann.alleles2[0], ht_ann.alleles2[1]),
        vep2=vep_ht[ht_ann.key].vep
    )
    ht_ann.write(f'gs://gnomad-tmp/compound_hets/{data_type}_ann2.ht', overwrite=True)
    ht_ann = hl.read_table(f'gs://gnomad-tmp/compound_hets/{data_type}_ann2.ht')
    ht_ann = ht_ann.key_by('locus1', 'alleles1')
    ht_ann = ht_ann.annotate(
        filters1=rf_ht[ht_ann.key].filters,
        freq1=freq_dict.map_values(lambda i: freq_ht[ht_ann.key].freq[i]),
        cpg1=methyation_ht[ht_ann.locus1].MEAN > 0.6,
        snv1=hl.is_snp(ht_ann.alleles1[0], ht_ann.alleles1[1]),
        vep1=vep_ht[ht_ann.key].vep
    )
    ht_ann = ht_ann.key_by('locus1', 'alleles1', 'locus2', 'alleles2')
    ht_ann.write(vp_ann_ht_path(*path_args), overwrite=args.overwrite)


def create_pbt_summary(data_type, path_args, args):

    pbt = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))

    # Compute counts, grouped by pop
    meta = get_gnomad_meta('exomes').select('pop')
    pbt = pbt.key_cols_by('s')
    pbt = pbt.annotate_cols(**meta[pbt.col_key])

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
    pbt = pbt.repartition(1000, shuffle=False)
    pbt.write(pbt_phase_count_ht_path(*path_args), overwrite=args.overwrite)


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

    et.write(pbt_trio_et_path(data_type, True, args.least_consequence, args.max_freq, args.chrom), overwrite=args.overwrite)


def main(args):

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]

    if args.create_adj_missing_mt:
        mt = get_gnomad_data(data_type).select_cols() if not args.pbt else hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
        mt = mt.select_rows()
        mt = mt.select_entries(
            GT=hl.or_missing(mt.GT.is_non_ref(), mt.GT),
            missing=hl.is_missing(mt.GT),
            adj=mt.adj
        ).select_cols().select_rows()

        if args.pbt:
            mt = mt.key_cols_by('s', trio_id=mt.source_trio.id)
            mt = extract_pbt_probands(mt, data_type)
            mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
            mt = mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)
        else:
            meta = get_gnomad_meta('exomes')
            mt = mt.filter_cols(meta[mt.col_key].high_quality)

        mt.write(adj_missing_mt_path(data_type, args.pbt), overwrite=args.overwrite)

    if args.create_vp_list:
        mt = hl.read_matrix_table(adj_missing_mt_path(data_type, args.pbt))

        if args.chrom:
            print(f"Selecting chrom {args.chrom}")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.chrom)])

        mt = filter_freq_and_csq(mt, data_type, args.max_freq, args.least_consequence)
        vp_ht = create_variant_pair_ht(mt, ['gene_id'])
        vp_ht.write(vp_list_ht_path(*path_args), overwrite=args.overwrite)

    if args.create_full_vp:
        create_full_vp(data_type, path_args, args)

    if args.create_vp_ann:
        create_vp_ann(path_args, args)

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
    parser.add_argument('--create_adj_missing_mt', help='Creates a gnomAD MT with only missing and adj fields.', action='store_true')
    parser.add_argument('--create_vp_list', help='Creates a HT containing all variant pairs but no other data.', action='store_true')
    parser.add_argument('--create_vp_ann', help='Creates a  HT with freq and methylation information for all variant pairs.', action='store_true')
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




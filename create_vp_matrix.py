from gnomad.utils.vep import CSQ_ORDER
import hail as hl
from resources import *
from gnomad_qc.v2.resources import annotations_ht_path, pbt_phased_trios_mt_path, get_gnomad_meta, fam_path, methylation_sites_ht_path, get_gnomad_data
from gnomad.resources.grch37 import lcr_intervals, decoy_intervals, seg_dup_intervals
import argparse
import logging
from typing import List
from chet_utils import vep_genes_expr

logger = logging.getLogger("create_vp_matrix")


def _get_ordered_vp_struct(v1: hl.expr.StructExpression, v2: hl.expr.StructExpression):
    return hl.if_else(
        v1.locus.position < v2.locus.position,
        hl.struct(v1=v1, v2=v2),
        hl.if_else(
            v1.locus.position == v2.locus.position,  # If positions are equal, sort on alt allele
            hl.if_else(
                v1.alleles[1] < v2.alleles[1],
                hl.struct(v1=v1, v2=v2),
                hl.struct(v1=v2, v2=v1)

            ),
            hl.struct(v1=v2, v2=v1)
        )

    )


def create_variant_pair_ht2(mt: hl.MatrixTable, row_groups: List[str]):
    """
    Create a variant-pair MatrixTable containing all variant-pairs that appear both in the same individual within the row group.
    E.g., if the row group is a gene, this creates a variant-pair MT with all variant-pairs in a given gene.

    Note that this produces an empty shell of an MT with no data in the rows, columns or entries.

    :param MatrixTable mt: Input MatrixTable
    :param list of str row_groups: Row annotations for delimiting variant-pairs territory
    :return: Variant-pair MT
    :rtype: MatrixTable
    """

    # mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = mt.select_cols().select_rows(*row_groups)
    mt = mt.select_entries(x=True)
    gmt = mt.group_rows_by(*row_groups).aggregate(
        vgt=hl.agg.filter(
            hl.is_defined(mt.x),
            hl.agg.collect(hl.struct(locus=mt.locus, alleles=mt.alleles))
        )
    )

    ht = gmt.select_rows(
        vps=hl.agg.explode(
            lambda x: hl.agg.collect_as_set(x),
            hl.range(0, hl.len(gmt.vgt))
                .flatmap(lambda i1: hl.range(i1 + 1, hl.len(gmt.vgt))
                         .map(lambda i2: _get_ordered_vp_struct(gmt.vgt[i1], gmt.vgt[i2])))
        )
    ).rows()

    ht = ht.explode(ht.vps)
    ht = ht.key_by(locus2=ht.vps.v2.locus, alleles2=ht.vps.v2.alleles, locus1=ht.vps.v1.locus, alleles1=ht.vps.v1.alleles)
    ht = ht.distinct()
    ht = ht.select()

    return ht


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

    # mt = mt.filter_entries(mt.GT.is_non_ref())
    et = mt.select_cols().select_rows(*row_groups).entries()
    et = et.group_by(*row_groups, *mt.col_key)._set_buffer_size(5).aggregate(
        vgt=(hl.agg.collect(hl.struct(locus=et.locus, alleles=et.alleles)))
    )

    et = et.annotate(
        vgt=hl.range(0, hl.len(et.vgt))
                      .flatmap(lambda i1: hl.range(i1+1, hl.len(et.vgt))
                               .map(lambda i2: _get_ordered_vp_struct(et.vgt[i1], et.vgt[i2])))
    )

    et = et.explode(et.vgt)
    et = et.key_by(locus2=et.vgt.v2.locus, alleles2=et.vgt.v2.alleles, locus1=et.vgt.v1.locus, alleles1=et.vgt.v1.alleles)
    et = et.select().distinct()

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

    mt = mt.select_rows(
        vep=vep_genes_expr(vep_ht[mt.row_key].vep, least_consequence),
        af=hl.float32(freq[mt.row_key].freq[0].AF)
    )

    mt = mt.filter_rows(hl.is_defined(mt.vep) & (hl.len(mt.vep) > 0) & (mt.af > 0) & (mt.af <= max_freq))
    mt = mt.explode_rows(mt.vep)
    mt = mt.rename({'vep': 'gene_id'})
    return mt


def get_counts_agg_expr(mt: hl.MatrixTable):
    return (
        hl.case(missing_false=True)
            # 0x
            .when(hl.is_missing(mt.GT1) & ~mt.missing1,
                  hl.case(missing_false=True)
                  .when(hl.is_missing(mt.GT2) & ~mt.missing2, [1, 0, 0, 0, 0, 0, 0, 0, 0])
                  .when(mt.GT2.is_het(), [0, 1, 0, 0, 0, 0, 0, 0, 0])
                  .when(mt.GT2.is_hom_var(), [0, 0, 1, 0, 0, 0, 0, 0, 0])
                  .default([0, 0, 0, 0, 0, 0, 0, 0, 0]))
            # 1x
            .when(mt.GT1.is_het(),
                  hl.case(missing_false=True)
                  .when(hl.is_missing(mt.GT2) & ~mt.missing2, [0, 0, 0, 1, 0, 0, 0, 0, 0])
                  .when(mt.GT2.is_het(), [0, 0, 0, 0, 1, 0, 0, 0, 0])
                  .when(mt.GT2.is_hom_var(), [0, 0, 0, 0, 0, 1, 0, 0, 0])
                  .default([0, 0, 0, 0, 0, 0, 0, 0, 0]))
            # 2x
            .when(mt.GT1.is_hom_var(),
                  hl.case(missing_false=True)
                  .when(hl.is_missing(mt.GT2) & ~mt.missing2, [0, 0, 0, 0, 0, 0, 1, 0, 0])
                  .when(mt.GT2.is_het(), [0, 0, 0, 0, 0, 0, 0, 1, 0])
                  .when(mt.GT2.is_hom_var(), [0, 0, 0, 0, 0, 0, 0, 0, 1])
                  .default([0, 0, 0, 0, 0, 0, 0, 0, 0]))
            .default([0, 0, 0, 0, 0, 0, 0, 0, 0])
    )


def create_full_vp(
        mt: hl.MatrixTable,
        vp_list_ht: hl.Table,
        data_type: str
):
    # TODO: This implementation was causing memory challenges.

    vp_list_ht = vp_list_ht.key_by('locus2', 'alleles2')
    vp_list_ht = vp_list_ht.select(locus1=vp_list_ht.locus1, alleles1=vp_list_ht.alleles1)
    vp_mt = mt.annotate_rows(v1=vp_list_ht.index(mt.row_key, all_matches=True))
    vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)
    vp_mt = vp_mt.rename({x: f'{x}2' for x in vp_mt.entry})

    vp_mt = vp_mt.explode_rows(vp_mt.v1)
    vp_mt = vp_mt.transmute_rows(**vp_mt.v1)
    vp_mt = vp_mt.checkpoint(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp0.mt', overwrite=True)

    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
    vp_mt = vp_mt.checkpoint(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp1.mt', overwrite=True)

    mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
    vp_mt = vp_mt.annotate_entries(**{f'{x}1': mt_joined[x] for x in mt.entry})
    vp_mt = vp_mt.checkpoint(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp2.mt', overwrite=True)
    vp_mt = vp_mt.repartition(10000, shuffle=True)
    vp_mt = vp_mt.checkpoint(f'gs://gnomad-tmp/compound_hets/{data_type}_vp_mt_tmp3.mt', overwrite=True)
    vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')

    return vp_mt


def create_vp_summary(mt: hl.MatrixTable) -> hl.Table:
    mt = mt.select_entries('adj1', 'adj2', gt_array=get_counts_agg_expr(mt))
    ht = mt.annotate_rows(
        gt_counts=hl.agg.group_by(
            mt.pop,
            hl.struct(
                raw=hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.gt_array),
                adj=hl.or_else(
                    hl.agg.filter(mt.adj1 & mt.adj2, hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.gt_array)),
                    [0, 0, 0, 0, 0, 0, 0, 0, 0] # In case there are no adj entries
                )
            )
        )
    ).rows()

    ht = ht.select(
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

    ht = ht.checkpoint(f'gs://gnomad-tmp/compound_hets/ht_sites_by_pop.ht', overwrite=True)
    ht = ht.key_by('locus1', 'alleles1', 'locus2', 'alleles2')
    return ht.repartition(1000, shuffle=False)


def create_vp_ann(
        vp_ht: hl.Table,
        data_type
) -> hl.Table:


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
    lcr_ht =lcr_intervals.ht()
    decoy_ht = decoy_intervals.ht()
    seg_dup_ht = seg_dup_intervals.ht()

    freq_meta = hl.eval(freq_ht.globals.freq_meta)
    freq_dict = {f['pop']: i for i, f in enumerate(freq_meta[:10]) if 'pop' in f}
    freq_dict['all'] = 0
    freq_dict = hl.literal(freq_dict)
    ht_ann = vp_ht.key_by('locus2', 'alleles2').select('locus1', 'alleles1')

    _freq_ht_indexed = freq_ht[ht_ann.key]
    ht_ann = ht_ann.annotate(
        filters2=rf_ht[ht_ann.key].filters,
        freq2=freq_dict.map_values(lambda i: _freq_ht_indexed.freq[i]),
        popmax2=_freq_ht_indexed.popmax,
        cpg2=methyation_ht[ht_ann.locus2].MEAN > 0.6,
        snv2=hl.is_snp(ht_ann.alleles2[0], ht_ann.alleles2[1]),
        vep2=vep_ht[ht_ann.key].vep,
        lcr2=hl.is_defined(lcr_ht[ht_ann.locus2]),
        decoy2=hl.is_defined(decoy_ht[ht_ann.locus2]),
        segdup2=hl.is_defined(seg_dup_ht[ht_ann.locus2])
    )
    ht_ann = ht_ann.checkpoint(f'gs://gnomad-tmp/compound_hets/{data_type}_ann2.ht', overwrite=True)
    ht_ann = ht_ann.key_by('locus1', 'alleles1')
    _freq_ht_indexed = freq_ht[ht_ann.key]
    ht_ann = ht_ann.annotate(
        filters1=rf_ht[ht_ann.key].filters,
        freq1=freq_dict.map_values(lambda i: _freq_ht_indexed.freq[i]),
        popmax1=_freq_ht_indexed.popmax,
        cpg1=methyation_ht[ht_ann.locus1].MEAN > 0.6,
        snv1=hl.is_snp(ht_ann.alleles1[0], ht_ann.alleles1[1]),
        vep1=vep_ht[ht_ann.key].vep,
        lcr1=hl.is_defined(lcr_ht[ht_ann.locus1]),
        decoy1=hl.is_defined(decoy_ht[ht_ann.locus1]),
        segdup1=hl.is_defined(seg_dup_ht[ht_ann.locus1])
    )
    return ht_ann.key_by('locus1', 'alleles1', 'locus2', 'alleles2')


def create_pbt_summary(data_type, path_args, args):

    pbt = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))

    # Compute counts, grouped by pop
    meta = get_gnomad_meta('exomes').select('pop')
    pbt = pbt.key_cols_by('s')
    pbt = pbt.annotate_cols(**meta[pbt.col_key])

    def get_trio_phase_expr(pbt) -> hl.expr.StructExpression:
        s_info = hl.agg.collect(
            hl.struct(
                same_hap=pbt.GT1[0] == pbt.GT2[0], # Works since het/het phased GTs only
                adj=pbt.trio_adj1 & pbt.trio_adj2,
                s=pbt.s
            )
        )

        def get_return_struct(s_info):
            ret = hl.struct(
                    same_hap_samples=s_info.filter(lambda x: x.same_hap).map(lambda x: x.s),
                    chet_samples=s_info.filter(lambda x: ~x.same_hap).map(lambda x: x.s),
                )
            return ret.annotate(
                n_same_hap=hl.len(ret.same_hap_samples),
                n_chet = hl.len(ret.chet_samples)
            )

        return hl.struct(
                adj=get_return_struct(s_info.filter(lambda x: x.adj)),
                raw=get_return_struct(s_info)
            )

    pbt = pbt.filter_entries(
        pbt.GT1.phased &
        pbt.GT2.phased &
        (pbt.GT1.ploidy == pbt.GT2.ploidy) &
        pbt.GT1.is_het() &
        pbt.GT2.is_het()
    )

    #TEST
    pbt = pbt.annotate_cols(pop=[pbt.pop, 'all'])
    pbt = pbt.annotate_rows(
        phase_by_pop=hl.agg.explode(
            lambda pop: hl.agg.group_by(
                pop,
                get_trio_phase_expr(pbt)
            ),
            pbt.pop
        )
    ).rows()

    # pbt = pbt.filter(pbt.phase_by_pop['all'].raw.n_same_hap + pbt.phase_by_pop['all'].raw.n_chet > 0) # I think that's not needed
    pbt = pbt.filter(hl.len(pbt.phase_by_pop)>0)
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
                    GT1=pbt.entry.GT1,
                    missing1=pbt.entry.missing1,
                    adj1=pbt.entry.adj1,
                    GT2=pbt.entry.GT1,
                    missing2=pbt.entry.missing1,
                    adj2=pbt.entry.adj1,
                    sex=pbt.sex,
                    pop=pbt.pop,
                    chet=hl.or_missing(
                        pbt.GT1.ploidy == pbt.GT2.ploidy,
                        hl.range(0, pbt.GT1.ploidy).any(lambda a: pbt.GT1[a] != pbt.GT2[a])
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
    tm = tm.checkpoint("gs://gnomad-tmp/pbt_vp_tm.tmp.mt", overwrite=True)

    # annotate variant-phase phase per trio / overall
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


def get_pbt_mt(data_type) -> hl.MatrixTable:
    mt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
    mt = mt.key_cols_by('s', trio_id=mt.source_trio.id)
    mt = extract_pbt_probands(mt, data_type)
    mt = mt.key_cols_by('s')
    return mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)


def main(args):

    hl.init(log="/tmp/hail_vp.log")

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]

    if args.create_vp_list:

        if args.pbt:
            mt = get_pbt_mt(data_type)
        else:
            mt = get_gnomad_data(data_type)
            mt = mt.filter_cols(mt.meta.high_quality)

        mt = mt.select_cols().select_rows()
        mt = mt.filter_entries(mt.GT.is_non_ref())
        mt = mt.select_entries()

        if not args.pbt:
            mt = mt.filter_cols(get_gnomad_meta('exomes')[mt.col_key].high_quality)

        if args.chrom:
            print(f"Selecting chrom {args.chrom}")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(args.chrom)])

        mt = filter_freq_and_csq(mt, data_type, args.max_freq, args.least_consequence)
        mt = mt.checkpoint('gs://gnomad-tmp/pre_vp_ht2.mt', overwrite=True)
        mt = mt.filter_rows(hl.is_defined(mt.gene_id))
        mt = mt.repartition(11000)
        mt = mt.checkpoint('gs://gnomad-tmp/pre_vp_ht_rep.mt', overwrite=True)

        if args.vp_list_by_chrom:
            chroms = [str(x) for x in range(1,23)] + ['X']
            for chrom in chroms:
                logger.info(f"Now writing VP list HT for chrom {chrom}")

                c_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chrom)])
                vp_ht = create_variant_pair_ht(c_mt, ['gene_id'])
                vp_ht.write(vp_list_ht_path(*path_args[:-1], chrom=chrom), overwrite=args.overwrite)

            chrom_hts = [hl.read_table(vp_list_ht_path(*path_args[:-1], chrom=chrom)) for chrom in chroms]
            vp_ht = chrom_hts[0].union(*chrom_hts[1:])
        else:
            vp_ht = create_variant_pair_ht(mt, ['gene_id'])

        vp_ht.write(vp_list_ht_path(*path_args[:-1]), overwrite=args.overwrite)

    if args.create_full_vp:
        if args.pbt:
            mt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
            mt = mt.key_cols_by('s', trio_id=mt.source_trio.id)
            mt = extract_pbt_probands(mt, data_type)
            mt = mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)
            mt = mt.select_entries(
                GT=hl.case()
                    .when(mt.PBT_GT.is_non_ref(), mt.PBT_GT)
                    .when(mt.GT.is_non_ref(), hl.call(mt.GT[0], mt.GT[1]))  # Unphase genotpes phased by GATK HC
                    .or_missing(),
                missing=hl.is_missing(mt.GT),
                adj=mt.adj,
                trio_adj=mt.trio_adj
            ).select_cols().select_rows()
        else:
            mt = get_gnomad_data(data_type)
            mt = mt.select_entries(
                GT=hl.or_missing(mt.GT.is_non_ref(), mt.GT),
                PID=mt.PID,
                missing=hl.is_missing(mt.GT),
                adj=mt.adj
            ).select_cols().select_rows()
            meta = get_gnomad_meta('exomes')
            mt = mt.filter_cols(meta[mt.col_key].high_quality)

        logger.info(f"Reading VP list from {vp_list_ht_path(*path_args)}")
        vp_mt = create_full_vp(
            mt,
            vp_list_ht=hl.read_table(vp_list_ht_path(*path_args)),
            data_type=data_type
        )
        vp_mt.write(full_mt_path(*path_args), overwrite=args.overwrite)

    if args.create_vp_ann:
        vp_ht = hl.read_matrix_table(full_mt_path(*path_args)).rows()
        ht_ann = create_vp_ann(
            vp_ht,
            data_type
        )
        ht_ann.write(vp_ann_ht_path(*path_args), overwrite=args.overwrite)

    if args.create_vp_summary:
        mt = hl.read_matrix_table(full_mt_path(data_type, False, args.least_consequence, args.max_freq, args.chrom))
        meta = get_gnomad_meta(data_type).select('pop', 'release')
        mt = mt.annotate_cols(**meta[mt.col_key])
        mt = mt.filter_cols(mt.release)

        if args.pbt:
            pbt_samples = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom)).cols().key_by('s')
            mt = mt.filter_cols(hl.is_missing(pbt_samples[mt.col_key]))

        ht = create_vp_summary(mt)
        ht.write(vp_count_ht_path(*path_args), overwrite=args.overwrite)

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
    parser.add_argument('--create_vp_list', help='Creates a HT containing all variant pairs but no other data.', action='store_true')
    parser.add_argument('--vp_list_by_chrom', help=f'If set, computes the VP HT by chrom first and then union them', action='store_true')
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
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chrom', help='Only run on given chromosome')

    args = parser.parse_args()
    main(args)




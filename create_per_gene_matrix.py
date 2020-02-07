from gnomad_hail import *
from gnomad_qc.v2.resources import get_gnomad_meta, get_gnomad_data, annotations_ht_path

import argparse
from resources import *

CHET_THRESHOLD = 0.5
SAME_HAP_THRESHOLD = 0.5

CSQ_CODES = {
        'lof': 0,
        'damaging_missense': 1,
        'missense_variant': 2,
        'synonymous_variant': 3
    }


def filter_to_chr20(tables: List[Union[hl.Table, hl.MatrixTable]]) ->  List[Union[hl.Table, hl.MatrixTable]]:
    return [hl.filter_intervals(t, [hl.parse_locus_interval('20')]) for t in tables]


def get_group_to_counts_expr(k: hl.expr.StructExpression, counts: hl.expr.DictExpression) -> hl.expr.ArrayExpression:
    return hl.range(1, k.snv - 1, step=-1).flatmap(
        lambda snv: hl.range(0, k.af_gt_0_001 + 1).flatmap(
            lambda af: hl.range(0, k.csq + 1).map(
                lambda csq: hl.struct(snv=hl.bool(snv), af_gt_0_001=hl.bool(af), csq=csq)
            )
        )
    ).map(
        lambda key: counts.get(key)
    )

def get_worst_gene_csq_code_expr(vep_expr: hl.expr.StructExpression) -> hl.expr.DictExpression:

    worst_gene_csq_expr = vep_expr.transcript_consequences.map(
        lambda ts: ts.select(
            'gene_id',
            'gene_symbol',
            csq=(
                hl.case(missing_false=True)
                .when(ts.lof == 'HC', CSQ_CODES['lof'])
                .when(ts.polyphen_prediction == 'probably_damaging', CSQ_CODES['damaging_missense'])
                .when(ts.consequence_terms.any(lambda x: x == 'missense_variant'), CSQ_CODES['missense_variant'])
                .when(ts.consequence_terms.all(lambda x: x == 'synonymous_variant'), CSQ_CODES['synonymous_variant'])
                .or_missing()
            )
        )
    )

    worst_gene_csq_expr = worst_gene_csq_expr.filter(lambda x: hl.is_defined(x.csq))
    worst_gene_csq_expr = worst_gene_csq_expr.group_by(lambda x: x.gene_id)
    worst_gene_csq_expr = worst_gene_csq_expr.map_values(
        lambda x: hl.sorted(x, key=lambda y: y.csq)[0]
    )

    return worst_gene_csq_expr


def compute_from_vp_mt(chr20: bool, overwrite: bool):

    def get_popmax_af(freq_dict):
        pops_to_exclude = hl.literal({'asj', 'fin', 'oth'})
        return hl.max(freq_dict.keys().filter(lambda x: ~pops_to_exclude.contains(x)).map(lambda x: freq_dict[x].AF))

    vp_mt = hl.read_matrix_table(full_mt_path('exomes'))
    ann_ht = hl.read_table(vp_ann_ht_path('exomes'))
    phase_ht = hl.read_table(phased_vp_count_ht_path('exomes'))

    if chr20:
        vp_mt, ann_ht, phase_ht = filter_to_chr20([vp_mt, ann_ht, phase_ht])

    vep1_expr = get_worst_gene_csq_code_expr(ann_ht.vep1)
    vep2_expr = get_worst_gene_csq_code_expr(ann_ht.vep2)
    ann_ht = ann_ht.select(
        'snv1',
        'snv2',
        popmax=hl.max(get_popmax_af(ann_ht.freq1), get_popmax_af(ann_ht.freq2)),
        filtered=(hl.len(ann_ht.filters1) > 0) | (hl.len(ann_ht.filters2) > 0),
        vep=vep1_expr.keys().filter(
            lambda k: vep2_expr.contains(k)
        ).map(
            lambda k: vep1_expr[k].annotate(
                csq=hl.if_else(
                    vep1_expr[k].csq < vep2_expr[k].csq,
                    vep1_expr[k].csq,
                    vep2_expr[k].csq
                )
            )
        )
    )

    meta = get_gnomad_meta('exomes')
    vp_mt = vp_mt.annotate_cols(
        pop=meta[vp_mt.col_key].pop
    )
    vp_mt = vp_mt.annotate_rows(
        **ann_ht[vp_mt.row_key],
        phase_info=phase_ht[vp_mt.row_key].phase_info
    )

    vp_mt = vp_mt.filter_rows(
        ~vp_mt.filtered &
        (vp_mt.popmax <= 0.01)
    )

    vp_mt = vp_mt.filter_entries(
        vp_mt.gt1.is_het() & vp_mt.gt2.is_het() & vp_mt.adj1 & vp_mt.adj2
    )

    vp_mt = vp_mt.select_entries(
        x=True
    )

    vp_mt = vp_mt.explode_rows('vep')
    vp_mt = vp_mt.transmute_rows(
        **vp_mt.vep
    )

    vp_mt.describe()

    vp_mt = vp_mt.group_rows_by(
        'gene_id',
        'gene_symbol'
    ).aggregate(
        counts=hl.agg.filter(
            vp_mt.x,
            hl.agg.group_by(
                hl.struct(
                    snv=vp_mt.snv1 & vp_mt.snv2,
                    af_gt_0_001=(vp_mt.popmax > 0.001),
                    csq=vp_mt.csq
                ),
                hl.struct(
                    n=hl.agg.count(),
                    max_p_chet=hl.agg.max(vp_mt.phase_info[vp_mt.pop].em.adj.p_chet)
                )
            )
        )
    )

    vp_mt = vp_mt.annotate_entries(
        counts=vp_mt.counts.keys().map(
            lambda k: (
                k,
                get_group_to_counts_expr(k, vp_mt.counts).fold(
                    lambda i, j: hl.struct(
                        n=i.n + j.n,
                        max_p_chet=hl.if_else(
                            i.max_p_chet > j.max_p_chet,
                            i.max_p_chet,
                            j.max_p_chet
                        )
                    ),
                    hl.struct(n=0, max_p_chet=-1.0)
                )
            )
        )
    )

    vp_mt.describe()

    gene_ht = vp_mt.annotate_rows(
        row_counts=hl.array(
            hl.agg.explode(
                lambda x:
                hl.agg.group_by(
                    x[0].annotate(
                        pop=vp_mt.pop,
                        csq=hl.literal(list(CSQ_CODES.keys()))[x[0].csq]
                    ),
                    hl.struct(
                        total=hl.agg.count_where(x[1].n > 0),
                        phased_c_het=hl.agg.count_where(x[1].max_p_chet > CHET_THRESHOLD),
                        phased_same_hap=hl.agg.count_where(x[1].max_p_chet < SAME_HAP_THRESHOLD)
                    )
                ),
                vp_mt.counts
            )
        )
    ).rows()

    gene_ht = gene_ht.explode('row_counts')
    gene_ht = gene_ht.select(
        **gene_ht.row_counts[0],
        **gene_ht.row_counts[1]
    )

    gene_ht.describe()

    # vp_mt = vp_mt.group_rows_by(
    #     'gene_id'
    # ).aggregate_rows(
    #     gene_symbol=hl.agg.take(vp_mt.gene_symbol, 1)[0]
    # ).aggregate(
    #     _lof_counts_0_01=get_entry_agg(vp_mt, hl.literal({'lof'}), 0.01),
    #     _lof_counts_0_001=get_entry_agg(vp_mt, hl.literal({'lof'}), 0.001),
    #     _damaging_missense_counts_0_01=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense'}), 0.01),
    #     _damaging_missense_counts_0_001=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense'}), 0.001),
    #     _missense_counts_0_01=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense', 'missense_variant'}), 0.01),
    #     _missense_counts_0_001=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense', 'missense_variant'}), 0.001),
    #     _synonymous_counts_0_01=get_entry_agg(vp_mt, hl.literal({'synonymous_variant'}), 0.01),
    #     _synonymous_counts_0_001=get_entry_agg(vp_mt, hl.literal({'synonymous_variant'}), 0.001)
    # ) #.checkpoint('gs://gnomad-tmp/vp_gene.mt', overwrite=overwrite) # TODO: Check whether this is actually harming performance
    #
    # gene_ht = vp_mt.annotate_rows(
    #     **{x[1:]: get_row_agg(vp_mt[x]) for x in vp_mt.entry if x != '_gene_symbol'}
    # ).rows()

    gene_ht = gene_ht.checkpoint('gs://gnomad-lfran/compound_hets/chet_per_gene.ht', overwrite=overwrite)

    gene_ht.flatten().export('gs://gnomad-lfran/compound_hets/chet_per_gene.tsv.gz')

def compute_from_full_mt(chr20: bool, overwrite: bool):

    mt = get_gnomad_data('exomes', adj=True, release_samples=True)
    freq_ht = hl.read_table(annotations_ht_path('exomes', 'frequencies'))
    vep_ht = hl.read_table(annotations_ht_path('exomes', 'vep'))
    rf_ht = hl.read_table(annotations_ht_path('exomes', 'rf'))

    if chr20:
        mt, freq_ht, vep_ht, rf_ht = filter_to_chr20([mt, freq_ht, vep_ht, rf_ht])

    vep_ht = vep_ht.annotate(
        vep=get_worst_gene_csq_code_expr(
            vep_ht.vep.annotate(
                transcript_consequences=vep_ht.vep.transcript_consequences.filter(
                    lambda tc: (tc.biotype == 'protein_coding')
                )
            )
        ).values()
    )
    mt = mt.annotate_rows(
        popmax=freq_ht[mt.row_key].popmax,
        vep=vep_ht[mt.row_key].vep,
        filters=rf_ht[mt.row_key].filters
    )
    mt = mt.filter_rows(
        (mt.popmax.AF <= 0.01) &
        (hl.len(mt.vep) > 0) &
        (hl.len(mt.filters) == 0)
    )

    mt = mt.explode_rows(mt.vep)
    mt = mt.transmute_rows(**mt.vep)

    mt = mt.group_rows_by(
        'gene_id'
    ).aggregate_rows(
        gene_symbol=hl.agg.take(mt.gene_symbol, 1)[0]
    ).aggregate(
        counts=hl.agg.group_by(
            hl.struct(
                snv=hl.is_snp(mt.alleles[0], mt.alleles[1]),
                af_gt_0_001=mt.popmax.AF > 0.001,
                csq=mt.csq
            ),
            hl.struct(
                n_hom=hl.agg.count_where(mt.GT.is_hom_var()),
                n_het=hl.agg.count_where(mt.GT.is_het())
            )
        )
    )

    mt.describe()

    mt = mt.annotate_entries(
        counts=mt.counts.keys().map(
            lambda k: (
                k,
                get_group_to_counts_expr(k, mt.counts).fold(
                    lambda i,j: hl.struct(
                        n_hom=i.n_hom+j.n_hom,
                        n_het=i.n_het+j.n_het
                    ),
                    hl.struct(n_hom=0, n_het=0)
                )
            )
        )
    )

    gene_ht = mt.annotate_rows(
        row_counts=hl.array(
            hl.agg.explode(
                lambda x:
                hl.agg.group_by(
                    x[0].annotate(
                        pop=mt.meta.pop,
                        csq=hl.literal(list(CSQ_CODES.keys()))[x[0].csq]
                    ),
                    hl.struct(
                        n_hom=hl.agg.count_where(x[1].n_hom > 0),
                        n_het=hl.agg.count_where(x[1].n_het == 1),
                        n_het_het=hl.agg.count_where(x[1].n_het > 1)
                    )
                ),
                mt.counts
            )
        )
    ).rows()

    gene_ht = gene_ht.explode('row_counts')
    gene_ht = gene_ht.select(
        'gene_symbol',
        **gene_ht.row_counts[0],
        **gene_ht.row_counts[1]
    )

    gene_ht = gene_ht.checkpoint(
        'gs://gnomad-lfran/compound_hets/het_and_hom_per_gene{}.ht'.format(
            '.chr20' if chr20 else ''
        ),
        overwrite=overwrite
    )

    gene_ht.flatten().export('gs://gnomad-lfran/compound_hets/het_and_hom_per_gene{}.tsv.gz'.format(
        '.chr20' if chr20 else ''
    ))


def main(args):
    if args.compute_from_vp_mt:
        compute_from_vp_mt(args.chr20, args.overwrite)
    if args.compute_from_full_mt:
        compute_from_full_mt(args.chr20, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compute_from_vp_mt', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--compute_from_full_mt', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chr20', help='Computes on chrom20 only', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    main(args)
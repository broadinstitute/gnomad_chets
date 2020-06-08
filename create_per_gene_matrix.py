from gnomad_qc.v2.resources import get_gnomad_meta, get_gnomad_data, annotations_ht_path
import hail as hl
from typing import List, Union
import argparse
from resources import *

CHET_THRESHOLD = 0.5
SAME_HAP_THRESHOLD = 0.5

CSQ_CODES = [
    'lof',
    'damaging_missense',
    'missense_variant',
    'synonymous_variant'
]


def filter_to_chr20(tables: List[Union[hl.Table, hl.MatrixTable]]) -> List[Union[hl.Table, hl.MatrixTable]]:
    return [hl.filter_intervals(t, [hl.parse_locus_interval('20')]) for t in tables]


def get_group_to_counts_expr(k: hl.expr.StructExpression, counts: hl.expr.DictExpression) -> hl.expr.ArrayExpression:
    return hl.range(1, k.snv - 1, step=-1).flatmap(
        lambda snv: hl.range(0, k.af_gt_0_001 + 1).flatmap(
            lambda af: hl.range(0, k.csq + 1).map(
                lambda csq: hl.struct(snv=hl.bool(snv), af_gt_0_001=hl.bool(af), csq=csq)
            )
        )
    ).filter(
        lambda key: counts.contains(key)
    ).map(
        lambda key: counts[key]
    )


def get_worst_gene_csq_code_expr(vep_expr: hl.expr.StructExpression) -> hl.expr.DictExpression:
    worst_gene_csq_expr = vep_expr.transcript_consequences.filter(
        lambda tc: tc.biotype == 'protein_coding'
    ).map(
        lambda ts: ts.select(
            'gene_id',
            'gene_symbol',
            csq=(
                hl.case(missing_false=True)
                    .when(ts.lof == 'HC', CSQ_CODES.index('lof'))
                    .when(ts.polyphen_prediction == 'probably_damaging', CSQ_CODES.index('damaging_missense'))
                    .when(ts.consequence_terms.any(lambda x: x == 'missense_variant'), CSQ_CODES.index('missense_variant'))
                    .when(ts.consequence_terms.all(lambda x: x == 'synonymous_variant'), CSQ_CODES.index('synonymous_variant'))
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
        pops_to_exclude = hl.literal({'asj', 'fin', 'oth', 'all'})  # TODO:  decide what to include/exclude here
        # pops_to_exclude = hl.literal({'oth'})
        return hl.max(
            freq_dict.keys().filter(
                lambda x: ~pops_to_exclude.contains(x) &
                          (freq_dict[x].AF > 0)
            ).map(lambda x: freq_dict[x].AF),
            filter_missing=False
        )

    meta = get_gnomad_meta('exomes')
    vp_mt = hl.read_matrix_table(full_mt_path('exomes'))
    vp_mt = vp_mt.filter_cols(meta[vp_mt.col_key].release)
    ann_ht = hl.read_table(vp_ann_ht_path('exomes'))
    phase_ht = hl.read_table(phased_vp_count_ht_path('exomes'))

    if chr20:
        vp_mt, ann_ht, phase_ht = filter_to_chr20([vp_mt, ann_ht, phase_ht])

    vep1_expr = get_worst_gene_csq_code_expr(ann_ht.vep1)
    vep2_expr = get_worst_gene_csq_code_expr(ann_ht.vep2)
    ann_ht = ann_ht.select(
        'snv1',
        'snv2',
        popmax=hl.max(get_popmax_af(ann_ht.freq1), get_popmax_af(ann_ht.freq2), filter_missing=False),
        filtered=(hl.len(ann_ht.filters1) > 0) | (hl.len(ann_ht.filters2) > 0),
        vep=vep1_expr.keys().filter(
            lambda k: vep2_expr.contains(k)
        ).map(
            lambda k: vep1_expr[k].annotate(
                csq=hl.max(vep1_expr[k].csq, vep2_expr[k].csq)
            )
        )
    )

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

    def get_grouped_phase_agg():
        return hl.agg.group_by(
            hl.case()
                .when(vp_mt.phase_info[vp_mt.pop].em.adj.p_chet > CHET_THRESHOLD, 1)
                .when(vp_mt.phase_info[vp_mt.pop].em.adj.p_chet < SAME_HAP_THRESHOLD, 2)
                .default(3)
            ,
            hl.agg.min(vp_mt.csq)
        )

    vp_mt = vp_mt.group_rows_by(
        'gene_id',
        'gene_symbol'
    ).aggregate(
        af_gt_0_001=hl.agg.filter(
            vp_mt.x,
            get_grouped_phase_agg()
        ),
        af_le_0_001=hl.agg.filter(
            (vp_mt.popmax <= 0.001) & vp_mt.x,
            get_grouped_phase_agg()
        )
    )

    vp_mt = vp_mt.checkpoint('gs://gnomad-tmp/test.mt', overwrite=True)

    gene_ht = vp_mt.annotate_rows(
        row_counts=hl.flatten([
            hl.array(
                hl.agg.group_by(
                    vp_mt.pop,
                    hl.struct(
                        csq=csq,
                        af=af,
                        # TODO: Review this
                        # These will only kept the worst csq -- now maybe it'd be better to keep either
                        # - the worst csq for chet or
                        # - the worst csq for both chet and same_hap
                        n_chet=hl.agg.count_where((vp_mt[af].get(1) == csq_i) & (vp_mt[af].get(2, 9) >= csq_i) & (vp_mt[af].get(3, 9) >= csq_i)),
                        n_same_hap=hl.agg.count_where((vp_mt[af].get(2) == csq_i) & (vp_mt[af].get(1, 9) > csq_i) & (vp_mt[af].get(1, 9) >= csq_i)),
                        n_unphased=hl.agg.count_where((vp_mt[af].get(3) == csq_i) & (vp_mt[af].get(1, 9) > csq_i) & (vp_mt[af].get(2, 9) > csq_i))
                    )
                )
            ).filter(
                lambda x: (x[1].n_chet > 0) | (x[1].n_same_hap > 0) | (x[1].n_unphased > 0)
            ).map(
                lambda x: x[1].annotate(
                    pop=x[0]
                )
            )
            for csq_i, csq in enumerate(CSQ_CODES)
            for af in ['af_gt_0_001', 'af_le_0_001']
        ])
    ).rows()

    gene_ht = gene_ht.explode('row_counts')
    gene_ht = gene_ht.select(
        **gene_ht.row_counts
    )

    gene_ht.describe()
    gene_ht = gene_ht.checkpoint(
        'gs://gnomad-lfran/compound_hets/chet_per_gene{}.ht'.format(
            '.chr20' if chr20 else ''
        ),
        overwrite=overwrite
    )

    gene_ht.flatten().export(
        'gs://gnomad-lfran/compound_hets/chet_per_gene{}.tsv.gz'.format(
            '.chr20' if chr20 else ''
        )
    )


def compute_from_full_mt(chr20: bool, overwrite: bool):
    mt = get_gnomad_data('exomes', adj=True, release_samples=True)
    mt = mt.filter_rows(mt.a_index == 1)  # TODO remove when VP_MT has been re-generated
    freq_ht = hl.read_table(annotations_ht_path('exomes', 'frequencies'))
    vep_ht = hl.read_table(annotations_ht_path('exomes', 'vep'))
    rf_ht = hl.read_table(annotations_ht_path('exomes', 'rf'))

    if chr20:
        mt, freq_ht, vep_ht, rf_ht = filter_to_chr20([mt, freq_ht, vep_ht, rf_ht])

    vep_ht = vep_ht.annotate(
        vep=get_worst_gene_csq_code_expr(vep_ht.vep).values()
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

    # # mt.describe()
    #
    # het_het_agg = hl.agg.take(mt.csq, 2, ordering=mt.csq)
    mt = mt.group_rows_by(
        'gene_id'
    ).aggregate_rows(
        gene_symbol=hl.agg.take(mt.gene_symbol, 1)[0]
    ).aggregate(
        counts=hl.agg.group_by(
            mt.popmax.AF > 0.001,
            hl.struct(
                hom_csq=hl.agg.filter(mt.GT.is_hom_var(), hl.agg.min(mt.csq)),
                het_csq=hl.agg.filter(mt.GT.is_het(), hl.agg.min(mt.csq)),
                het_het_csq=hl.sorted(
                    hl.array(
                        hl.agg.filter(mt.GT.is_het(), hl.agg.counter(mt.csq))
                    ),
                    key=lambda x: x[0]
                ).scan(
                    lambda i, j: (j[0], i[1] + j[1]),
                    (0, 0)
                ).find(
                    lambda x: x[1] > 1
                )[0]
                # het_het_csq=hl.agg.take(mt.csq, 2, ordering=mt.csq) -- this runs out of mem :(
                # het_het_csq=hl.agg.filter(mt.GT.is_het(), hl.or_missing(hl.len(het_het_agg)>1, het_het_agg[1]))
            )
        )
    )

    mt = mt.annotate_entries(
        counts=hl.struct(
            af_gt_0_001=hl.struct(
                hom_csq=hl.min(mt.counts.get(True).hom_csq, mt.counts.get(False).hom_csq),
                het_csq=hl.min(mt.counts.get(True).het_csq, mt.counts.get(False).het_csq),
                het_het_csq=hl.min(
                    mt.counts.get(True).het_het_csq,
                    mt.counts.get(False).het_het_csq,
                    hl.or_missing(
                        hl.is_defined(mt.counts.get(True).het_csq) & hl.is_defined(mt.counts.get(False).het_csq),
                        hl.max(mt.counts.get(True).het_csq, mt.counts.get(False).het_csq)
                    )
                ),
                # het_het_csq=hl.min(het_het_expr(mt.counts.get(True).het_het_csq), af_le_0_001_expr.het_het_csq)
            ),
            af_le_0_001=mt.counts.get(False)
        )
    )

    mt = mt.checkpoint('gs://gnomad-tmp/test_2.mt', overwrite=True)

    gene_ht = mt.annotate_rows(
        row_counts=hl.flatten([
            hl.array(
                hl.agg.group_by(
                    mt.meta.pop,
                    hl.struct(
                        csq=csq,
                        af=af,
                        n_hom=hl.agg.count_where(mt.counts[af].hom_csq == csq_i),
                        n_het=hl.agg.count_where(mt.counts[af].het_csq == csq_i),
                        n_het_het=hl.agg.count_where(mt.counts[af].het_het_csq == csq_i)
                    )
                )
            ).filter(
                lambda x: (x[1].n_het > 0) | (x[1].n_hom > 0) | (x[1].n_het_het > 0)
            ).map(
                lambda x: x[1].annotate(
                    pop=x[0]
                )
            )
            for csq_i, csq in enumerate(CSQ_CODES)
            for af in ['af_gt_0_001', 'af_le_0_001']
        ])
    ).rows()

    gene_ht = gene_ht.explode('row_counts')
    gene_ht = gene_ht.select(
        'gene_symbol',
        **gene_ht.row_counts
    )

    gene_ht.describe()

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

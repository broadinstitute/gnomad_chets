from gnomad_qc.v2.resources import get_gnomad_meta, get_gnomad_data, annotations_ht_path
from gnomad.resources.grch37.gnomad import public_release
import hail as hl
from typing import List, Union
import argparse
from resources import *

CHET_THRESHOLD = 0.505
SAME_HAP_THRESHOLD = 0.0164
ALLELE_FREQUENCY_CUTOFFS = [0.001] # 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.0015, 0.002

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
        lambda snv: hl.range(0, k.all + 1).flatmap(
            lambda af: hl.range(0, k.csq + 1).map(
                lambda csq: hl.struct(snv=hl.bool(snv), all=hl.bool(af), csq=csq)
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
        is_singleton_vp=(ann_ht.freq1['all'].AC < 2) & (ann_ht.freq2['all'].AC < 2),
        pop_af=hl.dict(
            ann_ht.freq1.key_set().intersection(ann_ht.freq2.key_set())
                .map(
                lambda pop: hl.tuple([pop, hl.max(ann_ht.freq1[pop].AF, ann_ht.freq2[pop].AF)])
            )
        ),
        popmax_af=hl.max(ann_ht.popmax1.AF, ann_ht.popmax2.AF, filter_missing=False),
        filtered=(hl.len(ann_ht.filters1) > 0) | (hl.len(ann_ht.filters2) > 0),
        vep=vep1_expr.keys().filter(
            lambda k: vep2_expr.contains(k)
        ).map(
            lambda k: vep1_expr[k].annotate(
                csq=hl.max(vep1_expr[k].csq, vep2_expr[k].csq)
            )
        )
    )
    ann_ht.describe()
    vp_mt = vp_mt.annotate_cols(
        pop=meta[vp_mt.col_key].pop
    )
    vp_mt = vp_mt.annotate_rows(
        **ann_ht[vp_mt.row_key],
        phase_info=phase_ht[vp_mt.row_key].phase_info
    )

    vp_mt = vp_mt.filter_rows(
        ~vp_mt.filtered
    )

    vp_mt = vp_mt.filter_entries(
        vp_mt.GT1.is_het() & vp_mt.GT2.is_het() & vp_mt.adj1 & vp_mt.adj2
    )

    vp_mt = vp_mt.select_entries(
        x=True
    )

    vp_mt = vp_mt.annotate_cols(
        pop=['all', vp_mt.pop]
    )
    vp_mt = vp_mt.explode_cols('pop')

    vp_mt = vp_mt.explode_rows('vep')
    vp_mt = vp_mt.transmute_rows(
        **vp_mt.vep
    )
    vp_mt.describe()
    def get_grouped_phase_agg():
        return hl.agg.group_by(
            hl.case()
                .when(~vp_mt.is_singleton_vp & (vp_mt.phase_info[vp_mt.pop].em.adj.p_chet > CHET_THRESHOLD), 1)
                .when(~vp_mt.is_singleton_vp & (vp_mt.phase_info[vp_mt.pop].em.adj.p_chet < SAME_HAP_THRESHOLD), 2)
                .default(3)
            ,
            hl.agg.min(vp_mt.csq)
        )

    vp_mt = vp_mt.group_rows_by(
        'gene_id',
        'gene_symbol'
    ).aggregate(
        all=hl.agg.filter(
            vp_mt.x &
            hl.if_else(
                vp_mt.pop == 'all',
                hl.is_defined(vp_mt.popmax_af) &
                (vp_mt.popmax_af <= MAX_FREQ),
                vp_mt.pop_af[vp_mt.pop] <= MAX_FREQ
            ),
            get_grouped_phase_agg()
        ),
        **{
            f"af_le_{af}": hl.agg.filter(
                hl.if_else(
                    vp_mt.pop == 'all',
                    hl.is_defined(vp_mt.popmax_af) &
                    (vp_mt.popmax_af <= af),
                    vp_mt.pop_af[vp_mt.pop] <= af
                )
                & vp_mt.x,
                get_grouped_phase_agg()
            )
            for af in ALLELE_FREQUENCY_CUTOFFS
        }
    )
    vp_mt.describe()
    vp_mt = vp_mt.checkpoint('gs://gnomad-tmp/compound_hets/chet_per_gene{}.2.mt'.format(
        '.chr20' if chr20 else ''
    ), overwrite=True)

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
                        n_worst_chet=hl.agg.count_where(vp_mt[af].get(1) == csq_i),
                        n_chet=hl.agg.count_where((vp_mt[af].get(1) == csq_i) & (vp_mt[af].get(2, 9) >= csq_i) & (vp_mt[af].get(3, 9) >= csq_i)),
                        n_same_hap=hl.agg.count_where((vp_mt[af].get(2) == csq_i) & (vp_mt[af].get(1, 9) > csq_i) & (vp_mt[af].get(3, 9) >= csq_i)),
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
            for af in ['all'] + [f'af_le_{af}' for af in ALLELE_FREQUENCY_CUTOFFS]
        ])
    ).rows()
    gene_ht.describe()
    gene_ht = gene_ht.explode('row_counts')
    gene_ht = gene_ht.select(
        **gene_ht.row_counts
    )

    gene_ht.describe()
    gene_ht = gene_ht.checkpoint(
        'gs://gnomad-sarah/compound_hets/chet_per_gene{}.ht'.format(
            '.chr20' if chr20 else ''
        ),
        overwrite=overwrite
    )
    gene_ht.describe()

    gene_ht.flatten().export(
        'gs://gnomad-sarah/compound_hets/chet_per_gene{}.tsv.gz'.format(
            '.chr20' if chr20 else ''
        )
    )


def compute_from_full_mt(chr20: bool, overwrite: bool):
    mt = get_gnomad_data('exomes', adj=True, release_samples=True)
    freq_ht = public_release('exomes').ht().select('freq','popmax')
    vep_ht = public_release('exomes').ht().select('vep')
    rf_ht = hl.read_table(annotations_ht_path('exomes', 'rf'))

    if chr20:
        mt, freq_ht, vep_ht, rf_ht = filter_to_chr20([mt, freq_ht, vep_ht, rf_ht])
    
    vep_ht = vep_ht.annotate(
        vep=get_worst_gene_csq_code_expr(vep_ht.vep).values()
    )

    freq_ht = freq_ht.select(
        freq=freq_ht.freq[:10],
        popmax=freq_ht.popmax
    )
    freq_meta = hl.eval(freq_ht.globals.freq_meta)
    freq_dict = {f['pop']: i for i, f in enumerate(freq_meta[:10]) if 'pop' in f}
    freq_dict['all'] = 0
    freq_dict = hl.literal(freq_dict)
    mt = mt.annotate_rows(
        **freq_ht[mt.row_key],
        vep=vep_ht[mt.row_key].vep,
        filters=rf_ht[mt.row_key].filters
    )
    mt = mt.filter_rows(
        (mt.freq[0].AF <= MAX_FREQ) &
        (hl.len(mt.vep) > 0) &
        (hl.len(mt.filters) == 0)
    )

    mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = mt.select_entries(
        is_het=mt.GT.is_het()
    )

    mt = mt.explode_rows(mt.vep)
    mt = mt.transmute_rows(**mt.vep)

    mt = mt.annotate_cols(
        pop=['all', mt.meta.pop]
    )

    mt.describe()

    mt = mt.explode_cols(mt.pop)
    mt.show(n_cols=10)
    
    mt = mt.group_rows_by(
        'gene_id'
    ).aggregate_rows(
        gene_symbol=hl.agg.take(mt.gene_symbol, 1)[0]
	).aggregate(
        all=hl.agg.filter(
                hl.if_else(
                    mt.pop == 'all',
                    hl.is_defined(mt.popmax[0]) & (mt.popmax[0].AF <= MAX_FREQ),
                    mt.freq[freq_dict[mt.pop]].AF <= MAX_FREQ
                ),
                hl.agg.group_by(
                    hl.if_else(
                        mt.pop == 'all',
                        mt.popmax[0].AF > 0,
                        mt.freq[freq_dict[mt.pop]].AF > 0
                    ),
                    hl.struct(
                        hom_csq=hl.agg.filter(~mt.is_het, hl.agg.min(mt.csq)),
                        het_csq=hl.agg.filter(mt.is_het, hl.agg.min(mt.csq)),
                        het_het_csq=hl.sorted(
                            hl.array(
                                hl.agg.filter(mt.is_het, hl.agg.counter(mt.csq))
                            ),
                            key=lambda x: x[0]
                        ).scan(
                            lambda i, j: (j[0], i[1] + j[1]),
                            (0,0)
                        ).find(
                            lambda x: x[1] > 1
                        )[0]
                    )
                )
            ),
    **{f"af_le_{af}": hl.agg.filter(
                hl.if_else(
                    mt.pop == 'all',
                    hl.is_defined(mt.popmax[0]) & (mt.popmax[0].AF <= MAX_FREQ),
                    mt.freq[freq_dict[mt.pop]].AF <= MAX_FREQ
                ),
                hl.agg.group_by(
                    hl.if_else(
                        mt.pop == 'all',
                        mt.popmax[0].AF <= af,
                        mt.freq[freq_dict[mt.pop]].AF <= af
                    ),
                    hl.struct(
                        hom_csq=hl.agg.filter(~mt.is_het, hl.agg.min(mt.csq)),
                        het_csq=hl.agg.filter(mt.is_het, hl.agg.min(mt.csq)),
                        het_het_csq=hl.sorted(
                            hl.array(
                                hl.agg.filter(mt.is_het, hl.agg.counter(mt.csq))
                            ),
                            key=lambda x: x[0]
                        ).scan(
                            lambda i, j: (j[0], i[1] + j[1]),
                            (0,0)
                        ).find(
                            lambda x: x[1] > 1
                        )[0]
                    )
                )
            )
        for af in ALLELE_FREQUENCY_CUTOFFS
    }
    )
    
    mt = mt.annotate_entries(
            all=hl.struct(
                hom_csq=hl.min(mt.all.get(True).hom_csq, mt.all.get(False).hom_csq),
                het_csq=hl.min(mt.all.get(True).het_csq, mt.all.get(False).het_csq),
                het_het_csq=hl.min(mt.all.get(True).het_het_csq, mt.all.get(False).het_het_csq)
            ),
            **{f"af_le_{af}":hl.struct(
                hom_csq=mt[f"af_le_{af}"].get(True).hom_csq,
                het_csq=mt[f"af_le_{af}"].get(True).het_csq,
                het_het_csq=mt[f"af_le_{af}"].get(True).het_het_csq
            )
            for af in ALLELE_FREQUENCY_CUTOFFS
            }
    )

    gene_ht = mt.annotate_rows(
        row_counts=hl.flatten([
            hl.array(
                hl.agg.group_by(
                    mt.pop,
                    hl.struct(
                        csq=csq,
                        af=af,
                        n_hom=hl.agg.count_where(mt[af].hom_csq == csq_i),
                        n_het=hl.agg.count_where(mt[af].het_csq == csq_i),
                        n_het_het=hl.agg.count_where(mt[af].het_het_csq == csq_i)
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
            for af in ['all'] + [f'af_le_{af}' for af in ALLELE_FREQUENCY_CUTOFFS]
        ])
    ).rows()

    gene_ht = gene_ht.explode('row_counts')
    gene_ht = gene_ht.select(
        'gene_symbol',
        **gene_ht.row_counts
    )

    gene_ht = gene_ht.checkpoint(
        'gs://gnomad-sarah/compound_hets/het_and_hom_per_gene{}.ht'.format(
            '.chr20' if chr20 else ''
        ),
        overwrite=overwrite
    )

    gene_ht.flatten().export('gs://gnomad-sarah/compound_hets/het_and_hom_per_gene{}.tsv.gz'.format(
        '.chr20' if chr20 else ''
    ))


def main(args):
    hl.init(log="/tmp/hail.log")
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

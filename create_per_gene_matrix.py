from gnomad_hail import *
from resources import LEAST_CONSEQUENCE, MAX_FREQ
from hail import MatrixTable

from resources import *

CHET_THRESHOLD = 0.5
SAME_HAP_THRESHOLD = 0.5

def get_worst_gene_csq_expr(vep_expr: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    return vep_expr.transcript_consequences.map(
        lambda ts: hl.tuple((
            ts.gene_id,
            ts.gene_symbol,
            (
                hl.case(missing_false=True)
                .when(ts.lof == 'HC', 'lof')
                .when((ts.polyphen_prediction == 'probably damaging'), 'damaging_missense')
                .when(ts.consequence_terms.any(lambda x: x == 'missense_variant'), 'missense_variant')
                .when(ts.consequence_terms.all(lambda x: x == 'synonymous_variant'), 'synonymous_variant')
                .or_missing()
            )
        ))
    ).filter(lambda x: hl.is_defined(x[2]))


def compute_from_vp_mt(overwrite: bool):

    def get_entry_agg(vp_mt: hl.MatrixTable, csqs: hl.expr.SetExpression, max_af: float) -> hl.expr.Expression:
        def get_popmax_af(freq_dict):
            return hl.max(freq_dict.keys().filter(lambda x: x != 'oth').map(lambda x: freq_dict[x].AF))

        return hl.agg.filter(
            csqs.contains(vp_mt.csq1) &
            csqs.contains(vp_mt.csq2) &
            (get_popmax_af(vp_mt.freq1) <= max_af) &
            (get_popmax_af(vp_mt.freq2) <= max_af) &
            vp_mt.is_het_het,
            hl.struct(
                n=hl.agg.count(),
                max_p_chet=hl.agg.max(vp_mt.phase_info[vp_mt.pop].em.adj.p_chet)
            )
        )

    def get_row_agg(entry_exp: hl.expr.StructExpression) -> hl.expr.StructExpression:
        return hl.struct(
            total=hl.agg.count_where(entry_exp.n > 0),
            phased_c_het=hl.agg.count_where(entry_exp.max_p_chet > CHET_THRESHOLD),
            phased_same_hap=hl.agg.count_where(entry_exp.max_p_chet < SAME_HAP_THRESHOLD)
        )

    vp_mt = hl.read_matrix_table(full_mt_path('exomes'))
    ann_ht = hl.read_table(vp_ann_ht_path('exomes'))
    ann_ht = ann_ht.select(
        'freq1',
        'freq2',
        filter1=hl.len(ann_ht.filters1) > 0,
        filter2=hl.len(ann_ht.filters2) > 0,
        vep1=get_worst_gene_csq_expr(ann_ht.vep1),
        vep2=get_worst_gene_csq_expr(ann_ht.vep2)
    )
    ann_ht = ann_ht.annotate(
        genes=hl.set(ann_ht.vep1.map(lambda x: x[0])).intersection(
            hl.set(ann_ht.vep2.map(lambda x: x[0]))
        )
    )

    phase_ht = hl.read_table(phased_vp_count_ht_path('exomes'))
    meta = get_gnomad_meta('exomes')
    vp_mt = vp_mt.select_entries(
        is_het_het=vp_mt.gt1.is_het() & vp_mt.gt2.is_het() & vp_mt.adj1 & vp_mt.adj2
    )
    vp_mt = vp_mt.annotate_cols(
        pop=meta[vp_mt.col_key].pop
    )
    vp_mt: MatrixTable = vp_mt.annotate_rows(
        **ann_ht[vp_mt.row_key],
        phase_info=phase_ht[vp_mt.row_key].phase_info
    )
    vp_mt = vp_mt.filter_rows(
        vp_mt.filter1 | vp_mt.filter2
    )

    vp_mt = vp_mt.explode_rows('genes')
    vp_mt = vp_mt.transmute_rows(
        gene_id=vp_mt.genes,
        gene_symbol=vp_mt.vep1.find(lambda x: x[0]==vp_mt.genes)[1],
        csq1=vp_mt.vep1.find(lambda x: x[0]==vp_mt.genes)[2],
        csq2=vp_mt.vep2.find(lambda x: x[0]==vp_mt.genes)[2]
    )

    vp_mt = vp_mt.group_rows_by(
        'gene_id'
    ).aggregate(
        _gene_symbol=hl.agg.take(vp_mt.gene_symbol, 1)[0],
        _lof_counts_0_01=get_entry_agg(vp_mt, hl.literal({'lof'}), 0.01),
        _lof_counts_0_001=get_entry_agg(vp_mt, hl.literal({'lof'}), 0.001),
        _damaging_missense_counts_0_01=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense'}), 0.01),
        _damaging_missense_counts_0_001=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense'}), 0.001),
        _missense_counts_0_01=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense', 'missense_variant'}), 0.01),
        _missense_counts_0_001=get_entry_agg(vp_mt, hl.literal({'lof', 'damaging_missense', 'missense_variant'}), 0.001),
        _synonymous_counts_0_01=get_entry_agg(vp_mt, hl.literal({'synonymous_variant'}), 0.01),
        _synonymous_counts_0_001=get_entry_agg(vp_mt, hl.literal({'synonymous_variant'}), 0.001)
    ).checkpoint('gs://gnomad-tmp/vp_gene.mt', overwrite=overwrite)

    gene_ht = vp_mt.annotate_rows(
        gene_symbol=hl.agg.take(vp_mt._gene_symbol, 1)[0],
        **{x[1:]: get_row_agg(vp_mt[x]) for x in vp_mt.entry if x != '_gene_symbol'}
    ).rows()

    gene_ht.flatten().export('gs://gnomad-lfran/compound_hets/chet_per_gene.tsv.gz')

def compute_from_full_mt(overwrite: bool):

    def get_entry_agg(mt: hl.MatrixTable, csqs: hl.expr.SetExpression, max_af: float) -> hl.expr.Expression:
        return hl.agg.filter(
            csqs.contains(mt.csq) &
            (mt.popmax.AF <= max_af),
            hl.struct(
                n_hom=hl.agg.count_where(mt.GT.is_hom_var()),
                n_het=hl.agg.count_where(mt.GT.is_het())
            )
        )

    def get_row_agg(entry_exp: hl.expr.StructExpression) -> hl.expr.StructExpression:
        return hl.struct(
            het=hl.agg.count_where(entry_exp.n_het == 1),
            het_het=hl.agg.count_where(entry_exp.n_het > 1),
            hom=hl.agg.count_where(entry_exp.n_hom > 0)
        )

    mt = get_gnomad_data('exomes', adj=True, release_samples=True)
    freq_ht = hl.read_table(annotations_ht_path('exomes', 'frequencies'))
    vep_ht = hl.read_table(annotations_ht_path('exomes', 'vep'))
    rf_ht = hl.read_table(annotations_ht_path('exomes', 'rf'))

    vep_ht = vep_ht.annotate(
        vep=get_worst_gene_csq_expr(
            vep_ht.vep.annotate(
                transcript_consequences=vep_ht.vep.transcript_consequences.filter(
                    lambda tc: (tc.biotype == 'protein_coding')
                )
            )
        )
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
    mt = mt.transmute_rows(
        gene_id=mt.vep[0],
        gene_symbol=mt.vep[1],
        csq=mt.vep[2]
    )

    mt = mt.group_rows_by(
        'gene_id'
    ).aggregate(
        _gene_symbol=hl.agg.take(mt.gene_symbol, 1)[0],
        _lof_counts_0_01=get_entry_agg(mt, hl.literal({'lof'}), 0.01),
        _lof_counts_0_001=get_entry_agg(mt, hl.literal({'lof'}), 0.001),
        _damaging_missense_counts_0_01=get_entry_agg(mt, hl.literal({'lof', 'damaging_missense'}), 0.01),
        _damaging_missense_counts_0_001=get_entry_agg(mt, hl.literal({'lof', 'damaging_missense'}), 0.001),
        _missense_counts_0_01=get_entry_agg(mt, hl.literal({'lof', 'damaging_missense', 'missense_variant'}), 0.01),
        _missense_counts_0_001=get_entry_agg(mt, hl.literal({'lof', 'damaging_missense', 'missense_variant'}), 0.001),
        _synonymous_counts_0_01=get_entry_agg(mt, hl.literal({'synonymous_variant'}), 0.01),
        _synonymous_counts_0_001=get_entry_agg(mt, hl.literal({'synonymous_variant'}), 0.001)
    ).checkpoint('gs://gnomad-tmp/gene.mt', overwrite=overwrite)

    gene_ht = mt.annotate_rows(
            gene_symbol=hl.agg.take(mt._gene_symbol, 1)[0],
            **{x[1:]: get_row_agg(mt[x]) for x in mt.entry if x != '_gene_symbol'}
        ).rows()

    gene_ht.flatten().export('gs://gnomad-lfran/compound_hets/het_and_hom_per_gene.tsv.gz')


def main(args):
    if args.compute_from_vp_mt:
        compute_from_vp_mt(args.overwrite)
    if args.compute_from_full_mt:
        compute_from_full_mt(args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compute_from_vp_mt', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--compute_from_full_mt', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    main(args)
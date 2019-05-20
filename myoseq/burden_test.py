from gnomad_hail import *

print("WARN: This script requires highmem machines!")

mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.mt')
meta = hl.read_table('gs://gnomad/projects/compound_hets/myoseq/sample_qc/MacArthur_LGMD_Callset_Jan2019.full_meta.ht')
pop_distance = hl.read_table('gs://gnomad-lfran/compound_hets/myoseq/sample_qc/myoseq_pop_distance_to_max_kde.ht')
variant_annotations_ht = hl.read_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.annotations.ht')
variant_annotations_ht.drop('was_split', 'a_index')
mt = mt.annotate_cols(
    **meta[mt.col_key],
    **pop_distance[mt.col_key],
)
mt = mt.annotate_rows(
    **variant_annotations_ht[mt.row_key]
)

# Filter samples failing QC
mt = mt.filter_cols((hl.len(mt.sample_filters) == 0))
counts = mt.aggregate_cols(hl.agg.counter(mt.is_case))
print(f'Found {counts[True]} cases and {counts[False]} controls for gene aggregation.')


# Filter sites failing QC, without any tx_annotation (i.e. without a protein-coding variant) or too common
mt = mt.filter_rows(
    (hl.len(mt.filters)==0) &
    hl.is_defined(mt.tx_annotation) &
    (hl.or_else(mt.gnomad_exomes_popmax.AF, hl.or_else(mt.gnomad_genomes_popmax.AF, 0.0)) < 0.05)
)

# Keep non-ref entries only
mt = mt.filter_entries(mt.GT.is_non_ref())

# Annotate genes and popmax
mt = mt.annotate_rows(
    gene=hl.set(mt.tx_annotation.map(
        lambda x: hl.struct(gene_symbol=x.symbol, gene_id=x.ensg)
    )),
    gnomad_popmax_af=hl.or_else(mt.gnomad_exomes_popmax.AF, mt.gnomad_genomes_popmax.AF)
)

# Aggregate by gene
mt = mt.explode_rows(mt.gene)
mt.write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_filtered_gene_exploded.mt', overwrite=True)

mt = hl.read_matrix_table('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_filtered_gene_exploded.mt')
mt = mt.group_rows_by(**mt.gene).aggregate(
    variants=hl.agg.collect(
        hl.struct(
            n_alt_alleles=mt.GT.n_alt_alleles(),
            adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD, haploid_adj_dp=5),
            gnomad_popmax_af=mt.gnomad_popmax_af,
            tx_annotation=mt.tx_annotation.filter(
                lambda x: (x.symbol == mt.gene.gene_symbol) & (x.ensg == mt.gene.gene_id)
            ).map(
                lambda x: x.select(
                    'csq',
                    'lof',
                    'lof_flag',
                    'Muscle_Skeletal',
                    'WholeBlood',
                    'mean_proportion',
                    'polyphen_prediction',
                    'sift_prediction'
                )
            )
        )
    )
).write('gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019_gene_by_sample.mt')

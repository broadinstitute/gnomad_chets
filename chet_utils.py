import hail as hl
from gnomad.utils.vep import CSQ_ORDER
from gnomad_qc.v2.resources import get_gnomad_meta, fam_path
from logging import getLogger

logger = getLogger('chet_utils')
BAD_THAI_TRIOS_PROJECT_ID = 'C978'


def vep_genes_expr(vep_expr: hl.expr.StructExpression, least_consequence: str) -> hl.expr.SetExpression:
    vep_consequences = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1]))
    return (
                hl.set(
                    vep_expr.transcript_consequences
                        .filter(
                        lambda tc: (tc.biotype == 'protein_coding') &
                                   (tc.consequence_terms.any(lambda c: vep_consequences.contains(c)))
                    )
                        .map(lambda x: x.gene_id)
                )
            )


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
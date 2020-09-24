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


def get_phased_gnomad_ht(
        ht: hl.Table,
        em: bool = True,
        lr: bool = True,
        shr: bool = True
) -> hl.Table:
    expr_fun = []

    if em:
        expr_fun.append(get_em_expressions)

    if lr:
        expr_fun.append(get_lr_expressions)

    if shr:
        expr_fun.append(get_single_het_expressions)

    if not expr_fun:
        raise (Exception("No expressions to annotate"))

    # Support for both exploded or dict versions of gt_counts
    # dict
    if isinstance(ht.gt_counts, hl.expr.DictExpression):
        ht = ht.select(
            phase_info=ht.gt_counts.map_values(
                lambda pop_count: hl.bind(
                    lambda x: hl.struct(
                        gt_counts=x,
                        **{
                            k: v for f in expr_fun for k, v in f(x).items()
                        }
                    ),
                    hl.struct(
                        raw=pop_count.raw.map(lambda y: hl.int32(y)),
                        adj=pop_count.adj.map(lambda z: hl.int32(z))
                    )
                )
            )
        )
    # exploded
    else:
        ht = ht.annotate(
            **{
                k: v for f in expr_fun for k, v in f(ht.gt_counts).items()
            }
        )

    return ht
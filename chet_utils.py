import hail as hl
from gnomad.utils.vep import CSQ_ORDER


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
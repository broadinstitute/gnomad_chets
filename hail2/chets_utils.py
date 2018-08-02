import hail as hl
from gnomad_hail.resources import annotations_ht_path
from typing import List


def add_minimal_vep(mt: hl.MatrixTable, data_type: str, consequences_to_keep: List[str]):
    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))

    mt = mt.annotate_globals(
        vep_consequences_to_keep=consequences_to_keep
        # polyphen_pred=POLYPHEN_PREDICTIONS,
        # sift_pred=SIFT_PREDICTIONS
    )

    mt = mt.annotate_rows(
        vep=(
            vep_ht[mt.row_key].vep.transcript_consequences
                .filter(
                lambda tc: (tc.biotype == 'protein_coding') &
                           (tc.allele_num == mt.a_index) &
                           ((hl.is_defined(tc.lof) & (tc.lof == "HC")) | tc.consequence_terms.any(lambda c: consequences_to_keep.contains(c)))
            ).map(
                lambda x: x.select(
                    'transcript_id',
                    lof=hl.or_missing(x.lof == "HC", hl.struct()),
                    # lof=hl.cond(hl.is_defined(x.lof), hl.struct(), hl.null(hl.tstruct())),
                    # polyphen_prediction=polyphen_pred_index[x.polyphen_prediction],
                    # sift_prediction=sift_pred_index[x.sift_prediction],
                    canonical=hl.or_missing(x.canonical == 1, hl.struct())
                )
            )
        )
    )

    return mt


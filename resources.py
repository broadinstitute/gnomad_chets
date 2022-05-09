from gnomad_qc.v2.resources import get_gnomad_data, get_gnomad_meta, pbt_phased_trios_mt_path
from gnomad_chets.chet_utils import extract_pbt_probands
import hail as hl

LEAST_CONSEQUENCE = '3_prime_UTR_variant'
MAX_FREQ = 0.05


def mini_mt_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'mt', 'mini_mt', pbt, least_consequence, max_freq, chrom)


def vp_list_ht_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'list', pbt, least_consequence, max_freq, chrom)


def vp_ann_ht_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'ann', pbt, least_consequence, max_freq, chrom)


def full_mt_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'mt', '', pbt, least_consequence, max_freq, chrom)


def vp_count_ht_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'counts', pbt, least_consequence, max_freq, chrom)


def phased_vp_count_ht_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None, release: bool = False):
    if release:
        return f"gs://gnomad/release/2.1.1/ht/{data_type}_phased_counts_{max_freq}_{least_consequence}_vp{f'_chrom{chrom}' if chrom else ''}.ht"
    else:
        return _chets_out_path(data_type, 'ht', 'phased_counts', pbt, least_consequence, max_freq, chrom)

def pbt_phase_count_ht_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'ht', 'pbt_phase_count', False, least_consequence, max_freq, chrom)


def pbt_trio_mt_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'mt', 'pbt_trio', False, least_consequence, max_freq, chrom)


def pbt_trio_et_path(data_type: str, pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'ht', 'pbt_trio', False, least_consequence, max_freq, chrom)


def pbt_comparison_full_mt_path(data_type: str, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'mt', 'pbt_comparison', False, least_consequence, max_freq, chrom)


# def pbt_comparison_vp_ann_ht_path(data_type: str, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
#     return _chets_out_path(data_type, 'ht', 'pbt_comparison_ann', False, least_consequence, max_freq, chrom)
#

def pbt_comparison_vp_count_ht_path(data_type: str, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'pbt_comparison_counts', False, least_consequence, max_freq, chrom)


def pbt_comparison_phased_vp_count_ht_path(data_type: str, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'pbt_comparison_phased_counts', False, least_consequence, max_freq, chrom)


def get_adj_missing_mt(data_type: str, pbt: bool) -> hl.MatrixTable:
    mt = get_gnomad_data(data_type).select_cols() if not pbt else hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))
    mt = mt.select_rows()
    mt = mt.select_entries(
        GT=hl.or_missing(mt.GT.is_non_ref(), mt.GT),
        missing=hl.is_missing(mt.GT),
        adj=mt.adj
    ).select_cols().select_rows()

    if pbt:
        mt = mt.key_cols_by('s', trio_id=mt.source_trio.id)
        mt = extract_pbt_probands(mt, data_type)
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        mt = mt.key_cols_by(s=mt.s, trio_id=mt.source_trio.id)
    else:
        meta = get_gnomad_meta('exomes')
        mt = mt.filter_cols(meta[mt.col_key].high_quality)

    return mt


def _chets_out_path(data_type: str, extension: str, stage: str = '', pbt: bool = False, least_consequence: str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return 'gs://gnomad{}/compound_hets/{}{}{}_{}_{}_vp{}.{}'.format(
        '-tmp/' if stage == 'mini_mt' else '/projects',
        data_type,
        '_pbt' if pbt else '',
        f'_{stage}' if stage else '',
        max_freq,
        least_consequence,
        f'_chrom{chrom}' if chrom else '',
        extension)

LEAST_CONSEQUENCE = '3_prime_UTR_variant'
MAX_FREQ = 0.02


def mini_mt_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'mini_mt', pbt, least_consequence, max_freq, chrom)


def vp_list_mt_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'list', pbt, least_consequence, max_freq, chrom)


def full_mt_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, '', pbt, least_consequence, max_freq, chrom)


def vp_count_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'counts', pbt, least_consequence, max_freq, chrom)

def pbt_phase_count_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'pbt_phase_count', False, least_consequence, max_freq, chrom)


def gnomad_adj_missing_path(data_type: str):
    return f'gs://gnomad/projects/compound_hets/gnomad_{data_type}_adj_missing.mt'

def _chets_out_path(data_type: str, stage: str = '', pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return 'gs://gnomad{}/compound_hets/{}{}{}_{}_{}_vp{}.{}'.format(
        '-tmp/' if stage == 'mini_mt' else '/projects',
        data_type,
        '_pbt' if pbt else '',
        f'_{stage}' if stage else '',
        max_freq,
        least_consequence,
        f'_chrom{chrom}' if chrom else '',
        'ht' if 'count' in stage else 'mt')
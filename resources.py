LEAST_CONSEQUENCE = '3_prime_UTR_variant'
MAX_FREQ = 0.02


def mini_mt_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'mt', 'mini_mt', pbt, least_consequence, max_freq, chrom)


def vp_list_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'list', pbt, least_consequence, max_freq, chrom)


def vp_ann_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'ann', pbt, least_consequence, max_freq, chrom)

def full_mt_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'mt', '', pbt, least_consequence, max_freq, chrom)


def vp_count_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'counts', pbt, least_consequence, max_freq, chrom)


def phased_vp_count_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return _chets_out_path(data_type, 'ht', 'phased_counts', pbt, least_consequence, max_freq, chrom)


def pbt_phase_count_ht_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'ht', 'pbt_phase_count', False, least_consequence, max_freq, chrom)


def pbt_trio_mt_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'mt', 'pbt_trio', False, least_consequence, max_freq, chrom)


def pbt_trio_et_path(data_type: str, pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    # Keeping pbt arg just so the signature mimics others
    return _chets_out_path(data_type, 'ht', 'pbt_trio', False, least_consequence, max_freq, chrom)


def adj_missing_mt_path(data_type: str, pbt: bool) -> str:
    return 'gs://gnomad/projects/compound_hets/gnomad_{}_adj_missing{}.mt'.format(
        data_type,
        ".trios_pbt_phased" if pbt else ""
    )


def _chets_out_path(data_type: str, extension: str, stage: str = '', pbt: bool = False, least_consequence : str = LEAST_CONSEQUENCE, max_freq: float = MAX_FREQ, chrom: str = None):
    return 'gs://gnomad{}/compound_hets/{}{}{}_{}_{}_vp{}.{}'.format(
        '-tmp/' if stage == 'mini_mt' else '/projects',
        data_type,
        '_pbt' if pbt else '',
        f'_{stage}' if stage else '',
        max_freq,
        least_consequence,
        f'_chrom{chrom}' if chrom else '',
        extension)
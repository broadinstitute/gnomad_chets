import hail as hl
import logging
import timeit
import argparse

from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_LEAST_CONSEQUENCE,
    DEFAULT_MAX_FREQ,
    DEFAULT_TMP_DIR,
    TEST_INTERVAL,
    get_phasing_resources,
)
from gnomad_chets.v4.utils import filter_for_testing

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("phase_gnomad")
logger.setLevel(logging.INFO)


def get_em_expr(gt_counts):
    gt_counts_int32 = gt_counts.map(lambda x: hl.int32(x))
    hap_counts = hl.experimental.haplotype_freq_em(gt_counts_int32)
    return hl.bind(
        lambda x: hl.struct(
            hap_counts=x,
            p_chet=(x[1] * x[2]) / (x[0] * x[3] + x[1] * x[2])
        ),
        hap_counts
    )


def get_phased_gnomad_ht(
        ht: hl.Table
) -> hl.Table:
    print(ht.describe())
    
    return dict(
        em=hl.struct(
            raw=get_em_expr(ht.gt_counts_raw),
            adj=get_em_expr(ht.gt_counts_adj),
        ),
        em_plus_one=hl.struct(
            raw=get_em_expr(ht.gt_counts_raw + [0, 0, 0, 0, 1, 0, 0, 0, 0]),
            adj=get_em_expr(ht.gt_counts_adj + [0, 0, 0, 0, 1, 0, 0, 0, 0]),
        )
    )

def main(args):
    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix or ""
    data_type = args.data_type
    test = args.test
    
    hl.init(
        log="/create_vp_matrix.log",
        tmp_dir=tmp_dir,
    )
    
    logger.info(
        f"""
        Running script with the following parameters:

            Data type: {data_type}
            Test: {test}
            Output postfix: {output_postfix}
            Overwrite: {overwrite}
            Tmp dir: {tmp_dir}
        """
    )
    
    resources = get_phasing_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir,
        output_postfix=output_postfix,
        overwrite=overwrite,
    )
    
    if args.phase:
        logger.info("Phasing variant pairs...")
        res=resources.phase
        
        #phase variant pairs
        ht=hl.read_table(args.file_to_phase)
        phased_dict=get_phased_gnomad_ht(ht)
        
        ht = ht.annotate(**dict(phased_dict))
                
        #write phased data
        ht=ht.write(res.phase.path,overwrite=overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tmp-dir",
        default=DEFAULT_TMP_DIR,
        help="Temporary directory for intermediate files.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to PCNT gene (chr21:46324141-46445769) for testing purposes.",
    )
    parser.add_argument(
        "--output-postfix",
        help=(
            'Postfix to append to output file names (e.g., "pcnt_test" for files like '
            "exomes.vp_list.pcnt_test.ht)."
        ),
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Whether to overwrite existing files."
    )
    parser.add_argument(
        "--data-type",
        default=DEFAULT_DATA_TYPE,
        choices=DATA_TYPE_CHOICES,
        help=(
            f'Data type to use. Must be one of {", ".join(DATA_TYPE_CHOICES)}. Default '
            f"is {DEFAULT_DATA_TYPE}.",
        ),
    )
    parser.add_argument(
        "--phase",
        action="store_true",
        help="Whether to phase variant pairs.",
    )
    
    parser.add_argument(
        "--file-to-phase",
        help="input file for phasing",
    )
        
    
    
    args = parser.parse_args()
    main(args)

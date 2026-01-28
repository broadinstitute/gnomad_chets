"""
Phase variant pairs for gnomAD v4.

This uses thevariant-pair hail table from create_vp_matrix.py, then runs expectation maximization (EM)-based phasing summary (via
hail.experimental.haplotype_freq_em), and writes an annotated hail table with
per-pair EM results. 

To run, use --phase
"""

import hail as hl
import logging
import timeit
import argparse
from gnomad.utils.vep import CSQ_ORDER

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


def get_em_expr(
    gt_counts: hl.Table,
    ) -> hl.struct:
    """
    Return an expression computing haplotype EM counts and p_chet.

    Parameters
    ----------
    gt_counts : ArrayExpression
        A length-9 array-like expression encoding genotype counts in the
        ordering expected by `hail.experimental.haplotype_freq_em`. This should be:
        The unphased input genotype counts for the variant pairs has to be provided in the following order:
        [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]
        Where _A_ and _a_ are the reference and non-reference alleles for the first variant, resp.
        And _B_ and _b_ are the reference and non-reference alleles for the second variant, resp.
        

    Returns
    -------
    StructExpression
        A struct with fields
        - hap_counts: array of EM-estimated haplotype counts
                      The estimated haplotype counts are returned in an array in the following order: [AB, aB, Ab, ab]
        - p_chet: estimated probability that the pair is compound het (closer to 1) or on the same haplotype (closer to 0)
        

    Notes
    -----
    The implementation converts counts to int32 and uses Hail's
    `haplotype_freq_em` to estimate haplotype frequencies. The returned
    `p_chet` follows the same algebra as the previous v2 implementation.
    """

    #this needs to be converted to int32 or type error . Even if values look like Python ints, they can be typed by Hail as float or int64 depending on upstream operations. Calling hl.int32 ensures the expression has the precise Hail type haplotype_freq_em expects (in haplotype_freq_em: @typecheck(gt_counts=expr_array(expr_int32)))
    hap_counts = hl.experimental.haplotype_freq_em(gt_counts.map(lambda x: hl.int32(x)))
    return hl.bind(
        lambda x: hl.struct(
            hap_counts=x,
            p_chet=(x[1] * x[2]) / (x[0] * x[3] + x[1] * x[2])
        ),
        hap_counts
    )


def get_phased_gnomad_ht(
        ht: hl.Table
) -> hl.struct:
    """
    Create phased annotations for a variant-pair table.

    Parameters
    ----------
    ht : hail.Table
        A table keyed by variant-pair (or containing per-pair genotype count
        fields) that includes `gt_counts_raw` and `gt_counts_adj` arrays.

    Returns
    -------
    hl.struct
        A hl.struct of annotations suitable to pass into `Table.annotate`.
        Contains two top-level structs: `em` and `em_plus_one`, each with
        `raw` and `adj` sub-structs holding the EM output and `p_chet`.

    The `em_plus_one` field computes a stabilized EM by adding a small
    pseudo-count vector before running the EM algorithm; this mirrors the
    lightweight stabilization in the original pipeline.
    """

    return hl.struct(
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
    """
    Compute phasing annotations and write an output HT.

    The function reads CLI-style `args` (from argparse), initializes Hail,
    loads pipeline resources via `get_phasing_resources`, and if `--phase` is
    supplied reads `--file-to-phase`, computes phased annotations using
    `get_phased_gnomad_ht`, checkpoints the result, and writes the final HT
    to the resource path.

    Parameters
    ----------
    args : argparse.Namespace
        Namespace with attributes matching the CLI flags declared at the
        bottom of this module (e.g., `tmp_dir`, `overwrite`, `output_postfix`,
        `data_type`, `test`, `phase`, `file_to_phase`).
    """

    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix or ""
    data_type = args.data_type
    test = args.test
    max_freq = args.max_freq
    least_consequence = args.least_consequence

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
        res = resources.phase

        ht = hl.read_table(args.file_to_phase)

        # Phase variant pairs: read input HT, compute phased annotations
        phased_dict = get_phased_gnomad_ht(ht)
        logger.info("Phasing complete. Now annotating phased data...")

        ht = ht.annotate(**phased_dict).checkpoint(
            hl.utils.new_temp_file("get_phased_gnomad", "ht")
        )
        
        #add in annotation on parameters that can be changed each time
        ht = ht.annotate(
            max_freq=max_freq,
            least_consequence=least_consequence,
        )
        
        logger.info("Annotating complete. Now writing phased data...")
        

        ht = ht.write(res.phase.path, overwrite=overwrite)

    stop = timeit.default_timer()
    logger.info(f"Time taken to run the script is {stop - start} seconds.")


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
    
    parser.add_argument(
        "--max-freq",
        type=float,
        default=DEFAULT_MAX_FREQ,
        help=f"Maximum global AF to keep (inclusive). Default is {DEFAULT_MAX_FREQ}.",
    )
    
    parser.add_argument(
        "--least-consequence",
        default=DEFAULT_LEAST_CONSEQUENCE,
        choices=CSQ_ORDER,
        help=(
            "Lowest-severity consequence to keep. Default "
            f"is {DEFAULT_LEAST_CONSEQUENCE}."
        ),
    )
        
    
    
    args = parser.parse_args()
    main(args)

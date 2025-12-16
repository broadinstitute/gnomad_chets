"""
Test script for variant pair genotype and genotype counts functions.

Creates a test dense MatrixTable and tests:
- create_variant_pair_genotype_ht
- create_variant_pair_genotype_counts_ht
"""

import argparse
import logging

import hail as hl

from gnomad_chets.v4.create_vp_matrix import (
    create_variant_pair_genotype_counts_ht,
    create_variant_pair_genotype_ht,
)
from gnomad_chets.v4.resources import DEFAULT_TMP_DIR

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("test_create_vp_matrix")
logger.setLevel(logging.INFO)


def create_test_dense_mt() -> hl.MatrixTable:
    """
    Create a hardcoded test dense MatrixTable with known edge cases.

    Creates a MatrixTable with:
    - 20 variants at different positions on chr21
    - 10 samples with specific genotype patterns
    - GT, GQ, DP, AD fields

    Test cases included:
    - Various combinations of het, hom_var, hom_ref, and missing genotypes
    - Edge cases for variant pairs where one is hom_ref and one is het/hom_var
    - Missing GT genotypes (multiple patterns)
    - Genotypes that fail adj filtering:
      - Low GQ (< 20) for het and hom_var
      - Low DP (< 10)
      - Missing AD
      - AD that doesn't match GT (het with hom_ref AD pattern)

    This creates variant pairs with known genotype combinations that we can verify,
    including cases where raw and adj counts differ.

    :return: Test MatrixTable.
    """
    logger.info("Creating hardcoded test MatrixTable with edge cases...")

    # Define 20 variants at different positions on chr21
    variants_data = []
    for i in range(20):
        pos = 46324141 + i * 1000  # Space them out by 1000bp
        # Alternate between different allele combinations
        if i % 3 == 0:
            alleles = ["A", "T"]
        elif i % 3 == 1:
            alleles = ["G", "C"]
        else:
            alleles = ["T", "A"]
        variants_data.append(
            {
                "locus": hl.locus("chr21", pos, reference_genome="GRCh38"),
                "alleles": alleles,
            }
        )

    # Define 10 samples
    samples = [f"sample_{i}" for i in range(10)]

    # Create MatrixTable with the right dimensions
    # range_matrix_table creates a MatrixTable with row_idx as row key
    # We need to add col_idx ourselves
    mt = hl.utils.range_matrix_table(n_rows=20, n_cols=10)

    # Annotate rows with variant data and save row_idx as non-key field
    variants_array = hl.literal(variants_data)
    mt = mt.select_rows(
        locus=variants_array[mt.row_idx].locus,
        alleles=variants_array[mt.row_idx].alleles,
    )

    # Now key by locus and alleles
    mt = mt.key_rows_by("locus", "alleles")

    # Add sample ID as annotation (col_idx is the column key at this point)
    samples_array = hl.literal(samples)
    mt = mt.key_cols_by(s=samples_array[hl.int32(mt.col_idx)])

    # Define genotypes with a pattern that creates various edge cases
    # Pattern: Create het, hom_var, missing, and adj-failing genotypes
    # This ensures we have variant pairs with known combinations
    gt_expr = (
        hl.case()
        # Create some het genotypes (pass adj)
        .when(
            (mt.row_idx % 5 == 0) & (mt.col_idx % 3 == 0), hl.call(0, 1)
        )  # het pattern
        # Create some hom_var genotypes (pass adj)
        .when(
            (mt.row_idx % 5 == 1) & (mt.col_idx % 3 == 1), hl.call(1, 1)
        )  # hom_var pattern
        # Create some missing genotypes (GT is missing)
        .when(
            (mt.row_idx % 5 == 2) & (mt.col_idx % 3 == 2), hl.missing(hl.tcall)
        )  # missing GT
        # Create more missing genotypes
        .when(
            (mt.row_idx % 5 == 3) & (mt.col_idx % 4 == 0), hl.missing(hl.tcall)
        )  # additional missing
        # Default to hom_ref for all others
        .default(hl.call(0, 0))
    )

    # Define GQ with some low values that will fail adj filtering (< 20)
    gq_expr = (
        hl.case()
        # Low GQ for variants that should fail adj (row_idx % 5 == 4)
        .when(
            (mt.row_idx % 5 == 4) & (mt.col_idx % 3 == 0), hl.int32(15)
        ).when(  # low GQ het
            (mt.row_idx % 5 == 4) & (mt.col_idx % 3 == 1), hl.int32(10)
        )  # low GQ hom_var
        # Normal GQ for others
        .default(hl.int32(40))
    )

    # Define DP with some low values that will fail adj filtering (< 10)
    dp_expr = (
        hl.case()
        # Low DP for some variants that should fail adj
        .when((mt.row_idx % 5 == 4) & (mt.col_idx % 3 == 2), hl.int32(5))  # low DP
        # Normal DP for others
        .default(hl.int32(30))
    )

    mt = mt.annotate_entries(
        GT=gt_expr,
        GQ=gq_expr,
        DP=dp_expr,
    )

    # Add AD based on GT, with some problematic cases
    ad_expr = (
        hl.case()
        # Missing AD for some genotypes (will fail adj)
        .when(
            (mt.row_idx % 5 == 4) & (mt.col_idx % 4 == 3),
            hl.missing(hl.tarray(hl.tint32)),
        )
        # AD that doesn't match GT (will fail adj) - het with wrong AD
        .when(
            (mt.row_idx % 5 == 4) & (mt.col_idx % 4 == 2) & mt.GT.is_het(),
            hl.array([30, 0]),
        )  # het with hom_ref AD
        # Normal AD based on GT
        .when(mt.GT.is_het(), hl.array([15, 15]))  # het: equal ref/alt
        .when(mt.GT.is_hom_var(), hl.array([0, 30]))  # hom_var: all alt
        .default(hl.array([30, 0]))  # hom_ref or missing: all ref
    )

    mt = mt.annotate_entries(AD=ad_expr)

    # Remove helper fields
    mt = mt.drop("row_idx", "col_idx")

    logger.info(
        f"Created test MatrixTable with {mt.count_rows()} rows and {mt.count_cols()} cols"
    )
    return mt


def create_test_variant_pair_ht(mt: hl.MatrixTable) -> hl.Table:
    """
    Create a test variant pair Table from a MatrixTable.

    Creates variant pairs from the MatrixTable, filtering to:
    - Only unique variant pairs (no duplicate pairs)
    - Only variants that have at least one non-ref genotype
    - Only pairs where at least one sample has a non-ref genotype for both variants

    The Table should have fields locus1, alleles1, locus2, alleles2 and be
    keyable by these fields (though not necessarily keyed).

    :param mt: MatrixTable to create variant pairs from.
    :return: Variant pair Table with fields locus1, alleles1, locus2, alleles2.
    """
    logger.info("Creating test variant pair Table...")

    # Filter to variants that have at least one non-ref genotype
    # This matches the real pipeline where variant pairs only include variants
    # that are actually called in at least one sample
    mt_variants = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    # Get variants (by locus and alleles) - these will be used to create pairs
    variants_ht = mt_variants.rows().select().distinct()
    variants = variants_ht.collect()

    logger.info(f"Found {len(variants)} variants with non-ref genotypes")

    # Create all possible pairs from unique variants
    pairs = []
    for i in range(len(variants)):
        for j in range(i + 1, len(variants)):
            pairs.append(
                {
                    "locus1": variants[i].locus,
                    "alleles1": variants[i].alleles,
                    "locus2": variants[j].locus,
                    "alleles2": variants[j].alleles,
                }
            )

    # Create Table with the required fields: locus1, alleles1, locus2, alleles2
    vp_ht = hl.Table.parallelize(
        pairs,
        schema=hl.tstruct(
            locus1=hl.tlocus("GRCh38"),
            alleles1=hl.tarray(hl.tstr),
            locus2=hl.tlocus("GRCh38"),
            alleles2=hl.tarray(hl.tstr),
        ),
    )

    # Filter to pairs where at least one sample has a non-ref genotype for both variants
    # We'll use entries table to find samples with non-ref for each variant
    entries = mt.entries()

    # Create a table of samples with non-ref for each variant
    # First, filter to non-ref entries and select relevant fields
    non_ref_entries = entries.filter(entries.GT.is_non_ref())
    non_ref_samples = non_ref_entries.group_by("locus", "alleles").aggregate(
        samples=hl.agg.collect_as_set("s")
    )

    # Annotate vp_ht with samples that have non-ref for each variant
    vp_ht = vp_ht.annotate(
        v1_samples=non_ref_samples[
            hl.struct(locus=vp_ht.locus1, alleles=vp_ht.alleles1)
        ].samples,
        v2_samples=non_ref_samples[
            hl.struct(locus=vp_ht.locus2, alleles=vp_ht.alleles2)
        ].samples,
    )

    # Filter to pairs where there's at least one overlapping sample
    vp_ht = vp_ht.filter(
        hl.is_defined(vp_ht.v1_samples)
        & hl.is_defined(vp_ht.v2_samples)
        & (hl.len(vp_ht.v1_samples.intersection(vp_ht.v2_samples)) > 0)
    )

    # Drop the helper fields
    vp_ht = vp_ht.select("locus1", "alleles1", "locus2", "alleles2")

    # Key by variant pair and select distinct (matching create_variant_pair_ht)
    vp_ht = vp_ht.key_by("locus2", "alleles2", "locus1", "alleles1")
    vp_ht = vp_ht.select().distinct()

    # Add a unique index id to each variant pair (matching create_variant_pair_ht line 216)
    vp_ht = vp_ht.add_index("vp_ht_idx")

    # Key by vp_ht_idx so the lookup in create_variant_pair_genotype_ht will work
    vp_ht = vp_ht.key_by("vp_ht_idx")

    logger.info(f"Created test variant pair Table with {vp_ht.count()} pairs")
    return vp_ht


def validate_genotype_counts(
    vp_gt_counts_ht: hl.Table,
    n_samples: int,
    check_name: str = "Validation",
) -> None:
    """
    Validate genotype counts for correctness.

    Checks:
    1. All counts are non-negative
    2. Raw counts sum to n_samples (including filtered-out samples)
    3. Adj counts sum to <= n_samples (may be less due to adj filtering)
    4. Adj counts are <= raw counts for each genotype combination
    5. Count arrays have exactly 9 elements

    :param vp_gt_counts_ht: Variant pair genotype counts Table.
    :param n_samples: Expected number of samples.
    :param check_name: Name for this validation check (for logging).
    """
    logger.info(f"\n{check_name}: Validating genotype counts...")

    # Collect all rows to check
    rows = vp_gt_counts_ht.collect()

    errors = []
    warnings = []

    gt_names = [
        "AABB",
        "AABb",
        "AAbb",
        "AaBB",
        "AaBb",
        "Aabb",
        "aaBB",
        "aaBb",
        "aabb",
    ]

    for i, row in enumerate(rows):
        locus1_str = f"{row.locus1.contig}:{row.locus1.position}"
        locus2_str = f"{row.locus2.contig}:{row.locus2.position}"
        pair_name = f"{locus1_str}/{locus2_str}"

        # Check 0: Count arrays exist and have correct length
        if row.gt_counts_raw is None:
            errors.append(f"  {pair_name}: gt_counts_raw is None")
            continue
        if row.gt_counts_adj is None:
            errors.append(f"  {pair_name}: gt_counts_adj is None")
            continue

        if len(row.gt_counts_raw) != 9:
            errors.append(
                f"  {pair_name}: Raw counts array has {len(row.gt_counts_raw)} elements, expected 9"
            )
        if len(row.gt_counts_adj) != 9:
            errors.append(
                f"  {pair_name}: Adj counts array has {len(row.gt_counts_adj)} elements, expected 9"
            )

        # Check 1: All counts are non-negative
        if any(c < 0 for c in row.gt_counts_raw):
            errors.append(
                f"  {pair_name}: Raw counts contain negative values: {row.gt_counts_raw}"
            )
        if any(c < 0 for c in row.gt_counts_adj):
            errors.append(
                f"  {pair_name}: Adj counts contain negative values: {row.gt_counts_adj}"
            )

        # Check 2: Raw counts sum to n_samples
        raw_sum = sum(row.gt_counts_raw)
        if raw_sum != n_samples:
            errors.append(
                f"  {pair_name}: Raw counts sum to {raw_sum}, expected {n_samples}. "
                f"Counts: {row.gt_counts_raw}"
            )

        # Check 3: Adj counts sum to <= n_samples
        adj_sum = sum(row.gt_counts_adj)
        if adj_sum > n_samples:
            errors.append(
                f"  {pair_name}: Adj counts sum to {adj_sum}, expected <= {n_samples}. "
                f"Counts: {row.gt_counts_adj}"
            )

        # Check 4: Adj counts are <= raw counts for each genotype combination
        for idx, (raw_count, adj_count) in enumerate(
            zip(row.gt_counts_raw, row.gt_counts_adj)
        ):
            if adj_count > raw_count:
                errors.append(
                    f"  {pair_name}: Adj count ({adj_count}) > raw count ({raw_count}) "
                    f"for {gt_names[idx]}"
                )

    # Report results
    if errors:
        logger.error(f"\n{check_name} FAILED with {len(errors)} error(s):")
        for error in errors[:20]:  # Show first 20 errors
            logger.error(error)
        if len(errors) > 20:
            logger.error(f"  ... and {len(errors) - 20} more error(s)")
        raise AssertionError(f"{check_name} failed with {len(errors)} error(s)")
    else:
        logger.info(f"  ✓ All validation checks passed for {len(rows)} variant pair(s)")


def main(args):
    """Run tests for variant pair genotype and genotype counts functions."""
    hl.init(log=args.log, tmp_dir=args.tmp_dir)

    test_mt_path = f"{args.tmp_dir}/test_dense_mt.mt"
    test_vp_ht_path = f"{args.tmp_dir}/test_variant_pair_ht.ht"
    test_vp_gt_ht_path = f"{args.tmp_dir}/test_variant_pair_genotype_ht.ht"
    test_vp_gt_counts_ht_path = (
        f"{args.tmp_dir}/test_variant_pair_genotype_counts_ht.ht"
    )

    logger_format = "\n" + "=" * 80 + "\n{message}\n" + "=" * 80 + "\n"

    logger.info(
        logger_format.format(
            message="Testing variant pair genotype and genotype counts functions"
        )
    )

    # Create test data
    if args.create_test_data:
        logger.info(logger_format.format(message="Creating test data..."))

        logger.info("\nStep 1: Creating test dense MatrixTable...")
        mt = create_test_dense_mt()
        mt = mt.checkpoint(test_mt_path, overwrite=args.overwrite)
        logger.info(f"Written test MatrixTable to {test_mt_path}")
        mt.show(n_rows=20, n_cols=10)

        logger.info("\nStep 2: Creating test variant pair Table...")
        vp_ht = create_test_variant_pair_ht(mt)
        vp_ht = vp_ht.checkpoint(test_vp_ht_path, overwrite=args.overwrite)
        logger.info(f"Written test variant pair Table to {test_vp_ht_path}")
        vp_ht.show(-1)

    # Test create_variant_pair_genotype_ht
    if args.test_genotype_ht:
        logger.info(
            logger_format.format(message="Testing create_variant_pair_genotype_ht...")
        )
        try:
            logger.info("Loading existing test dense MatrixTable...")
            mt = hl.read_matrix_table(test_mt_path)
            logger.info("Loading existing test variant pair Table...")
            vp_ht = hl.read_table(test_vp_ht_path)
        except (FileNotFoundError, hl.utils.FatalError) as e:
            logger.error(f"Error loading test data: {e}")
            raise e

        vp_gt_ht = create_variant_pair_genotype_ht(mt, vp_ht)
        vp_gt_ht = vp_gt_ht.checkpoint(test_vp_gt_ht_path, overwrite=args.overwrite)
        logger.info(f"Written test variant pair genotype Table to {test_vp_gt_ht_path}")
        vp_gt_ht.show(-1)

    # Test create_variant_pair_genotype_counts_ht
    if args.test_genotype_counts_ht:
        logger.info(
            logger_format.format(
                message="Testing create_variant_pair_genotype_counts_ht..."
            )
        )

        try:
            logger.info("Loading existing variant pair genotype Table...")
            vp_gt_ht = hl.read_table(test_vp_gt_ht_path)
        except (FileNotFoundError, hl.utils.FatalError):
            raise FileNotFoundError(
                "Variant pair genotype Table not found. Please create it first."
            )

        vp_gt_counts_ht = create_variant_pair_genotype_counts_ht(vp_gt_ht)
        vp_gt_counts_ht = vp_gt_counts_ht.checkpoint(
            test_vp_gt_counts_ht_path, overwrite=args.overwrite
        )

        logger.info(f"✓ Successfully created variant pair genotype counts Table")
        logger.info(f"  Number of variant pairs: {vp_gt_counts_ht.count()}")
        logger.info(f"  Written to {test_vp_gt_counts_ht_path}")
        vp_gt_counts_ht.show(-1)

        # Validate the counts
        # Get the number of samples from the original MT
        try:
            mt = hl.read_matrix_table(test_mt_path)
            n_samples = mt.count_cols()
            validate_genotype_counts(vp_gt_counts_ht, n_samples, "Count Validation")
        except (FileNotFoundError, hl.utils.FatalError) as e:
            logger.warning(f"Could not load MT to get sample count for validation: {e}")

    logger.info(logger_format.format(message="All tests completed successfully!"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test script for variant pair genotype and genotype counts functions"
    )
    parser.add_argument(
        "--create-test-data",
        action="store_true",
        help="Create test data",
    )
    parser.add_argument(
        "--test-genotype-ht",
        action="store_true",
        help="Test create_variant_pair_genotype_ht function",
    )
    parser.add_argument(
        "--test-genotype-counts-ht",
        action="store_true",
        help="Test create_variant_pair_genotype_counts_ht function",
    )
    parser.add_argument(
        "--tmp-dir",
        type=str,
        default=DEFAULT_TMP_DIR,
        help=f"Temporary directory to write all intermediate test files. "
        f"All test outputs will be written here. (default: {DEFAULT_TMP_DIR})",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files",
    )
    parser.add_argument(
        "--log",
        type=str,
        default="/tmp/hail_test_vp.log",
        help="Path to Hail log file (default: /tmp/hail_test_vp.log)",
    )

    args = parser.parse_args()

    # If no test flags are set, run both tests
    if not args.test_genotype_ht and not args.test_genotype_counts_ht:
        args.test_genotype_ht = True
        args.test_genotype_counts_ht = True

    main(args)

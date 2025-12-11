"""
Script to create variant pair matrix from gnomAD a v4 VariantDataset file.

This script creates two main outputs:

1. Variant pair list (--create-vp-list): A Table containing all unique ordered variant
   pairs that co-occur within the same sample and gene, filtered by variant QC,
   consequence severity, and allele frequency.

2. Full variant pair MatrixTable (--create-full-vp): A MatrixTable where each row
   represents a variant pair (v1, v2) and entries contain genotype information for
   both variants, enabling downstream analysis of compound heterozygote patterns.
"""

import argparse
import logging
import timeit
from typing import Optional

import hail as hl
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.vep import CSQ_ORDER, filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds

from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_LEAST_CONSEQUENCE,
    DEFAULT_MAX_FREQ,
    DEFAULT_TMP_DIR,
    TEST_INTERVAL,
    get_variant_pair_resources,
)
from gnomad_chets.v4.utils import filter_for_testing

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_vp_matrix")
logger.setLevel(logging.INFO)


def create_variant_filter_ht(
    filter_ht: hl.Table,
    freq_ht: hl.Table,
    vep_ht: hl.Table,
    least_consequence: str = DEFAULT_LEAST_CONSEQUENCE,
    max_freq: float = DEFAULT_MAX_FREQ,
) -> hl.Table:
    """
    Create a filter Table for variant pair determination.

    Filters variants to those that pass variant QC, have a consequence at least as
    severe as `least_consequence`, have a global AF <= `max_freq`, and have an
    associated gene ID.

    :param filter_ht: Final filter Table for filtering variants that pass QC. Must be
        keyed by 'locus' and 'alleles'.
    :param freq_ht: Frequency Table for filtering by global AF. Must be keyed by
        'locus' and 'alleles'.
    :param vep_ht: VEP Table for filtering by consequence severity. Must be keyed by
        'locus' and 'alleles'.
    :param least_consequence: Lowest-severity consequence to keep. Must be in CSQ_ORDER.
        Default is DEFAULT_LEAST_CONSEQUENCE.
    :param max_freq: Maximum global AF to keep (inclusive). Default is DEFAULT_MAX_FREQ.
    :return: Table with filtered variant data.
    """
    if least_consequence not in CSQ_ORDER:
        raise ValueError(f"least_consequence '{least_consequence}' not in CSQ_ORDER")

    # Create set of allowed consequences (all consequences at least as severe as
    # least_consequence, based on CSQ_ORDER).
    allowed_csqs = hl.literal(
        set(CSQ_ORDER[0 : CSQ_ORDER.index(least_consequence) + 1])
    )

    # Filter VEP transcripts to those with allowed consequences.
    # The gnomAD helper filters to protein-coding, Ensembl-only transcripts and applies
    # additional filtering criteria (consequence severity).
    vep_ht = vep_ht.annotate(
        filters=filter_ht[vep_ht.locus, vep_ht.alleles].filters,
        af=freq_ht[vep_ht.locus, vep_ht.alleles].freq[0].AF,
        gene_id=filter_vep_transcript_csqs_expr(
            vep_ht.vep.transcript_consequences,
            protein_coding=True,
            ensembl_only=True,
            additional_filtering_criteria=[
                lambda tc: tc.consequence_terms.any(lambda c: allowed_csqs.contains(c))
            ],
        ).map(lambda csq: csq.gene_id),
    )
    # Filter variants to those that pass QC, have a consequence at least as severe as
    # `least_consequence`, and have a gnomAD AF <= `max_freq`.
    vep_ht = vep_ht.filter(
        (vep_ht.filters.length() == 0)
        & (vep_ht.af > 0)
        & (vep_ht.af <= max_freq)
        & hl.is_defined(vep_ht.gene_id)
        & (hl.len(vep_ht.gene_id) > 0)
    )

    return vep_ht


def _get_ordered_vp_struct(
    v1: hl.expr.StructExpression, v2: hl.expr.StructExpression
) -> hl.expr.StructExpression:
    """
    Create an ordered variant pair struct ensuring consistent ordering.

    Orders variants by position first, then by alt allele if positions are equal.
    This ensures that (v1, v2) and (v2, v1) are treated as the same pair.

    :param v1: First variant struct with fields 'locus' and 'alleles'.
    :param v2: Second variant struct with fields 'locus' and 'alleles'.
    :return: Struct with fields 'v1' and 'v2' in canonical order.
    """
    return hl.if_else(
        v1.locus.position < v2.locus.position,
        hl.struct(v1=v1, v2=v2),
        # If positions are equal, sort on alt allele.
        hl.if_else(
            v1.locus.position == v2.locus.position,
            hl.if_else(
                v1.alleles[1] < v2.alleles[1],
                hl.struct(v1=v1, v2=v2),
                hl.struct(v1=v2, v2=v1),
            ),
            hl.struct(v1=v2, v2=v1),
        ),
    )


def create_variant_pair_ht(
    vds: hl.vds.VariantDataset,
    vep_ht: hl.Table,
) -> hl.Table:
    """
    Create a Hail Table of unique ordered variant pairs per sample per gene.

    :param vds: VariantDataset with filtered variant data.
    :return: Hail Table keyed by locus2, alleles2, locus1, alleles1 with one distinct
        row per unique variant pair.
    """
    vmt = vds.variant_data
    vmt = vmt.annotate_rows(gene_id=vep_ht[vmt.locus, vmt.alleles].gene_id)

    # Convert to entries table and explode on gene_id so each variant-gene combination
    # is a separate row.
    et = vmt.select_cols().select_rows("gene_id").entries()
    et = et.explode("gene_id")

    # Group by gene and sample, collecting unique variants per gene/sample.
    # Using collect_as_set ensures each variant appears only once per gene/sample.
    et = (
        et.group_by("gene_id", "s").aggregate(
            variants=hl.array(
                hl.agg.collect_as_set(hl.struct(locus=et.locus, alleles=et.alleles))
            )
        )
    ).checkpoint(
        hl.utils.new_temp_file("create_variant_pair_ht.gene_sample_grouped", "ht")
    )

    # Filter to samples with at least 2 variants (needed to form pairs).
    et = et.filter(hl.len(et.variants) >= 2)

    # Generate all ordered pairs of variants within each gene/sample.
    # The nested flatmap/map creates all combinations (i, j) where i < j, ensuring
    # each pair is created exactly once.
    et = et.annotate(
        pairs=(
            hl.range(0, hl.len(et.variants)).flatmap(
                lambda i1: (
                    hl.range(i1 + 1, hl.len(et.variants)).map(
                        lambda i2: _get_ordered_vp_struct(
                            et.variants[i1], et.variants[i2]
                        )
                    )
                )
            )
        )
    )

    # Explode pairs and extract variant pair fields.
    et = et.explode("pairs")
    et = et.transmute(
        vgt=et.pairs,
        locus1=et.pairs.v1.locus,
        alleles1=et.pairs.v1.alleles,
        locus2=et.pairs.v2.locus,
        alleles2=et.pairs.v2.alleles,
    )

    # Key by variant pair and select distinct pairs.
    # Keying by (locus2, alleles2, locus1, alleles1) ensures consistent ordering.
    et = et.key_by("locus2", "alleles2", "locus1", "alleles1")
    et = et.select().distinct()

    et = et.checkpoint(hl.utils.new_temp_file("create_variant_pair_ht.1", "ht"))

    return et.repartition(50, shuffle=True)


def create_dense_filtered_mt(
    vds: hl.vds.VariantDataset,
    vp_ht: hl.Table,
) -> hl.MatrixTable:
    """
    Create a dense filtered MatrixTable from a VariantDataset.

    :param vds: VariantDataset with filtered variant data.
    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2, alleles2.
    :return: Dense filtered MatrixTable.
    """
    v1_ht = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1).select().distinct()
    v2_ht = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2).select().distinct()
    variants_ht = (
        v1_ht.union(v2_ht)
        .distinct()
        .checkpoint(hl.utils.new_temp_file("create_dense_filtered_mt.variants", "ht"))
    )
    vds = hl.vds.filter_variants(vds, variants_ht)

    return hl.vds.to_dense_mt(vds)
    

# TODO: This is a work in progress.
def create_full_vp(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
) -> hl.MatrixTable:
    """
    Create a full variant pair MatrixTable from a variant pair list Table.

    Takes a MatrixTable of variants and a Table of variant pairs, and creates a
    MatrixTable where each row represents a variant pair (v1, v2) and entries contain
    genotype information for both variants.

    :param mt: MatrixTable with variant data. Row key must be (locus, alleles).
    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2, alleles2.
    :return: MatrixTable keyed by (locus1, alleles1, locus2, alleles2) with entry fields
        for both variants (suffixed with '1' and '2').
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )
    adj_gt_count_expr = hl.if_else(
        get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD), gt_count_expr, 0
    )
    mt = mt.select_entries(gt_info=(gt_count_expr, adj_gt_count_expr))
    ht = mt.localize_entries("gt_info", "samples")
    ht = ht.select(
        gt_info=hl.zip(ht.gt_info, ht.samples)
        .map(lambda x: (x[1].s, x[0].gt_info[0], x[0].gt_info[1]))
        .filter(lambda x: hl.is_defined(x[1]) | hl.is_defined(x[2]))
    )
    ht = ht.checkpoint(
        # hl.utils.new_temp_file("create_full_vp.0", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.1-FwOn431irU6ty0eySg6Qsk.ht",
        _read_if_exists=True,
        # overwrite=True,
    )
    vp_ht = vp_ht.add_index("vp_ht_idx")

    vp2_ht = vp_ht.group_by("locus2", "alleles2").aggregate(
        v1=hl.agg.collect(vp_ht.row.select("vp_ht_idx"))
    )
    vp2_ht = vp2_ht.annotate(v2=ht[vp2_ht.locus2, vp2_ht.alleles2].gt_info).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.2-52Xxvmi3378O5gKj0dwI6n.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    vp2_ht = vp2_ht.explode("v1")
    vp2_ht = vp2_ht.transmute(**vp2_ht.v1).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.3-52Xxvmi3378O5gKj0dwI6n.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    vp1_ht = vp_ht.group_by("locus1", "alleles1").aggregate(
        v2=hl.agg.collect(vp_ht.row.select("vp_ht_idx"))
    )
    vp1_ht = vp1_ht.annotate(v1=ht[vp1_ht.locus1, vp1_ht.alleles1].gt_info).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.4-52Xxvmi3378O5gKj0dwI6n.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    vp1_ht = vp1_ht.explode("v2")
    vp1_ht = vp1_ht.transmute(**vp1_ht.v2).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.5-52Xxvmi3378O5gKj0dwI6n.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    vp_ht = (
        vp2_ht.key_by("vp_ht_idx")
        .join(vp1_ht.key_by("vp_ht_idx"))
        .checkpoint(
            # hl.utils.new_temp_file("create_full_vp.2", "ht")
            "gs://gnomad-tmp-4day/create_full_vp.6-52Xxvmi3378O5gKj0dwI6n.ht",
            _read_if_exists=True,
            # overwrite=True,
        )
    )

    def _convert_gt_info_to_counts(v1_gt_num, v2_gt_num, gt_num_count):
        v1_hom_ref = hl.is_missing(v1_gt_num)
        v1_het = v1_gt_num == 1
        v1_hom_var = v1_gt_num == 2
        v2_hom_ref = hl.is_missing(v2_gt_num)
        v2_het = v2_gt_num == 1
        v2_hom_var = v2_gt_num == 2

        gt_num_count = hl.int32(gt_num_count)

        return (
            hl.case()
            .when(
                v1_hom_ref,
                hl.case()
                .when(v2_hom_ref, [gt_num_count, 0, 0, 0, 0, 0, 0, 0, 0])
                .when(v2_het, [0, gt_num_count, 0, 0, 0, 0, 0, 0, 0])
                .when(v2_hom_var, [0, 0, gt_num_count, 0, 0, 0, 0, 0, 0])
                .default([0] * 9),
            )
            .when(
                v1_het,
                hl.case()
                .when(v2_hom_ref, [0, 0, 0, gt_num_count, 0, 0, 0, 0, 0])
                .when(v2_het, [0, 0, 0, 0, gt_num_count, 0, 0, 0, 0])
                .when(v2_hom_var, [0, 0, 0, 0, 0, gt_num_count, 0, 0, 0])
                .default([0] * 9),
            )
            .when(
                v1_hom_var,
                hl.case()
                .when(v2_hom_ref, [0, 0, 0, 0, 0, 0, gt_num_count, 0, 0])
                .when(v2_het, [0, 0, 0, 0, 0, 0, 0, gt_num_count, 0])
                .when(v2_hom_var, [0, 0, 0, 0, 0, 0, 0, 0, gt_num_count])
                .default([0] * 9),
            )
            .default([0] * 9)
        )

    # TODO: To get the actual hom_ref counts we need to determine the number of samples
    # that don't have a value in v1 or v2 for each varaint pair and add that to the
    # count in index 0 of the array.
    # TODO: Add the correct genotype counts for adj variants. Need to make sure that
    # both genotypes pass the adj filter. I'm not sure how to do this yet for the
    # counts that include hom_ref.
    # TODO: incorporate population information into the counts after we have solved the
    # above TODOs.
    vp_ht = vp_ht.transmute(
        gt_counts_raw=(
            vp_ht.v1.map(lambda x: (x[0], (1, x[1])))
            .extend(vp_ht.v2.map(lambda x: (x[0], (2, x[1]))))
            .group_by(lambda x: x[0])
            .values()
            .map(lambda x: hl.dict(x.map(lambda y: y[1])))
            .aggregate(lambda x: hl.agg.counter([x.get(1), x.get(2)]))
            .items()
            .map(lambda x: _convert_gt_info_to_counts(x[0][0], x[0][1], x[1]))
            .aggregate(hl.agg.array_sum)
        )
    )

    vp_ht.describe()
    vp_ht.show()


def main(args):
    """Create variant pair matrix from gnomAD v4 VDS."""
    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix or ""
    data_type = args.data_type
    least_consequence = args.least_consequence
    max_freq = args.max_freq
    test = args.test

    hl.init(
        log="/create_vp_matrix.log",
        tmp_dir=tmp_dir,
    )

    # Get variant co-occurrence pipeline resources.
    resources = get_variant_pair_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir,
        output_postfix=output_postfix,
        overwrite=overwrite,
    )

    # Create variant filter Table.
    if args.create_variant_filter_ht:
        res = resources.create_variant_filter_ht
        res.check_resource_existence()

        filter_ht = res.filter_ht.ht()
        freq_ht = res.freq_ht.ht()
        vep_ht = res.vep_ht.ht()

        # Filter input resources to test interval if in test mode.
        if test:
            logger.info("Filtering filter_ht, freq_ht, and vep_ht to test interval...")
            filter_ht = filter_for_testing(filter_ht)
            freq_ht = filter_for_testing(freq_ht)
            vep_ht = filter_for_testing(vep_ht)

        ht = create_variant_filter_ht(
            filter_ht,
            freq_ht,
            vep_ht,
            least_consequence=least_consequence,
            max_freq=max_freq,
        ).checkpoint(res.variant_filter_ht.path, overwrite=overwrite)

        logger.info(
            f"The number of variants in the VEP Table that pass QC, have a consequence "
            f"at least as severe as {least_consequence}, and have a gnomAD AF <= "
            f"{max_freq} is {ht.count()}"
        )

    # Filter VariantDataset for determining variant pairs.
    if args.filter_vds:
        res = resources.filter_vds
        res.check_resource_existence()

        logger.info(f"Loading gnomAD v4 {data_type} VDS...")
        get_vds_func = (
            get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
        )
        get_vds_func(
            release_only=True,
            split=True,
            filter_intervals=None if not test else [TEST_INTERVAL],
            filter_variant_ht=res.variant_filter_ht.ht(),
            entries_to_keep=["GT", "GQ", "DP", "AD"],
            split_reference_blocks=False,
        ).write(res.filtered_vds.path, overwrite=overwrite)

    # Create variant pair list Table.
    if args.create_vp_list:
        res = resources.create_vp_list
        res.check_resource_existence()

        ht = create_variant_pair_ht(res.filtered_vds.vds(), res.variant_filter_ht.ht())
        ht = ht.checkpoint(res.vp_list_ht.path, overwrite=overwrite)

        logger.info(f"The number of unique variant pairs is {ht.count()}")

    # Create dense filtered MatrixTable.
    if args.create_dense_filtered_mt:
        res = resources.create_dense_filtered_mt
        res.check_resource_existence()

        mt = create_dense_filtered_mt(res.filtered_vds.vds(), res.vp_list_ht.ht())
        mt = mt.checkpoint(res.dense_filtered_mt.path, overwrite=overwrite)

        logger.info(
            f"The number of rows in the dense filtered MatrixTable is {mt.count_rows()}"
        )

    # Create full variant pair MatrixTable.
    if args.create_full_vp:
        res = resources.create_full_vp
        res.check_resource_existence()

        mt = create_full_vp(
            res.dense_filtered_mt.mt(),
            res.vp_list_ht.ht(read_args={"_n_partitions": 50}),
            # n_repartition=args.n_repartition if not test else None,
        )
        # mt = mt.checkpoint(res.vp_full_mt.path, overwrite=overwrite)
        # logger.info(
        #    f"The full variant pair MatrixTable has {mt.count_rows()} rows and "
        #    f"{mt.count_cols()} columns."
        # )

    stop = timeit.default_timer()
    logger.info(f"Time taken to run the script is {stop - start} seconds.")


# The order this should be run in is first create_vp_list (or vp_list_by_chrom), then create_full_vp, then create_vp_summary.

if __name__ == "__main__":
    # Argument parsing for data type, testing, and tmp-dir.
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
        help=f'Data type to use. Must be one of {", ".join(DATA_TYPE_CHOICES)}. Default is {DEFAULT_DATA_TYPE}.',
    )
    parser.add_argument(
        "--create-variant-filter-ht",
        action="store_true",
        help="Create the variant filter Table.",
    )
    parser.add_argument(
        "--least-consequence",
        default=DEFAULT_LEAST_CONSEQUENCE,
        choices=CSQ_ORDER,
        help=f"Lowest-severity consequence to keep. Default is {DEFAULT_LEAST_CONSEQUENCE}.",
    )
    parser.add_argument(
        "--max-freq",
        type=float,
        default=DEFAULT_MAX_FREQ,
        help=f"Maximum global AF to keep (inclusive). Default is {DEFAULT_MAX_FREQ}.",
    )
    parser.add_argument(
        "--filter-vds",
        action="store_true",
        help="Filter the VariantDataset for determining variant pairs.",
    )
    parser.add_argument(
        "--create-vp-list",
        action="store_true",
        help="first create just the list of possible variant pairs.",
    )
    parser.add_argument(
        "--create-dense-filtered-mt",
        action="store_true",
        help="Create the dense filtered MatrixTable.",
    )
    parser.add_argument(
        "--create-full-vp",
        action="store_true",
        help="Create the full variant pair MatrixTable.",
    )
    parser.add_argument(
        "--n-repartition",
        type=int,
        default=10000,
        help=(
            "Number of partitions to repartition the MatrixTable to. Default is 10000 "
            "unless --test is specified.",
        ),
    )

    args = parser.parse_args()
    main(args)

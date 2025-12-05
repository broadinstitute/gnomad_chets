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

# def create_full_vp(
#         mt: hl.MatrixTable,
#         vp_list_ht: hl.Table,
#         data_type: str,
#         tmp_dir: str
# ):
#     # TODO: This implementation was causing memory challenges.

#     vp_list_ht = vp_list_ht.key_by('locus2', 'alleles2')
#     vp_list_ht = vp_list_ht.select(locus1=vp_list_ht.locus1, alleles1=vp_list_ht.alleles1)
#     vp_mt = mt.annotate_rows(v1=vp_list_ht.index(mt.row_key, all_matches=True))
#     vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)
#     vp_mt = vp_mt.rename({x: f'{x}2' for x in vp_mt.entry})

#     vp_mt = vp_mt.explode_rows(vp_mt.v1)
#     vp_mt = vp_mt.transmute_rows(**vp_mt.v1)
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp0.mt', overwrite=True)

#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp1.mt', overwrite=True)

#     mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
#     vp_mt = vp_mt.annotate_entries(**{f'{x}1': mt_joined[x] for x in mt.entry})
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp2.mt', overwrite=True)
#     vp_mt = vp_mt.repartition(10000, shuffle=True)
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp3.mt', overwrite=True)
#     vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})
#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')

#     return vp_mt


# # An efficient implementation of creating a VariantDataset of variant pairs.
# def create_full_vp_vds_efficient(
#     vds_in: hl.vds.VariantDataset,
#     vp_list_ht: hl.Table,                # expected fields: locus1, alleles1, locus2, alleles2
#     tmp_dir: str,
#     entry_fields: Iterable[str] = None,  # None => keep all entry fields; otherwise list like ['GT','DP','GQ']
#     n_partitions: int = 2000,
#     flatten_pairs: bool = False          # if True: create GT1,GT2,... flattened fields (slower)
# ) -> hl.vds.VariantDataset:
#     """
#     Create a VDS of variant pairs efficiently from an input VDS.

#     Returns the VariantDataset object and writes it to out_path_vds.
#     """

#     # Work on variant_data (MatrixTable)
#     mt = vds_in.variant_data

#     # Choose which entry fields to carry. By default, use all.
#     if entry_fields is None:
#         entry_fields = list(mt.entry)

#     # 0) Minimal preselection on the original mt:
#     #    keep only row keys, row fields used by vp_list lookup (locus/alleles), and minimal entry fields.
#     #    This avoids carrying large row/col annotations through the pipeline.
#     mt_min = mt.select_rows().select_cols().select_entries(**{f: mt[f] for f in entry_fields})

#     # 1) Prepare vp_list keyed by locus2/alleles2
#     vp_list_ht = vp_list_ht.key_by('locus2', 'alleles2')
#     vp_list_ht = vp_list_ht.select(locus1=vp_list_ht.locus1, alleles1=vp_list_ht.alleles1)

#     # 2) Attach to mt_min an array of v1 partners (for each v2)
#     #    This requires indexing by mt_min.row_key (which is (locus, alleles))
#     vp_mt = mt_min.annotate_rows(v1 = vp_list_ht.index(mt_min.row_key, all_matches=True))

#     # 3) Keep only v2 rows that actually appear in pairs
#     vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)

#     # 4) Turn existing entry fields into a single nested struct v2 to avoid many columns
#     #    This is much more compact in Spark and Hail than duplicating many top-level fields.
#     vp_mt = vp_mt.annotate_entries(v2 = hl.struct(**{f: vp_mt[f] for f in entry_fields}))
#     # drop the original entry fields, keep only v2
#     vp_mt = vp_mt.select_entries('v2')

#     # 5) Explode v1 list so each row is one (v1, v2) pair (v2 = current row)
#     vp_mt = vp_mt.explode_rows(vp_mt.v1)
#     vp_mt = vp_mt.transmute_rows(**vp_mt.v1)  # brings locus1/alleles1 to top-level row fields

#     # checkpoint an intermediate MT to break lineage (helps memory)
#     cp0 = f'{tmp_dir}/vp_pairs_tmp0.mt'
#     vp_mt = vp_mt.checkpoint(cp0, overwrite=True)

#     # 6) Re-key by locus1/alleles1 so we can slice original mt_min to fetch v1 entries cheaply
#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')

#     # checkpoint again to persist the rekeyed state
#     cp1 = f'{tmp_dir}/vp_pairs_tmp1.mt'
#     vp_mt = vp_mt.checkpoint(cp1, overwrite=True)

#     # 7) Create a tiny mt that contains entry struct 'v1' for the original variants to slice from
#     mt_for_v1 = mt_min.annotate_entries(v1 = hl.struct(**{f: mt_min[f] for f in entry_fields}))
#     mt_for_v1 = mt_for_v1.select_entries('v1')

#     # 8) Slice mt_for_v1 by the vp_mt keys to get v1 entries aligned to vp_mt rows/cols.
#     #    This is a partition-local slice and is much cheaper than full join/repartition.
#     mt_joined = mt_for_v1[vp_mt.row_key, vp_mt.col_key]

#     # 9) Annotate vp_mt entries with v1 struct
#     vp_mt = vp_mt.annotate_entries(v1 = mt_joined.v1)

#     # At this point entry has two nested structs: v1 and v2

#     # checkpoint a compact paired MT
#     cp2 = f'{tmp_dir}/vp_pairs_tmp2.mt'
#     vp_mt = vp_mt.checkpoint(cp2, overwrite=True)

#     # 10) Repartition moderately (n_partitions) to make subsequent writing / VDS packaging efficient.
#     #     Avoid excessively large shuffles; n_partitions should match cluster size.
#     vp_mt = vp_mt.repartition(n_partitions, shuffle=True)

#     cp3 = f'{tmp_dir}/vp_pairs_tmp3.mt'
#     vp_mt = vp_mt.checkpoint(cp3, overwrite=True)

#     # 11) Rename locus/alleles -> locus2/alleles2 (current row fields represent the original v2)
#     vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})

#     # 12) Final canonical row-key ordering
#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')

#     # 13) Optionally flatten the nested structs into top-level fields like GT1, GT2, DP1, DP2, ...
#     if flatten_pairs:
#         # Flatten v1 and v2 structs into top-level entry fields with suffixes
#         # This increases width of the entry schema but may be necessary for downstream code.
#         def flatten_struct_to_entries(mt_in, struct_name, suffix):
#             # mt_in: MatrixTable with entry struct named struct_name (e.g., 'v1' or 'v2')
#             # returns a MatrixTable with additional top-level entry fields like GT1 = mt_in.entry.v1.GT
#             fields = list(entry_fields)
#             # annotate_entries with each flattened field
#             mt_out = mt_in.annotate_entries(**{f'{f}{suffix}': getattr(mt_in.entry[struct_name], f) for f in fields})
#             # we do not drop the nested struct because downstream code might expect it; you can drop it if desired:
#             mt_out = mt_out.select_entries(*[f'{f}{suffix}' for f in fields])
#             return mt_out

#         vp_mt = flatten_struct_to_entries(vp_mt, 'v1', '1')
#         vp_mt = flatten_struct_to_entries(vp_mt, 'v2', '2')

#         # optionally keep nested structs as well or drop them (we dropped them above)

#         # checkpoint after flatten
#         cp_flat = f'{tmp_dir}/vp_pairs_flattened.mt'
#         vp_mt = vp_mt.checkpoint(cp_flat, overwrite=True)

#     # 14) Build final VariantDataset with original sample_data (unchanged)
#     vds_out = hl.vds.VariantDataset(variant_data=vp_mt, sample_data=vds_in.sample_data)
#     return vds_out


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
    filter_ht: hl.Table,
    freq_ht: hl.Table,
    vep_ht: hl.Table,
    least_consequence: str = DEFAULT_LEAST_CONSEQUENCE,
    max_freq: float = DEFAULT_MAX_FREQ,
) -> hl.Table:
    """
    Create a Hail Table of unique ordered variant pairs per sample per gene.

    Filters variants to those that pass variant QC, have a consequence at least as 
    severe as `least_consequence`, and have a global AF <= `max_freq`. Then creates all 
    unique ordered variant pairs that co-occur within the same sample and gene.

    :param vds: VariantDataset with attribute `variant_data` (a MatrixTable).
    :param filter_ht: Final filter Table for filtering variants that pass QC. Must be
        keyed by 'locus' and 'alleles'.
    :param freq_ht: Frequency Table for filtering by global AF. Must be keyed by
        'locus' and 'alleles'.
    :param vep_ht: VEP Table for filtering by consequence severity. Must be keyed by
        'locus' and 'alleles'.
    :param least_consequence: Lowest-severity consequence to keep. Must be in CSQ_ORDER.
        Default is DEFAULT_LEAST_CONSEQUENCE.
    :param max_freq: Maximum global AF to keep (inclusive). Default is DEFAULT_MAX_FREQ.
    :return: Hail Table keyed by locus2, alleles2, locus1, alleles1 with one distinct
        row per unique variant pair.
    """
    if least_consequence not in CSQ_ORDER:
        raise ValueError(f"least_consequence '{least_consequence}' not in CSQ_ORDER")

    mt = vds.variant_data

    # Create set of allowed consequences (all consequences at least as severe as
    # least_consequence, based on CSQ_ORDER).
    allowed_csqs = hl.literal(
        set(CSQ_ORDER[0 : CSQ_ORDER.index(least_consequence) + 1])
    )

    # Filter VEP transcripts to those with allowed consequences.
    # The gnomAD helper filters to protein-coding, Ensembl-only transcripts and applies
    # additional filtering criteria (consequence severity).
    vep_ht = vep_ht.annotate(
        csqs=filter_vep_transcript_csqs_expr(
            vep_ht.vep.transcript_consequences,
            protein_coding=True,
            ensembl_only=True,
            additional_filtering_criteria=[
                lambda tc: tc.consequence_terms.any(lambda c: allowed_csqs.contains(c))
            ],
        )
    )
    # Keep only variants with at least one matching transcript.
    vep_ht = vep_ht.filter(hl.is_defined(vep_ht.csqs) & (hl.len(vep_ht.csqs) > 0))

    # Annotate rows with QC status, allele frequency, and gene IDs.
    mt = mt.annotate_rows(
        passes_qc=hl.is_defined(filter_ht[mt.locus, mt.alleles]),
        af=hl.float32(freq_ht[mt.locus, mt.alleles].freq[0].AF),
        gene_id=vep_ht[mt.locus, mt.alleles].csqs.gene_id,
    )

    # Filter variants to those that pass QC, have a consequence at least as severe as
    # `least_consequence`, and have a gnomAD AF <= `max_freq`.
    mt = mt.filter_rows(
        mt.passes_qc & hl.is_defined(mt.gene_id) & (mt.af > 0) & (mt.af <= max_freq)
    )

    # Convert to entries table and explode on gene_id so each variant-gene combination
    # is a separate row.
    et = mt.select_cols().select_rows("gene_id").entries()
    et = et.explode("gene_id")

    # Group by gene and sample, collecting unique variants per gene/sample.
    # Using collect_as_set ensures each variant appears only once per gene/sample.
    et = (
        et.group_by("gene_id", "s")
        ._set_buffer_size(20)
        .aggregate(
            variants=hl.array(
                hl.agg.collect_as_set(hl.struct(locus=et.locus, alleles=et.alleles))
            )
        )
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

    return et


def create_full_vp(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    n_repartition: Optional[int] = None,
) -> hl.MatrixTable:
    """
    Create a full variant pair MatrixTable from a variant pair list Table.

    Takes a MatrixTable of variants and a Table of variant pairs, and creates a
    MatrixTable where each row represents a variant pair (v1, v2) and entries contain
    genotype information for both variants.

    :param mt: MatrixTable with variant data. Row key must be (locus, alleles).
    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2, alleles2.
    :param n_repartition: Number of partitions to repartition the MatrixTable to. If
        None, will not repartition. Default is None.
    :return: MatrixTable keyed by (locus1, alleles1, locus2, alleles2) with entry fields
        for both variants (suffixed with '1' and '2').
    """
    # Prepare variant pair table for lookup by v2 (locus2, alleles2).
    vp_ht = vp_ht.key_by('locus2', 'alleles2')
    vp_ht = vp_ht.select(locus1=vp_ht.locus1, alleles1=vp_ht.alleles1)

    # Annotate each v2 row with all possible v1 partners.
    # all_matches=True returns an array since one v2 can pair with multiple v1s.
    vp_mt = mt.annotate_rows(v1=vp_ht.index(mt.row_key, all_matches=True))
    # Keep only rows that are part of at least one variant pair.
    vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)

    # Rename existing entry fields to v2 suffix (these represent the second variant).
    vp_mt = vp_mt.rename({x: f'{x}2' for x in vp_mt.entry})

    # Explode on v1 array to create one row per (v1, v2) pair.
    vp_mt = vp_mt.explode_rows(vp_mt.v1)
    # Move v1 fields (locus1, alleles1) to top-level row fields.
    vp_mt = vp_mt.transmute_rows(**vp_mt.v1)
    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.1", "mt"))

    # Re-key by v1 to enable efficient lookup of v1 entries.
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.2.keyed", "mt"))

    # Lookup v1 entries from original MatrixTable and annotate with v1 suffix.
    # This slice operation is efficient because vp_mt is keyed by (locus1, alleles1).
    mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
    vp_mt = vp_mt.annotate_entries(**{f'{x}1': mt_joined[x] for x in mt.entry})
    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.3.annotated", "mt"))

    # Repartition for efficient final write.
    if n_repartition is not None:
        vp_mt = vp_mt.repartition(n_repartition, shuffle=True)

    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.4.repartitioned", "mt"))

    # Rename current row key (locus, alleles) to v2 and set final row key.
    vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})
    vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')

    return vp_mt


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

    logger.info(f"Loading gnomAD v4 {data_type} VDS...")
    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )
    vds = get_vds_func(
        test=test,
        release_only=True,
        chrom=None if not test else TEST_INTERVAL.split(":")[0],
        split=True,
        filter_intervals=None if not test else [TEST_INTERVAL],
    )

    if test:
        logger.info("Filtered to PCNT gene (chr21:46324141-46445769) for testing...")

    # Create variant pair list Table.
    if args.create_vp_list:
        res = resources.create_vp_list
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

        ht = create_variant_pair_ht(
            vds,
            filter_ht,
            freq_ht,
            vep_ht,
            least_consequence=least_consequence,
            max_freq=max_freq,
        )
        ht = ht.checkpoint(res.vp_list_ht.path, overwrite=overwrite)

        logger.info(f"The number of unique variant pairs is {ht.count()}")

    # Create full variant pair MatrixTable.
    if args.create_full_vp:
        res = resources.create_full_vp
        res.check_resource_existence()
        
        # TODO: Will need to change with a densify added if necessary.
        mt = create_full_vp(
            vds.variant_data, 
            res.vp_list_ht.ht(), 
            n_repartition=args.n_repartition if not test else None,
        )
        mt = mt.checkpoint(res.vp_full_mt.path, overwrite=overwrite)
        logger.info(
            f"The full variant pair MatrixTable has {mt.count_rows()} rows and "
            f"{mt.count_cols()} columns."
        )

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
        "--create-vp-list",
        action="store_true",
        help="first create just the list of possible variant pairs.",
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

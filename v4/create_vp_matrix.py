"""
Script to create variant co-occurrence pipeline outputs from gnomAD v4 VariantDataset.

Pipeline steps (run in order):

1. Variant filter Table (--create-variant-filter-ht): Filters variants to those that
   pass QC, have a consequence at least as severe as the specified threshold, and have
   a global AF <= the specified maximum frequency.

2. Filtered VariantDataset (--filter-vds): Filters the gnomAD v4 VariantDataset to
   only include variants that pass the variant filter criteria.

3. Variant pair list Table (--create-variant-pair-list-ht): Creates a Table containing
   all unique ordered variant pairs that co-occur within the same sample and gene.

4. Dense filtered MatrixTable (--create-dense-filtered-mt): Creates a dense MatrixTable
   containing only the variants present in the variant pair list.

5. Variant pair genotype Table (--create-variant-pair-genotype-ht): Creates a Table
   with genotype information for both variants in each variant pair.

6. Variant pair genotype counts Table (--create-variant-pair-genotype-counts-ht):
   Creates a Table with genotype count arrays (raw and adj) for each variant pair,
   enabling downstream analysis of compound heterozygote patterns.
"""

import argparse
import logging
import timeit
from typing import Optional

import hail as hl
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.vep import CSQ_ORDER, filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds, add_meta

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


def create_variant_pair_mt(
    vds: hl.vds.VariantDataset,
    vep_ht: hl.Table,
) -> hl.Table:
    """
    Create a Hail Table of unique ordered variant pairs per sample per gene.

    :param vds: VariantDataset with filtered variant data.
    :param vep_ht: VEP Table with gene_id annotation. Must be keyed by 'locus' and
        'alleles'.
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
                hl.agg.collect_as_set(hl.struct(locus=et.locus, alleles=et.alleles,GT=et.GT))
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
        GT1=et.pairs.v1.GT,
        GT2=et.pairs.v2.GT
    )

    # Key by variant pair and select distinct pairs.
    # Keying by (locus2, alleles2, locus1, alleles1) ensures consistent ordering.
    et = et.key_by("locus2", "alleles2", "locus1", "alleles1")
    et = et.select("gene_id", "s", "GT1", "GT2").distinct()
    

    # Add a unique index id to each variant pair.
    et = et.add_index("vp_ht_idx")
    
    #convert to matrix table
    mt = et.to_matrix_table(
    row_key=["locus1", "alleles1", "locus2", "alleles2"],
    col_key=["s"],
    row_fields=["gene_id"]
)

    return mt


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
    return vds


def _encode_and_localize_genotypes(mt: hl.MatrixTable) -> hl.Table:
    """
    Encode genotypes for efficiency and localize to a Table.

    For most rare variants, most samples are hom_ref, so code them as missing to save
    space.

    Encodes genotypes as:

        - missing = hom_ref (space saving)
        - 0 = missing data (actual missing call)
        - 1 = het
        - 2 = hom_var

    For adj genotypes:

        - missing = hom_ref adj (space saving)
        - 0 = missing data or not adj (actual missing call)
        - 1 = het adj
        - 2 = hom_var adj

    Filters to keep only genotypes where the variant is called (reduces array size).

    :param mt: MatrixTable with variant data. Row key must be (locus, alleles).
    :return: Table with localized genotype info, filtered to called variants only.
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )

    # For adj genotypes, set to 0 (missing data) if the variant doesn't pass the adj filter.
    adj_gt_count_expr = hl.if_else(
        get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD), gt_count_expr, 0, missing_false=True
    )

    mt = mt.select_entries(gt_info=(gt_count_expr, adj_gt_count_expr))
    ht = mt.localize_entries("gt_info", "samples")

    # Store sample information: (sample_id, raw_gt_count, adj_gt_count).
    # Filter to keep only genotypes where the variant is called (reduces array size).
    ht = ht.select(
        gt_info=hl.enumerate(ht.gt_info)
        .map(lambda x: (x[0], x[1].gt_info[0], x[1].gt_info[1]))
        .filter(lambda x: hl.is_defined(x[1]) | hl.is_defined(x[2]))
    )
    ht = ht.cache()

    return ht


def _prepare_variant_pair_index(vp_ht: hl.Table) -> hl.Table:
    """
    Prepare variant pair Table for genotype annotation.

    Creates separate entries for each variant in the pair (v1 and v2), then unions them.
    This helps with performance issues when annotating genotype info for both variants.

    :param vp_ht: Variant pair Table with fields locus1, alleles1, locus2, alleles2.
    :return: Unioned Table with index field and variant pair indicator (vp=1 or vp=2).
    """
    # Capture number of partitions before modifying the table.
    n_partitions = vp_ht.n_partitions()

    vp1_ht = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1)
    vp1_ht = vp1_ht.select("vp_ht_idx", vp=1)

    vp2_ht = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2)
    vp2_ht = vp2_ht.select("vp_ht_idx", vp=2)

    # Repartition to make the partitions more even to avoid performance issues when
    # annotating with genotype info.
    vp_union_ht = (
        vp1_ht.union(vp2_ht)
        .collect_by_key()
        .repartition(n_partitions, shuffle=True)
        .cache()
    )

    return vp_union_ht


def _annotate_variant_pairs_with_genotypes(
    vp_union_ht: hl.Table,
    ht: hl.Table,
) -> hl.Table:
    """
    Annotate variant pairs with genotype information and group by variant pair index.

    Takes a unioned variant pair Table and annotates each variant with its genotype
    info, then groups by variant pair index to collect genotype info for both variants
    in each pair.

    :param vp_union_ht: Unioned variant pair Table with index field and variant pair
        indicator.
    :param ht: Localized entries Table with genotype info.
    :return: Variant pair Table with genotype info grouped by variant pair index.
    """
    vp_union_ht = vp_union_ht.annotate(
        gt_info=ht[vp_union_ht.locus, vp_union_ht.alleles].gt_info
    )
    vp_union_ht = vp_union_ht.explode("values")
    vp_union_ht = vp_union_ht.transmute(**vp_union_ht.values)

    vp_union_ht = vp_union_ht.annotate(
        gt_info=ht[vp_union_ht.locus, vp_union_ht.alleles].gt_info
    )
    vp_union_ht = vp_union_ht.group_by("vp_ht_idx").aggregate(
        gt_info=hl.agg.collect((vp_union_ht.vp, vp_union_ht.gt_info))
    )

    # Store sample information in the variant pair Table.
    vp_union_ht = vp_union_ht.annotate_globals(samples=ht.index_globals().samples)

    return vp_union_ht


def create_variant_pair_genotype_ht(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
) -> hl.Table:
    """
    Create a variant pair genotype Table from a MatrixTable and variant pair list.

    Encodes genotypes for efficiency, prepares the variant pair Table for annotation,
    and annotates each variant pair with genotype information for both variants.

    :param mt: MatrixTable with variant data. Row key must be (locus, alleles).
    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2, alleles2.
    :return: Variant pair Table with genotype info for both variants in each pair.
    """
    # Encode genotypes and localize to Table.
    ht = _encode_and_localize_genotypes(mt)

    # Add index to variant pair Table so both variants can be annotated with genotype
    # info separately and then joined together by the common index. This helps with
    # performance issues observed when trying to annotate both variants with genotype
    # info simultaneously.
    vp_ht = vp_ht.key_by("locus1", "alleles1", "locus2", "alleles2")
    vp_ht = vp_ht.add_index("vp_ht_idx").key_by("vp_ht_idx").cache()

    # Prepare variant pair Table for genotype annotation.
    vp_union_ht = _prepare_variant_pair_index(vp_ht)

    # Annotate variant pairs with genotype information.
    vp_union_ht = _annotate_variant_pairs_with_genotypes(vp_union_ht, ht)

    # Add back the original variant pair fields.
    vp_ht = vp_union_ht.annotate(**vp_ht[vp_union_ht.vp_ht_idx])

    return vp_ht


def _convert_gt_info_to_counts(
    gt_counts_dict: hl.expr.DictExpression,
    n_samples_filtered_out_expr: hl.expr.Int32Expression,
) -> hl.expr.ArrayExpression:
    """
    Convert genotype count dictionary to genotype count array.

    Takes a dictionary from `hl.agg.counter([v1_gt, v2_gt])` where keys are arrays of
    [v1_gt, v2_gt] and values are counts, and converts it to a 9-element array.

    Genotype encoding:

        - 0 = hom_ref
        - 1 = het
        - 2 = hom_var

    Genotype count array indices represent:

        0: AABB (hom_ref/hom_ref)
        1: AABb (hom_ref/het)
        2: AAbb (hom_ref/hom_var)
        3: AaBB (het/hom_ref)
        4: AaBb (het/het)
        5: Aabb (het/hom_var)
        6: aaBB (hom_var/hom_ref)
        7: aaBb (hom_var/het)
        8: aabb (hom_var/hom_var)

    The array index is calculated as: v1_genotype * 3 + v2_genotype. This maps to:

        0 * 3 + 0 = 0 (AABB)
        0 * 3 + 1 = 1 (AABb)
        0 * 3 + 2 = 2 (AAbb)
        1 * 3 + 0 = 3 (AaBB)
        1 * 3 + 1 = 4 (AaBb)
        1 * 3 + 2 = 5 (Aabb)
        2 * 3 + 0 = 6 (aaBB)
        2 * 3 + 1 = 7 (aaBb)
        2 * 3 + 2 = 8 (aabb)

    The `n_samples_filtered_out_expr` parameter is used to add hom_ref/hom_ref counts
    for filtered-out samples. These samples are filtered out during encoding because
    both variants are hom_ref, which we code as missing for space saving. Both raw and
    adj counts need `n_samples_filtered_out_expr` added to get the correct total count.
    For raw counts, the [0, 0] entry in the dictionary may already contain some counts
    from samples where both variants are hom_ref (if any were not filtered out). For
    adj counts, the [0, 0] entry in the dictionary will always be 0 because all samples
    where both variants are hom_ref are filtered out during encoding, and samples where
    variants are hom_ref but fail adj filtering are represented as 0 (missing data)
    rather than hom_ref.

    :param gt_counts_dict: Dictionary from counter aggregation with keys as tuples
        (v1_gt, v2_gt) and values as counts.
    :param n_samples_filtered_out_expr: Number of samples filtered out (missing both
        variants) to add to hom_ref/hom_ref count (index 0).
    :return: Array of 9 integers representing counts for each genotype combination.
    """
    # Create array by looking up counts for each index in the dictionary.
    # Build arrays of indices and corresponding dict keys for readability.
    indices = hl.range(0, 9)  # [0, 1, 2, 3, 4, 5, 6, 7, 8]
    dict_keys = hl.array(
        [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    )

    # Zip indices and keys, then look up counts in dictionary.
    return hl.zip(indices, dict_keys).map(
        lambda x: hl.if_else(
            x[0] == 0,
            # Index 0: hom_ref/hom_ref - add filtered-out samples.
            gt_counts_dict.get(x[1], 0) + n_samples_filtered_out_expr,
            gt_counts_dict.get(x[1], 0),
        )
    )


def _calculate_genotype_counts(
    per_sample_gt_expr: hl.expr.ArrayExpression,
    n_samples_filtered_out_expr: hl.expr.Int32Expression,
    gt_field: str,
) -> hl.expr.ArrayExpression:
    """
    Calculate genotype counts for variant pairs from per-sample genotype expressions.

    This function processes both raw and adj genotype counts using the same logic:

        1. Extracts genotype values for each variant in the pair.
        2. Filters to keep samples where both variants are defined (not missing),
           which means they are either hom_ref (missing, will be remapped to 0), het (1),
           or hom_var (2), or missing data (0).
        3. Remaps missing (hom_ref) to 0, keeping missing data (0) as 0.
        4. Counts (v1_gt, v2_gt) combinations.
        5. Converts to array format using _convert_gt_info_to_counts.
        6. Adds hom_ref/hom_ref counts for filtered-out samples.

    Before remapping, genotype values are:

        - missing = hom_ref (space saving)
        - 0 = missing data (actual missing call or failed adj for adj_gt)
        - 1 = het
        - 2 = hom_var

    After remapping (before calling _convert_gt_info_to_counts), genotype values are:

        - 0 = hom_ref (remapped from missing) or missing data (already 0)
        - 1 = het
        - 2 = hom_var

    This applies to both raw and adj genotypes after remapping.

    :param per_sample_gt_expr: Array of per-sample genotype info tuples (sample_id,
        v1_gt_info, v2_gt_info).
    :param n_samples_filtered_out_expr: Number of samples filtered out (missing both
        variants).
    :param gt_field: Field name to extract from gt_info ("raw_gt" or "adj_gt").
    :return: Array of genotype counts [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb,
        aabb].
    """
    # Extract genotype counts for each variant in the pair.
    # per_sample_gt_expr is an array of dicts where keys are 1 (v1) and 2 (v2),
    # and values are structs with fields: s, vp, raw_gt, adj_gt.
    # Before remapping: missing=hom_ref, 0=missing data, 1=het, 2=hom_var.
    gt_expr = per_sample_gt_expr.map(lambda x: [x.get(1)[gt_field], x.get(2)[gt_field]])

    # Keep samples where both variants are defined (not both missing data).
    # Before remapping: missing=hom_ref, 0=missing data, 1=het, 2=hom_var.
    # Filter out cases where both are 0 (missing data), but keep cases where one or both
    # are missing (hom_ref) or have valid genotypes (1 or 2).
    gt_expr = gt_expr.filter(
        lambda x: (
            (hl.is_missing(x[0]) | (x[0] != 0)) & (hl.is_missing(x[1]) | (x[1] != 0))
        )
    )
    # Remap missing (hom_ref) to 0. Missing data (already 0) stays as 0.
    # After this, all values are 0 (hom_ref or missing data), 1 (het), or 2 (hom_var).
    # Note: After remapping, we cannot distinguish hom_ref from missing data (both are 0),
    # but this is acceptable since we've already filtered out cases where both are missing data.
    gt_expr = gt_expr.map(lambda x: x.map(lambda y: hl.or_else(y, 0)))

    # Count (v1_gt, v2_gt) combinations and convert to array of counts.
    # At this point, all genotypes are 0 (hom_ref or missing data), 1 (het), or 2 (hom_var).
    # Structure: [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]
    gt_expr = hl.or_missing(
        gt_expr.length() > 0,
        _convert_gt_info_to_counts(
            gt_expr.aggregate(hl.agg.counter),
            n_samples_filtered_out_expr,
        ),
    )

    return gt_expr


def create_variant_pair_genotype_counts_ht(ht: hl.Table) -> hl.Table:
    """
    Create a variant pair genotype counts Table from a variant pair genotype Table.

    Calculates raw and adj genotype counts for each variant pair, grouping samples by
    their genotype combinations and converting to count arrays.

    :param ht: Variant pair genotype Table with gt_info field containing per-sample
        genotype information for both variants.
    :return: Variant pair Table with gt_counts_raw and gt_counts_adj fields containing
        genotype count arrays.
    """

    # Calculate genotype counts for raw and adj genotypes.
    # Structure: [variant pair index, (sample ID, raw GT, adj GT)].
    ht = ht.annotate(
        gt_info=ht.gt_info.flatmap(
            lambda x: x[1].map(
                lambda y: hl.struct(s=y[0], vp=x[0], raw_gt=y[1], adj_gt=y[2])
            )
        )
    )

    # Calculate the total hom_ref/hom_ref that have been filtered out.
    # Samples missing both variants (filtered out) are hom_ref/hom_ref.
    n_samples = ht.samples.length()
    n_samples_with_data_expr = hl.set(ht.gt_info.map(lambda x: x.s)).length()
    n_samples_filtered_out_expr = n_samples - n_samples_with_data_expr

    # Group by sample and create a dict of variant IDs to their genotype info.
    # Per sample array of variants in the pair in the form of:
    # {1: Struct(vp=1, s=sample ID, raw GT, adj GT),
    # 2: Struct(vp=2, s=sample ID, raw GT, adj GT)}
    ht = ht.annotate(
        gt_info=(
            ht.gt_info.group_by(lambda x: x.s)
            .values()
            .map(lambda x: hl.dict(x.map(lambda y: (y.vp, y))))
        ),
        n_samples_filtered_out_expr=n_samples_filtered_out_expr,
    )

    # Calculate raw and adj genotype counts for each variant pair.
    ht = ht.select(
        "locus1",
        "alleles1",
        "locus2",
        "alleles2",
        **{
            f"gt_counts_{n}": _calculate_genotype_counts(
                ht.gt_info, ht.n_samples_filtered_out_expr, gt_field=f"{n}_gt"
            )
            for n in ["raw", "adj"]
        },
    ).cache()

    return ht.key_by("locus1", "alleles1", "locus2", "alleles2")

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
        log=tmp_dir+"/create_vp_matrix.log",
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
            Least consequence: {least_consequence}
            Max freq: {max_freq}
        """
    )

    # Get variant co-occurrence pipeline resources.
    resources = get_variant_pair_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir,
        output_postfix=output_postfix,
        overwrite=overwrite,
    )

    if args.create_variant_filter_ht:
        logger.info("Creating variant filter Table...")
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
            f"Number of variants in the VEP Table that pass QC, have a consequence "
            f"at least as severe as {least_consequence}, and have a gnomAD AF <= "
            f"{max_freq}: {ht.count()}"
        )

    if args.filter_vds:
        logger.info(f"Filtering gnomAD v4 {data_type} VDS...")
        res = resources.filter_vds
        res.check_resource_existence()

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
        logger.info("The filtered VDS has been written...")

    if args.create_variant_pair_list_mt:
        logger.info("Creating variant pair list MatrixTable...")
        res = resources.create_variant_pair_list_mt
        res.check_resource_existence()

        mt = create_variant_pair_mt(res.filtered_vds.vds(), res.variant_filter_ht.ht())
        mt = mt.checkpoint(res.vp_list_mt.path, overwrite=overwrite)
        logger.info(
            "The variant pair list MatrixTable has been written...\n"
            f"The number of unique variant pairs is {mt.count_rows()}"
        )  
        
    if args.create_dense_filtered_mt:
        logger.info("Creating dense filtered MatrixTable...")
        res = resources.create_dense_filtered_mt
        res.check_resource_existence()

        mt = create_dense_filtered_mt(res.filtered_vds.vds(), res.vp_list_ht.ht())
        mt = mt.checkpoint(res.dense_filtered_mt.path, overwrite=overwrite) #FAILED HERE
        logger.info(
            "The dense filtered MatrixTable has been written...\n"
            f"The number of rows in the dense filtered MatrixTable is {mt.count_rows()}"
        )

    # TODO: Decide if we need this step, or if it should just be merged with the counts
    # step. The output of this step will be very large, so if we don't need to actually
    # write it out, that might be better.
    if args.create_variant_pair_genotype_ht:
        logger.info("Creating variant pair genotype Table...")
        res = resources.create_variant_pair_genotype_ht
        res.check_resource_existence()

        ht = create_variant_pair_genotype_ht(
            res.dense_filtered_mt.mt(),
            res.vp_list_ht.ht(),  # (read_args={"_n_partitions": 50}),
            # n_repartition=args.n_repartition if not test else None,
        )
        ht.write(res.vp_gt_ht.path, overwrite=overwrite)
        logger.info("The variant pair genotype Table has been written...")

    # TODO: Add population-specific counts.
    if args.create_variant_pair_genotype_counts_ht:
        logger.info("Creating variant pair genotype counts Table...")
        res = resources.create_variant_pair_genotype_counts_ht
        res.check_resource_existence()

        ht = create_variant_pair_genotype_counts_ht(res.vp_gt_ht.ht())
        ht = ht.checkpoint(res.vp_gt_counts_ht.path, overwrite=overwrite)
        logger.info("The variant pair genotype counts Table has been written...")

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
        "--create-variant-filter-ht",
        action="store_true",
        help="Create the variant filter Table.",
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
        "--create-variant-pair-list-mt",
        action="store_true",
        help="first create just the list of possible variant pairs.",
    )
    parser.add_argument(
        "--create-dense-filtered-mt",
        action="store_true",
        help="Create the dense filtered MatrixTable.",
    )
    parser.add_argument(
        "--create-variant-pair-genotype-ht",
        action="store_true",
        help="Create the full variant pair MatrixTable.",
    )
    parser.add_argument(
        "--create-variant-pair-genotype-counts-ht",
        action="store_true",
        help="Create the variant pair genotype counts Table.",
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

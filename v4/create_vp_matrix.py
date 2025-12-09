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
from gnomad_qc.v4.resources.basics import (get_gnomad_v4_genomes_vds,
                                           get_gnomad_v4_vds)

from gnomad_chets.v4.resources import (DATA_TYPE_CHOICES, DEFAULT_DATA_TYPE,
                                       DEFAULT_LEAST_CONSEQUENCE,
                                       DEFAULT_MAX_FREQ, DEFAULT_TMP_DIR,
                                       TEST_INTERVAL,
                                       get_variant_pair_resources)
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


def create_full_vp_old(
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
    vp_ht = vp_ht.key_by("locus2", "alleles2")
    vp_ht = vp_ht.select(locus1=vp_ht.locus1, alleles1=vp_ht.alleles1)

    # Annotate each v2 row with all possible v1 partners.
    # all_matches=True returns an array since one v2 can pair with multiple v1s.
    vp_mt = mt.annotate_rows(v1=vp_ht.index(mt.row_key, all_matches=True))
    # Keep only rows that are part of at least one variant pair.
    vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)

    # Rename existing entry fields to v2 suffix (these represent the second variant).
    vp_mt = vp_mt.rename({x: f"{x}2" for x in vp_mt.entry})

    # Explode on v1 array to create one row per (v1, v2) pair.
    vp_mt = vp_mt.explode_rows(vp_mt.v1)
    # Move v1 fields (locus1, alleles1) to top-level row fields.
    vp_mt = vp_mt.transmute_rows(**vp_mt.v1)
    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.1", "mt"))

    # Re-key by v1 to enable efficient lookup of v1 entries.
    vp_mt = vp_mt.key_rows_by("locus1", "alleles1")
    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.2.keyed", "mt"))

    # Lookup v1 entries from original MatrixTable and annotate with v1 suffix.
    # This slice operation is efficient because vp_mt is keyed by (locus1, alleles1).
    mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
    vp_mt = vp_mt.annotate_entries(**{f"{x}1": mt_joined[x] for x in mt.entry})
    vp_mt = vp_mt.checkpoint(hl.utils.new_temp_file("create_full_vp.3.annotated", "mt"))

    # Repartition for efficient final write.
    if n_repartition is not None:
        vp_mt = vp_mt.repartition(n_repartition, shuffle=True)

    vp_mt = vp_mt.checkpoint(
        hl.utils.new_temp_file("create_full_vp.4.repartitioned", "mt")
    )

    # Rename current row key (locus, alleles) to v2 and set final row key.
    vp_mt = vp_mt.rename({"locus": "locus2", "alleles": "alleles2"})
    vp_mt = vp_mt.key_rows_by("locus1", "alleles1", "locus2", "alleles2")

    return vp_mt


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
    mt = mt.select_entries(
        GT=hl.or_missing(mt.GT.is_non_ref(), mt.GT),
        missing=hl.is_missing(mt.GT),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    ht = mt.localize_entries("entry_structs", "samples")
    ht = ht.annotate_globals(samples=ht.index_globals().samples.map(lambda x: x.s))
    ht = ht.checkpoint(
        # hl.utils.new_temp_file("create_full_vp.0", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.0-FwOn431irU6ty0eySg6Qsk.ht",
        _read_if_exists=True,
        # overwrite=True,
    )

    vp_ht = vp_ht.add_index("idx").checkpoint(
        # hl.utils.new_temp_file("create_full_vp.1", "ht"),
        "gs://gnomad-tmp-4day/create_full_vp.1-eQ8IM7Aq61Cm31LpTPhMKF.ht",
        # _read_if_exists=True,
        overwrite=True,
    )

    vp1_ht = vp_ht.group_by("locus1", "alleles1").aggregate(
        v2=hl.agg.collect(
            hl.struct(
                idx=vp_ht.idx,
                locus2=vp_ht.locus2,
                alleles2=vp_ht.alleles2,
            )
        )
    )
    vp1_ht = vp1_ht.annotate(
        v1=hl.zip(
            ht.index_globals().samples, ht[vp1_ht.locus1, vp1_ht.alleles1].entry_structs
        ).map(lambda x: x[1].annotate(s=x[0]))
    ).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.2-52Xxvmi3378O5gKj0dwI6n.ht",
        # _read_if_exists=True,
        overwrite=True,
    )

    vp2_ht = vp_ht.group_by("locus2", "alleles2").aggregate(
        v1=hl.agg.collect(
            hl.struct(
                idx=vp_ht.idx,
                locus1=vp_ht.locus1,
                alleles1=vp_ht.alleles1,
            )
        )
    )
    vp2_ht = vp2_ht.annotate(
        v2=hl.zip(
            ht.index_globals().samples, ht[vp2_ht.locus2, vp2_ht.alleles2].entry_structs
        ).map(lambda x: x[1].annotate(s=x[0]))
    ).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.3-52Xxvmi3378O5gKj0dwI6n.ht",
        # _read_if_exists=True,
        overwrite=True,
    )

    vp1_ht = vp1_ht.explode("v2").checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.4-52Xxvmi3378O5gKj0dwI6n.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    vp2_ht = vp2_ht.explode("v1").checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.5-52Xxvmi3378O5gKj0dwI6n.ht",
        # _read_if_exists=True,
        overwrite=True,
    )

    vp1_ht = (
        vp1_ht.annotate(**vp1_ht.v2)
        .key_by("idx")
        .select("v1")
        .checkpoint(
            # hl.utils.new_temp_file("create_full_vp.2", "ht")
            "gs://gnomad-tmp-4day/create_full_vp.6-52Xxvmi3378O5gKj0dwI6n.ht",
            # _read_if_exists=True,
            overwrite=True,
        )
    )
    vp1_ht.describe()
    vp1_ht.show()
    vp2_ht = vp2_ht.annotate(**vp2_ht.v1).key_by("idx")
    vp2_ht = vp2_ht.select("locus2", "alleles2", "locus1", "alleles1", "v2").checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.7-52Xxvmi3378O5gKj0dwI6n.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    vp2_ht.describe()
    vp2_ht.show()

    ht = vp1_ht.join(vp2_ht).checkpoint(
        # hl.utils.new_temp_file("create_full_vp.2", "ht")
        "gs://gnomad-tmp-4day/create_full_vp.8-52Xxvmi3378O5gKj0dwI6n.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    ht.describe()
    ht.show()

    # vp_ht = vp_ht.annotate_globals(
    #    samples=ht.index_globals().columns.map(lambda x: x.s)
    # ).checkpoint(hl.utils.new_temp_file("create_full_vp.1", "ht"))
    # vp_ht = vp_ht.annotate(
    #    v1=ht[vp_ht.locus1, vp_ht.alleles1].entry_structs,
    # ).checkpoint(hl.utils.new_temp_file("create_full_vp.2", "ht"))
    # vp_ht = vp_ht.annotate(
    #    v2=ht[vp_ht.locus2, vp_ht.alleles2].entry_structs,
    # ).checkpoint(hl.utils.new_temp_file("create_full_vp.3", "ht"))
    # vp_ht = vp_ht.select(
    #    entry_structs=hl.zip(vp_ht.samples, vp_ht.v1, vp_ht.v2).map(
    #        lambda x: hl.struct(s=x[0], v1=x[1], v2=x[2])
    #    )
    # ).checkpoint(hl.utils.new_temp_file("create_full_vp.4", "ht"))
    # vp_ht.describe()
    # vp_ht.show()


def get_counts_agg_expr(mt: hl.MatrixTable):
    return (
        hl.case(missing_false=True)
        # 0x
        .when(
            hl.is_missing(mt.GT1) & ~mt.missing1,
            hl.case(missing_false=True)
            .when(hl.is_missing(mt.GT2) & ~mt.missing2, [1, 0, 0, 0, 0, 0, 0, 0, 0])
            .when(mt.GT2.is_het(), [0, 1, 0, 0, 0, 0, 0, 0, 0])
            .when(mt.GT2.is_hom_var(), [0, 0, 1, 0, 0, 0, 0, 0, 0])
            .default([0, 0, 0, 0, 0, 0, 0, 0, 0]),
        )
        # 1x
        .when(
            mt.GT1.is_het(),
            hl.case(missing_false=True)
            .when(hl.is_missing(mt.GT2) & ~mt.missing2, [0, 0, 0, 1, 0, 0, 0, 0, 0])
            .when(mt.GT2.is_het(), [0, 0, 0, 0, 1, 0, 0, 0, 0])
            .when(mt.GT2.is_hom_var(), [0, 0, 0, 0, 0, 1, 0, 0, 0])
            .default([0, 0, 0, 0, 0, 0, 0, 0, 0]),
        )
        # 2x
        .when(
            mt.GT1.is_hom_var(),
            hl.case(missing_false=True)
            .when(hl.is_missing(mt.GT2) & ~mt.missing2, [0, 0, 0, 0, 0, 0, 1, 0, 0])
            .when(mt.GT2.is_het(), [0, 0, 0, 0, 0, 0, 0, 1, 0])
            .when(mt.GT2.is_hom_var(), [0, 0, 0, 0, 0, 0, 0, 0, 1])
            .default([0, 0, 0, 0, 0, 0, 0, 0, 0]),
        ).default([0, 0, 0, 0, 0, 0, 0, 0, 0])
    )


def create_vp_summary(mt: hl.MatrixTable, tmp_dir: str) -> hl.Table:
    mt = mt.select_entries("adj1", "adj2", gt_array=get_counts_agg_expr(mt))
    ht = mt.annotate_rows(
        gt_counts=hl.agg.group_by(
            mt.pop,
            hl.struct(
                raw=hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.gt_array),
                adj=hl.or_else(
                    hl.agg.filter(
                        mt.adj1 & mt.adj2,
                        hl.agg.array_agg(lambda x: hl.agg.sum(x), mt.gt_array),
                    ),
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],  # In case there are no adj entries
                ),
            ),
        )
    ).rows()

    ht = ht.select(
        gt_counts=hl.bind(
            lambda x: hl.dict(
                hl.zip(ht.gt_counts.keys(), ht.gt_counts.values()).append(
                    (
                        "all",
                        hl.fold(
                            lambda i, j: hl.struct(
                                raw=i.raw + j.raw, adj=i.adj + j.adj
                            ),
                            x[0],
                            x[1:],
                        ),
                    )
                )
            ),
            ht.gt_counts.values(),
        )
    )

    ht = ht.checkpoint(f"{tmp_dir}/ht_sites_by_pop.ht", overwrite=True)
    ht = ht.key_by("locus1", "alleles1", "locus2", "alleles2")
    return ht.repartition(1000, shuffle=False)


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

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

Use --backend batch to run on Hail Query-on-Batch instead of Spark (local or
Dataproc). Requires hailctl auth login and hailctl config set batch/remote_tmpdir,
batch/billing_project, and query/backend batch (or pass --backend batch).
See https://hail.is/docs/0.2/cloud/query_on_batch.html.
"""

import argparse
import logging
import os
import tempfile
import timeit
from typing import Optional

from packaging import version


import hail as hl
from gnomad.resources.grch38.reference_data import clinvar, gencode
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.filtering import filter_to_clinvar_pathogenic
from gnomad.utils.vep import CSQ_ORDER, filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.annotations import get_insilico_predictors
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds

from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_LEAST_CONSEQUENCE,
    DEFAULT_MAX_FREQ,
    DEFAULT_EXON_DOWNSTREAM_PADDING,
    DEFAULT_EXON_UPSTREAM_PADDING,
    DEFAULT_MIN_PANGOLIN,
    DEFAULT_MIN_SPLICE_AI,
    DEFAULT_TMP_DIR,
    TEST_INTERVALS,
    _get_output_postfix,
    get_variant_filter_ht,
    get_variant_pair_resources,
)
from gnomad_chets.v4.utils import (
    calculate_partitions_by_size,
    compute_v2_independent_set,
    filter_for_testing,
)

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
        v1 <= v2,
        hl.struct(v1=v1, v2=v2),
        hl.struct(v1=v2, v2=v1),
    )


def create_variant_pair_ht(
    mt: hl.MatrixTable,
    vep_ht: hl.Table,
) -> hl.Table:
    """
    Create a Hail Table of unique ordered variant pairs per sample per gene.

    :param mt: MatrixTable with filtered variant data.
    :param vep_ht: VEP Table with gene_id annotation. Must be keyed by 'locus' and
        'alleles'.
    :return: Hail Table keyed by locus2, alleles2, locus1, alleles1 with one distinct
        row per unique variant pair.
    """
    n_partitions = vep_ht.n_partitions()
    mt = mt.add_row_index("variant_idx")
    variant_index_ht = mt.rows().key_by("variant_idx").cache()

    mt = mt.annotate_rows(gene_id=vep_ht[mt.locus, mt.alleles].gene_id)

    # Note: We do it this way because a row grouping by gene_id results in an 
    # aggregation with one gene per partition, and this can lead to memory issues.
    # Convert to entries table and explode on gene_id so each variant-gene combination
    # is a separate row.
    ht = mt.select_cols().select_rows("variant_idx", "gene_id").entries()
    ht = ht.filter(ht.GT.is_non_ref())
    ht = ht.explode("gene_id")

    # Group by gene and sample, collecting unique variants per gene/sample.
    # Using collect_as_set ensures each variant appears only once per gene/sample.
    ht = ht.group_by("gene_id", "s").aggregate(
        variants=hl.array(hl.agg.collect_as_set(ht.variant_idx))
    )

    # Filter to samples with at least 2 variants (needed to form pairs).
    ht = ht.filter(ht.variants.length() >= 2)
    ht = ht.checkpoint(
        hl.utils.new_temp_file("create_variant_pair_ht.gene_sample_grouped", "ht")
    )
    #ht = hl.read_table(
    #    "gs://gnomad-tmp-4day/create_variant_pair_ht.gene_sample_grouped-67uCBVffY8RCWc4MZYjc0A.ht"
    #)

    # Generate all ordered pairs of variants within each gene/sample.
    # The nested flatmap/map creates all combinations (i, j) where i < j, ensuring
    # each pair is created exactly once.
    ht = ht.annotate(
        pairs=(
            hl.range(0, hl.len(ht.variants)).flatmap(
                lambda i1: (
                    hl.range(i1 + 1, hl.len(ht.variants)).map(
                        lambda i2: _get_ordered_vp_struct(
                            ht.variants[i1], ht.variants[i2]
                        )
                    )
                )
            )
        )
    )

    # Explode pairs.
    ht = ht.explode("pairs")

    # Key by variant pair and select distinct pairs.
    # Use new shuffle method for apply models to prevent shuffle errors.
    hl._set_flags(use_new_shuffle="1")
    ht = ht.group_by(v1=ht.pairs.v1, v2=ht.pairs.v2).aggregate(
        gene_id=hl.agg.collect_as_set(ht.gene_id)
    )
    # Restore partition count; group_by shuffle often coalesces to few partitions
    # (e.g. spark.sql.shuffle.partitions=24), which would carry through to the
    # written variant pair table and downstream steps.
    ht = ht.repartition(n_partitions, shuffle=True).cache()
    hl._set_flags(use_new_shuffle=None)

    # Add a unique index id to each variant pair.
    # This is used later so both variants can be annotated with genotype
    # info separately and then joined together by the common index. This helps with
    # performance issues observed when trying to annotate both variants with genotype
    # info simultaneously.
    ht = ht.add_index("vp_ht_idx").key_by("vp_ht_idx")

    variant_index_keyed_v1 = variant_index_ht[ht.v1]
    variant_index_keyed_v2 = variant_index_ht[ht.v2]
    ht = ht.select(
        "gene_id",
        locus1=variant_index_keyed_v1.locus,
        alleles1=variant_index_keyed_v1.alleles,
        locus2=variant_index_keyed_v2.locus,
        alleles2=variant_index_keyed_v2.alleles,
    )

    return ht


def create_variant_pair_filter_ht(vp_ht: hl.Table) -> hl.Table:
    """
    Create a filter Table for variant pairs (unique variants appearing in any pair).

    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2, alleles2.
    :return: Filter Table keyed by locus, alleles.
    """
    v1_ht = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1).select().distinct()
    v2_ht = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2).select().distinct()
    n_partitions = vp_ht.n_partitions()
    ht = (
        v1_ht.union(v2_ht)
        .distinct()
        .repartition(n_partitions, shuffle=True)
        .checkpoint(hl.utils.new_temp_file("create_dense_filtered_mt.variants", "ht"))
    )
    return ht


def _create_var_idx_ht(mt: hl.MatrixTable) -> hl.Table:
    """
    Assign a unique int64 ``var_idx`` to each variant in ``mt``.

    Integer keys are dramatically faster for joins than ``(locus, alleles)`` struct
    keys and let us co-partition the encoded genotype Table and variant pair Table
    on the same int64 space.

    :param mt: MatrixTable whose row key is ``(locus, alleles)``.
    :return: Table keyed by ``(locus, alleles)`` with a ``var_idx`` field.
    """
    return mt.rows().add_index("var_idx").select("var_idx")


def _encode_genotype_sets_by_var_idx(
    mt: hl.MatrixTable,
    var_idx_ht: hl.Table,
) -> hl.Table:
    """
    Encode genotypes as per-variant sample sets keyed by ``var_idx``.

    For each variant, produces sets of sample indices by genotype category for
    both the raw and adj call, plus an ``all_samples`` set used to derive
    hom-ref counts by inclusion/exclusion downstream.

    Raw/adj genotypes are encoded as:

        - missing (None) = hom-ref (space saving)
        - 0 = missing data (no GT call, or failed adj for ``adj_gt``)
        - 1 = het
        - 2 = hom-var

    Output schema:

        ----------------------------------------
        Global fields:
            'samples': array<struct { s: str }>
        ----------------------------------------
        Row fields:
            'v_idx': int64
            'all_samples': set<int32>     # smaller of data vs complement set
            'all_samples_is_complement': bool
            'n_with_data': int32          # count of samples with non-hom-ref data
            'raw_het': set<int32>
            'raw_hv': set<int32>
            'adj_het': set<int32>
            'adj_hv': set<int32>
        ----------------------------------------
        Key: ['v_idx']
        ----------------------------------------

    :param mt: MatrixTable with variant data. Row key must be ``(locus, alleles)``.
    :param var_idx_ht: Var-idx Table from :func:`_create_var_idx_ht` keyed by
        ``(locus, alleles)`` with a ``var_idx`` field.
    :return: Table keyed by ``v_idx`` with per-variant sample sets.
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )
    adj_gt_count_expr = hl.if_else(
        get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD), gt_count_expr, 0, missing_false=True
    )

    mt = mt.select_entries(raw_gt=gt_count_expr, adj_gt=adj_gt_count_expr)
    ht = mt.localize_entries("_entries", "samples")

    # Build per-variant sample-index sets for each genotype category.
    # To save space, all_samples stores the smaller of "samples with data"
    # vs "samples without data" (complement). A boolean flag and n_with_data
    # let downstream code resolve which form is stored without touching the
    # set itself.
    gt = hl.enumerate(ht._entries)
    with_data = gt.filter(
        lambda x: hl.is_defined(x[1].raw_gt) | hl.is_defined(x[1].adj_gt)
    )
    without_data = gt.filter(
        lambda x: hl.is_missing(x[1].raw_gt) & hl.is_missing(x[1].adj_gt)
    )
    n_with = with_data.length()
    use_complement = n_with > without_data.length()
    ht = ht.select(
        all_samples=hl.if_else(
            use_complement,
            hl.set(without_data.map(lambda x: x[0])),
            hl.set(with_data.map(lambda x: x[0])),
        ),
        all_samples_is_complement=use_complement,
        n_with_data=hl.int32(n_with),
        raw_het=hl.set(
            with_data.filter(lambda x: x[1].raw_gt == 1).map(lambda x: x[0])
        ),
        raw_hv=hl.set(
            with_data.filter(lambda x: x[1].raw_gt == 2).map(lambda x: x[0])
        ),
        adj_het=hl.set(
            with_data.filter(lambda x: x[1].adj_gt == 1).map(lambda x: x[0])
        ),
        adj_hv=hl.set(
            with_data.filter(lambda x: x[1].adj_gt == 2).map(lambda x: x[0])
        ),
    )

    # Rekey by var_idx and drop the locus/alleles fields to shrink row size.
    ht = ht.annotate(v_idx=var_idx_ht[ht.locus, ht.alleles].var_idx)
    return ht.key_by("v_idx").drop("locus", "alleles")



def _count_from_sets(
    v1_het: hl.expr.SetExpression,
    v1_hv: hl.expr.SetExpression,
    v1_all: hl.expr.SetExpression,
    v1_n: hl.expr.Int32Expression,
    v1_is_complement: hl.expr.BooleanExpression,
    v2_het: hl.expr.SetExpression,
    v2_hv: hl.expr.SetExpression,
    v2_all: hl.expr.SetExpression,
    v2_n: hl.expr.Int32Expression,
    v2_is_complement: hl.expr.BooleanExpression,
    n_samples: hl.expr.Int32Expression,
) -> hl.expr.ArrayExpression:
    """
    Compute 9-element genotype count array from per-variant sample sets.

    ``v_all`` stores either the "samples with data" set (positive) or its
    complement ("samples without data"), indicated by ``v_is_complement``.
    ``v_n`` is always the count of samples *with data* regardless of which
    form is stored.

    All branching on the complement flag operates on ints (set lengths and
    intersection lengths), never on sets themselves.

    Count array layout: ``[AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]``
    where A/a = v1 ref/alt, B/b = v2 ref/alt.

    :param v1_het: Sample set for v1 het genotype.
    :param v1_hv: Sample set for v1 hom-var genotype.
    :param v1_all: Stored set for v1 (positive or complement).
    :param v1_n: Number of samples with non-hom-ref data at v1.
    :param v1_is_complement: Whether ``v1_all`` is the complement form.
    :param v2_het: Sample set for v2 het genotype.
    :param v2_hv: Sample set for v2 hom-var genotype.
    :param v2_all: Stored set for v2 (positive or complement).
    :param v2_n: Number of samples with non-hom-ref data at v2.
    :param v2_is_complement: Whether ``v2_all`` is the complement form.
    :param n_samples: Total number of samples in the cohort.
    :return: 9-element genotype count array.
    """
    # Compute |data1 ∩ data2| from the stored intersection.
    # Let S = stored intersection length, |A| = v1_all.length(), etc.
    #   both positive:  |d1 ∩ d2| = S
    #   both complement: |d1 ∩ d2| = n - |A| - |B| + S
    #   v1 pos, v2 comp: |d1 ∩ d2| = |A| - S  (pos ∩ comp = pos \ comp)
    #   v1 comp, v2 pos: |d1 ∩ d2| = |B| - S
    stored_isect = v1_all.intersection(v2_all).length()
    v1_sz = v1_all.length()
    v2_sz = v2_all.length()
    v1v2_data_overlap = (
        hl.case()
        .when(~v1_is_complement & ~v2_is_complement, stored_isect)
        .when(v1_is_complement & v2_is_complement,
              n_samples - v1_sz - v2_sz + stored_isect)
        .when(~v1_is_complement & v2_is_complement, v1_sz - stored_isect)
        .default(v2_sz - stored_isect)
    )

    # AABB: both truly hom-ref.
    hom_ref_both = n_samples - v1_n - v2_n + v1v2_data_overlap

    # Non-ref x non-ref bins (het/hv sets are always positive).
    het_het = v1_het.intersection(v2_het).length()
    het_hv = v1_het.intersection(v2_hv).length()
    hv_het = v1_hv.intersection(v2_het).length()
    hv_hv = v1_hv.intersection(v2_hv).length()

    # Edge bins: |{v_gt at this variant, hom-ref at other}|.
    # = |v_gt| - |v_gt ∩ other_data|.
    # |v_gt ∩ other_data| depends on whether other_all is complement:
    #   positive:   |v_gt ∩ other_all|
    #   complement: |v_gt| - |v_gt ∩ other_all|  (= |v_gt \ other_comp|)
    def _in_other_data(v_gt, other_all, other_is_comp):
        raw_isect = v_gt.intersection(other_all).length()
        return hl.if_else(other_is_comp, v_gt.length() - raw_isect, raw_isect)

    v1_het_in_v2 = _in_other_data(v1_het, v2_all, v2_is_complement)
    v1_hv_in_v2 = _in_other_data(v1_hv, v2_all, v2_is_complement)
    v2_het_in_v1 = _in_other_data(v2_het, v1_all, v1_is_complement)
    v2_hv_in_v1 = _in_other_data(v2_hv, v1_all, v1_is_complement)

    return hl.array([
        hom_ref_both,                       # AABB
        v2_het.length() - v2_het_in_v1,     # AABb
        v2_hv.length() - v2_hv_in_v1,       # AAbb
        v1_het.length() - v1_het_in_v2,     # AaBB
        het_het,                             # AaBb
        het_hv,                              # Aabb
        v1_hv.length() - v1_hv_in_v2,       # aaBB
        hv_het,                              # aaBb
        hv_hv,                               # aabb
    ])


def create_variant_pair_genotype_counts(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Compute variant pair genotype counts directly from MT and pair list.

    Pipeline:

        1. Assign an int64 ``var_idx`` to each variant in ``mt``.
        2. Encode genotypes as per-variant sample sets (het/hv/all for raw and
           adj), keyed by ``v_idx``.
        3. Rekey ``vp_ht`` by ``v1_idx`` with ``v2_idx`` field (drops locus/alleles
           for the inner join; they're re-attached at the end for output).
        4. Compute per-variant pair-participation counts to weight partition
           boundaries.
        5. Partition both tables on identical size-balanced intervals over
           ``v_idx``, so the v1 lookup is a partition-local zip-join.
        6. Look up v1 and v2 sets for each pair and compute 9-element count
           arrays (raw and adj) inline via set intersections.

    :param mt: Dense filtered MatrixTable with variant data. Row key must be
        ``(locus, alleles)``.
    :param vp_ht: Table of variant pairs with fields ``locus1, alleles1, locus2,
        alleles2`` (other fields are ignored).
    :param max_join_partitions: Upper bound on the number of partitions used for
        the co-partitioned join.
    :return: Variant pair Table keyed by ``(locus1, alleles1, locus2, alleles2)``
        with ``gt_counts_raw`` and ``gt_counts_adj`` fields.
    """
    # (1) var_idx table: int64 index per variant.
    #var_idx_path = hl.utils.new_temp_file("var_idx", "ht")
    #_create_var_idx_ht(mt).write(var_idx_path, overwrite=True)
    var_idx_path = "gs://gnomad-tmp-4day/var_idx-2sCNgNj3LQQW4sNMtIUdxq.ht"
    var_idx_ht = hl.read_table(var_idx_path)

    # (2) Encoded genotype sets keyed by v_idx.
    #encoded_path = hl.utils.new_temp_file("encoded_gt_sets_by_var_idx", "ht")
    #_encode_genotype_sets_by_var_idx(mt, var_idx_ht).write(
    #    encoded_path, overwrite=True
    #)
    encoded_path = "gs://gnomad-tmp-4day/encoded_gt_sets_by_var_idx-unutLSbNBFBbr3Ham082aT.ht"
    encoded_gt_ht = hl.read_table(encoded_path)

    # (3) Rekey vp_ht by (v1_idx, v2_idx).
    #vp_idx_path = hl.utils.new_temp_file("vp_by_idx", "ht")
    #vp_ht = vp_ht.annotate(
    #    v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
    #    v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    #)
    #vp_ht = vp_ht.key_by("v1_idx", "v2_idx").select(
    #    "locus1", "alleles1", "locus2", "alleles2"
    #)
    #vp_ht.write(vp_idx_path, overwrite=True)
    vp_idx_path = "gs://gnomad-tmp-4day/vp_by_v1_idx-WX06VNXUP7fWwCjjW9CB11.ht"
    vp_ht = hl.read_table(vp_idx_path)

    # (4) Size-balanced partitions computed from encoded_gt_ht (keyed by
    # v_idx), weighted by how many pairs each variant appears in as v1.
    # Using encoded_gt_ht's single-field key means the same intervals can
    # be applied to both tables: encoded_gt_ht directly, and vp_ht via
    # prefix key match on v1_idx. This makes the v1 lookup a partition-local
    # zip-join with no shuffle.
    #v1_counts = (
    #    vp_ht.key_by()
    #    .select("v1_idx")
    #    .group_by("v1_idx")
    #    .aggregate(_n_pairs=hl.int32(hl.agg.count()))
    #)
    #v1_counts = v1_counts.key_by(v_idx=v1_counts.v1_idx).drop("v1_idx")
    #v1_counts = v1_counts.checkpoint(
    #    hl.utils.new_temp_file("v1_pair_counts", "ht")
    #)
    v1_counts = hl.read_table("gs://gnomad-tmp-30day/v1_pair_counts-NOTDhKG71rLQG5M8Ysac4X.ht")
    #n_join_partitions = min(encoded_gt_ht.n_partitions() * 3, max_join_partitions)
    n_join_partitions = 1000
    partition_intervals = calculate_partitions_by_size(
        encoded_gt_ht,
        n_join_partitions,
        size_field="n_with_data",
        weight_ht=v1_counts,
        weight_field="_n_pairs",
    )

    # Re-read both tables with the same v_idx intervals. vp_ht's (v1_idx,
    # v2_idx) key prefix-matches the v_idx intervals on v1_idx, so Hail
    # co-partitions them and the v1 annotate is a partition-local zip-join.
    encoded_gt_ht = hl.read_table(encoded_path, _intervals=partition_intervals)
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())
    #vp_ht = hl.read_table(vp_idx_path, _intervals=partition_intervals)

    # (5) Both lookups + counts in a single expression. The v1 lookup is a
    # partition-local zip-join (no shuffle) since vp_ht and encoded_gt_ht
    # are co-partitioned on v1_idx/v_idx. The v2 lookup shuffles only the
    # ~10 GB encoded_gt_ht (not the 5M-row vp_ht), and the set
    # intersections produce ints — no large data is ever materialized on
    # vp_ht rows.
    vp_ht = hl.read_table(vp_idx_path, _intervals=partition_intervals)
    v1 = encoded_gt_ht[vp_ht.v1_idx]
    v2 = encoded_gt_ht[vp_ht.v2_idx]
    vp_ht = vp_ht.select(
        "locus1",
        "alleles1",
        "locus2",
        "alleles2",
        gt_counts_raw=_count_from_sets(
            v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
        gt_counts_adj=_count_from_sets(
            v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
    )

    return vp_ht.key_by("locus1", "alleles1", "locus2", "alleles2")


def create_variant_pair_genotype_counts_v2(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Compute variant pair genotype counts with two shuffle-free zip-joins.

    Pipeline:

        1. Assign int64 ``var_idx`` to each variant.
        2. Encode genotypes as per-variant sample sets keyed by ``v_idx``.
        3. Collect v2 indices per v1 into an array, keyed by ``v1_idx``.
        4. Partition by v1 weight, co-partition with encoded_gt_ht, annotate v1 sets
           via zip-join (no shuffle).
        5. Add random key, explode v2s, rekey by ``(rand_key, v2_idx)``.
        6. Re-partition by v2 weight, co-partition with encoded_gt_ht, annotate v2
           sets via zip-join (no shuffle).
        7. Compute 9-element count arrays inline.

    :param mt: Dense filtered MatrixTable.
    :param vp_ht: Variant pair list Table.
    :param max_join_partitions: Upper bound on partition count.
    :return: Variant pair Table with ``gt_counts_raw`` and ``gt_counts_adj``.
    """
    # --- (1) var_idx ---
    #var_idx_path = hl.utils.new_temp_file("var_idx", "ht")
    #_create_var_idx_ht(mt).write(var_idx_path, overwrite=True)
    var_idx_path = "gs://gnomad-tmp-4day/var_idx-2sCNgNj3LQQW4sNMtIUdxq.ht"
    var_idx_ht = hl.read_table(var_idx_path)

    # --- (2) Encoded genotype sets keyed by v_idx ---
    #encoded_path = hl.utils.new_temp_file("encoded_gt_sets_by_var_idx", "ht")
    #_encode_genotype_sets_by_var_idx(mt, var_idx_ht).write(
    #    encoded_path, overwrite=True
    #)
    encoded_path = "gs://gnomad-tmp-4day/encoded_gt_sets_by_var_idx-unutLSbNBFBbr3Ham082aT.ht"
    encoded_gt_ht = hl.read_table(encoded_path)
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())

    # --- (3) Collect v2s per v1 ---
    # Rekey vp_ht to (v1_idx) and collect v2_idx + locus/alleles into arrays.
    # = vp_ht.annotate(
    #    v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
    #    v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    #)
    #vp_collected_path = hl.utils.new_temp_file("vp_collected_by_v1", "ht")
    vp_collected_path = "gs://gnomad-tmp-30day/vp_collected_by_v1-p56ZkYzIWbWj2efYZDUfD8.ht"
    #vp_collected = vp_ht.group_by("v1_idx").aggregate(
    #    pairs=hl.agg.collect(
    #        hl.struct(
    #            v2_idx=vp_ht.v2_idx,
    #            locus1=vp_ht.locus1,
    #            alleles1=vp_ht.alleles1,
    #            locus2=vp_ht.locus2,
    #            alleles2=vp_ht.alleles2,
    #        )
    #    )
    #)
    #vp_collected = vp_collected.key_by(v_idx=vp_collected.v1_idx).drop("v1_idx")
    #vp_collected.write(vp_collected_path, overwrite=True)
    vp_collected = hl.read_table(vp_collected_path)

    # --- (4) Partition by v1 weight, annotate v1 sets via zip-join ---
    n_join_partitions = min(encoded_gt_ht.n_partitions() * 3, max_join_partitions)
    v1_intervals_path = hl.utils.new_temp_file("v1_partition_intervals", "he")
    v1_partition_intervals = calculate_partitions_by_size(
        encoded_gt_ht,
        n_join_partitions,
        size_field="n_with_data",
        weight_ht=vp_collected,
        weight_field="pairs",
    )
    hl.experimental.write_expression(
        hl.literal(v1_partition_intervals), v1_intervals_path
    )
    # To reuse: v1_partition_intervals = hl.eval(
    #     hl.experimental.read_expression(v1_intervals_path)
    # )
    encoded_gt_ht_v1 = hl.read_table(encoded_path, _intervals=v1_partition_intervals)
    vp_collected = hl.read_table(vp_collected_path, _intervals=v1_partition_intervals)

    # Annotate v1 sets — zip-join, no shuffle.
    #vp_with_v1_path = hl.utils.new_temp_file("vp_with_v1_sets", "ht")
    vp_with_v1_path = "gs://gnomad-tmp-30day/vp_with_v1_sets-aEQoMrAZMLlOWlz9yBTHN7.ht"
    #vp_collected = vp_collected.annotate(v1=encoded_gt_ht_v1[vp_collected.v_idx])
    #vp_collected.write(vp_with_v1_path, overwrite=True)
    vp_collected = hl.read_table(vp_with_v1_path)

    # --- (5) Explode pairs, look up v2 sets inline, compute counts ---
    # Explode pairs so each row is one (v1, v2) pair with v1 sets attached.
    # The v2 lookup is an indexed read against the full encoded_gt_ht (~10 GB
    # shuffle). Since v1 sets are already on each row, only the small gt
    # table shuffles — not the 5M-row vp table. No rekey needed.
    vp_exploded = vp_collected.explode("pairs")
    vp_exploded = vp_exploded.transmute(
        rand_key=hl.rand_int64(100000000),
        v2_idx=vp_exploded.pairs.v2_idx,
        locus1=vp_exploded.pairs.locus1,
        alleles1=vp_exploded.pairs.alleles1,
        locus2=vp_exploded.pairs.locus2,
        alleles2=vp_exploded.pairs.alleles2,
    ).checkpoint(hl.utils.new_temp_file("vp_exploded", "ht"))

    ###encoded_gt_ht_v2 = hl.read_table(encoded_path)
    ###v1 = vp_exploded.v1
    ###v2 = encoded_gt_ht_v2[vp_exploded.v2_idx]
    ###vp_exploded = vp_exploded.select(
    ###    "locus1",
    ###    "alleles1",
    ###    "locus2",
    ###    "alleles2",
    ###    gt_counts_raw=_count_from_sets(
    ###        v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
    ###        v1.all_samples_is_complement,
    ###        v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
    ###        v2.all_samples_is_complement,
    ###        n_samples,
    ###    ),
    ###    gt_counts_adj=_count_from_sets(
    ###        v1.all_samples_is_complement,
    ###        v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
    ###        v2.all_samples_is_complement,
    ###        n_samples,
    ###    ),
    ###)

    return vp_exploded #.key_by("locus1", "alleles1", "locus2", "alleles2")


def create_variant_pair_genotype_counts_v3(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Compute variant pair genotype counts with V1/V2 role assignment.

    Assigns each variant globally to V1 (zip-join, free) or V2 (shuffle)
    via a greedy independent set algorithm that minimizes total V2 weight.
    Then:

        1. Assign var_idx, encode genotype sets (reuse from earlier runs).
        2. Collect variants + edges to the driver, run greedy V2 assignment.
        3. Restructure pair table: V2 variant always on v2 side. For V1-V1
           pairs, the smaller variant goes on v2 side.
        4. Build a small V2-only genotype table for the shuffle lookup.
        5. Collect v2s per v1, co-partition with full encoded_gt_ht, annotate
           v1 sets via zip-join (free).
        6. Explode, look up v2 sets against the small V2 table (small shuffle),
           compute counts inline.

    :param mt: Dense filtered MatrixTable.
    :param vp_ht: Variant pair list Table.
    :param max_join_partitions: Upper bound on partition count.
    :return: Variant pair Table with ``gt_counts_raw`` and ``gt_counts_adj``.
    """
    # --- (1) var_idx + encoded genotype sets ---
    ##var_idx_path = hl.utils.new_temp_file("var_idx", "ht")
    ##_create_var_idx_ht(mt).write(var_idx_path, overwrite=True)
    var_idx_path = "gs://gnomad-tmp-30day/var_idx-d1OBqRwTN80wMYwyLGoFhz.ht"
    var_idx_ht = hl.read_table(var_idx_path)

    ##encoded_path = hl.utils.new_temp_file("encoded_gt_sets_by_var_idx", "ht")
    ##_encode_genotype_sets_by_var_idx(mt, var_idx_ht).write(
    ##    encoded_path, overwrite=True
    ##)
    encoded_path = "gs://gnomad-tmp-30day/encoded_gt_sets_by_var_idx-QzikI9hqzQhpoxO56Gv0QF.ht"
    encoded_gt_ht = hl.read_table(encoded_path)
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())

    # --- (2) V2 assignment on driver ---
    # Collect variant weights.
    logger.info("Collecting variants and edges for V2 assignment...")
    ##variant_data = encoded_gt_ht.select("n_with_data").collect()
    ##variants = [(row.v_idx, row.n_with_data) for row in variant_data]

    # Rekey vp_ht to get v1_idx/v2_idx for edge collection.
    ##vp_ht = vp_ht.annotate(
    ##    _v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
    ##    _v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    ##)
    ##edge_data = vp_ht.key_by().select("_v1_idx", "_v2_idx").collect()
    ##edges = [(row._v1_idx, row._v2_idx) for row in edge_data]

    ##v2_set = compute_v2_independent_set(variants, edges)
    ##logger.info("V2 set total n_with_data: %d", sum(n for v, n in variants if v in v2_set))

    # --- (3) Restructure pair table ---
    # For each pair, decide which side is v1/v2:
    #   - If one variant is V2 and the other is V1: V2 goes to v2 side
    #   - If both are V1: smaller n_with_data goes to v2 side
    #   - Both V2: shouldn't happen (independent set), but handle gracefully
    ##v2_broadcast = hl.literal(v2_set, dtype=hl.tset(hl.tint64))
    ##n_with_data_map = hl.literal(
    ##    {v: n for v, n in variants}, dtype=hl.tdict(hl.tint64, hl.tint32)
    ##)

    ##vp_ht = vp_ht.annotate(
    ##    _v1_is_v2=v2_broadcast.contains(vp_ht._v1_idx),
    ##    _v2_is_v2=v2_broadcast.contains(vp_ht._v2_idx),
    ##)
    # Swap when: v1 is in V2 and v2 is not, OR both are V1 and v1 has smaller
    # n_with_data (we want the smaller one on v2 side for V1-V1 pairs).
    ##swap = (
    ##    (vp_ht._v1_is_v2 & ~vp_ht._v2_is_v2)
    ##    | (
    ##        ~vp_ht._v1_is_v2
    ##        & ~vp_ht._v2_is_v2
    ##        & (n_with_data_map.get(vp_ht._v1_idx) < n_with_data_map.get(vp_ht._v2_idx))
    ##    )
    ##)
    ##vp_ht = vp_ht.select(
    ##    v1_idx=hl.if_else(swap, vp_ht._v2_idx, vp_ht._v1_idx),
    ##    v2_idx=hl.if_else(swap, vp_ht._v1_idx, vp_ht._v2_idx),
    ##    locus1=hl.if_else(swap, vp_ht.locus2, vp_ht.locus1),
    ##    alleles1=hl.if_else(swap, vp_ht.alleles2, vp_ht.alleles1),
    ##    locus2=hl.if_else(swap, vp_ht.locus1, vp_ht.locus2),
    ##    alleles2=hl.if_else(swap, vp_ht.alleles1, vp_ht.alleles2),
    ##)

    # --- (4) Build small genotype table for all v2-side variants ---
    # This includes V2 set variants AND V1 variants that appear on the v2
    # side of V1-V1 pairs. Much smaller than the full encoded_gt_ht since
    # the heavy/common variants are always on the v1 side.
    # Repartition with calculate_partitions_by_size to avoid hot partitions
    # in the v2 shuffle lookup.
    #v2_variants = vp_ht.key_by(v_idx=vp_ht.v2_idx).select().distinct()
    #v2_gt_ht = encoded_gt_ht.semi_join(v2_variants)
    #v2_gt_unpartitioned_path = hl.utils.new_temp_file("v2_side_gt_unpartitioned", "ht")
    #v2_gt_ht.write(v2_gt_unpartitioned_path, overwrite=True)
    v2_gt_unpartitioned_path = "gs://gnomad-tmp-30day/v2_side_encoded_gt_sets-yH2UQj4jmgwqHNmZ4tf4wA.ht"
    v2_gt_ht = hl.read_table(v2_gt_unpartitioned_path)

    #v2_partition_intervals = calculate_partitions_by_size(
    #    v2_gt_ht,
    #    min(v2_gt_ht.n_partitions() * 3, max_join_partitions),
    #    size_field="n_with_data",
    #)
    #v2_gt_path = hl.utils.new_temp_file("v2_side_encoded_gt_sets_repartitioned", "ht")
    v2_gt_path = "gs://gnomad-tmp-30day/v2_side_encoded_gt_sets_repartitioned-LG5lWceRJwjy0fC3hhpsnJ.ht"
    #hl.read_table(v2_gt_unpartitioned_path, _intervals=v2_partition_intervals).write(
    #    v2_gt_path, overwrite=True
    #)
    v2_gt_ht = hl.read_table(v2_gt_path)
    logger.info(
        "V2-side genotype table: %d variants (full table: %d)",
        v2_gt_ht.count(),
        encoded_gt_ht.count(),
    )

    # --- (5) Collect v2s per v1, co-partition, annotate v1 via zip-join ---
    ##vp_collected_path = hl.utils.new_temp_file("vp_collected_by_v1_v3", "ht")
    vp_collected_path = "gs://gnomad-tmp-30day/vp_collected_by_v1_v3-wrjT4wUDhJIjwCJ7O44UKS.ht"
    ##vp_collected = vp_ht.group_by(vp_ht.v1_idx).aggregate(
    ##    pairs=hl.agg.collect(
    ##        hl.struct(
    ##            v2_idx=vp_ht.v2_idx,
    ##            locus1=vp_ht.locus1,
    ##            alleles1=vp_ht.alleles1,
    ##            locus2=vp_ht.locus2,
    ##            alleles2=vp_ht.alleles2,
    ##        )
    ##    )
    ##)
    ##vp_collected = vp_collected.key_by(v_idx=vp_collected.v1_idx).drop("v1_idx")
    ##vp_pre_chunk_path = hl.utils.new_temp_file("vp_collected_pre_chunk_v3", "ht")
    vp_pre_chunk_path = "gs://gnomad-tmp-30day/vp_collected_pre_chunk_v3-kDwZW2vDGQFn5h4YYBz7pR.ht"
    ##vp_collected = vp_collected.checkpoint(vp_pre_chunk_path)
    vp_collected = hl.read_table(vp_pre_chunk_path)

    # Compute partition intervals BEFORE chunking, so the weight reflects
    # the total pair count per v1 (not just the first chunk).
    n_join_partitions = min(encoded_gt_ht.n_partitions() * 3, max_join_partitions)
    v1_intervals_path = hl.utils.new_temp_file("v1_partition_intervals_v3", "he")
    v1_partition_intervals = calculate_partitions_by_size(
        encoded_gt_ht,
        n_join_partitions,
        size_field="n_with_data",
        weight_ht=vp_collected,
        weight_field="pairs",
    )
    hl.experimental.write_expression(
        hl.literal(v1_partition_intervals), v1_intervals_path
    )

    # Split large pairs arrays into chunks to prevent any single row from
    # producing a huge partition after explode. Each chunk becomes its own
    # row with the same v_idx key.
    ##PAIR_BATCH_SIZE = 1000
    ##vp_collected = vp_collected.annotate(
    ##    pairs=hl.range(
    ##        0, hl.len(vp_collected.pairs), PAIR_BATCH_SIZE
    ##    ).map(lambda i: vp_collected.pairs[i:i + PAIR_BATCH_SIZE])
    ##)
    ##vp_collected = vp_collected.explode("pairs")

    ##vp_collected.write(vp_collected_path, overwrite=True)
    vp_collected = hl.read_table(vp_collected_path)

    encoded_gt_ht_v1 = hl.read_table(encoded_path, _intervals=v1_partition_intervals)
    vp_collected = hl.read_table(vp_collected_path, _intervals=v1_partition_intervals)

    ##vp_with_v1_path = hl.utils.new_temp_file("vp_with_v1_sets_v3", "ht")
    vp_with_v1_path = "gs://gnomad-tmp-30day/vp_with_v1_sets_v3-q3mqBBmgS1J7AmuxPkpCND.ht"
    ##vp_collected = vp_collected.annotate(v1=encoded_gt_ht_v1[vp_collected.v_idx])
    ##vp_collected.write(vp_with_v1_path, overwrite=True)
    vp_collected = hl.read_table(vp_with_v1_path)

    # --- (6) Explode, look up v2 from small table, compute counts ---
    vp_exploded = vp_collected.explode("pairs")
    vp_exploded = vp_exploded.transmute(
        v2_idx=vp_exploded.pairs.v2_idx,
        locus1=vp_exploded.pairs.locus1,
        alleles1=vp_exploded.pairs.alleles1,
        locus2=vp_exploded.pairs.locus2,
        alleles2=vp_exploded.pairs.alleles2,
    )

    # v2 lookup against the small v2-side table. Force a broadcast join so
    # Spark sends the 6 GB table to every executor's memory instead of
    # shuffling to disk (which overflows 40 GB secondary worker disks).
    hl.utils.java.Env.spark_session().conf.set(
        "spark.sql.autoBroadcastJoinThreshold",
        str(8 * 1024 * 1024 * 1024),  # 8 GB
    )
    v2_gt_ht = hl.read_table(v2_gt_path)
    v2 = v2_gt_ht[vp_exploded.v2_idx]
    v1 = vp_exploded.v1
    vp_exploded = vp_exploded.select(
        "locus1",
        "alleles1",
        "locus2",
        "alleles2",
        gt_counts_raw=_count_from_sets(
            v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
        gt_counts_adj=_count_from_sets(
            v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
    )

    return vp_exploded.key_by("locus1", "alleles1", "locus2", "alleles2")


def _compute_counts_for_subset(
    vp_subset: hl.Table,
    encoded_gt_ht: hl.Table,
    n_samples: hl.expr.Int32Expression,
    label: str,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Compute genotype counts for a subset of variant pairs using gene_idx
    co-partitioning (steps 4-6).

    Both v1 and v2 lookups are partition-local zip-joins against a gt table
    keyed by ``(gene_idx, v_idx)``.

    :param vp_subset: Restructured pair Table with gene_idx, v1_idx, v2_idx,
        locus1, alleles1, locus2, alleles2.
    :param encoded_gt_ht: Encoded genotype sets keyed by v_idx.
    :param n_samples: Total sample count expression.
    :param label: Label for temp file naming (e.g., "light" or "heavy").
    :param max_join_partitions: Upper bound on partition count.
    :return: Counts Table with gt_counts_raw and gt_counts_adj.
    """
    # --- Build gt table keyed by (gene_idx, v_idx) ---
    unkeyed_sub = vp_subset.key_by()
    v1_genes = unkeyed_sub.select(
        gene_idx=unkeyed_sub.gene_idx, v_idx=unkeyed_sub.v1_idx
    )
    v2_genes = unkeyed_sub.select(
        gene_idx=unkeyed_sub.gene_idx, v_idx=unkeyed_sub.v2_idx
    )
    variant_genes = (
        v1_genes.union(v2_genes).key_by("gene_idx", "v_idx").distinct()
    )
    variant_genes = variant_genes.key_by("v_idx")

    gt_by_gene_path = hl.utils.new_temp_file(
        f"gt_by_gene_idx_v_idx_{label}", "ht"
    )
    gt_by_gene = variant_genes.annotate(
        **encoded_gt_ht[variant_genes.v_idx]
    )
    gt_by_gene = gt_by_gene.key_by("gene_idx", "v_idx")
    gt_by_gene.write(gt_by_gene_path, overwrite=True)
    gt_by_gene = hl.read_table(gt_by_gene_path)
    logger.info("%s GT by (gene_idx, v_idx): %d rows", label, gt_by_gene.count())

    # --- Collect pairs per (gene_idx, v1), co-partition, annotate v1 ---
    vp_collected_path = hl.utils.new_temp_file(
        f"vp_collected_by_gene_v1_{label}", "ht"
    )
    vp_collected = vp_subset.group_by("gene_idx", "v1_idx").aggregate(
        pairs=hl.agg.collect(
            hl.struct(
                v2_idx=vp_subset.v2_idx,
                locus1=vp_subset.locus1,
                alleles1=vp_subset.alleles1,
                locus2=vp_subset.locus2,
                alleles2=vp_subset.alleles2,
            )
        )
    )
    vp_collected = vp_collected.key_by(
        gene_idx=vp_collected.gene_idx,
        v_idx=vp_collected.v1_idx,
    ).drop("v1_idx")

    PAIR_BATCH_SIZE = 1000
    vp_collected = vp_collected.annotate(
        pairs=hl.range(
            0, hl.len(vp_collected.pairs), PAIR_BATCH_SIZE
        ).map(lambda i: vp_collected.pairs[i:i + PAIR_BATCH_SIZE])
    )
    vp_collected = vp_collected.explode("pairs")
    vp_collected.write(vp_collected_path, overwrite=True)
    vp_collected = hl.read_table(vp_collected_path)

    n_join_partitions = min(
        gt_by_gene.n_partitions() * 3, max_join_partitions
    )
    partition_intervals = calculate_partitions_by_size(
        gt_by_gene,
        n_join_partitions,
        size_field="n_with_data",
        weight_ht=vp_collected,
        weight_field="pairs",
    )

    gt_by_gene_v1 = hl.read_table(
        gt_by_gene_path, _intervals=partition_intervals
    )
    vp_collected = hl.read_table(
        vp_collected_path, _intervals=partition_intervals
    )

    vp_with_v1_path = hl.utils.new_temp_file(
        f"vp_with_v1_sets_{label}", "ht"
    )
    vp_collected = vp_collected.annotate(
        v1=gt_by_gene_v1[vp_collected.gene_idx, vp_collected.v_idx]
    )
    vp_collected.write(vp_with_v1_path, overwrite=True)
    vp_collected = hl.read_table(vp_with_v1_path)

    # --- Explode, rekey by (gene_idx, v2_idx), v2 zip-join, counts ---
    vp_exploded = vp_collected.explode("pairs")
    vp_exploded = vp_exploded.transmute(
        v2_idx=vp_exploded.pairs.v2_idx,
        locus1=vp_exploded.pairs.locus1,
        alleles1=vp_exploded.pairs.alleles1,
        locus2=vp_exploded.pairs.locus2,
        alleles2=vp_exploded.pairs.alleles2,
    )

    vp_by_v2_path = hl.utils.new_temp_file(
        f"vp_by_gene_v2_{label}", "ht"
    )
    vp_exploded = vp_exploded.key_by("gene_idx", "v2_idx")
    vp_exploded.write(vp_by_v2_path, overwrite=True)

    gt_by_gene_v2 = hl.read_table(
        gt_by_gene_path, _intervals=partition_intervals
    )
    vp_exploded = hl.read_table(
        vp_by_v2_path, _intervals=partition_intervals
    )

    v1 = vp_exploded.v1
    v2 = gt_by_gene_v2[vp_exploded.gene_idx, vp_exploded.v2_idx]
    return vp_exploded.select(
        "locus1",
        "alleles1",
        "locus2",
        "alleles2",
        gt_counts_raw=_count_from_sets(
            v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
        gt_counts_adj=_count_from_sets(
            v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
    ).key_by("locus1", "alleles1", "locus2", "alleles2")


def prepare_genotype_count_inputs(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    output_dir: str,
    all_samples_len_threshold: int = 50000,
) -> None:
    """
    Prepare all intermediate tables for genotype count computation.

    Encodes genotypes, assigns V1/V2 roles, computes gene_idx, restructures
    the pair table, and splits into light/heavy subsets. Writes all outputs
    to ``output_dir`` so subsequent steps can read them on different clusters.

    :param mt: Dense filtered MatrixTable.
    :param vp_ht: Variant pair list Table with gene_id field.
    :param output_dir: GCS directory for intermediate outputs.
    :param all_samples_len_threshold: Threshold for splitting light/heavy.
    """
    # --- (1) var_idx + encoded genotype sets ---
    var_idx_path = f"{output_dir}/var_idx.ht"
    _create_var_idx_ht(mt).write(var_idx_path, overwrite=True)
    var_idx_ht = hl.read_table(var_idx_path)

    encoded_path = f"{output_dir}/encoded_gt_sets_by_var_idx.ht"
    _encode_genotype_sets_by_var_idx(mt, var_idx_ht).write(
        encoded_path, overwrite=True
    )
    encoded_gt_ht = hl.read_table(encoded_path)

    # --- (2) V2 assignment on driver ---
    logger.info("Collecting variants and edges for V2 assignment...")
    variant_data = encoded_gt_ht.select("n_with_data").collect()
    variants = [(row.v_idx, row.n_with_data) for row in variant_data]

    vp_ht = vp_ht.annotate(
        _v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        _v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    edge_data = vp_ht.key_by().select("_v1_idx", "_v2_idx").collect()
    edges = [(row._v1_idx, row._v2_idx) for row in edge_data]

    v2_set = compute_v2_independent_set(variants, edges)
    logger.info("V2 set total n_with_data: %d", sum(n for v, n in variants if v in v2_set))

    # --- (3) Restructure pair table with gene_idx ---
    vp_ht = vp_ht.annotate(_gene_id=hl.array(vp_ht.gene_id)[0])

    gene_idx_ht = vp_ht.key_by().select("_gene_id", "_v1_idx", "_v2_idx")
    gene_idx_ht = gene_idx_ht.annotate(
        _min_v=hl.min(gene_idx_ht._v1_idx, gene_idx_ht._v2_idx)
    )
    gene_idx_ht = gene_idx_ht.group_by("_gene_id").aggregate(
        gene_idx=hl.agg.min(gene_idx_ht._min_v)
    )
    gene_idx_ht = gene_idx_ht.key_by("_gene_id")

    v2_broadcast = hl.literal(v2_set, dtype=hl.tset(hl.tint64))
    n_with_data_map = hl.literal(
        {v: n for v, n in variants}, dtype=hl.tdict(hl.tint64, hl.tint32)
    )
    vp_ht = vp_ht.annotate(
        _v1_is_v2=v2_broadcast.contains(vp_ht._v1_idx),
        _v2_is_v2=v2_broadcast.contains(vp_ht._v2_idx),
        _gene_idx=gene_idx_ht[vp_ht._gene_id].gene_idx,
    )
    swap = (
        (vp_ht._v1_is_v2 & ~vp_ht._v2_is_v2)
        | (
            ~vp_ht._v1_is_v2
            & ~vp_ht._v2_is_v2
            & (n_with_data_map.get(vp_ht._v1_idx) < n_with_data_map.get(vp_ht._v2_idx))
        )
    )
    vp_ht = vp_ht.select(
        gene_idx=vp_ht._gene_idx,
        v1_idx=hl.if_else(swap, vp_ht._v2_idx, vp_ht._v1_idx),
        v2_idx=hl.if_else(swap, vp_ht._v1_idx, vp_ht._v2_idx),
        locus1=hl.if_else(swap, vp_ht.locus2, vp_ht.locus1),
        alleles1=hl.if_else(swap, vp_ht.alleles2, vp_ht.alleles1),
        locus2=hl.if_else(swap, vp_ht.locus1, vp_ht.locus2),
        alleles2=hl.if_else(swap, vp_ht.alleles1, vp_ht.alleles2),
    )

    # --- Split into light and heavy ---
    heavy_variants = encoded_gt_ht.filter(
        encoded_gt_ht.all_samples.length() >= all_samples_len_threshold
    ).select().key_by("v_idx")

    vp_ht = vp_ht.annotate(
        _is_heavy=(
            hl.is_defined(heavy_variants[vp_ht.v1_idx])
            | hl.is_defined(heavy_variants[vp_ht.v2_idx])
        )
    )

    vp_light = vp_ht.filter(~vp_ht._is_heavy).drop("_is_heavy")
    vp_heavy = vp_ht.filter(vp_ht._is_heavy).drop("_is_heavy")

    vp_light.write(f"{output_dir}/vp_light.ht", overwrite=True)
    vp_heavy.write(f"{output_dir}/vp_heavy.ht", overwrite=True)

    light_count = hl.read_table(f"{output_dir}/vp_light.ht").count()
    heavy_count = hl.read_table(f"{output_dir}/vp_heavy.ht").count()
    logger.info(
        "Split pairs: %d light, %d heavy (threshold: all_samples.length() >= %d)",
        light_count,
        heavy_count,
        all_samples_len_threshold,
    )


def compute_counts_for_split(
    output_dir: str,
    split: str,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Compute genotype counts for a light or heavy split.

    Reads intermediate tables written by :func:`prepare_genotype_count_inputs`
    and runs the gene_idx co-partitioned computation (steps 4-6).

    :param output_dir: GCS directory containing intermediate outputs.
    :param split: Either ``"light"`` or ``"heavy"``.
    :param max_join_partitions: Upper bound on partition count.
    :return: Counts Table with gt_counts_raw and gt_counts_adj.
    """
    encoded_gt_ht = hl.read_table(f"{output_dir}/encoded_gt_sets_by_var_idx.ht")
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())
    vp_subset = hl.read_table(f"{output_dir}/vp_{split}.ht")

    logger.info("Computing counts for %s split (%d pairs)...", split, vp_subset.count())
    result = _compute_counts_for_subset(
        vp_subset, encoded_gt_ht, n_samples, split, max_join_partitions
    )
    return result.key_by("locus1", "alleles1", "locus2", "alleles2")


###############################################################################
# Clean pipeline (split strategy)
###############################################################################

DEFAULT_SPLIT_THRESHOLD = 37000
"""Default n_with_data threshold for splitting light/heavy variants (~P90)."""


###############################################################################
# Per-sample approach (use_new_shuffle to bypass local disk limits)
###############################################################################


def _create_gt_info_ht(mt: hl.MatrixTable) -> hl.Table:
    """
    Localize the dense MT into per-variant gt_info arrays.

    For each variant, ``gt_info`` is an array of ``(sample_idx, raw_gt,
    adj_gt)`` tuples for every sample with non-hom-ref data (same encoding
    as :func:`_encode_genotype_sets_by_var_idx`).

    :param mt: Dense filtered MatrixTable.
    :return: Table keyed by ``(locus, alleles)`` with ``gt_info`` field.
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )
    adj_gt_count_expr = hl.if_else(
        get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
        gt_count_expr, 0, missing_false=True,
    )
    mt = mt.select_entries(raw_gt=gt_count_expr, adj_gt=adj_gt_count_expr)
    ht = mt.localize_entries("_entries", "samples")

    gt = hl.enumerate(ht._entries)
    with_data = gt.filter(
        lambda x: hl.is_defined(x[1].raw_gt) | hl.is_defined(x[1].adj_gt)
    )
    return ht.select(
        gt_info=with_data.map(
            lambda x: (hl.int32(x[0]), x[1].raw_gt, x[1].adj_gt)
        ),
    )


def _build_variant_pair_map(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    output_dir: str,
) -> hl.Table:
    """
    Build a per-variant pair map: for each variant, its pairs and positions.

    Returns a Table keyed by ``var_idx`` with a ``vps`` field:
    ``array<(other_var_idx, vp_ht_idx, position)>`` where position 1 means
    this variant is v1 in the pair and 2 means v2.

    :param vp_ht: Variant pair list Table (keyed by vp_ht_idx).
    :param var_idx_ht: Var-idx Table keyed by ``(locus, alleles)``.
    :param output_dir: Directory for the intermediate checkpoint.
    :return: Per-variant pair map Table keyed by ``var_idx``.
    """
    vp_ht = vp_ht.annotate(
        var1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        var2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    vp_agg = vp_ht.explode("gene_id")
    vp_agg = vp_agg.group_by("gene_id").aggregate(
        v_by_vp=hl.agg.explode(
            lambda x: hl.agg.collect(x),
            [
                (vp_agg.var1_idx, vp_agg.var2_idx, vp_agg.vp_ht_idx, 1),
                (vp_agg.var2_idx, vp_agg.var1_idx, vp_agg.vp_ht_idx, 2),
            ],
        ).group_by(lambda x: x[0]),
    )
    agg_path = f"{output_dir}/vp_agg_by_gene.ht"
    vp_agg = vp_agg.checkpoint(agg_path, overwrite=True)

    # Flatten: one row per variant with its pair list.
    vp_map = vp_agg.select(v_by_vp=hl.array(vp_agg.v_by_vp)).explode("v_by_vp")
    vp_map = vp_map.select(
        var_idx=vp_map.v_by_vp[0],
        vps=vp_map.v_by_vp[1].map(lambda x: (x[1], x[2], x[3])),
    ).key_by("var_idx")
    return vp_map


def compute_genotype_counts_per_sample(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    variant_filter_ht: hl.Table,
    variant_filter_path: str,
    output_dir: str,
) -> hl.Table:
    """
    Compute genotype counts via per-sample-per-gene grouping.

    Instead of per-pair set intersections (which shuffle the large encoded
    GT table), this approach transposes the data:

    1. Localize the dense MT into per-variant ``gt_info`` arrays.
    2. Build a per-variant pair map (which pairs each variant is in).
    3. Annotate each variant with gene IDs, var_idx, and its pair map.
    4. Explode by gene, then by sample.
    5. Group by ``(sample_idx, gene_id)`` using ``use_new_shuffle``
       (shuffles to GCS, bypassing local disk limits).
    6. For each (sample, gene) group, iterate over each variant's pair
       list, look up the other variant's genotype, and emit per-pair
       contributions.
    7. Aggregate contributions by pair to get the 9-element count arrays.

    :param mt: Dense filtered MatrixTable.
    :param vp_ht: Variant pair list Table.
    :param variant_filter_ht: Variant filter Table with ``gene_id`` field.
    :param output_dir: Directory for intermediate checkpoints.
    :return: Counts Table with ``gt_counts_raw`` and ``gt_counts_adj``.
    """
    def _read_or_compute(path, compute_fn, step_name):
        """Read existing table or compute and write it."""
        try:
            ht = hl.read_table(path)
            logger.info("%s: reusing existing %s", step_name, path)
            return ht
        except Exception:
            logger.info("%s: computing...", step_name)
            ht = compute_fn()
            ht.write(path, overwrite=True)
            return hl.read_table(path)

    # --- Step 1: gt_info encoding ---
    gt_info_path = f"{output_dir}/gt_info.ht"
    gt_info_ht = _read_or_compute(
        gt_info_path, lambda: _create_gt_info_ht(mt), "Step 1 (gt_info)",
    )
    n_samples = hl.eval(hl.len(gt_info_ht.index_globals().samples))

    # --- Step 2: var_idx ---
    var_idx_path = f"{output_dir}/var_idx.ht"
    var_idx_ht = _read_or_compute(
        var_idx_path, lambda: _create_var_idx_ht(mt), "Step 2 (var_idx)",
    )

    # --- Step 3: variant pair map ---
    flat_path = f"{output_dir}/vp_map_flat.ht"
    vp_map = _read_or_compute(
        flat_path,
        lambda: _build_variant_pair_map(vp_ht, var_idx_ht, output_dir),
        "Step 3 (pair map)",
    )

    # --- Step 3b: gene_idx (min var_idx per gene) ---
    gene_idx_path = f"{output_dir}/gene_idx.ht"
    def _compute_gene_idx():
        vp_with_idx = vp_ht.annotate(
            var1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
            var2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
            _gene_id=hl.array(vp_ht.gene_id)[0],
        )
        t = vp_with_idx.key_by().select("_gene_id", "var1_idx", "var2_idx")
        t = t.annotate(_min_v=hl.min(t.var1_idx, t.var2_idx))
        return t.group_by("_gene_id").aggregate(
            gene_idx=hl.agg.min(t._min_v),
        ).key_by("_gene_id")
    gene_idx_ht = _read_or_compute(
        gene_idx_path, _compute_gene_idx, "Step 3b (gene_idx)",
    )

    # ===================================================================
    # From here on, gene_idx is the primary partition key. All tables
    # are keyed (gene_idx, ...) so that per-gene aggregations are
    # prefix-key scans with NO shuffle.
    # ===================================================================

    # Materialize gene_idx as a literal dict for fast lookups (no join).
    gene_idx_map = hl.dict(
        gene_idx_ht.aggregate(
            hl.agg.group_by(
                gene_idx_ht._gene_id,
                hl.agg.take(gene_idx_ht.gene_idx, 1).first(),
            )
        )
    )

    # Co-partition var_idx and variant_filter with gt_info for zip-joins.
    gt_info_intervals = gt_info_ht._calculate_new_partitions(
        gt_info_ht.n_partitions()
    )
    var_idx_ht = hl.read_table(var_idx_path, _intervals=gt_info_intervals)
    variant_filter_ht = hl.read_table(
        variant_filter_path, _intervals=gt_info_intervals,
    )

    # Build gene_ht: per-variant gene_idx array, var_idx, and vps.
    # All joins are zip-joins (co-partitioned) or broadcasts (vp_map).
    gene_ht = variant_filter_ht.select(
        gene_id=hl.set(variant_filter_ht.gene_id)
    )
    var_idx_expr = var_idx_ht[gene_ht.locus, gene_ht.alleles].var_idx
    gene_ht = gene_ht.transmute(
        gene_idx=gene_ht.gene_id.map(lambda x: gene_idx_map.get(x)),
        var_idx=var_idx_expr,
        vps=vp_map[var_idx_expr].vps,
    )
    gene_ht = gene_ht.filter(hl.is_defined(gene_ht.var_idx)).cache()

    # --- Step 4: Annotate gt_info + explode by gene_idx ---
    # Zip-join gene_ht onto gt_info_ht (co-partitioned), then explode
    # gene_idx and key by (gene_idx, var_idx).
    by_gene_path = f"{output_dir}/v5_gt_info_by_gene.ht"
    def _annotate_and_explode():
        return (
            gt_info_ht
            .annotate(**gene_ht[gt_info_ht.key])
            .explode("gene_idx")
            .key_by("gene_idx", "var_idx")
        )
    by_gene_ht = _read_or_compute(
        by_gene_path, _annotate_and_explode, "Step 4 (annotate + explode by gene)",
    )

    # --- Step 5: Repartition by size on (gene_idx, var_idx) key ---
    repart_path = f"{output_dir}/v5_gt_info_by_gene_repart.ht"
    def _repartition():
        n_parts = max(400, min(by_gene_ht.n_partitions() * 4, 10000))
        intervals = calculate_partitions_by_size(
            by_gene_ht, n_parts, size_field="gt_info", weight_field="vps",
        )
        return hl.read_table(by_gene_path, _intervals=intervals)
    by_gene_ht_repart = _read_or_compute(
        repart_path, _repartition, "Step 5 (repartition by size)",
    )

    # --- Step 6: Explode gt_info, group by (gene_idx, sample_idx) ---
    # The group_by key starts with gene_idx, which is already the first
    # partition key → data stays gene-local.
    # use_new_shuffle must be set BEFORE _read_or_compute because Hail
    # is lazy — the group_by executes during the write, not during
    # expression construction.
    grouped_path = f"{output_dir}/v5_per_sample_by_gene.ht"
    def _explode_and_group():
        ht = by_gene_ht_repart.explode("gt_info")
        ht = ht.transmute(
            sample_idx=ht.gt_info[0],
            raw_gt=ht.gt_info[1],
            adj_gt=ht.gt_info[2],
        )
        return ht.group_by(ht.gene_idx, ht.sample_idx).aggregate(
            gt_info=hl.agg.collect_as_set(
                hl.struct(
                    var_idx=ht.var_idx,
                    vps=ht.vps,
                    raw_gt=ht.raw_gt,
                    adj_gt=ht.adj_gt,
                )
            )
        )
    hl._set_flags(use_new_shuffle="1")
    per_sample = _read_or_compute(
        grouped_path, _explode_and_group, "Step 6 (explode + group by gene, sample)",
    )
    hl._set_flags(use_new_shuffle=None)

    # --- Step 6b: Repartition per-sample by gt_info size ---
    # use_new_shuffle produces many partitions. Repartition to a
    # reasonable count balanced by gt_info size.
    repart_grouped_path = f"{output_dir}/v5_per_sample_by_gene_repart.ht"
    def _repartition_per_sample():
        n_parts = max(200, min(1000, per_sample.n_partitions() // 20))
        logger.info("Step 6b: repartitioning to %d partitions", n_parts)
        intervals = calculate_partitions_by_size(
            per_sample, n_parts, size_field="gt_info",
        )
        return hl.read_table(grouped_path, _intervals=intervals)
    per_sample = _read_or_compute(
        repart_grouped_path, _repartition_per_sample,
        "Step 6b (repartition per-sample)",
    )

    # --- Step 7: Compute contribs per sample, key by (gene_idx, sample) ---
    # Flatten vps into (pair_id, position, gt) entries, group_by pair_id,
    # then extract v1/v2 genotypes. No dict lookup or double-counting filter.
    contribs_path = f"{output_dir}/v5_per_sample_contribs.ht"
    def _compute_contribs():
        # Flatten: for each variant, emit one entry per pair it's in.
        # has_data = True when either raw_gt or adj_gt is defined (matches
        # the shared all_samples definition in the set-based approach).
        pair_entries = hl.array(per_sample.gt_info).flatmap(lambda v:
            v.vps.map(lambda p: hl.struct(
                pair_id=p[1], position=p[2],
                raw_gt=v.raw_gt, adj_gt=v.adj_gt,
                has_data=hl.is_defined(v.raw_gt) | hl.is_defined(v.adj_gt),
            ))
        )
        # Group by pair_id → 1 or 2 entries per pair.
        by_pair = pair_entries.group_by(lambda x: x.pair_id)

        def _bin(v1_gt, v2_gt, v1_has_data, v2_has_data):
            """Compute bin index matching _count_from_sets behavior.

            A sample is "in data" when has_data is True (either raw or adj
            is defined). Samples "in data" but with gt==0 or gt==None are
            uncounted (bin=-1). Samples NOT "in data" at a variant are
            hom_ref (category 0) for that variant.
            """
            # If either side has data but its gt is 0 (missing) → uncounted.
            skip = (
                (v1_has_data & ((v1_gt == 0) | hl.is_missing(v1_gt)))
                | (v2_has_data & ((v2_gt == 0) | hl.is_missing(v2_gt)))
            )
            # Category: 0=hom_ref (not in data), 1=het, 2=hom_var.
            v1_cat = hl.if_else(v1_has_data, hl.or_else(v1_gt, 0), 0)
            v2_cat = hl.if_else(v2_has_data, hl.or_else(v2_gt, 0), 0)
            return hl.if_else(skip, -1, v1_cat * 3 + v2_cat)

        # For each pair, find v1 (position==1) and v2 (position==2).
        contribs = hl.array(by_pair).map(lambda kv: hl.bind(
            lambda v1, v2: hl.struct(
                pair_id=kv[0],
                has_overlap=hl.is_defined(v1) & hl.is_defined(v2),
                raw_bin=_bin(
                    hl.or_else(v1.raw_gt, -1),
                    hl.or_else(v2.raw_gt, -1),
                    hl.or_else(v1.has_data, False),
                    hl.or_else(v2.has_data, False),
                ),
                adj_bin=_bin(
                    hl.or_else(v1.adj_gt, -1),
                    hl.or_else(v2.adj_gt, -1),
                    hl.or_else(v1.has_data, False),
                    hl.or_else(v2.has_data, False),
                ),
            ),
            kv[1].find(lambda e: e.position == 1),
            kv[1].find(lambda e: e.position == 2),
        ))

        return per_sample.select(contribs=contribs)
    per_sample_contribs = _read_or_compute(
        contribs_path, _compute_contribs, "Step 7 (compute contribs)",
    )

    # --- Step 8: Explode contribs, group by pair_id ---
    # Group by pair_id only (not gene_idx) to deduplicate pairs that
    # appear in multiple genes. Produces one row per unique pair.
    logger.info("Step 8: Aggregating by pair_id...")
    exploded = per_sample_contribs.explode("contribs")
    c = exploded.contribs

    pair_counts = exploded.group_by(
        pair_id=c.pair_id,
    ).aggregate(
        counts=hl.agg.counter(
            hl.tuple([c.has_overlap, c.raw_bin, c.adj_bin])
        ),
    )

    # --- Step 9: Unpack counter, join with pair table, compute AABB ---
    logger.info("Step 9: Assembling final output...")

    # Helper: sum counts from the counter dict where a key field matches.
    # counts is dict<tuple(bool, int32, int32), int64>.
    cts = pair_counts.counts
    overlap = hl.sum(
        hl.array(cts).filter(lambda kv: kv[0][0]).map(lambda kv: kv[1])
    )

    def _bin_count(cts_expr, bin_idx, pos):
        """Sum counts where tuple position *pos* (1=raw, 2=adj) == bin_idx."""
        return hl.sum(
            hl.array(cts_expr)
            .filter(lambda kv: kv[0][pos] == bin_idx)
            .map(lambda kv: kv[1])
        )

    pair_counts = pair_counts.transmute(
        overlap=overlap,
        **{f"raw_{i}": _bin_count(cts, i, 1) for i in range(1, 9)},
        **{f"adj_{i}": _bin_count(cts, i, 2) for i in range(1, 9)},
    )

    # Join with pair table for locus/alleles and encoded GT for n_with_data.
    encoded_gt_ht = hl.read_table(f"{output_dir}/encoded_gt_sets_by_var_idx.ht")
    vp_keyed = vp_ht.annotate(
        v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    ).key_by("vp_ht_idx")

    pair_counts = pair_counts.annotate(vp=vp_keyed[pair_counts.pair_id])
    pair_counts = pair_counts.annotate(
        v1_n=encoded_gt_ht[pair_counts.vp.v1_idx].n_with_data,
        v2_n=encoded_gt_ht[pair_counts.vp.v2_idx].n_with_data,
    )
    AABB = (
        hl.int64(n_samples) - pair_counts.v1_n - pair_counts.v2_n
        + pair_counts.overlap
    )

    result = pair_counts.select(
        locus1=pair_counts.vp.locus1,
        alleles1=pair_counts.vp.alleles1,
        locus2=pair_counts.vp.locus2,
        alleles2=pair_counts.vp.alleles2,
        gt_counts_raw=hl.array(
            [AABB] + [pair_counts[f"raw_{i}"] for i in range(1, 9)]
        ),
        gt_counts_adj=hl.array(
            [AABB] + [pair_counts[f"adj_{i}"] for i in range(1, 9)]
        ),
    )
    return result.key_by("locus1", "alleles1", "locus2", "alleles2")


def encode_genotypes(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    output_dir: str,
) -> None:
    """
    Step A: Encode genotypes and write the encoded GT table.

    Creates ``var_idx.ht`` and ``encoded_gt_sets_by_var_idx.ht`` in *output_dir*.
    Run once; reuse for light and heavy compute steps with any threshold.

    :param mt: Dense filtered MatrixTable.
    :param vp_ht: Variant pair list Table (unused directly, but mt must
        contain exactly the variants in vp_ht).
    :param output_dir: GCS directory for outputs.
    """
    var_idx_path = f"{output_dir}/var_idx.ht"
    logger.info("Writing var_idx to %s", var_idx_path)
    _create_var_idx_ht(mt).write(var_idx_path, overwrite=True)
    var_idx_ht = hl.read_table(var_idx_path)

    encoded_path = f"{output_dir}/encoded_gt_sets_by_var_idx.ht"
    logger.info("Encoding genotypes to %s", encoded_path)
    _encode_genotype_sets_by_var_idx(mt, var_idx_ht).write(
        encoded_path, overwrite=True
    )
    encoded_gt_ht = hl.read_table(encoded_path)
    logger.info(
        "Encoded %d variants. n_with_data distribution: min=%d, median~=%d, max=%d",
        encoded_gt_ht.count(),
        encoded_gt_ht.aggregate(hl.agg.min(encoded_gt_ht.n_with_data)),
        encoded_gt_ht.aggregate(
            hl.agg.approx_quantiles(encoded_gt_ht.n_with_data, [0.5])[0]
        ),
        encoded_gt_ht.aggregate(hl.agg.max(encoded_gt_ht.n_with_data)),
    )


def compute_counts_light(
    output_dir: str,
    vp_ht: hl.Table,
    split_threshold: int = DEFAULT_SPLIT_THRESHOLD,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Step B: Compute genotype counts for the light split.

    Reads the encoded GT table from *output_dir*, classifies variants by
    ``n_with_data < split_threshold``, filters pairs to those where BOTH
    variants are light, builds a small GT table containing only those
    variants, and computes counts via a simple v1 zip-join + v2 shuffle.

    The v2 shuffle is safe on standard 40 GB-disk clusters because the
    light GT table contains only small-set variants (~1-2 GB total).

    :param output_dir: GCS directory written by :func:`encode_genotypes`.
    :param vp_ht: Variant pair list Table.
    :param split_threshold: ``n_with_data`` cutoff for heavy variants.
    :param max_join_partitions: Upper bound on partition count.
    :return: Counts Table with gt_counts_raw and gt_counts_adj.
    """
    var_idx_ht = hl.read_table(f"{output_dir}/var_idx.ht")
    encoded_gt_ht = hl.read_table(f"{output_dir}/encoded_gt_sets_by_var_idx.ht")
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())

    # Classify variants.
    heavy_variants = encoded_gt_ht.filter(
        encoded_gt_ht.n_with_data >= split_threshold
    ).select().key_by("v_idx")

    # Add v_idx to pair table and filter to light pairs.
    vp_ht = vp_ht.annotate(
        v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    vp_light = vp_ht.filter(
        ~hl.is_defined(heavy_variants[vp_ht.v1_idx])
        & ~hl.is_defined(heavy_variants[vp_ht.v2_idx])
    )
    vp_light = vp_light.key_by("v1_idx", "v2_idx").select(
        "locus1", "alleles1", "locus2", "alleles2"
    )
    logger.info("Light pairs: %d", vp_light.count())

    # Build small GT table for light variants only.
    light_v1 = vp_light.key_by(v_idx=vp_light.v1_idx).select().distinct()
    light_v2 = vp_light.key_by(v_idx=vp_light.v2_idx).select().distinct()
    light_variants = light_v1.union(light_v2).distinct()
    gt_light = encoded_gt_ht.semi_join(light_variants)

    gt_light_path = hl.utils.new_temp_file("gt_light", "ht")
    gt_light.write(gt_light_path, overwrite=True)
    gt_light = hl.read_table(gt_light_path)
    logger.info("Light GT table: %d variants", gt_light.count())

    # Write vp_light so we can re-read with partition intervals.
    vp_light_path = hl.utils.new_temp_file("vp_light", "ht")
    vp_light.write(vp_light_path, overwrite=True)

    # Co-partition on v1_idx for the v1 zip-join.
    n_parts = min(gt_light.n_partitions() * 3, max_join_partitions)
    partition_intervals = calculate_partitions_by_size(
        gt_light, n_parts, size_field="n_with_data"
    )
    gt_light = hl.read_table(gt_light_path, _intervals=partition_intervals)
    vp_light = hl.read_table(vp_light_path, _intervals=partition_intervals)

    # v1 zip-join (free), v2 shuffle of small gt_light (~1-2 GB).
    v1 = gt_light[vp_light.v1_idx]
    v2 = gt_light[vp_light.v2_idx]
    result = vp_light.select(
        "locus1", "alleles1", "locus2", "alleles2",
        gt_counts_raw=_count_from_sets(
            v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
        gt_counts_adj=_count_from_sets(
            v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
            v1.all_samples_is_complement,
            v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
            v2.all_samples_is_complement,
            n_samples,
        ),
    )
    return result.key_by("locus1", "alleles1", "locus2", "alleles2")


def compute_counts_heavy(
    output_dir: str,
    vp_ht: hl.Table,
    split_threshold: int = DEFAULT_SPLIT_THRESHOLD,
    max_join_partitions: int = 10000,
) -> hl.Table:
    """
    Step C: Compute genotype counts for the heavy split.

    Reads the encoded GT table from *output_dir*, classifies variants,
    filters pairs to those where at least one variant is heavy, and
    computes counts using gene_idx co-partitioning (both v1 and v2 lookups
    are zip-joins). Requires a cluster with SSDs or larger disks for the
    rekey shuffle in the v2 step.

    :param output_dir: GCS directory written by :func:`encode_genotypes`.
    :param vp_ht: Variant pair list Table (with gene_id field).
    :param split_threshold: ``n_with_data`` cutoff for heavy variants.
    :param max_join_partitions: Upper bound on partition count.
    :return: Counts Table with gt_counts_raw and gt_counts_adj.
    """
    var_idx_ht = hl.read_table(f"{output_dir}/var_idx.ht")
    encoded_gt_ht = hl.read_table(f"{output_dir}/encoded_gt_sets_by_var_idx.ht")
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())

    # Classify variants.
    heavy_variants = encoded_gt_ht.filter(
        encoded_gt_ht.n_with_data >= split_threshold
    ).select().key_by("v_idx")

    # Add v_idx to pair table and filter to heavy pairs.
    vp_ht = vp_ht.annotate(
        v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].var_idx,
        v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].var_idx,
    )
    vp_heavy = vp_ht.filter(
        hl.is_defined(heavy_variants[vp_ht.v1_idx])
        | hl.is_defined(heavy_variants[vp_ht.v2_idx])
    )
    logger.info("Heavy pairs: %d", vp_heavy.count())

    # --- V2 role assignment ---
    logger.info("Collecting variants and edges for V2 assignment...")
    variant_data = encoded_gt_ht.select("n_with_data").collect()
    variants = [(row.v_idx, row.n_with_data) for row in variant_data]

    edge_data = vp_heavy.key_by().select("v1_idx", "v2_idx").collect()
    edges = [(row.v1_idx, row.v2_idx) for row in edge_data]

    v2_set = compute_v2_independent_set(variants, edges)
    logger.info("V2 set: %d variants", len(v2_set))

    v2_broadcast = hl.literal(v2_set, dtype=hl.tset(hl.tint64))
    n_with_data_map = hl.literal(
        {v: n for v, n in variants}, dtype=hl.tdict(hl.tint64, hl.tint32)
    )

    # --- gene_idx ---
    vp_heavy = vp_heavy.annotate(_gene_id=hl.array(vp_heavy.gene_id)[0])
    gene_idx_ht = vp_heavy.key_by().select("_gene_id", "v1_idx", "v2_idx")
    gene_idx_ht = gene_idx_ht.annotate(
        _min_v=hl.min(gene_idx_ht.v1_idx, gene_idx_ht.v2_idx)
    )
    gene_idx_ht = gene_idx_ht.group_by("_gene_id").aggregate(
        gene_idx=hl.agg.min(gene_idx_ht._min_v)
    )
    gene_idx_ht = gene_idx_ht.key_by("_gene_id")

    # --- Restructure: apply V2 swap + gene_idx ---
    vp_heavy = vp_heavy.annotate(
        _v1_is_v2=v2_broadcast.contains(vp_heavy.v1_idx),
        _v2_is_v2=v2_broadcast.contains(vp_heavy.v2_idx),
        gene_idx=gene_idx_ht[vp_heavy._gene_id].gene_idx,
    )
    swap = (
        (vp_heavy._v1_is_v2 & ~vp_heavy._v2_is_v2)
        | (
            ~vp_heavy._v1_is_v2
            & ~vp_heavy._v2_is_v2
            & (
                n_with_data_map.get(vp_heavy.v1_idx)
                < n_with_data_map.get(vp_heavy.v2_idx)
            )
        )
    )
    vp_heavy = vp_heavy.select(
        gene_idx=vp_heavy.gene_idx,
        v1_idx=hl.if_else(swap, vp_heavy.v2_idx, vp_heavy.v1_idx),
        v2_idx=hl.if_else(swap, vp_heavy.v1_idx, vp_heavy.v2_idx),
        locus1=hl.if_else(swap, vp_heavy.locus2, vp_heavy.locus1),
        alleles1=hl.if_else(swap, vp_heavy.alleles2, vp_heavy.alleles1),
        locus2=hl.if_else(swap, vp_heavy.locus1, vp_heavy.locus2),
        alleles2=hl.if_else(swap, vp_heavy.alleles1, vp_heavy.alleles2),
    )

    # Use _compute_counts_for_subset for the gene_idx co-partitioned path.
    result = _compute_counts_for_subset(
        vp_heavy, encoded_gt_ht, n_samples, "heavy", max_join_partitions
    )
    return result.key_by("locus1", "alleles1", "locus2", "alleles2")


def main(args):
    """Create variant pair matrix from gnomAD v4 VDS."""
    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix or ""
    data_type = args.data_type
    least_consequence = args.least_consequence
    max_freq = args.max_freq
    test = args.test or bool(args.gene)
    test_intervals = (
        {args.gene: TEST_INTERVALS[args.gene]} if args.gene else TEST_INTERVALS
    )
    
    # Get current Hail version
    hail_version = hl.version().split('-')[0]  # Remove git hash suffix
    current_version = version.parse(hail_version)
    threshold_version = version.parse("0.2.120")
    
    if current_version > threshold_version:
        logger.warning(
            f"WARNING: Using Hail version {hl.__version__} which is greater than 0.2.120. "
            f"This will cause issues in create_variant_pair_ht, please use Hail version 0.2.120 or lower."
        )    

    hl.init(
        log=os.path.join(tempfile.gettempdir(), "create_vp_matrix.log"),
        tmp_dir=tmp_dir,
        backend=args.backend,
    )

    logger.info(
        f"""
        Running script with the following parameters:

            Data type: {data_type}
            Backend: {args.backend}
            Test: {test}
            Gene: {args.gene or 'all test intervals'}
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
    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
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
            filter_ht = filter_for_testing(filter_ht, test_intervals)
            freq_ht = filter_for_testing(freq_ht, test_intervals)
            vep_ht = filter_for_testing(vep_ht, test_intervals)

        ht = create_variant_filter_ht(
            filter_ht,
            freq_ht,
            vep_ht,
            least_consequence=least_consequence,
            max_freq=max_freq,
        ).checkpoint(res.variant_filter_ht.path, overwrite=overwrite)
        logger.info("Number of variants in the variant filter Table: %d", ht.count())

    if args.filter_vmt:
        logger.info(f"Filtering gnomAD v4 {data_type} variant data MatrixTable...")
        #if current_version > threshold_version:
            #raise ValueError(
            #    "Hail version 0.2.120 or lower is required handle the vds filtering "
            #    "correctly."
            #)
        res = resources.filter_vmt
        res.check_resource_existence()

        vds = get_vds_func(
            release_only=True,
            split=True,
            filter_intervals=None if not test else list(test_intervals.values()),
            filter_variant_ht=res.variant_filter_ht.ht(),
            entries_to_keep=["GT"],
            split_reference_blocks=False,
        )
        vds.variant_data.write(res.filtered_vmt.path, overwrite=overwrite)
        logger.info("The filtered VDS has been written...")

    if args.create_variant_pair_list_ht:
        logger.info("Creating variant pair list Table...")
        res = resources.create_variant_pair_list_ht
        res.check_resource_existence()

        ht = create_variant_pair_ht(res.filtered_vmt.mt(), res.variant_filter_ht.ht())
        ht = ht.checkpoint(res.vp_list_ht.path, overwrite=overwrite)
        logger.info(
            "The variant pair list Table has been written...\n"
            f"The number of unique variant pairs is {ht.count()}"
        )

    if args.create_dense_filtered_mt:
        logger.info("Creating dense filtered MatrixTable...")
        #if current_version > threshold_version:
        #    raise ValueError(
        #        "Hail version 0.2.120 or lower is required handle the vds filtering "
        #        "correctly."
        #    )
        res = resources.create_dense_filtered_mt
        res.check_resource_existence()

        ht = create_variant_pair_filter_ht(res.vp_list_ht.ht())
        vds = get_vds_func(
            release_only=True,
            split=True,
            filter_intervals=None if not test else list(test_intervals.values()),
            filter_variant_ht=ht,
            entries_to_keep=["GT", "GQ", "DP", "AD"],
            split_reference_blocks=False,
        )
        mt = hl.vds.to_dense_mt(vds)
        mt = mt.checkpoint(res.dense_filtered_mt.path, overwrite=overwrite)
        logger.info(
            "The dense filtered MatrixTable has been written...\n"
            f"The number of rows in the dense filtered MatrixTable is {mt.count_rows()}"
        )

    # --- Genotype count steps (4 phases, can run on different clusters) ---
    count_output_dir = f"{tmp_dir}/genotype_count_intermediates{_get_output_postfix(output_postfix, test)}"
    split_threshold = args.split_threshold

    if args.encode_genotypes:
        logger.info("Encoding genotypes...")
        res = resources.create_variant_pair_genotype_counts_ht
        res.check_resource_existence()

        encode_genotypes(
            res.dense_filtered_mt.mt(),
            res.vp_list_ht.ht(),
            output_dir=count_output_dir,
        )
        logger.info("Encoded genotypes written to %s", count_output_dir)

    if args.compute_counts_light:
        logger.info("Computing counts for light split (threshold=%d)...", split_threshold)
        res = resources.create_variant_pair_genotype_counts_ht

        ht = compute_counts_light(
            count_output_dir, res.vp_list_ht.ht(), split_threshold
        )
        ht.write(f"{count_output_dir}/counts_light.ht", overwrite=overwrite)
        logger.info("Light counts written.")

    if args.compute_counts_heavy:
        logger.info("Computing counts for heavy split (threshold=%d)...", split_threshold)
        res = resources.create_variant_pair_genotype_counts_ht

        ht = compute_counts_heavy(
            count_output_dir, res.vp_list_ht.ht(), split_threshold
        )
        ht.write(f"{count_output_dir}/counts_heavy.ht", overwrite=overwrite)
        logger.info("Heavy counts written.")

    if args.compute_counts_per_sample:
        logger.info("Computing counts via per-sample grouping...")
        res = resources.create_variant_pair_genotype_counts_ht
        res.check_resource_existence()

        vf_resource = get_variant_filter_ht(
            data_type=data_type, test=test, tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )

        ht = compute_genotype_counts_per_sample(
            res.dense_filtered_mt.mt(),
            res.vp_list_ht.ht(),
            vf_resource.ht(),
            variant_filter_path=vf_resource.path,
            output_dir=count_output_dir,
        )
        ht.write(res.vp_gt_counts_ht.path, overwrite=overwrite)
        logger.info("Per-sample counts written to %s", res.vp_gt_counts_ht.path)

    if args.combine_counts:
        logger.info("Combining light + heavy counts...")
        res = resources.create_variant_pair_genotype_counts_ht

        light_path = f"{count_output_dir}/counts_light.ht"
        heavy_path = f"{count_output_dir}/counts_heavy.ht"
        tables = []
        for path in [light_path, heavy_path]:
            try:
                tables.append(hl.read_table(path))
            except Exception:
                logger.info("Skipping %s (not found).", path)
        if tables:
            ht = tables[0] if len(tables) == 1 else tables[0].union(tables[1])
            ht = ht.key_by("locus1", "alleles1", "locus2", "alleles2")
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
        "--backend",
        default="spark",
        choices=("spark", "batch"),
        help=(
            "Hail Query backend: 'spark' (default, uses local/Dataproc Spark) or "
            "'batch' (Hail Query-on-Batch). Use 'batch' to run on Hail Batch instead "
            "of locally or on Dataproc."
        ),
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to test intervals (all genes in TEST_INTERVALS) for testing.",
    )
    parser.add_argument(
        "--gene",
        choices=list(TEST_INTERVALS),
        help=(
            "Run on a single gene; uses that gene's interval from TEST_INTERVALS and "
            "implies --test."
        ),
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
        "--include-region-variants",
        action="store_true",
        help=(
            "Include variants in retained exon regions (Laura's noncanonical "
            "exons) in the variant filter Table, regardless of VEP consequence "
            "severity. Used with --create-variant-filter-ht."
        ),
    )
    parser.add_argument(
        "--exon-upstream-padding",
        type=int,
        default=DEFAULT_EXON_UPSTREAM_PADDING,
        help=(
            "Number of bp to pad before each exon start for "
            "--include-region-variants. Default is "
            f"{DEFAULT_EXON_UPSTREAM_PADDING}."
        ),
    )
    parser.add_argument(
        "--exon-downstream-padding",
        type=int,
        default=DEFAULT_EXON_DOWNSTREAM_PADDING,
        help=(
            "Number of bp to pad after each exon end for "
            "--include-region-variants. Default is "
            f"{DEFAULT_EXON_DOWNSTREAM_PADDING}."
        ),
    )
    parser.add_argument(
        "--include-noncoding-pathogenic",
        action="store_true",
        help=(
            "Include noncoding variants that are ClinVar P/LP, have spliceAI "
            "> --min-splice-ai, or pangolin > --min-pangolin in the variant "
            "filter Table. Used with --create-variant-filter-ht."
        ),
    )
    parser.add_argument(
        "--min-splice-ai",
        type=float,
        default=DEFAULT_MIN_SPLICE_AI,
        help=(
            f"Minimum spliceAI delta score for --include-noncoding-pathogenic. "
            f"Default is {DEFAULT_MIN_SPLICE_AI}."
        ),
    )
    parser.add_argument(
        "--min-pangolin",
        type=float,
        default=DEFAULT_MIN_PANGOLIN,
        help=(
            f"Minimum pangolin delta score for --include-noncoding-pathogenic. "
            f"Default is {DEFAULT_MIN_PANGOLIN}."
        ),
    )
    parser.add_argument(
        "--filter-vmt",
        action="store_true",
        help="Filter the MatrixTable for determining variant pairs.",
    )
    parser.add_argument(
        "--create-variant-pair-list-ht",
        action="store_true",
        help="first create just the list of possible variant pairs.",
    )
    parser.add_argument(
        "--create-dense-filtered-mt",
        action="store_true",
        help="Create the dense filtered MatrixTable.",
    )
    parser.add_argument(
        "--encode-genotypes",
        action="store_true",
        help=(
            "Step A: Encode genotypes from the dense MT into per-variant "
            "sample sets. Run once; reuse for light/heavy with any threshold."
        ),
    )
    parser.add_argument(
        "--compute-counts-light",
        action="store_true",
        help=(
            "Step B: Compute genotype counts for light pairs (both variants "
            "below --split-threshold). Safe on standard 40 GB-disk clusters."
        ),
    )
    parser.add_argument(
        "--compute-counts-heavy",
        action="store_true",
        help=(
            "Step C: Compute genotype counts for heavy pairs (at least one "
            "variant above --split-threshold). Requires SSD or large-disk cluster."
        ),
    )
    parser.add_argument(
        "--compute-counts-per-sample",
        action="store_true",
        help=(
            "Compute genotype counts via per-sample-per-gene grouping. "
            "Uses use_new_shuffle to bypass local disk limits. "
            "Exact counts, works on standard clusters."
        ),
    )
    parser.add_argument(
        "--combine-counts",
        action="store_true",
        help="Step D: Union light + heavy counts into the final output Table.",
    )
    parser.add_argument(
        "--split-threshold",
        type=int,
        default=DEFAULT_SPLIT_THRESHOLD,
        help=(
            f"n_with_data cutoff for classifying variants as heavy. "
            f"Default is {DEFAULT_SPLIT_THRESHOLD} (~P90 for test genes)."
        ),
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

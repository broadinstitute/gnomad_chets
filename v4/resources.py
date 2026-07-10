"""Script containing variant co-occurrence related resources."""

from typing import Optional

import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    TableResource,
    VariantDatasetResource,
)
from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad.resources.grch38.gnomad import GEN_ANC_GROUPS, all_sites_an
from gnomad.resources.grch38.reference_data import clinvar
from gnomad_qc.v4.resources.annotations import (
    get_freq,
    get_insilico_predictors,
    get_vep,
)
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.variant_qc import final_filter

########################################################################################
### Constants
########################################################################################

DEFAULT_TMP_DIR = "gs://gnomad-tmp-30day"
"""Default temporary directory for variant co-occurrence pipeline output files."""

VARIANT_COOCCURRENCE_ROOT = "gs://gnomad/v4.1/variant_cooccurrence"
"""Official output root directory for variant co-occurrence pipeline output files."""

TEST_INTERVALS = {
    "PCNT": "chr21:46324141-46445769",
    "COL18A1": "chr21:45405123-45513721",
    "AHNAK2": "chr14:104937244-104978374",
    "TTN": "chr2:178525989-178830802",
    "FLG": "chr1:152302165-152325239",
    "OBSCN": "chr1:228208044-228378876",
    # "HRNR": "chr1:152212076-152224193", # Flagged low MQ
    # "NBPF10": "chr1:146064711-146229000", # Flagged low coverage and low MQ
    "PLEC": "chr8:143915153-143976734",
    # "PDE4DIP": "chr1:148808181-149048286", # Flagged low MQ
    # "FCGBP": "chr19:39863323-39934626", # Flagged low coverage and low MQ
    "NEB": "chr2:151485336-151734487",
    "LAMA5": "chr20:62307955-62367312",
    "SYNE1": "chr6:152121687-152637801",
    "CAPN3": "chr15:42359498-42412949",
    "SGCA": "chr17:50164214-50175928",
    # "DMD": "chrX:31097677-33339609", # chrX — hemizygous males encode as hom-var, breaking HWE and the in-trans test; needs sex-aware filter
    "DYSF": "chr2:71453561-71686763",
    "ANO5": "chr11:21782659-22283567",
    "MUC16": "chr19:8848844-8981342",
    "CHD2": "chr15:92900189-93027996",
    "ABCA4": "chr1:93992834-94121148",
    "ACO2": "chr22:41447830-41529273",
    "ATM": "chr11:108223044-108369102",
    "PKD1": "chr16:2088708-2135898",
    "PRF1": "chr10:70597348-70602759",
    "HFE": "chr6:26087281-26098343",
    "SERPINA1": "chr14:94376747-94390693",
    "NPHS2": "chr1:179550539-179575952",
}
"""Test intervals for genes used in testing mode."""

DATA_TYPE_CHOICES = ["exomes", "genomes"]
"""Valid data type choices for variant co-occurrence pipeline."""

DEFAULT_DATA_TYPE = "exomes"
"""Default data type for variant co-occurrence pipeline."""

GLOBAL_POP = "all"
"""Label for the full-cohort (global) stratum in per-population outputs.

The per-pop genotype-count / EM outputs are keyed by genetic-ancestry group
plus this ``"all"`` entry, which is the whole cohort (equals the flat global
``gt_counts_*`` / ``em`` fields, and the sum over the specific groups when
every sample carries a group label)."""


def get_pops(data_type: str = DEFAULT_DATA_TYPE) -> list:
    """Genetic-ancestry groups for per-population stratification, ``all`` first.

    Wraps ``gnomad.resources.grch38.gnomad.GEN_ANC_GROUPS["v4"][data_type]``
    (afr, amr, asj, eas, fin, mid, nfe, remaining, sas for exomes) and
    prepends :data:`GLOBAL_POP`.
    """
    return [GLOBAL_POP] + list(GEN_ANC_GROUPS["v4"][data_type])


def get_sample_pop_ht(data_type: str = DEFAULT_DATA_TYPE) -> hl.Table:
    """``s`` → genetic-ancestry group (``pop``) from the v4 sample-QC meta HT.

    Centralizes the ``population_inference.pop`` field path so the count and
    trio scripts attach the same per-sample pop label the v4 freq pipeline uses.
    """
    ht = meta(data_type=data_type).ht()
    return ht.select(pop=ht.population_inference.pop)

DEFAULT_MAX_FREQ = 0.05
"""Default maximum global AF to keep (inclusive)."""

DEFAULT_LEAST_CONSEQUENCE = "3_prime_UTR_variant"
"""Default lowest-severity consequence to keep."""

DEFAULT_EXON_UPSTREAM_PADDING = 3
"""Default padding in bp before each exon start (acceptor side, -1 to -3)."""

DEFAULT_EXON_DOWNSTREAM_PADDING = 8
"""Default padding in bp after each exon end (donor side, +1 to +8)."""

DEFAULT_MIN_SPLICE_AI = 0.2
"""Default minimum spliceAI delta score for noncoding pathogenic variant inclusion."""

DEFAULT_MIN_PANGOLIN = 0.14
"""Default minimum pangolin delta score for noncoding pathogenic variant inclusion."""

########################################################################################
### In-trans-OE candidate / intronic-padding defaults
########################################################################################

DEFAULT_IN_TRANS_OE_MAX_AF = 0.2
"""Default upper AF bound (inclusive) for in-trans-OE candidates.

Deliberately above the standard pipeline's 5% cap so that genuinely
common-but-suspect variants enter the analysis. The implicit lower bound
is 0 (any positive AF qualifies)."""

DEFAULT_IN_TRANS_OE_ACCEPTOR_PADDING = 50
"""Default acceptor-side padding (bp) for the in-trans-OE intronic-padding
source. Covers the branch-point region (lariat adenosine, typically -18 to
-40 from the splice site)."""

DEFAULT_IN_TRANS_OE_DONOR_PADDING = 15
"""Default donor-side padding (bp) for the in-trans-OE intronic-padding
source. Covers cryptic 5' splice signals slightly beyond VEP's built-in +8
splice region."""

########################################################################################
### Source-tag constants for the orthogonal --include-* flags
########################################################################################
# A variant in the filter HT can carry multiple tags if it qualifies under
# more than one path; the union is built when sources are merged in
# create_variant_filter_ht.

# ClinVar release pinned for the sites HT (see preprocess_sites_ht). Stamped
# onto the sites HT globals as ``clinvar_version`` for provenance.
CLINVAR_VERSION = "20250504"

SOURCE_CLINVAR_PLP = "clinvar_plp"
SOURCE_CLINVAR_BLB = "clinvar_blb"
SOURCE_CLINVAR_VUS = "clinvar_vus"
SOURCE_HC_LOF = "hc_lof"
SOURCE_SPLICE_PATH = "splice_path"
SOURCE_IN_TRANS_OE_CANDIDATE = "in_trans_oe_candidate"
SOURCE_IN_TRANS_OE_INTRONIC_PADDING = "in_trans_oe_intronic_padding"

IN_TRANS_OE_CANDIDATE_SOURCES = frozenset(
    [SOURCE_IN_TRANS_OE_CANDIDATE, SOURCE_IN_TRANS_OE_INTRONIC_PADDING]
)
"""Variants tagged ONLY with these sources are considered "OE-candidate-only"
— pairs where both sides are OE-candidate-only get dropped by the pair
filter (avoids candidate × candidate explosion)."""

IN_TRANS_OE_PARTNER_SOURCES = frozenset(
    [
        SOURCE_CLINVAR_PLP,
        SOURCE_CLINVAR_BLB,
        SOURCE_CLINVAR_VUS,
        SOURCE_HC_LOF,
        SOURCE_SPLICE_PATH,
    ]
)
"""Sources that mark a variant as a "partner" for the in-trans-OE analysis.
Kept for reference / inspection scripts; the pair filter uses
IN_TRANS_OE_CANDIDATE_SOURCES for the OE-only test."""

########################################################################################
### Sites HT schema field names
########################################################################################
# Field names emitted by ``assemble_sites_ht``. The variant-filter builder
# consults the schema to gate the per-source ``include_*`` flags.

SITES_FIELD_SPLICEAI = "spliceai_ds_max"
SITES_FIELD_PANGOLIN = "pangolin_largest_ds"
SITES_FIELD_CLINVAR = "clinvar"

CLINVAR_CATEGORY_FIELD_FMT = "is_{category}"
"""Per-category membership-bool field name template inside the
``SITES_FIELD_CLINVAR`` struct on the sites HT (e.g. ``is_plp``,
``is_blb``, ``is_vus`` — one per :data:`CLINVAR_CATEGORIES`). Precomputed
by :func:`gnomad_chets.v4.create_vp_list.assemble_sites_ht` via
:func:`gnomad_chets.v4.utils.clinvar_category_match_expr`. The sites
HT's ``clinvar`` struct also carries ``GENEINFO`` for later
VEP-symbol cross-referencing in
:func:`gnomad_chets.v4.create_vp_list._get_clinvar_gene_id_expr`."""

########################################################################################
### ClinVar significance categories
########################################################################################

CLINVAR_CATEGORY_PLP = "plp"
CLINVAR_CATEGORY_BLB = "blb"
CLINVAR_CATEGORY_VUS = "vus"
CLINVAR_CATEGORIES = (
    CLINVAR_CATEGORY_PLP,
    CLINVAR_CATEGORY_BLB,
    CLINVAR_CATEGORY_VUS,
)
CLINVAR_CATEGORY_SOURCE_TAG = {
    CLINVAR_CATEGORY_PLP: SOURCE_CLINVAR_PLP,
    CLINVAR_CATEGORY_BLB: SOURCE_CLINVAR_BLB,
    CLINVAR_CATEGORY_VUS: SOURCE_CLINVAR_VUS,
}
CLINVAR_CATEGORY_PARTNER_SET = {
    CLINVAR_CATEGORY_PLP: "CLINVAR_PLP",
    CLINVAR_CATEGORY_BLB: "CLINVAR_BLB",
    CLINVAR_CATEGORY_VUS: "CLINVAR_VUS",
}
"""Maps each ClinVar significance category to its corresponding
``partner_set`` label (the PARTNER_SET_* constants in v4.in_trans_oe).
Hardcoded string values to avoid a back-import; kept in sync by convention."""


########################################################################################
### Create Variant Co-occurrence Matrix Resource Functions
########################################################################################
def _get_output_dir(test: bool, tmp_dir: Optional[str]) -> str:
    """
    Determine the output directory based on test and tmp_dir parameters.

    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files. If None and test is False,
        uses VARIANT_COOCCURRENCE_ROOT. If None and test is True, uses DEFAULT_TMP_DIR.
    :return: Output directory path.
    """
    if tmp_dir is not None:
        return tmp_dir
    elif test:
        return DEFAULT_TMP_DIR
    else:
        return VARIANT_COOCCURRENCE_ROOT


def _get_output_postfix(output_postfix: Optional[str], test: bool) -> str:
    """
    Determine the output postfix based on output_postfix and test parameters.

    :param output_postfix: Postfix to append to output file names. If None and test is
        True, uses ".pcnt_test". If None and test is False, uses empty string.
    :param test: Whether to use a test postfix.
    :return: Output postfix string (with leading dot if not empty).
    """
    if output_postfix is not None:
        return f".{output_postfix}"
    elif test:
        return ".pcnt_test"
    else:
        return ""


def _get_resource_path(
    data_type: str,
    resource_name: str,
    extension: str,
    test: bool,
    tmp_dir: Optional[str],
    output_postfix: Optional[str],
) -> str:
    """
    Generate the full path for a variant co-occurrence resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param resource_name: Name of resource (e.g., 'vp_list' or 'vp_full').
    :param extension: File extension (e.g., '.ht' or '.vds').
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Full path string for the resource.
    """
    output_dir = _get_output_dir(test, tmp_dir)
    postfix = _get_output_postfix(output_postfix, test)
    return f"{output_dir}/{data_type}.{resource_name}{postfix}{extension}"


def get_sites_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get assembled per-variant sites Table resource.

    The sites HT is the join of every per-variant annotation source that
    downstream pipeline steps might need: ``filters`` / ``af`` / ``an`` /
    ``an_pct`` / full VEP / SpliceAI + Pangolin scores / ClinVar info.
    Built by :func:`gnomad_chets.v4.create_vp_list.assemble_sites_ht`
    and consumed by :func:`create_variant_filter_ht`.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Sites Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="sites",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_variant_filter_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get variant filter Table resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Variant filter Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="variant_filter",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_filtered_vmt(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> MatrixTableResource:
    """
    Get filtered MatrixTable resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Filtered MatrixTable resource.
    """
    return MatrixTableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="filtered_vmt",
            extension=".mt",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_variant_pair_list_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get variant pair list Table resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Variant pair list Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="variant_pairs",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_variant_pair_genotype_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get full variant pair genotype Table resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Full variant pair genotype Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="variant_pairs.genotypes",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_variant_pair_genotype_counts_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get variant pair genotype counts Table resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Variant pair genotype counts Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="variant_pairs.genotype_counts",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_variant_size_info_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the variant size-info Table resource.

    One row per variant in the encoded GT table, with diagnostic columns
    used to decide heavy/light routing and visualize the distribution:

    * ``v_idx`` (key) — variant index from var_idx.ht
    * ``locus``, ``alleles`` — joined from var_idx.ht
    * ``gene_id`` — joined from the variant filter HT
    * ``_contribution`` — degree(v) × payload(v) bytes
    * ``_cum_before`` — cumulative contribution of variants ranked higher
      in the descending-by-contribution order (= sum strictly preceding
      this row)
    * ``_bytes`` — partner-load metric: max(v1_load, v2_load) computed by
      the cross-role augmentation pass
    * ``is_heavy`` — pulled into the heavy set by the greedy budget at
      build time (does not include downstream cutoff filtering)
    * ``split_count`` — ``ceil(max(contribution, partner_load) / TARGET)``

    Downstream consumers (``--compute-counts-light`` and
    ``--compute-counts-heavy``) read this HT and apply a runtime contribution
    cutoff to redo the heavy/light decision without rebuilding the HT.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="variant_size_info",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_excluded_genes_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the excluded-genes Table resource.

    One row per ``gene_id`` to be **fully excluded** from both the light
    and heavy counting steps — typically the chr19/exome heavy hitters
    (MUC16, RYR1, FBN3, FCGBP, COL5A3) that we defer to per-gene jobs.

    Schema: keyed by ``gene_id`` (str), with an optional ``reason`` field.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="excluded_genes",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_phase(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
   Get phased variant pair Table resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Full variant pair genotype Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="phased",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_annotated_phase(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get annotated phased variant pair Table resource.

    Output of :func:`gnomad_chets.v4.phase_gnomad.annotate_phased_ht_with_sites`
    — the phased HT with ``v1_ann`` and ``v2_ann`` per-variant structs
    (source tags, gene_id, an_pct, AC/AN/AF, most_severe_consequence,
    gene_symbols, SpliceAI, Pangolin, ClinVar) and a pair-level
    ``shared_gene_ids`` field.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="phased.annotated",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


# Per-variant "only filters" HT covering all v4 exome variants (not restricted
# to AC>0 in release) — the row-set anchor for the sites HT. AC0-in-release
# variants that are otherwise QC-PASS are retained so trio-side annotations
# can find them.
FINAL_FILTER_ONLY_FILTERS_PATH = (
    "gs://gnomad/v4.1/variant_qc/exomes/"
    "gnomad.exomes.v4.1.final_filter.all_variants.only_filters.ht"
)

# Canonical variant-QC final-filter HT. Its ``filters`` field is byte-identical
# to the v4.1.1 public release ``filters`` (the release is populated from it),
# and differs from FINAL_FILTER_ONLY_FILTERS_PATH *only* in the InbreedingCoeff
# token: the two were built from different freq HTs, so InbreedingCoeff was
# recomputed (a genuine value change, not float noise). We carry both into the
# sites HT so the InbreedingCoeff filter can be resolved downstream rather than
# baked in here. Covers fewer variants than only_filters (release-QC set only).
FINAL_FILTER_PATH = (
    "gs://gnomad/v4.1/variant_qc/exomes/"
    "gnomad.exomes.v4.1.final_filter.ht"
)


def get_pbt_trio_matrix(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> MatrixTableResource:
    """
    Get the phase-by-transmission (PBT) trio MatrixTable resource.

    One column per complete trio, with ``proband_entry`` / ``father_entry`` /
    ``mother_entry`` structs carrying ``GT``, ``adj`` and the transmission-
    phased ``PBT_GT``. Produced by ``--create-pbt-trio-matrix``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: PBT trio MatrixTable resource.
    """
    return MatrixTableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="pbt_trio_matrix",
            extension=".mt",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_pbt_mt(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> MatrixTableResource:
    """
    Get the exploded per-sample PBT MatrixTable resource.

    The PBT trio matrix exploded back into a sample MatrixTable (one column
    per trio membership): entries carry ``GT``, ``adj``, ``PBT_GT`` and the
    trio-level ``trio_adj`` flag, plus a ``source_trio`` column struct.
    Produced by ``--explode-pbt``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Exploded PBT MatrixTable resource.
    """
    return MatrixTableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="pbt",
            extension=".mt",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_pbt_multi_families_mt(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> MatrixTableResource:
    """
    Get the multi-offspring-family consensus PBT MatrixTable resource.

    For samples appearing in more than one trio with consistent parents,
    holds the consensus phased call (``consensus_gt``), the
    ``phase_concordance`` of the phased votes, and a ``discordant_gts`` flag.
    Produced by ``--phase-multi-families``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Multi-offspring-family consensus PBT MatrixTable resource.
    """
    return MatrixTableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="pbt_multi_families",
            extension=".mt",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_trio_variant_pair_list_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the trio variant-pair list Table resource.

    Unique variant pairs co-carried by trio probands (derived from the
    exploded PBT MatrixTable via ``create_variant_pair_ht``). This is the
    pair universe the trio-truth and gnomAD-comparison counts share.
    Produced by ``--derive-trio-vps``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Trio variant-pair list Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="trio_variant_pairs",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_trio_phase_counts_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the trio phase-count Table resource (trio "truth").

    One row per trio variant pair with ``raw`` / ``adj`` structs holding
    ``n_same_hap`` and ``n_chet`` — the count of probands phased cis vs trans
    by transmission. Produced by ``--call-trio-chet``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Trio phase-count Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="trio_phase_counts",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_gnomad_no_pbt_counts_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the gnomAD-minus-PBT genotype-counts Table resource.

    ``gt_counts_raw`` / ``gt_counts_adj`` for the trio variant pairs computed
    over release samples with all PBT trio members removed (so the trios' own
    genotypes don't leak into the population estimate). Produced by
    ``--gnomad-counts-no-pbt``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: gnomAD-minus-PBT genotype-counts Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="gnomad_no_pbt.genotype_counts",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_gnomad_no_pbt_phased_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the EM-phased gnomAD-minus-PBT Table resource.

    The gnomAD-minus-PBT counts annotated with ``em`` / ``em_plus_one``
    haplotype-EM phase estimates (via
    :func:`gnomad_chets.v4.phase_gnomad.get_phased_gnomad_ht`). Produced by
    ``--phase-gnomad-counts``.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: EM-phased gnomAD-minus-PBT Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="gnomad_no_pbt.phased",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_trio_comparison_ht(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> TableResource:
    """
    Get the trio-vs-gnomAD phase comparison Table resource.

    Joins the trio truth (cis/trans calls) with the gnomAD-minus-PBT EM
    ``p_chet`` for each pair, plus distance / gene annotations. Produced by
    ``--export-comparison`` (which also writes a flattened ``.tsv`` alongside).

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Trio-vs-gnomAD comparison Table resource.
    """
    return TableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="trio_comparison",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )

########################################################################################
### Pipeline Resource Collections
########################################################################################


def get_variant_pair_resources(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
    overwrite: bool = False,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant co-occurrence pipeline.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use test resources.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the variant
        co-occurrence pipeline.
    """
    # Initialize variant co-occurrence pipeline resource collection.
    vp_pipeline = PipelineResourceCollection(
        pipeline_name="variant_cooccurrence",
        overwrite=overwrite,
    )

    # Create resource collection for assembling the per-variant sites HT.
    # This joins every annotation source (filters, freq, vep, an_pct,
    # splice predictors, clinvar) into one wide HT that downstream steps
    # consume.
    preprocess_sites_ht = PipelineStepResourceCollection(
        "--preprocess-sites-ht",
        input_resources={
            "v4 QC + reference resources": {
                "filter_ht": TableResource(FINAL_FILTER_ONLY_FILTERS_PATH),
                "release_filter_ht": TableResource(FINAL_FILTER_PATH),
                "freq_ht": get_freq(data_type=data_type),
                "vep_ht": get_vep(data_type=data_type),
                "an_ht": all_sites_an(data_type=data_type),
                "spliceai_ht": get_insilico_predictors("spliceai"),
                "pangolin_ht": get_insilico_predictors("pangolin"),
                "clinvar_ht": clinvar.versions[CLINVAR_VERSION],
            },
        },
        output_resources={
            "sites_ht": get_sites_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Create resource collection for creating variant filter Table.
    create_variant_filter_ht = PipelineStepResourceCollection(
        "--create-variant-filter-ht",
        pipeline_input_steps=[preprocess_sites_ht],
        output_resources={
            "variant_filter_ht": get_variant_filter_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Create resource collection for filtering VariantDataset.
    filter_vmt = PipelineStepResourceCollection(
        "--filter-vmt",
        pipeline_input_steps=[create_variant_filter_ht],
        output_resources={
            "filtered_vmt": get_filtered_vmt(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Create resource collection for creating variant co-occurrence list.
    create_vp_list = PipelineStepResourceCollection(
        "--create-variant-pair-list-ht",
        pipeline_input_steps=[create_variant_filter_ht, filter_vmt],
        output_resources={
            "vp_list_ht": get_variant_pair_list_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Create resource collection for creating variant pair genotype counts Table
    # directly from the variant pair list (--encode-genotypes densifies the
    # pair-list variants out of the VDS transiently; no persisted dense MT).
    create_vp_gt_counts_ht = PipelineStepResourceCollection(
        "--create-variant-pair-genotype-counts-ht",
        pipeline_input_steps=[create_vp_list],
        output_resources={
            "vp_gt_counts_ht": get_variant_pair_genotype_counts_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Variant size-info Table (one row per encoded variant with the
    # diagnostic columns used to redo heavy/light routing at runtime).
    build_variant_size_info_ht = PipelineStepResourceCollection(
        "--build-variant-size-info",
        pipeline_input_steps=[create_vp_list],
        output_resources={
            "variant_size_info_ht": get_variant_size_info_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Excluded-genes Table — gene_ids dropped from both light and heavy.
    build_excluded_genes_ht = PipelineStepResourceCollection(
        "--exclude-gene-ids",
        pipeline_input_steps=[],
        output_resources={
            "excluded_genes_ht": get_excluded_genes_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Add all steps to the variant co-occurrence pipeline resource collection.
    vp_pipeline.add_steps(
        {
            "preprocess_sites_ht": preprocess_sites_ht,
            "create_variant_filter_ht": create_variant_filter_ht,
            "filter_vmt": filter_vmt,
            "create_variant_pair_list_ht": create_vp_list,
            "create_variant_pair_genotype_counts_ht": create_vp_gt_counts_ht,
            "build_variant_size_info_ht": build_variant_size_info_ht,
            "build_excluded_genes_ht": build_excluded_genes_ht,
        }
    )

    return vp_pipeline


def get_phasing_resources(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
    overwrite: bool = False,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the phasing pipeline.
    """
    # Initialize variant co-occurrence pipeline resource collection.
    phasing_pipeline = PipelineResourceCollection(
        pipeline_name="phasing",
        overwrite=overwrite,
    )
    

    # Create resource collection for creating variant pair genotype counts Table.
    create_phase = PipelineStepResourceCollection(
        "--phase",
        output_resources={
            "phase": get_phase(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Attach per-variant sites annotations (source tags, VEP, SpliceAI,
    # Pangolin, freq, ClinVar) onto the phased pair HT.
    annotate_phased_ht = PipelineStepResourceCollection(
        "--annotate-phased-ht",
        input_resources={
            "v4 phased + sites + variant-filter HTs": {
                "phased": get_phase(
                    data_type=data_type,
                    test=test,
                    tmp_dir=tmp_dir,
                    output_postfix=output_postfix,
                ),
                "sites_ht": get_sites_ht(
                    data_type=data_type,
                    test=test,
                    tmp_dir=tmp_dir,
                    output_postfix=output_postfix,
                ),
                "variant_filter_ht": get_variant_filter_ht(
                    data_type=data_type,
                    test=test,
                    tmp_dir=tmp_dir,
                    output_postfix=output_postfix,
                ),
            },
        },
        output_resources={
            "annotated_phase": get_annotated_phase(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Add all steps to the variant co-occurrence pipeline resource collection.
    phasing_pipeline.add_steps(
        {
            "phase": create_phase,
            "annotate_phased_ht": annotate_phased_ht,
        }
    )

    return phasing_pipeline


def get_trio_phasing_resources(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
    overwrite: bool = False,
    trio_set: str = "pedigree",
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for the trio-phasing pipeline.

    Trio side: ``--create-pbt-trio-matrix`` builds the PBT trio MatrixTable
    from the gnomAD VDS + finalized pedigree, ``--explode-pbt`` explodes it
    into a per-sample MatrixTable, ``--phase-multi-families`` derives the
    multi-offspring-family consensus, ``--derive-trio-vps`` builds the trio
    variant-pair list, and ``--call-trio-chet`` counts cis/trans probands per
    pair (trio truth).

    Comparison side: ``--gnomad-counts-no-pbt`` computes gnomAD genotype counts
    over release samples with PBT members removed, ``--phase-gnomad-counts``
    EM-phases them, and ``--export-comparison`` joins trio truth vs gnomAD
    ``p_chet``.

    The ``trio_set`` value is appended to every trio-phasing output path so
    that the ``pedigree`` (all trios) and ``trios`` (one per family) runs never
    overwrite each other. The shared variant-filter HT input keeps the plain
    ``output_postfix`` (it is trio-set independent, built by ``create_vp_list``).

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use test resources.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :param overwrite: Whether to overwrite resources if they exist.
    :param trio_set: Which pedigree resource the run uses — ``"pedigree"`` (all
        trios, multiple per family) or ``"trios"`` (one per family). Folded into
        the trio-phasing output paths.
    :return: PipelineResourceCollection for the trio-phasing pipeline.
    """
    trio_pipeline = PipelineResourceCollection(
        pipeline_name="trio_phasing",
        overwrite=overwrite,
    )

    # Qualify every trio output with the trio set so pedigree/trios runs are
    # isolated; the variant-filter input stays on the plain postfix.
    trio_postfix = (
        f"{output_postfix}.{trio_set}" if output_postfix is not None else trio_set
    )

    # Build the PBT-phased trio MatrixTable from the gnomAD VDS + pedigree.
    create_pbt_trio_matrix = PipelineStepResourceCollection(
        "--create-pbt-trio-matrix",
        output_resources={
            "pbt_trio_matrix": get_pbt_trio_matrix(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # Explode the trio matrix back into a per-sample MatrixTable.
    explode_pbt = PipelineStepResourceCollection(
        "--explode-pbt",
        pipeline_input_steps=[create_pbt_trio_matrix],
        output_resources={
            "pbt_mt": get_pbt_mt(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # Consensus phase across families with multiple offspring.
    phase_multi_families = PipelineStepResourceCollection(
        "--phase-multi-families",
        pipeline_input_steps=[explode_pbt],
        output_resources={
            "pbt_multi_families_mt": get_pbt_multi_families_mt(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # Derive the trio variant-pair list from the exploded PBT MT. Also reads
    # the variant filter HT (built by create_vp_list --create-variant-filter-ht
    # for the same postfix) for gene_id / an_pct.
    derive_trio_vps = PipelineStepResourceCollection(
        "--derive-trio-vps",
        pipeline_input_steps=[explode_pbt],
        add_input_resources={
            "create_vp_list.py --create-variant-filter-ht": {
                "variant_filter_ht": get_variant_filter_ht(
                    data_type=data_type,
                    test=test,
                    tmp_dir=tmp_dir,
                    output_postfix=output_postfix,
                )
            }
        },
        output_resources={
            "trio_vp_list_ht": get_trio_variant_pair_list_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # Count cis/trans probands per pair (trio truth).
    call_trio_chet = PipelineStepResourceCollection(
        "--call-trio-chet",
        pipeline_input_steps=[explode_pbt, derive_trio_vps],
        output_resources={
            "trio_phase_counts_ht": get_trio_phase_counts_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # gnomAD genotype counts on the trio VPs, release samples minus PBT members.
    # Inputs include explode_pbt so the reuse path can read the densified PBT MT
    # (pbt_mt) to compute the PBT∩release subtraction without re-densifying.
    gnomad_counts_no_pbt = PipelineStepResourceCollection(
        "--gnomad-counts-no-pbt",
        pipeline_input_steps=[derive_trio_vps, explode_pbt],
        output_resources={
            "gnomad_no_pbt_counts_ht": get_gnomad_no_pbt_counts_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # EM-phase the gnomAD-minus-PBT counts.
    phase_gnomad_counts = PipelineStepResourceCollection(
        "--phase-gnomad-counts",
        pipeline_input_steps=[gnomad_counts_no_pbt],
        output_resources={
            "gnomad_no_pbt_phased_ht": get_gnomad_no_pbt_phased_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    # Join trio truth vs gnomAD EM phase and export.
    export_comparison = PipelineStepResourceCollection(
        "--export-comparison",
        pipeline_input_steps=[call_trio_chet, phase_gnomad_counts],
        output_resources={
            "trio_comparison_ht": get_trio_comparison_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=trio_postfix,
            )
        },
    )

    trio_pipeline.add_steps(
        {
            "create_pbt_trio_matrix": create_pbt_trio_matrix,
            "explode_pbt": explode_pbt,
            "phase_multi_families": phase_multi_families,
            "derive_trio_vps": derive_trio_vps,
            "call_trio_chet": call_trio_chet,
            "gnomad_counts_no_pbt": gnomad_counts_no_pbt,
            "phase_gnomad_counts": phase_gnomad_counts,
            "export_comparison": export_comparison,
        }
    )

    return trio_pipeline
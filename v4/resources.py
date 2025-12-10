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
from gnomad_qc.v4.resources.annotations import get_freq, get_vep
from gnomad_qc.v4.resources.variant_qc import final_filter

########################################################################################
### Constants
########################################################################################

DEFAULT_TMP_DIR = "gs://gnomad-tmp-4day"
"""Default temporary directory for variant co-occurrence pipeline output files."""

VARIANT_COOCCURRENCE_ROOT = "gs://gnomad/v4.1/variant_cooccurrence"
"""Official output root directory for variant co-occurrence pipeline output files."""

TEST_INTERVAL = "chr21:46324141-46445769"
"""Test interval for PCNT gene used in testing mode."""

DATA_TYPE_CHOICES = ["exomes", "genomes"]
"""Valid data type choices for variant co-occurrence pipeline."""

DEFAULT_DATA_TYPE = "exomes"
"""Default data type for variant co-occurrence pipeline."""

DEFAULT_MAX_FREQ = 0.05
"""Default maximum global AF to keep (inclusive)."""

DEFAULT_LEAST_CONSEQUENCE = "3_prime_UTR_variant"
"""Default lowest-severity consequence to keep."""


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


def get_filtered_vds(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> VariantDatasetResource:
    """
    Get filtered VariantDataset resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Filtered VariantDataset resource.
    """
    return VariantDatasetResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="filtered_vds",
            extension=".vds",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_vp_list_ht(
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
            resource_name="vp_list",
            extension=".ht",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_filtered_dense_mt(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> MatrixTableResource:
    """
    Get filtered dense MatrixTable resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Filtered dense MatrixTable resource.
    """
    return MatrixTableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="filtered.dense",
            extension=".mt",
            test=test,
            tmp_dir=tmp_dir,
            output_postfix=output_postfix,
        )
    )


def get_vp_full_mt(
    data_type: str = DEFAULT_DATA_TYPE,
    test: bool = False,
    tmp_dir: Optional[str] = None,
    output_postfix: Optional[str] = None,
) -> MatrixTableResource:
    """
    Get full variant pair MatrixTable resource.

    :param data_type: Data type to use. Must be one of 'exomes' or 'genomes'.
    :param test: Whether to use a tmp path for testing.
    :param tmp_dir: Temporary directory for output files.
    :param output_postfix: Postfix to append to output file names.
    :return: Full variant pair MatrixTable resource.
    """
    return MatrixTableResource(
        _get_resource_path(
            data_type=data_type,
            resource_name="vp_full",
            extension=".mt",
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

    # Create resource collection for creating variant filter Table.
    create_variant_filter_ht = PipelineStepResourceCollection(
        "--create-variant-filter-ht",
        input_resources={
            "v4 QC resources": {
                "filter_ht": final_filter(data_type=data_type),
                "freq_ht": get_freq(data_type=data_type),
                "vep_ht": get_vep(data_type=data_type),
            },
        },
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
    filter_vds = PipelineStepResourceCollection(
        "--filter-vds",
        pipeline_input_steps=[create_variant_filter_ht],
        output_resources={
            "filtered_vds": get_filtered_vds(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Create resource collection for creating variant co-occurrence list.
    create_vp_list = PipelineStepResourceCollection(
        "--create-vp-list",
        pipeline_input_steps=[create_variant_filter_ht, filter_vds],
        output_resources={
            "vp_list_ht": get_vp_list_ht(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    create_dense_filtered_mt = PipelineStepResourceCollection(
        "--create-dense-filtered-mt",
        pipeline_input_steps=[filter_vds, create_vp_list],
        output_resources={
            "dense_filtered_mt": get_filtered_dense_mt(
                data_type=data_type,
                test=test,
                tmp_dir=tmp_dir,
                output_postfix=output_postfix,
            )
        },
    )

    # Create resource collection for creating full variant co-occurrence VDS.
    create_full_vp = PipelineStepResourceCollection(
        "--create-full-vp",
        pipeline_input_steps=[create_dense_filtered_mt, create_vp_list],
        output_resources={
            "vp_full_mt": get_vp_full_mt(
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
            "create_variant_filter_ht": create_variant_filter_ht,
            "filter_vds": filter_vds,
            "create_vp_list": create_vp_list,
            "create_dense_filtered_mt": create_dense_filtered_mt,
            "create_full_vp": create_full_vp,
        }
    )

    return vp_pipeline

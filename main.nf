#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Parameters ---
// New, descriptive parameters for launching the pipeline.
params.gcs_tmp_dir = null
params.gcs_data_path = null

// Original parameters that the script blocks expect. We map the new params to them.
params.cluster   = null
params.project   = null
params.tmp_dir   = params.gcs_tmp_dir
params.data_path = params.gcs_data_path


log.info """
--- Hail Dataproc Pipeline (Success Flag Pattern) ---
Cluster     : ${params.cluster}
Project     : ${params.project}
Temp Dir    : ${params.tmp_dir}
Data Path   : ${params.data_path}
-----------------------------------------------------
"""

// --- Workflow ---
workflow {
    // --- This section robustly prepares all inputs ---
    def py_module_files = [ file('resources.py'), file('chet_utils.py') ]
    def main_script_file = file('create_vp_matrix.py')

    def initial_config_tuple = tuple(
        params.cluster,
        params.project,
        main_script_file,
        py_module_files
    )
    ch_initial_config = Channel.of( initial_config_tuple )


    // --- This section correctly chains the processes using their success flags ---
    create_vp_list( ch_initial_config )
    create_full_vp( create_vp_list.out.vp_list_success )
    create_vp_summary( create_full_vp.out.full_vp_success )

    // You can view the final success flag to know the pipeline is done.
    create_vp_summary.out.summary_success.view { "✅ Pipeline finished successfully. Final flag: ${it}" }
}

// --- Processes ---

process create_vp_list {
    executor 'local'

    input:
    // This receives the initial configuration.
    tuple val(cluster), val(project), path(main_script), val(pyfiles)

    output:
    // The output forwards the config AND adds the new success file for the next step.
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path("vp_list.success"), emit: vp_list_success

    script:
    // This script block IS NOT CHANGED, except for adding `touch` at the end.
    """
    # 1. Run the original hailctl command.
    hailctl dataproc submit ${cluster} ${main_script} \\
        --exomes \\
        --least_consequence 3_prime_UTR_variant \\
        --max_freq 0.05 \\
        --chrom 21 \\
        --create_vp_list \\
        --tmp_dir ${params.tmp_dir} \\
        --gnomad_data_path ${params.data_path} \\
        --project ${project} \\
        --testing \\
        --pyfiles ${pyfiles.join(',')}

    # 2. If the command above succeeds, create a flag file. This is our "proof" of completion.
    touch vp_list.success
    """
}

process create_full_vp {
    executor 'local'

    input:
    // This requires the output from the previous step, including its success flag.
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path(previous_success_flag)

    output:
    // This creates and outputs its OWN success flag, while forwarding the config.
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path("full_vp.success"), emit: full_vp_success

    script:
    // This script block IS IDENTICAL to your original request, plus `touch`.
    """
    hailctl dataproc submit ${cluster} ${main_script} \\
        --exomes \\
        --least_consequence 3_prime_UTR_variant \\
        --max_freq 0.05 \\
        --chrom 21 \\
        --create_full_vp \\
        --tmp_dir ${params.tmp_dir} \\
        --gnomad_data_path ${params.data_path} \\
        --project ${project} \\
        --testing \\
        --pyfiles ${pyfiles.join(',')}

    touch full_vp.success
    """
}

process create_vp_summary {
    executor 'local'

    input:
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path(previous_success_flag)

    output:
    // This is the final step, so it only needs to output its success flag.
    path "summary.success", emit: summary_success

    script:
    // This script block IS IDENTICAL to your original request, plus `touch`.
    """
    hailctl dataproc submit ${cluster} ${main_script} \\
        --exomes \\
        --least_consequence 3_prime_UTR_variant \\
        --max_freq 0.05 \\
        --chrom 21 \\
        --create_vp_summary \\
        --tmp_dir ${params.tmp_dir} \\
        --gnomad_data_path ${params.data_path} \\
        --project ${project} \\
        --testing \\
        --pyfiles ${pyfiles.join(',')}

    touch summary.success
    """
}
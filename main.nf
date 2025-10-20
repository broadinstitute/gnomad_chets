#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
This script is a refactored version that follows Nextflow best practices.
*/

// --- Parameters ---
params.cluster   = null
params.project   = null
params.tmp_dir   = null
params.data_path = null
params.outdir    = './results'

log.info """
--- Hail Dataproc Pipeline ---
Cluster     : ${params.cluster}
Project     : ${params.project}
Temp Dir    : ${params.tmp_dir}
Data Path   : ${params.data_path}
Output Dir  : ${params.outdir}
------------------------------
"""

// --- Workflow ---
workflow {
    ch_config = Channel.of(
        tuple(params.cluster, params.project, params.tmp_dir, params.data_path)
    )

    // *** CHANGE 1: Create a channel for the main Python script ***
    ch_main_script = Channel.fromPath('create_vp_matrix.py')

    ch_pyfiles = Channel.fromPath('{resources,chet_utils}.py').collect()

    // *** CHANGE 2: Pass all three channels to the process ***
    create_vp_list( ch_config, ch_main_script, ch_pyfiles )

    create_vp_list.out.vp_list_ht.view { ht_path -> "Successfully created Hail Table at: ${ht_path}" }
}

// --- Process ---
process create_vp_list {
    executor 'local'
    publishDir "${params.outdir}/vp_list", mode: 'copy'

    // *** CHANGE 3: Declare the new input channel for the main script ***
    input:
    tuple val(cluster), val(project), val(tmp_dir), val(data_path)
    path main_script
    path pyfiles

    output:
    path "exomes_v2_list_PCNT.ht", type: 'dir', emit: vp_list_ht

    script:
    // *** CHANGE 4: Use the main_script variable in the command ***
    """
    hailctl dataproc submit ${cluster} ${main_script} \\
        --exomes \\
        --least_consequence 3_prime_UTR_variant \\
        --max_freq 0.05 \\
        --chrom 21 \\
        --create_vp_list \\
        --tmp_dir ${tmp_dir} \\
        --gnomad_data_path ${data_path} \\
        --project ${project} \\
        --overwrite \\
        --pyfiles ${pyfiles.join(',')}
    """
}
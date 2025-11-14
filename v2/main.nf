#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Parameters ---
params.gcs_tmp_dir = null
params.gcs_data_path = null
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
    // --- Define all file inputs at the top level ---
    def main_py_modules = [ file('resources.py'), file('chet_utils.py') ]
    def main_script_file = file('create_vp_matrix.py')

    // <<< FIX #1: Define the files needed for the phasing step separately.
    def phasing_script_file = file('compute_gnomad_phase.py')
    def phasing_py_modules = [ file('resources.py'), file('chet_utils.py'), file('phasing.py') ]


    // --- Create Channels for process inputs ---
    def initial_config_tuple = tuple( params.cluster, params.project, main_script_file, main_py_modules )
    ch_initial_config = Channel.of( initial_config_tuple )

    // <<< FIX #2: Create a separate channel for the phasing step's files.
    ch_phasing_files = Channel.of( tuple(phasing_script_file, phasing_py_modules) )


    // --- Chain the processes using their success flags ---
    create_vp_list( ch_initial_config )
    create_full_vp( create_vp_list.out.vp_list_success )
    create_vp_summary( create_full_vp.out.full_vp_success )

    // <<< FIX #3: Call the final process by providing BOTH the success channel from the
    // previous step AND the new channel containing its specific files.
    create_phased_ht(
        create_vp_summary.out.summary_success,
        ch_phasing_files
    )

    create_phased_ht.out.phasing_success.view { "✅ Pipeline finished successfully. Final flag: ${it}" }
}

// --- Processes ---

process create_vp_list {
    executor 'local'
    input:
    tuple val(cluster), val(project), path(main_script), val(pyfiles)
    output:
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path("vp_list.success"), emit: vp_list_success
    script:
    """
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
    touch vp_list.success
    """
}

process create_full_vp {
    executor 'local'
    input:
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path(previous_success_flag)
    output:
    tuple val(cluster), val(project), path(main_script), val(pyfiles), path("full_vp.success"), emit: full_vp_success
    script:
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
    // This part from the previous fix remains correct.
    tuple val(cluster), val(project), path("summary.success"), emit: summary_success
    script:
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

process create_phased_ht {
    tag "Phasing (chr21)"
    executor 'local'

    // <<< FIX #4: Define a composite input. It now takes two tuples:
    // 1. The success signal from the previous step.
    // 2. The specific script and pyfiles it needs to run.
    input:
    tuple val(cluster), val(project), path(previous_success_flag)
    tuple path(phasing_script), val(phasing_pyfiles)

    output:
    path "phasing.success", emit: phasing_success

    script:
    // This block is now much cleaner. It uses the input variables directly.
    // The `phasing_script` and `phasing_pyfiles` variables are now populated
    // by Nextflow from the channel inputs, and the files are guaranteed to exist.
    def infile = "${params.tmp_dir}/exomes_gnomad_v2_PCNT_summary.ht"
    def outfile = "${params.tmp_dir}/exomes_gnomad_v2_PCNT_phased.ht"
    """
    echo "Starting phasing step..."
    hailctl dataproc submit ${cluster} ${phasing_script} \\
        --exome \\
        --least_consequence 3_prime_UTR_variant \\
        --max_freq 0.05 \\
        --chrom 21 \\
        --create_phased_vp_summary \\
        --no_lr \\
        --no_shr \\
        --infile ${infile} \\
        --outfile ${outfile} \\
        --pyfiles ${phasing_pyfiles.join(',')}

    touch phasing.success
    """
}
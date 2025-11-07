cluster="rungar"
project="rungar-sandbox-d7e7"

DATA_TYPES=["exomes"]
VERSIONS=["v2"]
GENES=["CBR3"]

rule all:
    input:
        expand("/Users/raungar/Documents/SamochaLab/CompoundHet/Logs/3.{version}_{gene}_summary.txt",
               version=VERSIONS, gene=GENES)

rule create_vp_list:
    params:
        mt="gs://rungar-sandbox-tmp-4day/gnomad_{version}_{gene}.mt",
        tmp_dir="gs://rungar-sandbox-tmp-month",
        least_consequence="3_prime_UTR_variant",
        max_freq=0.05,
        chrom=21
    output:
        "/Users/raungar/Documents/SamochaLab/CompoundHet/Logs/1.{version}_{gene}_list.txt"
    shell:
        """
        hailctl dataproc submit {cluster} create_vp_matrix.py \
            --exome \
            --least_consequence {params.least_consequence} \
            --max_freq {params.max_freq} \
            --chrom {params.chrom} \
            --create_vp_list \
            --tmp_dir {params.tmp_dir} \
            --gnomad_data_path {params.mt} \
            --project {project} \
            --testing \
            --pyfiles resources.py,chet_utils.py
        touch {output}
        """

rule create_full_vp:
    input:
        "/Users/raungar/Documents/SamochaLab/CompoundHet/Logs/1.{version}_{gene}_list.txt"
    params:
        infile="gs://rungar-sandbox-tmp-4day/gnomad_{version}_{gene}.mt",
        tmp_dir="gs://rungar-sandbox-tmp-month",
        least_consequence="3_prime_UTR_variant",
        max_freq=0.05,
        chrom=21,
        name="gnomad_{version}_{gene}"
    output:
        "/Users/raungar/Documents/SamochaLab/CompoundHet/Logs/2.{version}_{gene}_full.txt"
    shell:
        """
        hailctl dataproc submit {cluster} create_vp_matrix.py \
            --exomes \
            --least_consequence {params.least_consequence} \
            --max_freq {params.max_freq} \
            --chrom {params.chrom} \
            --create_full_vp \
            --tmp_dir {params.tmp_dir} \
            --gnomad_data_path {params.infile} \
            --project {project} \
            --name {params.name} \
            --testing \
            --pyfiles resources.py,chet_utils.py 
        touch {output}
        """

rule create_vp_summary:
    input:
        "/Users/raungar/Documents/SamochaLab/CompoundHet/Logs/2.{version}_{gene}_full.txt"
    params:
        infile="gs://rungar-sandbox-tmp-4day/gnomad_{version}_{gene}.mt",
        tmp_dir="gs://rungar-sandbox-tmp-month",
        least_consequence="3_prime_UTR_variant",
        max_freq=0.05,
        chrom=21,
        name="gnomad_{version}_{gene}"
    output:
        "/Users/raungar/Documents/SamochaLab/CompoundHet/Logs/3.{version}_{gene}_summary.txt"
    shell:
        """
        hailctl dataproc submit {cluster} create_vp_matrix.py \
            --exomes \
            --least_consequence {params.least_consequence} \
            --max_freq {params.max_freq} \
            --chrom {params.chrom} \
            --create_vp_summary \
            --tmp_dir {params.tmp_dir} \
            --gnomad_data_path {params.infile} \
            --project {project} \
            --name {params.name} \
            --testing \
            --pyfiles resources.py,chet_utils.py 
        touch {output}
        """  
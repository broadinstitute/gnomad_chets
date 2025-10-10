
cluster="rungar"

rule all:
    input: storage.gcs("gs://rungar-sandbox-tmp-month/test/test.txt")

rule write_gcs:
    output: storage.gcs("gs://rungar-sandbox-tmp-month/test/test.txt")
    shell:
        """
        hailctl dataproc submit {cluster} test.py > {output}
        """


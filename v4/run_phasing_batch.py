"""Phase a region with Eagle, Beagle or SHAPEIT5 via Hail Batch.

Generalizes :mod:`gnomad_chets.v4.run_eagle_batch` to the three statistical
phasers so the same exported VCF can be run through each and scored the same
way by :mod:`gnomad_chets.v4.phase_from_vcf`. All three come from one image
(``docker/phasing_tools/Dockerfile``) sharing one set of GRCh38 genetic maps,
so a difference in the output is a difference between tools rather than
between environments.

None of the three is given a reference panel: each phases the cohort against
itself, which is what the production local-ancestry pipeline does and what
makes sense when the cohort supplies its own haplotypes.

SHAPEIT5 is run in two stages when ``--shapeit-rare`` is set: ``phase_common``
builds a scaffold from variants above ``--shapeit-common-maf``, then
``phase_rare`` places the remainder onto it. That two-stage design is the whole
reason SHAPEIT5 is interesting here -- it is the only one of the three built
for the rare end of the spectrum, which is where compound-heterozygote
candidates live.

Usage::

    python v4/run_phasing_batch.py --tool beagle \\
      --vcf gs://gnomad-tmp-30day/pcnt_region.genomes.hq10k.vcf.bgz \\
      --contig chr21 \\
      --out-prefix gs://gnomad-tmp-30day/phasing/pcnt.hq10k.beagle
"""

import argparse
import logging

import hailtop.batch as hb

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_phasing_batch")
logger.setLevel(logging.INFO)

PHASING_IMAGE = "gcr.io/broad-mpg-gnomad/phasing_tools:5.1.1"


def eagle_command(j, vcf, contig: str, cpu: int) -> str:
    """Eagle v2.4.1, non-reference mode.

    :param j: Batch job, whose ``ofile`` receives the phased VCF.
    :param vcf: Localized input resource group.
    :param contig: Contig name, unused by Eagle but kept for a uniform signature.
    :param cpu: Threads.
    :return: Shell command.
    """
    return f"""
    $EAGLE_BIN \\
        --geneticMapFile $EAGLE_MAP \\
        --numThreads {cpu} \\
        --outPrefix {j.prefix} \\
        --vcfOutFormat z \\
        --vcf {vcf['vcf.gz']}
    mv {j.prefix}.vcf.gz {j.ofile['vcf.gz']}
    """


def beagle_command(j, vcf, contig: str, cpu: int) -> str:
    """Beagle 5.5, phasing only.

    ``impute=false`` matters: Beagle imputes ungenotyped markers by default, and
    imputed genotypes are not observations. Leaving it on would add sites the
    other two tools never emit and make the pair sets incomparable.

    :param j: Batch job.
    :param vcf: Localized input resource group.
    :param contig: Contig, used to select the genetic map.
    :param cpu: Threads.
    :return: Shell command.
    """
    return f"""
    java -Xmx$(( {cpu} * 6 ))g -jar $BEAGLE_JAR \\
        gt={vcf['vcf.gz']} \\
        map=$BEAGLE_MAP_DIR/{contig}.GRCh38.map \\
        out={j.prefix} \\
        impute=false \\
        nthreads={cpu}
    mv {j.prefix}.vcf.gz {j.ofile['vcf.gz']}
    """


def shapeit_command(j, vcf, contig: str, cpu: int, common_maf: float, rare: bool) -> str:
    """SHAPEIT5, one or two stages.

    :param j: Batch job.
    :param vcf: Localized input resource group.
    :param contig: Contig, used to select the genetic map and as the region.
    :param cpu: Threads.
    :param common_maf: MAF floor for the common-variant scaffold.
    :param rare: Run ``phase_rare`` onto the scaffold as a second stage.
    :return: Shell command.
    """
    common = f"""
    SHAPEIT5_phase_common \\
        --input {vcf['vcf.gz']} \\
        --map $SHAPEIT_MAP_DIR/{contig}.b38.gmap \\
        --region {contig} \\
        --filter-maf {common_maf} \\
        --thread {cpu} \\
        --output scaffold.bcf
    bcftools index -f scaffold.bcf
    """
    if not rare:
        return common + f"""
    bcftools view -Oz -o {j.ofile['vcf.gz']} scaffold.bcf
    """
    return common + f"""
    SHAPEIT5_phase_rare \\
        --input-plain {vcf['vcf.gz']} \\
        --scaffold scaffold.bcf \\
        --map $SHAPEIT_MAP_DIR/{contig}.b38.gmap \\
        --input-region {contig} \\
        --scaffold-region {contig} \\
        --thread {cpu} \\
        --output phased.bcf
    bcftools view -Oz -o {j.ofile['vcf.gz']} phased.bcf
    """


def add_phasing_job(
    batch: hb.Batch,
    vcf: hb.ResourceGroup,
    tool: str,
    contig: str,
    name: str,
    cpu: int,
    memory: str,
    storage: str,
    image: str,
    shapeit_common_maf: float,
    shapeit_rare: bool,
) -> hb.batch.job.Job:
    """Add the phasing job for one tool.

    :param batch: Batch to add to.
    :param vcf: Localized input, with its tabix index alongside.
    :param tool: ``eagle``, ``beagle`` or ``shapeit``.
    :param contig: Contig being phased.
    :param name: Job name.
    :param cpu: Cores, also passed to the tool as its thread count.
    :param memory: Batch memory tier.
    :param storage: Disk request.
    :param image: Docker image carrying all three tools.
    :param shapeit_common_maf: SHAPEIT5 scaffold MAF floor.
    :param shapeit_rare: Run SHAPEIT5's rare stage.
    :return: The job, whose ``ofile`` holds the phased VCF.
    """
    j = batch.new_job(name=name)
    j.image(image)
    j.cpu(cpu)
    j.memory(memory)
    j.storage(storage)
    # Non-preemptible: a multi-hour phasing run restarted from scratch on
    # preemption costs more than the spot discount saves.
    j.spot(False)
    j.declare_resource_group(ofile={"vcf.gz": "{root}.vcf.gz"})

    if tool == "eagle":
        cmd = eagle_command(j, vcf, contig, cpu)
    elif tool == "beagle":
        cmd = beagle_command(j, vcf, contig, cpu)
    else:
        cmd = shapeit_command(j, vcf, contig, cpu, shapeit_common_maf, shapeit_rare)

    j.command("set -euxo pipefail")
    j.command(cmd)
    return j


def main(args):
    """Build and run the Batch."""
    backend = hb.ServiceBackend(
        billing_project=args.billing_project, remote_tmpdir=args.remote_tmpdir
    )
    batch = hb.Batch(backend=backend, name=args.batch_name or f"{args.tool}-phasing")

    vcf = batch.read_input_group(
        **{"vcf.gz": args.vcf, "vcf.gz.tbi": f"{args.vcf}.tbi"}
    )

    j = add_phasing_job(
        batch,
        vcf,
        tool=args.tool,
        contig=args.contig,
        name=args.batch_name or f"{args.tool}-phasing",
        cpu=args.cpu,
        memory=args.memory,
        storage=args.storage,
        image=args.image,
        shapeit_common_maf=args.shapeit_common_maf,
        shapeit_rare=args.shapeit_rare,
    )
    batch.write_output(j.ofile, args.out_prefix)
    logger.info("Phased VCF will be written to %s.vcf.gz", args.out_prefix)

    batch.run(wait=args.wait)
    backend.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tool", required=True, choices=["eagle", "beagle", "shapeit"]
    )
    parser.add_argument("--vcf", required=True, help="Input VCF (.bgz/.vcf.gz).")
    parser.add_argument(
        "--contig", required=True, help="Contig, e.g. chr21. Selects the genetic map."
    )
    parser.add_argument(
        "--out-prefix", required=True, help="Output prefix; '.vcf.gz' is appended."
    )
    parser.add_argument("--cpu", type=int, default=16, help="Cores and tool threads.")
    parser.add_argument("--memory", default="highmem", help="Batch memory tier.")
    parser.add_argument("--storage", default="100G", help="Disk request.")
    parser.add_argument("--image", default=PHASING_IMAGE, help="Docker image.")
    parser.add_argument(
        "--shapeit-common-maf",
        type=float,
        default=0.001,
        help="MAF floor for the SHAPEIT5 common-variant scaffold.",
    )
    parser.add_argument(
        "--shapeit-rare",
        action="store_true",
        help="Run SHAPEIT5's phase_rare stage onto the scaffold.",
    )
    parser.add_argument("--billing-project", default="gnomad-production")
    parser.add_argument(
        "--remote-tmpdir",
        default="gs://gnomad-tmp-30day/batch-an-tests-julia/batch-tmp",
    )
    parser.add_argument("--batch-name")
    parser.add_argument(
        "--wait", action="store_true", help="Block until the batch finishes."
    )
    main(parser.parse_args())

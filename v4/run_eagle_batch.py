"""Run Eagle on a pre-exported VCF via Hail Batch.

Phases a region ourselves so a statistical phaser can be scored on genome data,
where it has the common-variant scaffold that exome capture denies it. The
companion export step is :mod:`gnomad_chets.v4.export_for_phasing`; the phased
output is consumed by :mod:`gnomad_chets.v4.phase_from_vcf`.

The job definition mirrors the production local-ancestry pipeline
(`broadinstitute/gnomad_local_ancestry`, ``batch/lai_batch_pipeline.py``): same
image, same Eagle invocation, same per-chromosome shape. Reusing it means the
phasing matches how gnomAD has actually run Eagle rather than a fresh guess at
the parameters.

Note there is no ``--vcfRef``. Eagle is run in non-reference mode, phasing the
cohort against itself, which is what the LAI pipeline does and what makes sense
when the cohort is already large enough to supply its own haplotypes.

Runtime scales roughly linearly in samples. The 2020 pilot took 3.25 h for
6,044 samples and 162k variants on 2 threads; ~136k samples and ~266k variants
on 16 threads is expected to take considerably longer, so the job is run
non-preemptible.

Usage::

    python v4/run_eagle_batch.py \\
      --vcf gs://gnomad-tmp-30day/pcnt_region.genomes.hq.vcf.bgz \\
      --out-prefix gs://gnomad-tmp-30day/eagle/pcnt_region.genomes.hq.phased
"""

import argparse
import logging

import hailtop.batch as hb

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_eagle_batch")
logger.setLevel(logging.INFO)

EAGLE_IMAGE = "gcr.io/broad-mpg-gnomad/lai_phasing:latest"
# Shipped inside the image alongside the binary.
GENETIC_MAP = "Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"


def add_eagle_job(
    batch: hb.Batch,
    vcf: hb.ResourceGroup,
    name: str,
    cpu: int,
    memory: str,
    storage: str,
    image: str,
) -> hb.batch.job.Job:
    """Add the Eagle phasing job.

    :param batch: Batch to add to.
    :param vcf: Localized input, with the tabix index alongside it.
    :param name: Job name.
    :param cpu: Cores, also passed to Eagle as ``--numThreads``.
    :param memory: Batch memory tier, e.g. ``highmem``.
    :param storage: Disk request; must hold the input plus the phased output.
    :param image: Docker image carrying the Eagle binary and genetic map.
    :return: The job, whose ``ofile`` resource group holds the phased VCF.
    """
    j = batch.new_job(name=name)
    j.image(image)
    j.cpu(cpu)
    j.memory(memory)
    j.storage(storage)
    # Non-preemptible: a multi-hour phasing run restarted from scratch on
    # preemption would cost more than the spot discount saves.
    j.spot(False)
    j.declare_resource_group(ofile={"vcf.gz": "{root}.vcf.gz"})

    j.command(f"""
        Eagle_v2.4.1/eagle \\
            --geneticMapFile {GENETIC_MAP} \\
            --numThreads {cpu} \\
            --outPrefix {j.ofile} \\
            --vcfOutFormat z \\
            --vcf {vcf['vcf.gz']}
        """)
    return j


def main(args):
    """Build and run the Batch."""
    backend = hb.ServiceBackend(
        billing_project=args.billing_project, remote_tmpdir=args.remote_tmpdir
    )
    batch = hb.Batch(backend=backend, name=args.batch_name)

    # Read the index in alongside the VCF so Eagle finds it under the expected
    # sibling name once localized.
    vcf = batch.read_input_group(
        **{"vcf.gz": args.vcf, "vcf.gz.tbi": f"{args.vcf}.tbi"}
    )

    j = add_eagle_job(
        batch,
        vcf,
        name=args.batch_name,
        cpu=args.cpu,
        memory=args.memory,
        storage=args.storage,
        image=args.image,
    )
    batch.write_output(j.ofile, args.out_prefix)
    logger.info("Phased VCF will be written to %s.vcf.gz", args.out_prefix)

    batch.run(wait=args.wait)
    backend.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input VCF (.bgz/.vcf.gz).")
    parser.add_argument(
        "--out-prefix", required=True, help="Output prefix; '.vcf.gz' is appended."
    )
    parser.add_argument("--cpu", type=int, default=16, help="Cores and Eagle threads.")
    parser.add_argument("--memory", default="highmem", help="Batch memory tier.")
    parser.add_argument("--storage", default="100G", help="Disk request.")
    parser.add_argument("--image", default=EAGLE_IMAGE, help="Docker image.")
    parser.add_argument("--billing-project", default="gnomad-production")
    parser.add_argument(
        "--remote-tmpdir",
        default="gs://gnomad-tmp-30day/batch-an-tests-julia/batch-tmp",
    )
    parser.add_argument("--batch-name", default="eagle-pcnt-region")
    parser.add_argument(
        "--wait", action="store_true", help="Block until the batch finishes."
    )
    main(parser.parse_args())

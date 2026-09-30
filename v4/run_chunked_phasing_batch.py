"""Phase a chromosome in overlapping chunks, one Batch, chunks in parallel.

Beagle phases a whole chromosome in one process, carrying haplotype state
across its internal sliding windows. At cohort scale that process runs out of
heap: full chr21 at 120,279 samples died with ``OutOfMemoryError: Java heap
space`` in its second window after 18h52m, at ``window=15`` on the largest heap
Hail Batch can give a container (104 GB, so ``-Xmx96g``).

Neither of the in-process levers fixes that. ``nthreads`` is shared-memory, so
16 threads work out of one heap and peak memory is unchanged; ``window`` was
already lowered from Beagle's default of 40 to the 15 All of Us used, and the
failure came *after* a window at that setting had succeeded.

What does fix it is a separate JVM per region. The chunk size here is not a
guess: ``s1-r2`` phased 120,279 samples over 2,095,026 markers standalone in
the same 96 GB heap at ``window=15`` ($18.61, 14h50m), while the *smaller*
1,959,476-marker window inside the whole-chromosome run OOMed. A standalone
~2M-marker chunk is a configuration that has already run.

Chunks also run concurrently, so wall-clock is one chunk rather than their sum.

Chunks overlap. For pair-level work the overlap means a pair can be scored from
whichever chunk contains both of its variants, which avoids ligation entirely --
pairs are within-gene and so far shorter than the overlap. ``bcftools concat
--ligate`` (in the image) remains available if a single chromosome-wide phased
VCF is wanted for something else.

The split is a single job: each phasing job then localizes only its own chunk
instead of all of them re-subsetting the full input.

Usage::

    python v4/run_chunked_phasing_batch.py \\
      --vcf gs://gnomad-tmp-30day/chr21_full/s1.vcf.bgz --contig chr21 \\
      --chunk chr21:5010000-16500000 --chunk chr21:13500000-27500000 \\
      --chunk chr21:24500000-37500000 --chunk chr21:34500000-46709983 \\
      --beagle-window 15 \\
      --out-prefix gs://gnomad-tmp-30day/chr21_full/phased/s1.beagle.win15
"""

import argparse
import logging

import hailtop.batch as hb

from gnomad_chets.v4.run_phasing_batch import PHASING_IMAGE, add_phasing_job

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_chunked_phasing_batch")
logger.setLevel(logging.INFO)


def add_split_job(batch, vcf, chunks, image, cpu, storage):
    """Cut the input VCF into one indexed VCF per chunk.

    Done once rather than inside each phasing job, so the full input is
    localized a single time. Record counts are echoed per chunk because the
    chunk boundaries are chosen to keep each one near a marker count already
    known to phase, and that only holds if the counts come out as intended.

    Cut with ``tabix | bgzip`` rather than ``bcftools view -Oz``. bcftools
    parses each record into BCF and re-serializes it, which at 120,279 samples
    per line means ~2e11 genotype fields per chunk: the first attempt took over
    two hours on one chunk and would have spent ~8h splitting before any
    phasing began. tabix does no parsing at all -- it seeks with the index and
    passes lines through -- so the only real work left is recompression, which
    ``bgzip -@`` threads.

    Chunks are independent, so they run concurrently within this one job rather
    than as separate jobs, which would each have to localize the whole input.

    :param batch: Batch to add to.
    :param vcf: Localized input resource group (``vcf.gz`` + ``vcf.gz.tbi``).
    :param chunks: Region strings, e.g. ``chr21:5010000-16500000``.
    :param image: Docker image carrying bcftools.
    :param cpu: Cores.
    :param storage: Disk request; must hold the input plus every chunk.
    :return: The job, with one resource group ``c{i}`` per chunk.
    """
    j = batch.new_job(name="split")
    j.image(image)
    j.cpu(cpu)
    j.memory("standard")
    j.storage(storage)
    j.spot(False)

    threads = max(1, cpu // len(chunks))
    block = ["pids=()"]
    for i, region in enumerate(chunks):
        name = f"c{i}"
        j.declare_resource_group(
            **{name: {"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"}}
        )
        out = j[name]
        block.append(
            f"( tabix -h {vcf['vcf.gz']} {region} "
            f"| bgzip -@ {threads} -c > {out['vcf.gz']} "
            f"&& tabix -p vcf {out['vcf.gz']} ) & pids+=($!)"
        )
    # wait on each pid rather than a bare `wait`, so a failed chunk fails the
    # job instead of being reported as success.
    block.append('for p in "${pids[@]}"; do wait "$p"; done')
    for i, region in enumerate(chunks):
        block.append(
            f"echo \"CHUNK {i} {region} records: "
            f"$(bcftools index -n {j[f'c{i}']['vcf.gz']})\""
        )

    j.command("set -euxo pipefail")
    j.command("\n".join(block))
    return j


def main(args):
    """Build and run the Batch."""
    backend = hb.ServiceBackend(
        billing_project=args.billing_project, remote_tmpdir=args.remote_tmpdir
    )
    batch = hb.Batch(backend=backend, name=args.batch_name or "chunked-phasing")

    vcf = batch.read_input_group(
        **{"vcf.gz": args.vcf, "vcf.gz.tbi": f"{args.vcf}.tbi"}
    )
    split = add_split_job(
        batch, vcf, args.chunk, args.image, args.split_cpu, args.split_storage
    )

    for i, region in enumerate(args.chunk):
        j = add_phasing_job(
            batch,
            split[f"c{i}"],
            tool=args.tool,
            contig=args.contig,
            name=f"{args.tool}-chunk{i}",
            cpu=args.cpu,
            memory=args.memory,
            storage=args.storage,
            image=args.image,
            shapeit_common_maf=args.shapeit_common_maf,
            shapeit_rare=args.shapeit_rare,
            beagle_window=args.beagle_window,
        )
        out = f"{args.out_prefix}.chunk{i}"
        batch.write_output(j.ofile, out)
        logger.info("chunk %d (%s) -> %s.vcf.gz", i, region, out)

    batch.run(wait=args.wait)
    backend.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input VCF (.bgz/.vcf.gz).")
    parser.add_argument("--contig", required=True, help="Contig, selects the map.")
    parser.add_argument(
        "--chunk",
        action="append",
        required=True,
        help=(
            "Region to phase as its own process, repeatable. Overlap adjacent "
            "chunks: the overlap is what lets a pair be scored from a chunk "
            "containing both its variants."
        ),
    )
    parser.add_argument("--out-prefix", required=True, help="'.chunk{i}.vcf.gz' appended.")
    parser.add_argument("--tool", default="beagle", choices=["eagle", "beagle", "shapeit"])
    parser.add_argument("--cpu", type=int, default=16, help="Cores per phasing job.")
    parser.add_argument("--beagle-window", type=float, help="Beagle window in cM.")
    parser.add_argument("--memory", default="highmem", help="Batch memory tier.")
    parser.add_argument("--storage", default="100G", help="Disk per phasing job.")
    parser.add_argument(
        "--split-cpu",
        type=int,
        default=16,
        help="Cores for the split; divided across the concurrent chunk cuts.",
    )
    parser.add_argument(
        "--split-storage",
        default="400G",
        help="Disk for the split; holds the input plus every chunk.",
    )
    parser.add_argument("--image", default=PHASING_IMAGE, help="Docker image.")
    parser.add_argument("--shapeit-common-maf", type=float, default=0.001)
    parser.add_argument("--shapeit-rare", action="store_true")
    parser.add_argument("--billing-project", default="gnomad-production")
    parser.add_argument("--remote-tmpdir", default="gs://gnomad-tmp-30day/batch")
    parser.add_argument("--batch-name", help="Batch name.")
    parser.add_argument("--wait", action="store_true")
    main(parser.parse_args())

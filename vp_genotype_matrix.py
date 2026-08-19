"""
Compute genotype count matrices for a user-supplied list of variant pairs.

Unlike the full variant co-occurrence pipelines in ``v2/create_vp_matrix.py`` and
``v4/create_vp_matrix.py`` (which *discover* co-occurring variant pairs genome-wide),
this script takes a pre-built list of specific variant pairs you already care about and
annotates each pair with genotype counts pulled from gnomAD.

Input: a Hail Table of variant pairs with ``locus1``, ``alleles1``, ``locus2``,
``alleles2`` row fields (the same schema produced by ``create_variant_pair_ht`` in
either pipeline). The table's actual key can be anything -- e.g. the CAPN3 variant
pair list this was built against is keyed by ``locus, alleles, gene, gene_id`` with
``locus1``/``alleles1`` duplicating ``locus``/``alleles``. Any additional columns
already on the table (``gene``, ``gene_id``, a pair ID, notes, etc.) are preserved
as-is.

Output: the same table with 18 new columns appended, one per two-variant genotype
class, for both raw and high-quality (``adj``) filtered genotypes::

    raw_AABB, raw_AABb, raw_AAbb, raw_AaBB, raw_AaBb, raw_Aabb, raw_aaBB, raw_aaBb, raw_aabb
    adj_AABB, adj_AABb, adj_AAbb, adj_AaBB, adj_AaBb, adj_Aabb, adj_aaBB, adj_aaBb, adj_aabb

where "A/a" refers to variant 1 (hom-ref/het/hom-var -> AA/Aa/aa) and "B/b" refers to
variant 2, e.g. ``AaBb`` = het at variant 1 AND het at variant 2 (the classic
compound-het-candidate genotype).

Genome build support
---------------------
``--genome-build grch37`` (default) pulls genotypes from the gnomAD v2 exomes/genomes
callset (matching ``gnomad_chets/v2``), which already carries a precomputed ``adj``
entry annotation.

``--genome-build grch38`` pulls genotypes from the gnomAD v4 VDS (matching
``gnomad_chets/v4``), computing ``adj`` from GQ/DP/AD via
``gnomad.utils.annotations.get_adj_expr``, exactly as ``v4/create_vp_matrix.py`` does
in its ``--filter-vds`` / ``--create-dense-filtered-mt`` steps. This path requires
access to the gnomAD v4 VDS resources (``gnomad_qc.v4``) and has not yet been run
end-to-end here -- double check output on a small test pair list before trusting it at
scale.

The core genotype-encoding/counting logic (``_encode_and_localize_genotypes``,
``create_variant_pair_genotype_ht``, ``create_variant_pair_genotype_counts_ht``, etc.)
is ported from ``v4/create_vp_matrix.py``, with the one change needed to make it
genome-build-agnostic: ``adj`` is now passed in as an argument rather than always being
computed from GQ/DP/AD, since v2 and v4 gnomAD data expose it differently.

Performance
-----------
For a small, specific pair list (a handful of variants in one or two genes), the
dominant cost is *not* the genotype counting itself -- it's naively reading the
full, genome-wide gnomAD MatrixTable/VDS just to pull out a few rows, which forces
Hail to scan every partition. Two fixes, mirrored from the partition-pruning /
shuffle-tuning approach used in the ``jg/v4-pipeline`` branch's
``create_vp_list.py``:

- ``get_covering_intervals()`` computes one small locus interval per *exact*
  variant in the input pair list (not one bounding-box interval per contig --
  a bounding box would still scan every partition between two variants that
  happen to be far apart), and both loaders apply these intervals
  (``hl.filter_intervals`` for v2, ``filter_intervals=`` for the v4 VDS
  reader) *before* the exact row-key filter, so Hail only scans partitions
  that could possibly contain one of these exact variants. Always safe to
  apply -- adjust ``--interval-padding`` (default 100bp) only if you want a
  tighter/looser margin; the exact filter downstream does the real precision
  filtering regardless.
- The real shuffles in the genotype-counting logic (pair-list repartitioning,
  ``_prepare_variant_pair_index``'s ``collect_by_key`` + ``repartition``, and
  the ``group_by("vp_ht_idx")`` aggregation) try the experimental
  ``use_new_shuffle="1")`` implementation first -- a meaningful speedup, and
  it works fine on some inputs (e.g. a 377,485-pair FBN1 run) -- but
  automatically fall back to the default shuffle implementation if that hits
  a known Hail-internal lowering bug seen on at least one larger input (a
  575,365-pair DYSF run; ``TableReaderWithExtraUID`` / "requirement failed"
  on Hail 0.2.135). See ``_with_new_shuffle_fallback``.
- Every major stage also checkpoints to a deterministic, resumable path (see
  ``_checkpoint`` / ``run_tag``) rather than a random one-off temp file, so a
  killed/restarted run against the same input doesn't redo already-finished
  work.

Example usage
-------------
.. code-block:: bash

    python vp_genotype_matrix.py \\
        --variant-pair-list-ht gs://my-bucket/my_variant_pairs.ht \\
        --genome-build grch37 \\
        --data-type exomes \\
        --output gs://my-bucket/my_variant_pairs.with_genotypes.ht \\
        --output-format ht

    # Or export a flattened TSV instead of a Hail Table:
    python vp_genotype_matrix.py \\
        --variant-pair-list-ht gs://my-bucket/my_variant_pairs.ht \\
        --output gs://my-bucket/my_variant_pairs.with_genotypes.tsv.bgz \\
        --output-format tsv
"""

import argparse
import http.client
import logging
import os
from typing import Callable, Dict, List, Optional

import hail as hl
from gnomad.utils.file_utils import file_exists

# Workaround for a Hail/py4j local-backend quirk, not an interval-count problem:
# hl.filter_intervals() round-trips through Hail's local Spark backend (an HTTP
# call to a py4j gateway on localhost) to evaluate the intervals literal. With a
# few thousand exact per-variant intervals, that response can exceed Python's
# http.client header-line-length limit (1MB), which surfaces as
# `http.client.LineTooLong: got more than 1048576 bytes when reading header
# line` -- not an error about the intervals themselves. Raising the limit here
# avoids that without changing interval granularity at all (still exactly one
# tight interval per variant, no merging).
http.client._MAXLINE = 64 * 1024 * 1024

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vp_genotype_matrix")
logger.setLevel(logging.INFO)

DEFAULT_TMP_DIR = "gs://rungar-sandbox-tmp-4day/"
"""Default temporary directory for intermediate checkpoints."""

GENOME_BUILD_REFERENCE = {"grch37": "GRCh37", "grch38": "GRCh38"}
"""Map from --genome-build CLI choice to Hail reference genome name."""

MIN_SHUFFLE_PARTITIONS = 32
"""
Absolute floor on the partition count used for this pipeline's shuffles.

`_prepare_variant_pair_index` (ported from v4/create_vp_matrix.py) originally reused
the *input* vp_ht's own partition count as the repartition target for its one real
shuffle. That's fine for the genome-wide discovery pipeline it was written for, where
vp_ht is derived from a genome-wide MatrixTable and naturally has many partitions. For
this script's use case -- a user-supplied pair list Hail Table -- that partition count
has nothing to do with how much data is actually in it (could be 1, could be
1,000, regardless of row count), which would force the union/collect_by_key/
repartition step (and everything downstream that depends on it) onto too few
partitions -- killing parallelism for a "network shuffle" stage regardless of cluster
size. Symptom: the job appears to hang for a very long time on a "wrote table with N
rows in 1 partition" / "Ordering unsorted dataset with network shuffle" step.

This floor is now a last-resort backstop inside `_prepare_variant_pair_index` itself;
the primary fix is `add_genotype_matrix` explicitly repartitioning the pair list by
actual row count (see ROWS_PER_SHUFFLE_PARTITION / MAX_SHUFFLE_PARTITIONS below)
*before* it ever reaches `_prepare_variant_pair_index`, so vp1_ht/vp2_ht -- the tables
that actually feed collect_by_key()'s shuffle -- inherit a sane partition count too,
not just the final .repartition() call's target.
"""

ROWS_PER_SHUFFLE_PARTITION = 5_000
"""
Target number of pair-instance rows per partition for this pipeline's shuffle steps.

`add_genotype_matrix` uses this to size the repartition of the pair list before
`_prepare_variant_pair_index` unions each pair into two rows (one per variant), so the
count used here is doubled relative to the raw pair count. A fixed partition-count
floor (MIN_SHUFFLE_PARTITIONS) is fine for a small pair list, but doesn't scale up for
a much larger one (e.g. a genome-wide gene like FBN1 with many more pairs) -- more
pairs should mean more partitions, not the same 32. 5,000 rows/partition is a
starting point sized for how lightweight each row is at this stage (just a locus/
alleles key + a small index/flag); tune up or down if partitions still end up too
small (task-scheduling overhead dominates) or too large (a few slow "straggler"
partitions dominate the wall-clock time) for a given run.
"""

MAX_SHUFFLE_PARTITIONS = 2_000
"""
Cap on the partition count from the ROWS_PER_SHUFFLE_PARTITION calculation above, so
an extremely large pair list doesn't get scaled into tens of thousands of tiny
partitions (each with its own non-trivial scheduling/task overhead) on a cluster that
can't usefully run that many concurrent tasks anyway.
"""

GENOTYPE_ROWS_PER_SHUFFLE_PARTITION = ROWS_PER_SHUFFLE_PARTITION
"""
Target rows/partition for the repartition in `_annotate_variant_pairs_with_genotypes`
right before group_by("vp_ht_idx").

This used to be a smaller value than ROWS_PER_SHUFFLE_PARTITION (1,000 vs. 5,000), on
the theory that since each row here also carries a joined-on `gt_info` field (a
per-sample genotype array), more/smaller partitions would balance total *bytes* per
partition better than sizing by row count alone. In practice, on a small cluster
without the (crash-prone, since removed -- see the use_new_shuffle notes below) new
shuffle implementation, that produced 5x more partitions than the other shuffle
stages in the same pipeline needed for a comparable row count, and the extra
scheduling overhead outweighed the benefit -- observed as this specific stage
freezing (0 actively running tasks, not even retries) while every other stage using
ROWS_PER_SHUFFLE_PARTITION finished quickly. Matching that same target here instead
uses a partition count already known empirically to work on this pipeline.
"""

SKEW_THRESHOLD_PAIRS = 500
"""
Partner count above which a variant is treated as a skewed "hub" key in
`_prepare_variant_pair_index` and gets its rows salted across SKEW_SALT_BUCKETS
sub-groups instead of collected into one.

Repartitioning by row count (ROWS_PER_SHUFFLE_PARTITION) balances partitions when
pairs are spread roughly evenly across variants, but it can't fix a genuinely uneven
*pairing design* -- e.g. a small set of known pathogenic variants each deliberately
cross-paired against every other candidate variant in a gene, so they end up with far
more partners than a typical variant. collect_by_key() puts every one of a variant's
partners into a single row, and the .explode() that consumes it downstream then dumps
all of them onto whichever single partition/task that one row landed on -- a straggler
that no amount of extra partitioning elsewhere helps, since the skew lives inside one
row, not across rows. 500 is set well above a typical pair list's per-variant partner
count (tens to low hundreds) so ordinary variants are never salted -- salting is *only*
worth its own overhead (an extra count-and-join pass over the full pair-instance
table) for genuinely lopsided keys.
"""

SKEW_SALT_BUCKETS = 32
"""
Number of sub-groups a skewed key (see SKEW_THRESHOLD_PAIRS) gets split across.

A hub variant with, say, 3,500 partners gets its rows randomly assigned one of 32
salt values before collect_by_key() groups by (locus, alleles, salt) instead of just
(locus, alleles) -- splitting one ~3,500-row group into ~32 groups of ~110 rows each,
which land across many partitions instead of one. Every downstream step re-collapses
by vp_ht_idx (the pair index) regardless of which salt bucket a row passed through, so
this only changes how the work is spread across partitions, not the result.
"""

GENOTYPE_CLASSES = [
    "AABB",
    "AABb",
    "AAbb",
    "AaBB",
    "AaBb",
    "Aabb",
    "aaBB",
    "aaBb",
    "aabb",
]
"""
Names for the 9 two-variant genotype classes, in the same order as the count arrays
produced by `_convert_gt_info_to_counts` (index = v1_genotype * 3 + v2_genotype, where
0/1/2 = hom-ref/het/hom-var). "A/a" = variant 1, "B/b" = variant 2.
"""


########################################################################################
### Checkpoint/resume helpers
########################################################################################
def _run_tag_from_path(path: str) -> str:
    """
    Derive a filesystem-safe checkpoint tag from an input Table path.

    e.g. "gs://bucket/dir/FBN1_variant_pairs.ht" -> "FBN1_variant_pairs". Used so
    reruns against the same input pair list land on the same checkpoint paths.

    :param path: Path to a Hail Table (or similar).
    :return: Basename with a trailing .ht/.mt stripped, if present.
    """
    name = path.rstrip("/").split("/")[-1]
    for suffix in (".ht", ".mt"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def _with_new_shuffle_fallback(compute_fn: Callable[[], hl.Table]) -> hl.Table:
    """
    Run a Hail computation with the experimental new-shuffle implementation enabled,
    falling back to the default (slower but stable) implementation only if that
    specific implementation crashes.

    `use_new_shuffle="1"` is a meaningful speedup for the shuffle-heavy operations in
    this pipeline (collect_by_key, repartition, group_by), and works fine at some
    scales (e.g. it completed a 377,485-pair FBN1 run without issue) -- but at least
    one larger/differently-shaped input (a 575,365-pair DYSF run) hit a Hail-internal
    lowering bug in it (`TableReaderWithExtraUID` / "requirement failed" in
    LowerAndExecuteShuffles/repartitionNoShuffle, on Hail 0.2.135), reproduced
    regardless of partition count or the shuffled table's key shape. Since it's
    apparently input-dependent rather than always broken, disabling it everywhere
    traded a lot of speed for safety on the (more common) inputs where it works fine.
    This tries the fast path first and only pays for the slow path on the specific
    inputs that need it.

    :param compute_fn: Zero-argument callable that builds *and forces execution of*
        (e.g. via `.cache()`/`.checkpoint()`) a Hail Table. Must be safe to call
        twice -- every call site here operates on an already-checkpointed or
        otherwise cheap-to-rebuild Table, so a retry only redoes the one shuffle in
        question, not any expensive upstream work.
    :return: The materialized Table.
    """
    hl._set_flags(use_new_shuffle="1")
    try:
        result = compute_fn()
    except hl.utils.java.FatalError as e:
        if "TableReaderWithExtraUID" not in str(e):
            raise
        logger.warning(
            "use_new_shuffle hit the known Hail-internal lowering bug "
            "(TableReaderWithExtraUID) on this operation; retrying with the "
            "default shuffle implementation instead."
        )
        hl._set_flags(use_new_shuffle=None)
        result = compute_fn()
    finally:
        hl._set_flags(use_new_shuffle=None)
    return result


def _checkpoint(
    ht: hl.Table,
    tmp_dir: str,
    run_tag: Optional[str],
    stage: str,
    resume: bool = True,
    try_new_shuffle: bool = False,
) -> hl.Table:
    """
    Checkpoint `ht`, reading back an existing checkpoint instead of recomputing it
    if one is already there.

    This is the resumability mechanism for this pipeline: every major stage
    (interval-filtered variants, encoded genotypes, the pair index, the salted/
    collect_by_key union) writes to a *deterministic* path (derived from `run_tag`)
    instead of the random, one-off temp files `hl.utils.new_temp_file` produces.
    Re-running `add_genotype_matrix` for the same input (same `run_tag`) after a kill
    -- e.g. because a later stage hit a slow/stuck shuffle -- reads each already-
    completed stage straight back from disk instead of redoing interval pruning,
    dense MatrixTable filtering, repartitioning, and skew salting all over again, and
    picks up compute right where it left off.

    :param ht: Table to checkpoint.
    :param tmp_dir: Base temporary directory (same one passed to add_genotype_matrix).
    :param run_tag: Stable identifier for this input/run (see `_run_tag_from_path`),
        or None to fall back to the old behavior -- a random, non-resumable temp file
        (used when these functions are called directly without a run_tag, e.g. from
        a notebook or test).
    :param stage: Short name for this checkpoint stage (e.g. "variants", "vp_union")
        -- combined with `run_tag` to form the checkpoint path. Must be unique per
        call site.
    :param resume: If False, always recompute and overwrite even if a checkpoint
        already exists at this path (still writes one, so a later call *with*
        resume=True can pick it up). Use this to force one stage to be redone after
        changing something upstream of it (e.g. --interval-padding) without having
        to manually delete checkpoint files.
    :param try_new_shuffle: Whether the computation feeding into this checkpoint
        includes a real shuffle worth trying under `use_new_shuffle` first (see
        `_with_new_shuffle_fallback`). Only set this at checkpoints that directly
        follow a repartition/collect_by_key/group_by -- it's the checkpoint's own
        `.checkpoint()`/`.cache()` call that actually forces execution of whatever
        lazy shuffle came before it, which is where a crash (if any) would surface.
    :return: The checkpointed (or, if reused, re-read) Table.
    """
    if run_tag is None:
        write = lambda: ht.checkpoint(hl.utils.new_temp_file("vp_genotype_matrix", "ht"))
        return _with_new_shuffle_fallback(write) if try_new_shuffle else write()

    path = os.path.join(
        tmp_dir, "vp_genotype_matrix_checkpoints", f"{run_tag}.{stage}.ht"
    )
    if resume and file_exists(path):
        logger.info("Reusing existing checkpoint for stage '%s': %s", stage, path)
        return hl.read_table(path)

    logger.info("Checkpointing stage '%s' to %s", stage, path)
    write = lambda: ht.checkpoint(path, overwrite=True)
    return _with_new_shuffle_fallback(write) if try_new_shuffle else write()


########################################################################################
### Build-agnostic genotype counting logic (ported from v4/create_vp_matrix.py)
########################################################################################
def _encode_and_localize_genotypes(
    mt: hl.MatrixTable,
    adj_expr: hl.expr.BooleanExpression,
    tmp_dir: str = DEFAULT_TMP_DIR,
    run_tag: Optional[str] = None,
    resume: bool = True,
) -> hl.Table:
    """
    Encode genotypes for efficiency and localize to a Table.

    For most rare variants, most samples are hom_ref, so code them as missing to save
    space.

    Encodes genotypes as:

        - missing = hom_ref (space saving)
        - 0 = missing data (actual missing call)
        - 1 = het
        - 2 = hom_var

    For adj genotypes:

        - missing = hom_ref adj (space saving)
        - 0 = missing data or not adj (actual missing call)
        - 1 = het adj
        - 2 = hom_var adj

    Filters to keep only genotypes where the variant is called (reduces array size).

    :param mt: MatrixTable with variant data. Row key must be (locus, alleles). Must
        have a 'GT' entry field.
    :param adj_expr: Boolean expression indicating whether each entry passes the
        high-quality ("adj") genotype filter. Passed in explicitly (rather than always
        computed from GQ/DP/AD) so this function works whether 'adj' is already
        precomputed (e.g. gnomAD v2) or needs to be derived (e.g. gnomAD v4).
    :param tmp_dir: Base temporary directory, used for checkpointing (see
        `_checkpoint`).
    :param run_tag: Stable identifier for this input/run, for resumable
        checkpointing (see `_checkpoint`); None disables resuming.
    :param resume: See `_checkpoint`.
    :return: Table with localized genotype info, filtered to called variants only.
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )

    # For adj genotypes, set to 0 (missing data) if the variant doesn't pass the adj
    # filter.
    adj_gt_count_expr = hl.if_else(
        adj_expr, gt_count_expr, 0, missing_false=True
    )

    mt = mt.select_entries(gt_info=(gt_count_expr, adj_gt_count_expr))
    ht = mt.localize_entries("gt_info", "samples")

    # Store sample information: (sample_id, raw_gt_count, adj_gt_count).
    # Filter to keep only genotypes where the variant is called (reduces array size).
    ht = ht.select(
        gt_info=hl.enumerate(ht.gt_info)
        .map(lambda x: (x[0], x[1].gt_info[0], x[1].gt_info[1]))
        .filter(lambda x: hl.is_defined(x[1]) | hl.is_defined(x[2]))
    )
    ht = _checkpoint(ht, tmp_dir, run_tag, "encoded_genotypes", resume=resume)

    return ht


def _prepare_variant_pair_index(
    vp_ht: hl.Table,
    tmp_dir: str = DEFAULT_TMP_DIR,
    run_tag: Optional[str] = None,
    resume: bool = True,
) -> hl.Table:
    """
    Prepare variant pair Table for genotype annotation.

    Creates separate entries for each variant in the pair (v1 and v2), then unions
    them. This helps with performance issues when annotating genotype info for both
    variants.

    Any variant with more than SKEW_THRESHOLD_PAIRS partners has its rows salted
    across SKEW_SALT_BUCKETS sub-groups before collect_by_key() (see that constant's
    docstring) -- purely an execution-strategy detail to avoid a data-skew straggler
    task; doesn't change which rows come out the other end.

    :param vp_ht: Variant pair Table with fields locus1, alleles1, locus2, alleles2.
    :param tmp_dir: Base temporary directory, used for checkpointing (see
        `_checkpoint`).
    :param run_tag: Stable identifier for this input/run, for resumable
        checkpointing (see `_checkpoint`); None disables resuming.
    :param resume: See `_checkpoint`.
    :return: Unioned Table with index field and variant pair indicator (vp=1 or vp=2).
    """
    # Backstop only -- add_genotype_matrix already repartitions vp_ht by actual row
    # count before calling this function (see ROWS_PER_SHUFFLE_PARTITION docstring),
    # so this is a no-op in the normal path. Kept here in case this function is ever
    # called directly with an un-repartitioned Table.
    n_partitions = max(vp_ht.n_partitions(), MIN_SHUFFLE_PARTITIONS)

    vp1_ht = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1)
    vp1_ht = vp1_ht.select("vp_ht_idx", vp=1)

    vp2_ht = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2)
    vp2_ht = vp2_ht.select("vp_ht_idx", vp=2)

    vp_all_ht = vp1_ht.union(vp2_ht)

    # Skew mitigation (see SKEW_THRESHOLD_PAIRS / SKEW_SALT_BUCKETS docstrings): find
    # variants with an unusually large number of partners and split *only* those
    # across multiple sub-groups before collect_by_key(), so their eventual
    # .explode() downstream spreads across many partitions instead of dumping
    # thousands of rows onto one. This costs one extra count-and-join pass over
    # vp_all_ht (itself a shuffle, but over lightweight rows with no large embedded
    # arrays -- cheap relative to what it's preventing), and is a no-op (constant
    # salt=0, same single group as before) for the vast majority of ordinary,
    # non-skewed variants.
    key_counts_ht = _checkpoint(
        vp_all_ht.group_by("locus", "alleles").aggregate(n=hl.agg.count()),
        tmp_dir, run_tag, "key_counts", resume=resume,
    )

    n_skewed = key_counts_ht.aggregate(
        hl.agg.count_where(key_counts_ht.n > SKEW_THRESHOLD_PAIRS)
    )
    if n_skewed:
        logger.info(
            "Found %d variant(s) with more than %d partners; salting their rows "
            "across %d sub-groups to avoid a data-skew straggler task.",
            n_skewed, SKEW_THRESHOLD_PAIRS, SKEW_SALT_BUCKETS,
        )

    vp_all_ht = vp_all_ht.annotate(
        _salt=hl.if_else(
            key_counts_ht[vp_all_ht.locus, vp_all_ht.alleles].n > SKEW_THRESHOLD_PAIRS,
            hl.rand_int32(SKEW_SALT_BUCKETS),
            0,
        )
    )
    vp_all_ht = vp_all_ht.key_by("locus", "alleles", "_salt")

    # _salt only needs to exist as part of the key long enough for collect_by_key()
    # above to actually split the hot groups -- nothing downstream reads it. Re-keying
    # to just (locus, alleles) -- a strict prefix of the current sort order, so this is
    # a free, shuffle-free operation -- and dropping _salt puts every operation after
    # this point back on the exact key shape the rest of this pipeline was built for.
    vp_all_ht = vp_all_ht.select("vp_ht_idx", "vp")

    # This collect_by_key + repartition is a real shuffle, and the checkpoint below
    # tries the faster use_new_shuffle="1" implementation first -- it's a meaningful
    # speedup here and works fine at some scales (e.g. a 377,485-pair FBN1 run) -- but
    # falls back to the default implementation if that hits a known Hail-internal
    # lowering bug seen on at least one larger input (a 575,365-pair DYSF run;
    # `TableReaderWithExtraUID` / "requirement failed", reproduced regardless of
    # whether _salt is still part of the key, so it's not about key shape). See
    # `_with_new_shuffle_fallback`.
    vp_union_ht = vp_all_ht.collect_by_key().key_by("locus", "alleles").drop("_salt")
    vp_union_ht = vp_union_ht.repartition(n_partitions, shuffle=True)
    vp_union_ht = _checkpoint(
        vp_union_ht, tmp_dir, run_tag, "vp_union", resume=resume, try_new_shuffle=True
    )

    return vp_union_ht


def _annotate_variant_pairs_with_genotypes(
    vp_union_ht: hl.Table,
    ht: hl.Table,
    n_pairs: Optional[int] = None,
) -> hl.Table:
    """
    Annotate variant pairs with genotype information and group by variant pair index.

    :param vp_union_ht: Unioned variant pair Table with index field and variant pair
        indicator.
    :param ht: Localized entries Table with genotype info.
    :param n_pairs: Row count of the original (pre-union) pair list, if already known
        by the caller -- lets the repartition before group_by("vp_ht_idx") below be
        sized without an extra `.count()` action, since exploding "values" always
        restores exactly 2 * n_pairs rows regardless of how collect_by_key grouped or
        salted them upstream. If omitted, counts vp_union_ht directly instead (a real
        action, since by this point it carries the joined gt_info payload).
    :return: Variant pair Table with genotype info grouped by variant pair index.
    """
    vp_union_ht = vp_union_ht.annotate(
        gt_info=ht[vp_union_ht.locus, vp_union_ht.alleles].gt_info
    )
    vp_union_ht = vp_union_ht.explode("values")
    vp_union_ht = vp_union_ht.transmute(**vp_union_ht.values)

    vp_union_ht = vp_union_ht.annotate(
        gt_info=ht[vp_union_ht.locus, vp_union_ht.alleles].gt_info
    )

    # Same right-sizing fix as add_genotype_matrix's pair_key_ht repartition and
    # _prepare_variant_pair_index's collect_by_key repartition: the
    # group_by("vp_ht_idx").aggregate() below is another real shuffle, and letting it
    # inherit whatever partition count fell out of collect_by_key/explode above --
    # rather than sizing it to the actual row count -- risks the same kind of
    # straggler-task imbalance seen upstream (observed in practice: a handful of
    # stuck tasks on this exact step even after the collect_by_key skew fix).
    n_rows = 2 * n_pairs if n_pairs is not None else vp_union_ht.count()
    target_partitions = min(
        MAX_SHUFFLE_PARTITIONS,
        max(MIN_SHUFFLE_PARTITIONS, -(-n_rows // GENOTYPE_ROWS_PER_SHUFFLE_PARTITION)),
    )

    # Both operations below (the explicit repartition() and
    # group_by("vp_ht_idx").aggregate(), since grouping by key is itself a shuffle)
    # are lazy and don't actually execute until the "genotype_ht" checkpoint back in
    # add_genotype_matrix forces them -- that's where use_new_shuffle-with-fallback is
    # applied (see _with_new_shuffle_fallback), covering both shuffles at once.
    vp_union_ht = vp_union_ht.repartition(target_partitions, shuffle=True)

    vp_union_ht = vp_union_ht.group_by("vp_ht_idx").aggregate(
        gt_info=hl.agg.collect((vp_union_ht.vp, vp_union_ht.gt_info))
    )

    vp_union_ht = vp_union_ht.annotate_globals(samples=ht.index_globals().samples)

    return vp_union_ht


def create_variant_pair_genotype_ht(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    adj_expr: hl.expr.BooleanExpression,
    n_pairs: Optional[int] = None,
    tmp_dir: str = DEFAULT_TMP_DIR,
    run_tag: Optional[str] = None,
    resume: bool = True,
) -> hl.Table:
    """
    Create a variant pair genotype Table from a MatrixTable and variant pair list.

    :param mt: MatrixTable with variant data. Row key must be (locus, alleles).
    :param vp_ht: Table of variant pairs with fields locus1, alleles1, locus2,
        alleles2. Should contain *only* these 4 key fields (no extra columns) --
        callers should strip other columns off before calling this, since this
        function adds its own 'vp_ht_idx' field.
    :param adj_expr: Boolean expression on `mt` indicating high-quality genotypes.
    :param n_pairs: Row count of `vp_ht`, if the caller already knows it (e.g.
        add_genotype_matrix computes this anyway to size an earlier repartition).
        Passed through to `_annotate_variant_pairs_with_genotypes` so it can size its
        own repartition without an extra `.count()` action; if omitted, that function
        falls back to counting itself.
    :param tmp_dir: Base temporary directory, used for checkpointing (see
        `_checkpoint`).
    :param run_tag: Stable identifier for this input/run, for resumable
        checkpointing (see `_checkpoint`); None disables resuming.
    :param resume: See `_checkpoint`.
    :return: Variant pair Table with genotype info for both variants in each pair.
    """
    ht = _encode_and_localize_genotypes(
        mt, adj_expr, tmp_dir=tmp_dir, run_tag=run_tag, resume=resume
    )

    vp_ht = vp_ht.key_by("locus1", "alleles1", "locus2", "alleles2")
    vp_ht = vp_ht.add_index("vp_ht_idx").key_by("vp_ht_idx")
    vp_ht = _checkpoint(vp_ht, tmp_dir, run_tag, "pair_index", resume=resume)

    vp_union_ht = _prepare_variant_pair_index(
        vp_ht, tmp_dir=tmp_dir, run_tag=run_tag, resume=resume
    )
    vp_union_ht = _annotate_variant_pairs_with_genotypes(vp_union_ht, ht, n_pairs=n_pairs)

    vp_ht = vp_union_ht.annotate(**vp_ht[vp_union_ht.vp_ht_idx])

    return vp_ht


def _convert_gt_info_to_counts(
    gt_counts_dict: hl.expr.DictExpression,
    n_samples_filtered_out_expr: hl.expr.Int32Expression,
) -> hl.expr.ArrayExpression:
    """
    Convert genotype count dictionary to a 9-element genotype count array.

    Array index = v1_genotype * 3 + v2_genotype (0=hom-ref, 1=het, 2=hom-var), matching
    the order of GENOTYPE_CLASSES: [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb,
    aabb].

    :param gt_counts_dict: Dictionary from counter aggregation with keys as [v1_gt,
        v2_gt] and values as counts.
    :param n_samples_filtered_out_expr: Number of samples filtered out (missing both
        variants) to add to the hom-ref/hom-ref count (index 0).
    :return: Array of 9 integers representing counts for each genotype combination.
    """
    indices = hl.range(0, 9)
    dict_keys = hl.array(
        [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    )

    return hl.zip(indices, dict_keys).map(
        lambda x: hl.if_else(
            x[0] == 0,
            gt_counts_dict.get(x[1], 0) + n_samples_filtered_out_expr,
            gt_counts_dict.get(x[1], 0),
        )
    )


def _calculate_genotype_counts(
    per_sample_gt_expr: hl.expr.ArrayExpression,
    n_samples_filtered_out_expr: hl.expr.Int32Expression,
    gt_field: str,
) -> hl.expr.ArrayExpression:
    """
    Calculate genotype counts for variant pairs from per-sample genotype expressions.

    :param per_sample_gt_expr: Array of per-sample genotype info tuples (sample_id,
        v1_gt_info, v2_gt_info).
    :param n_samples_filtered_out_expr: Number of samples filtered out (missing both
        variants).
    :param gt_field: Field name to extract from gt_info ("raw_gt" or "adj_gt").
    :return: Array of genotype counts [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb,
        aabb].
    """
    gt_expr = per_sample_gt_expr.map(lambda x: [x.get(1)[gt_field], x.get(2)[gt_field]])

    gt_expr = gt_expr.filter(
        lambda x: (
            (hl.is_missing(x[0]) | (x[0] != 0)) & (hl.is_missing(x[1]) | (x[1] != 0))
        )
    )
    gt_expr = gt_expr.map(lambda x: x.map(lambda y: hl.or_else(y, 0)))

    # No or_missing guard on an empty gt_expr: v2 falls back to a zero array here
    # (`hl.or_else(hl.agg.filter(adj1 & adj2, ...), [0] * 9)` in
    # v2/create_vp_matrix.py's create_vp_summary), so returning NA instead would make
    # these pairs look absent rather than uncounted. Counting an empty array is
    # well-defined -- hl.agg.counter over zero elements is an empty dict, leaving
    # [n_samples_filtered_out, 0, ...]. That still agrees with v2: a sample is only
    # in n_samples_filtered_out if it is an adj-passing hom-ref at *both* variants,
    # which is exactly what v2's aggregation would have put in the AABB cell, so
    # whenever it is non-zero v2's fallback would not have fired either.
    gt_expr = _convert_gt_info_to_counts(
        gt_expr.aggregate(hl.agg.counter),
        n_samples_filtered_out_expr,
    )

    return gt_expr


def create_variant_pair_genotype_counts_ht(ht: hl.Table) -> hl.Table:
    """
    Create a variant pair genotype counts Table from a variant pair genotype Table.

    :param ht: Variant pair genotype Table with gt_info field containing per-sample
        genotype information for both variants.
    :return: Variant pair Table keyed by locus1/alleles1/locus2/alleles2 with
        gt_counts_raw and gt_counts_adj fields (each a 9-element array in
        GENOTYPE_CLASSES order).
    """
    ht = ht.annotate(
        gt_info=ht.gt_info.flatmap(
            lambda x: x[1].map(
                lambda y: hl.struct(s=y[0], vp=x[0], raw_gt=y[1], adj_gt=y[2])
            )
        )
    )

    n_samples = ht.samples.length()
    n_samples_with_data_expr = hl.set(ht.gt_info.map(lambda x: x.s)).length()
    n_samples_filtered_out_expr = n_samples - n_samples_with_data_expr

    ht = ht.annotate(
        gt_info=(
            ht.gt_info.group_by(lambda x: x.s)
            .values()
            .map(lambda x: hl.dict(x.map(lambda y: (y.vp, y))))
        ),
        n_samples_filtered_out_expr=n_samples_filtered_out_expr,
    )

    ht = ht.select(
        "locus1",
        "alleles1",
        "locus2",
        "alleles2",
        **{
            f"gt_counts_{n}": _calculate_genotype_counts(
                ht.gt_info, ht.n_samples_filtered_out_expr, gt_field=f"{n}_gt"
            )
            for n in ["raw", "adj"]
        },
    ).cache()

    return ht.key_by("locus1", "alleles1", "locus2", "alleles2")


########################################################################################
### Partition pruning
########################################################################################
def get_covering_intervals(
    variants_ht: hl.Table, padding: int = 100
) -> List[hl.utils.Interval]:
    """
    Compute one small locus interval per *exact* variant in a Table.

    A handful of variants filtered directly out of a genome-wide gnomAD
    MatrixTable/VDS (``mt.filter_rows(hl.is_defined(variants_ht[mt.row_key]))``
    alone) still forces Hail to scan every partition, since row-key-set filters
    aren't pushed down to skip partitions. Restricting to covering intervals
    *first* (``hl.filter_intervals`` / ``filter_intervals=``) lets Hail prune
    partitions before the exact-match filter runs -- the same partition-pruning
    idea as the ``--chr`` / ``--test-chrom`` sharding in the ``jg/v4-pipeline``
    branch's ``create_vp_list.py``, just derived automatically from the input
    pair list instead of a CLI flag.

    Deliberately one tight interval *per variant* rather than one bounding-box
    interval per contig: a single covering interval (min position to max
    position) reads every partition in between even when the variants
    themselves are sparse across that span (e.g. two variants far apart on the
    same chromosome, or a gene with large intronic gaps between the exons
    actually carrying pairs). Per-variant intervals only ever touch partitions
    that could contain one of these exact variants -- at least as fast as a
    bounding box, and often much faster when variants aren't tightly clustered.
    Only becomes a concern for pair lists spanning many thousands of variants,
    where the interval list itself gets large; not a factor at the scale this
    script is meant for (a handful of genes at a time).

    :param variants_ht: Table keyed by (locus, alleles) (e.g. from
        `get_variants_ht`).
    :param padding: Small bp padding added on each side of each variant's exact
        position (covers multi-bp indels / off-by-one edge cases -- the
        downstream exact row-key filter still does the real precision
        filtering, so this only needs to be "safely non-zero", not tight).
        Default 100bp.
    :return: List of Hail locus intervals, one per unique variant locus in
        `variants_ht`.
    """
    reference_genome = variants_ht.locus.dtype.reference_genome.name
    loci = variants_ht.aggregate(hl.agg.collect_as_set(variants_ht.locus))

    return [
        hl.parse_locus_interval(
            f"{locus.contig}:{max(1, locus.position - padding)}-{locus.position + padding}",
            reference_genome=reference_genome,
        )
        for locus in loci
    ]


########################################################################################
### Build-specific genotype data loaders
###
### Each loader takes (data_type, variants_ht, intervals) and must return a dense
### hl.MatrixTable, row-keyed by (locus, alleles) and filtered down to (at least) the
### rows in `variants_ht`, with a 'GT' entry field and an 'adj' boolean entry field.
########################################################################################
def _get_dense_mt_grch37(
    data_type: str, variants_ht: hl.Table, intervals: Optional[List[hl.utils.Interval]]
) -> hl.MatrixTable:
    """
    Load gnomAD v2 (GRCh37) genotype data, filtered to the given variants.

    Relies on the gnomAD v2 release MatrixTable already carrying a precomputed 'adj'
    entry annotation (used as-is elsewhere in gnomad_chets/v2, e.g.
    ``resources.get_adj_missing_mt``).

    :param data_type: One of 'exomes' or 'genomes'.
    :param variants_ht: Table of variants (keyed by locus, alleles) to filter to.
    :param intervals: Covering intervals (see `get_covering_intervals`) applied before
        the exact row-key filter, to prune partitions. Skipped if None/empty.
    :return: Dense MatrixTable with 'GT' and 'adj' entry fields, release samples only.
    """
    from gnomad_qc.v2.resources import get_gnomad_data, get_gnomad_meta

    logger.info("Loading gnomAD v2 %s MatrixTable...", data_type)
    mt = get_gnomad_data(data_type)
    if intervals:
        logger.info("Pruning to %d covering interval(s) before exact filter...", len(intervals))
        mt = hl.filter_intervals(mt, intervals)
    # Release samples only, matching the sample set the published v2 counts were
    # computed over (`create_vp_summary`'s caller in v2/create_vp_matrix.py filters
    # to `meta.release` before counting). `high_quality` is the *pair-discovery*
    # sample set -- it's ~28k samples larger (153,927 vs 125,748 in v2 exomes) since
    # it still includes related and non-releasable individuals, so counting over it
    # inflates every genotype cell and won't reconcile against the published table.
    # The grch38 loader below is release-only for the same reason.
    meta = get_gnomad_meta(data_type)
    mt = mt.filter_cols(meta[mt.col_key].release)
    mt = mt.filter_rows(hl.is_defined(variants_ht[mt.row_key]))

    return mt


def _get_dense_mt_grch38(
    data_type: str, variants_ht: hl.Table, intervals: Optional[List[hl.utils.Interval]]
) -> hl.MatrixTable:
    """
    Load gnomAD v4 (GRCh38) genotype data, filtered to the given variants.

    Mirrors the ``--filter-vmt`` step in the ``jg/v4-pipeline`` branch's
    ``create_vp_list.py`` (and the older ``--filter-vds`` in
    ``v4/create_vp_matrix.py``): pulls only entries needed to compute 'adj'
    (GQ/DP/AD), restricts to covering intervals up front so partition scanning is
    bounded (same idea as that branch's ``--chr``/``--test-chrom`` sharding, just
    computed from the pair list instead of a CLI flag), filters the VDS down to the
    requested variants, and densifies.

    NOTE: this path requires gnomad_qc v4 / VDS access and has not yet been validated
    end-to-end -- sanity check output on a small test pair list first.

    :param data_type: One of 'exomes' or 'genomes'.
    :param variants_ht: Table of variants (keyed by locus, alleles) to filter to.
    :param intervals: Covering intervals (see `get_covering_intervals`) passed as
        `filter_intervals` to bound partition scanning. Skipped if None/empty.
    :return: Dense MatrixTable with 'GT' and 'adj' entry fields, release samples only.
    """
    from gnomad.utils.annotations import get_adj_expr
    from gnomad_qc.v4.resources.basics import (
        get_gnomad_v4_genomes_vds,
        get_gnomad_v4_vds,
    )

    logger.info("Loading gnomAD v4 %s VDS...", data_type)
    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )
    if intervals:
        logger.info("Pruning to %d covering interval(s) before exact filter...", len(intervals))
    vds = get_vds_func(
        release_only=True,
        split=True,
        filter_intervals=intervals or None,
        filter_variant_ht=variants_ht,
        entries_to_keep=["GT", "GQ", "DP", "AD"],
        split_reference_blocks=False,
    )
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD))

    return mt


BUILD_LOADERS: Dict[
    str, Callable[[str, hl.Table, Optional[List[hl.utils.Interval]]], hl.MatrixTable]
] = {
    "grch37": _get_dense_mt_grch37,
    "grch38": _get_dense_mt_grch38,
}


########################################################################################
### Orchestration
########################################################################################
def get_variants_ht(vp_ht: hl.Table) -> hl.Table:
    """
    Get the set of unique variants referenced by a variant pair Table.

    :param vp_ht: Variant pair Table with fields locus1, alleles1, locus2, alleles2.
    :return: Distinct Table of variants, keyed by (locus, alleles).
    """
    v1_ht = vp_ht.key_by(locus=vp_ht.locus1, alleles=vp_ht.alleles1).select().distinct()
    v2_ht = vp_ht.key_by(locus=vp_ht.locus2, alleles=vp_ht.alleles2).select().distinct()

    return v1_ht.union(v2_ht).distinct()


def _is_canonical_pair_order(
    locus1: hl.expr.LocusExpression,
    alleles1: hl.expr.ArrayExpression,
    locus2: hl.expr.LocusExpression,
    alleles2: hl.expr.ArrayExpression,
) -> hl.expr.BooleanExpression:
    """
    Check whether a pair is ordered the way the published gnomAD v2 tables order pairs.

    Mirrors ``_get_ordered_vp_struct`` in v2/create_vp_matrix.py: sort on locus
    position, tie-broken on the alt allele. Note that this compares position only, not
    contig -- v2 only ever formed pairs within a single gene, so contig was constant.

    Used purely to warn (see `add_genotype_matrix`); nothing here reorders pairs.

    :param locus1: Locus of the first variant in the pair.
    :param alleles1: Alleles of the first variant in the pair.
    :param locus2: Locus of the second variant in the pair.
    :param alleles2: Alleles of the second variant in the pair.
    :return: Boolean expression, True if the pair is in v2 canonical order.
    """
    return hl.if_else(
        locus1.position == locus2.position,
        alleles1[1] <= alleles2[1],
        locus1.position < locus2.position,
    )


def add_genotype_matrix(
    vp_ht: hl.Table,
    genome_build: str,
    data_type: str,
    tmp_dir: str = DEFAULT_TMP_DIR,
    interval_padding: int = 100,
    run_tag: Optional[str] = None,
    resume: bool = True,
) -> hl.Table:
    """
    Annotate a variant pair Table with per-pair genotype counts.

    :param vp_ht: Input variant pair Table. Must contain locus1/alleles1/locus2/
        alleles2 fields (as row fields -- the table's actual key can be anything,
        e.g. locus/alleles/gene/gene_id). Any additional columns are preserved.
    :param genome_build: One of 'grch37' or 'grch38'.
    :param data_type: One of 'exomes' or 'genomes'.
    :param tmp_dir: Temporary directory for intermediate checkpoints.
    :param interval_padding: Bp padding per variant for the partition-pruning step
        (see `get_covering_intervals`) -- one small interval per exact variant, not
        one bounding box per contig. Set to a negative value to disable interval
        pruning entirely and always fall back to the exact row-key filter alone.
    :param run_tag: Stable identifier for this input/run (see `_run_tag_from_path`)
        -- every major stage checkpoints to a deterministic path derived from this
        tag, so re-running with the same tag after a kill resumes from the last
        completed stage instead of starting over. None (the default for direct
        calls) disables resuming and falls back to random, one-off temp files, as
        before. `main()` derives this automatically from --variant-pair-list-ht.
    :param resume: If False, ignore any existing checkpoints and recompute/overwrite
        every stage (still writes fresh checkpoints, so a later run with resume=True
        can pick them up). Use to force a clean re-run after changing something
        upstream (e.g. --interval-padding) without deleting checkpoint files by hand.
    :return: `vp_ht` with 18 additional columns: raw_<GENOTYPE_CLASS> and
        adj_<GENOTYPE_CLASS> for each of the 9 genotype classes.
    """
    if genome_build not in BUILD_LOADERS:
        raise ValueError(
            f"Unknown genome_build '{genome_build}'. Must be one of "
            f"{list(BUILD_LOADERS)}."
        )

    orig_key = list(vp_ht.key)
    if not {"locus1", "alleles1", "locus2", "alleles2"}.issubset(vp_ht.row):
        raise ValueError(
            "vp_ht must contain locus1, alleles1, locus2, and alleles2 fields."
        )

    variants_ht = _checkpoint(
        get_variants_ht(vp_ht), tmp_dir, run_tag, "variants", resume=resume
    )
    n_variants = variants_ht.count()
    logger.info(
        "Filtering gnomAD %s %s data to %s variants...", genome_build, data_type, n_variants
    )

    intervals = get_covering_intervals(variants_ht, padding=interval_padding) if interval_padding >= 0 else None
    if intervals:
        logger.info(
            "Computed %d covering interval(s) (padding=%dbp) to prune partitions "
            "before the exact variant filter: %s",
            len(intervals), interval_padding, intervals,
        )

    mt = BUILD_LOADERS[genome_build](data_type, variants_ht, intervals)

    # Only keep the pair key for internal processing so the helper functions above
    # (which add their own fields, e.g. 'vp_ht_idx') don't clash with any extra
    # columns the user's input table already has. Note: vp_ht's *actual* key may be
    # something else entirely (e.g. locus/alleles/gene/gene_id) -- re-key by the pair
    # fields first, then select() with no args to drop everything else, matching the
    # pattern used in v4/create_vp_matrix.py's create_dense_filtered_mt().
    pair_key_ht = vp_ht.key_by("locus1", "alleles1", "locus2", "alleles2").select()

    # Right-size partitions here, once, based on actual pair count -- rather than
    # letting _prepare_variant_pair_index's MIN_SHUFFLE_PARTITIONS floor be the only
    # thing standing between this pipeline and a single-partition shuffle. This
    # matters because that floor only controls the *target* of its final
    # .repartition() call; vp1_ht/vp2_ht (the tables that actually feed
    # collect_by_key()'s own shuffle) are built directly from pair_key_ht via
    # key_by()/select() with no repartition of their own, so they inherit whatever
    # partition count pair_key_ht has *here*. Doing it once up front, sized to the
    # real row count (doubled, since _prepare_variant_pair_index unions each pair
    # into 2 rows), fixes that for every downstream shuffle in this function, not
    # just the last one. See ROWS_PER_SHUFFLE_PARTITION / MAX_SHUFFLE_PARTITIONS.
    #
    # The canonical-order check rides along in that same count action rather than
    # costing a pass of its own. It only warns: counts are computed for whatever
    # orientation the caller supplied and are correct for it, and reordering here
    # would desync the pair fields from any extra columns the caller carried
    # alongside them (e.g. a table keyed by locus/alleles/gene where locus1
    # duplicates locus).
    n_pairs, n_noncanonical = pair_key_ht.aggregate(
        (
            hl.agg.count(),
            hl.agg.count_where(
                ~_is_canonical_pair_order(
                    pair_key_ht.locus1,
                    pair_key_ht.alleles1,
                    pair_key_ht.locus2,
                    pair_key_ht.alleles2,
                )
            ),
        )
    )
    if n_noncanonical:
        logger.warning(
            "%d of %d input pairs are not in gnomAD v2 canonical order (sorted on "
            "locus position, tie-broken on alt allele). Their counts are correct for "
            "the orientation given, but will not line up with the published v2 "
            "co-occurrence table: the pair key will not join, and matching pairs "
            "unordered leaves the genotype cells transposed (AABb<->AaBB, "
            "AAbb<->aaBB, Aabb<->aaBb; AABB, AaBb and aabb are unaffected). Swap "
            "locus1/alleles1 with locus2/alleles2 on those rows before comparing.",
            n_noncanonical, n_pairs,
        )

    target_partitions = min(
        MAX_SHUFFLE_PARTITIONS,
        max(MIN_SHUFFLE_PARTITIONS, -(-2 * n_pairs // ROWS_PER_SHUFFLE_PARTITION)),
    )
    pair_key_ht = pair_key_ht.repartition(target_partitions, shuffle=True)
    # Checkpoint immediately (rather than leaving this repartition lazy, fused into
    # whatever executes it far downstream) so a shuffle crash here is caught and
    # retried right at the source -- see _with_new_shuffle_fallback -- and so a
    # resumed run doesn't redo this shuffle either.
    pair_key_ht = _checkpoint(
        pair_key_ht, tmp_dir, run_tag, "pair_key_repartitioned",
        resume=resume, try_new_shuffle=True,
    )
    logger.info(
        "Repartitioned %d variant pairs into %d partitions ahead of the "
        "genotype-counting shuffle.",
        n_pairs, target_partitions,
    )

    gt_ht = create_variant_pair_genotype_ht(
        mt, pair_key_ht, mt.adj, n_pairs=n_pairs,
        tmp_dir=tmp_dir, run_tag=run_tag, resume=resume,
    )
    gt_ht = _checkpoint(
        gt_ht, tmp_dir, run_tag, "genotype_ht", resume=resume, try_new_shuffle=True
    )
    counts_ht = create_variant_pair_genotype_counts_ht(gt_ht)

    counts_ht = counts_ht.select(
        **{
            f"raw_{name}": counts_ht.gt_counts_raw[i]
            for i, name in enumerate(GENOTYPE_CLASSES)
        },
        **{
            f"adj_{name}": counts_ht.gt_counts_adj[i]
            for i, name in enumerate(GENOTYPE_CLASSES)
        },
    )

    out_ht = vp_ht.annotate(
        **counts_ht[vp_ht.locus1, vp_ht.alleles1, vp_ht.locus2, vp_ht.alleles2]
    )

    return out_ht.key_by(*orig_key)


def _flatten_for_export(ht: hl.Table) -> hl.Table:
    """
    Convert locus/allele fields to plain strings so the Table can be exported to text.

    :param ht: Table with locus1/alleles1/locus2/alleles2 fields.
    :return: Table with those fields converted to strings (flattened otherwise).
    """
    ht = ht.annotate(
        locus1=hl.str(ht.locus1),
        alleles1=hl.delimit(ht.alleles1, ","),
        locus2=hl.str(ht.locus2),
        alleles2=hl.delimit(ht.alleles2, ","),
    )
    return ht.flatten()


def main(args):
    """Annotate a variant pair list with genotype counts from gnomAD."""
    reference_genome = GENOME_BUILD_REFERENCE[args.genome_build]
    hl.init(
        log="/vp_genotype_matrix.log",
        tmp_dir=args.tmp_dir,
        default_reference=reference_genome,
    )

    logger.info(
        "Running with genome_build=%s, data_type=%s, input=%s, output=%s",
        args.genome_build,
        args.data_type,
        args.variant_pair_list_ht,
        args.output,
    )

    run_tag = args.checkpoint_tag or _run_tag_from_path(args.variant_pair_list_ht)
    logger.info(
        "Checkpoint tag: %s (pass --checkpoint-tag to override, or --no-resume to "
        "ignore existing checkpoints for this run)",
        run_tag,
    )

    vp_ht = hl.read_table(args.variant_pair_list_ht)
    out_ht = add_genotype_matrix(
        vp_ht,
        genome_build=args.genome_build,
        data_type=args.data_type,
        tmp_dir=args.tmp_dir,
        interval_padding=args.interval_padding,
        run_tag=run_tag,
        resume=not args.no_resume,
    )

    if args.output_format == "ht":
        out_ht = out_ht.checkpoint(args.output, overwrite=args.overwrite)
        logger.info("Wrote %s variant pairs with genotype counts to %s", out_ht.count(), args.output)
    else:
        _flatten_for_export(out_ht).export(args.output)
        logger.info("Exported variant pairs with genotype counts to %s", args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant-pair-list-ht",
        required=True,
        help=(
            "Path to input Hail Table of variant pairs. Must have locus1/alleles1/"
            "locus2/alleles2 row fields (any actual key is fine). Other columns are "
            "preserved in the output."
        ),
    )
    parser.add_argument(
        "--genome-build",
        choices=list(BUILD_LOADERS),
        default="grch37",
        help="Genome build of the input variant pair list. Default is grch37.",
    )
    parser.add_argument(
        "--data-type",
        choices=["exomes", "genomes"],
        default="exomes",
        help="gnomAD data type to pull genotypes from. Default is exomes.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path. A Hail Table path if --output-format ht, otherwise a text file path.",
    )
    parser.add_argument(
        "--output-format",
        choices=["ht", "tsv"],
        default="ht",
        help="Output format. Default is ht.",
    )
    parser.add_argument(
        "--tmp-dir",
        default=DEFAULT_TMP_DIR,
        help="Temporary directory for intermediate checkpoints.",
    )
    parser.add_argument(
        "--checkpoint-tag",
        default=None,
        help=(
            "Identifier used to name resumable checkpoints under --tmp-dir. "
            "Defaults to the input --variant-pair-list-ht's filename, so re-running "
            "against the same input (e.g. after a kill) automatically resumes from "
            "the last completed stage instead of starting over. Pass a different "
            "tag to force a fresh run without touching existing checkpoint files."
        ),
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help=(
            "Ignore any existing checkpoints for this --checkpoint-tag and "
            "recompute every stage from scratch (still writes fresh checkpoints "
            "for a later resumable run)."
        ),
    )
    parser.add_argument(
        "--interval-padding",
        type=int,
        default=100,
        help=(
            "Bp padding per variant for the partition-pruning step applied before "
            "the exact variant filter (see get_covering_intervals) -- one small "
            "interval per exact variant, not one bounding box per contig, so "
            "scattered variants don't force reads of everything in between. This "
            "is the main speed lever for small pair lists -- without it, Hail "
            "scans every partition of the full gnomAD MatrixTable/VDS just to "
            "pull out a handful of variants. Default 100. Pass a negative value "
            "to disable and always fall back to the plain row-key filter."
        ),
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Whether to overwrite existing output."
    )

    args = parser.parse_args()
    main(args)

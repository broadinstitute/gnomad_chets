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
access to the gnomAD v4 VDS resources (``gnomad_qc.v4``).

Known gnomAD v4 differences (READ THIS before using ``--genome-build grch38``)
--------------------------------------------------------------------------------
The counts this script produces are **not** identical to gnomAD v4 release
frequencies, in two known ways. Neither affects ``--genome-build grch37``.

1. **The high-AB het -> hom-alt correction is NOT applied.** GATK versions before
   4.1.4.1 mis-called some true hom-alt genotypes as hets with a high allele
   balance. gnomAD v4's released frequencies correct this (``gnomad_qc``
   ``generate_freq.py``): a call is reclassified het -> hom-var when it is an
   adj-passing het-ref call with ``AD[1]/DP > 0.9``, is not a true het-non-ref, the
   sample is not already on the fixed hom-alt model, and the variant's adj AF is
   above 1%. The v4 release applies that to *both* raw and adj strata.

   This script does not, so for any variant meeting those conditions its hom-var
   counts run **low** and its het counts run **high** relative to gnomAD v4 --
   i.e. the ``Aabb`` / ``aaBb`` / ``aabb`` cells are undercounted and ``AaBb`` is
   overcounted. Note the AF > 1% gate: the correction never fires for rare
   variants, so a pair list of rare candidates (the usual case here) is unaffected.
   Sanity-check any pair with a common endpoint against the release before
   reporting it.

   To enable it, the dense MatrixTable has to carry three extra fields, at which
   point the borrowed encoder applies the correction on its own:

   - ``af``: per-variant adj AF, from the release freq HT (``get_freq().freq[0].AF``)
   - ``fixed_homalt_model``: per-sample, from ``meta().project_meta``
   - ``_het_non_ref``: per-entry ``LGT.is_het_non_ref()``, which must be captured
     **before** the multi-allelic split -- ``gnomad_qc``'s loader splits internally,
     so this needs a custom split (see ``_split_variant_data_keeping_phase`` in
     ``v4/compute_vp_counts.py``).

2. **Sex-ploidy adjustment is not applied**, so chrX/chrY results are not
   v4-consistent: hemizygous calls are not converted to haploid and XY hets on
   non-PAR X are not dropped. A strict no-op on autosomes, so autosomal pair lists
   are unaffected.

The core genotype-encoding/counting logic is borrowed from
``v4/compute_vp_counts.py`` on the ``jg/v4-pipeline`` branch, with ``adj`` passed in as
an argument rather than always computed from GQ/DP/AD, since v2 and v4 gnomAD data
expose it differently. Each variant is encoded once into sets of sample indices, and
the 9 genotype cells per pair come out of set algebra over those sets -- see the
"Build-agnostic genotype counting logic" section header for what was and wasn't taken
from upstream, and why.

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
- Counting no longer shuffles per-pair genotype payloads at all. Each variant
  is encoded once into sample-index sets; a pair row carries only two integer
  indices and is counted by indexed lookup, so cost scales with the number of
  *variants* and their carrier counts rather than with the number of pairs.
  This is what makes gnomAD v4 (730,947 release samples) tractable -- the
  previous design replicated a variant's full per-sample array once per pair
  and did not finish on a single gene. ``_checkpoint`` retains the
  ``use_new_shuffle`` fast-path-with-fallback (see
  ``_with_new_shuffle_fallback``) for any shuffle Hail still chooses to
  schedule.
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

ENCODED_ROWS_PER_PARTITION = 25
"""
Target variants per partition for the encoded genotype Table.

The encoded Table inherits its partitioning from the dense MatrixTable, which for a
narrow interval is a handful of partitions regardless of payload -- on the FKRP test,
688 MiB of sample-index sets in 3 partitions, so the per-pair join could only read it
3 ways. Repartitioning after encoding fixes the read side of that.

Sizing by row count is crude here, because row sizes span five orders of magnitude
(on FKRP, a median of 10 stored sample indices against a maximum of 561k), so
partitions end up very uneven in bytes. It is still far better than the 3 partitions
that fall out otherwise. Balancing by actual set volume is what
`v4/compute_vp_counts.py` does with its variant-size-info + light/heavy split, which
is not ported here.
"""

MAX_ENCODED_PARTITIONS = 1_000
"""
Cap on the partition count from ENCODED_ROWS_PER_PARTITION, so a very large variant
list doesn't produce tens of thousands of tiny partitions whose scheduling overhead
outweighs the parallelism.
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
produced by `_count_from_sets` (index = v1_genotype * 3 + v2_genotype, where
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
### Build-agnostic genotype counting logic
###
### Borrowed from v4/compute_vp_counts.py on branch jg/v4-pipeline (ff17082).
###
### The original implementation here localized per-sample genotype arrays per variant,
### joined those arrays onto *both* endpoints of every pair, and then shuffled on
### group_by("vp_ht_idx"). That replicates a variant's per-sample array once per pair
### it participates in, so the shuffle grows as O(n_pairs x n_samples). Against gnomAD
### v2 (125,748 samples) that was merely slow; against v4 release (730,947) a single
### ~30kb gene moves tens of GB and does not finish. The encoding also retained every
### uncovered sample, because a no-call gets raw_gt=0, which is *defined* -- at a v4
### exome site that is often 300k+ samples of pure padding per variant.
###
### The replacement encodes each *variant* once into sets of sample indices, keys pairs
### by (v1_idx, v2_idx), and derives all 9 genotype cells by set algebra -- including
### the hom-ref/hom-ref cell, which falls out of n_samples by inclusion/exclusion
### rather than by enumerating hom-ref samples. Per-variant data is joined by index and
### never replicated per pair, so cost scales with the number of variants and their
### carrier counts rather than with the number of pairs.
###
### Borrowed, in the order they appear below:
###   _encode_genotype_sets_by_var_idx <- _encode_genotype_sets_by_var_idx
###                                       (adapted: takes this script's `adj_expr`
###                                       parameter instead of the upstream
###                                       `use_precomputed_adj` flag; assigns v_idx
###                                       itself rather than joining a var_idx Table,
###                                       see below; phase sidecar and high-AB
###                                       correction dropped)
###   _count_from_sets                 <- _count_from_sets (verbatim)
###   _project_count_fields            <- _project_count_fields (minus phased_het)
###   _drop_pairs_missing_v_idx        <- _drop_pairs_missing_v_idx (verbatim)
###   count_pairs_via_index            <- count_all_pairs_via_index + the non-pop,
###                                       non-phase half of _count_pairs_via_index,
###                                       merged into one function
###
### Upstream's ``_create_var_idx_ht`` (``mt.rows().add_index("var_idx")``) is
### deliberately *not* used. It forces its own pass over the MatrixTable, and here `mt`
### is a lazy VDS-densify, so materializing it densified the VDS a second time -- on the
### FKRP test that was 57 minutes to produce a 0 MiB, 2,390-row table, a third of total
### runtime. The encoder assigns ``v_idx`` with ``add_index`` on its own localized rows
### instead (identical values: same rows, same order), keeps ``locus``/``alleles``, and
### the ``(locus, alleles) -> v_idx`` lookup is projected back off the checkpoint. One
### densify. Upstream doesn't have this problem because it encodes from an already
### materialized dense MT.
###
### Deliberately NOT borrowed:
###   - The v4 high-AB het -> hom-alt correction. gnomAD v4's released frequencies
###     reclassify some high-allele-balance hets as hom-var (a GATK <4.1.4.1 artifact).
###     Upstream applies it only when the dense MT carries `af`, `fixed_homalt_model`
###     and `_het_non_ref`; neither loader in this script supplies any of those, so
###     upstream would skip it too, and it is omitted rather than carried as dead code.
###     Consequence: for variants with adj AF > 1% *and* allele balance > 0.9, hom-var
###     counts here run lower (and het counts higher) than gnomAD v4 release
###     frequencies. It is gated on AF > 1%, so it never fires for the rare variants
###     this script is usually pointed at. Enabling it means joining the release freq
###     HT and sample meta onto the dense MT first.
###   - Sex-ploidy adjustment for hemizygous chrX/chrY calls. Upstream applies it at
###     densify; it is a strict no-op on autosomes. chrX/chrY output from this script
###     is therefore not v4-release-consistent.
###   - Per-population stratification, physical (PGT/PID) phase counts, PBT-sample
###     subtraction, and the light/heavy partition split -- all scale or feature
###     concerns for the genome-wide pipeline, none of them needed for a user-supplied
###     pair list.
###
### Counting semantics are unchanged from the original implementation in this file.
### The per-sample genotype classification, the adj gating, and the handling of
### no-call samples all agree cell-for-cell; see the equivalence test noted in the
### commit message. The one visible change is that the 18 output columns are now
### int32 rather than int64 (the counts are bounded by the sample count).
########################################################################################
def _encode_genotype_sets_by_var_idx(
    mt: hl.MatrixTable,
    adj_expr: hl.expr.BooleanExpression,
) -> hl.Table:
    """
    Encode genotypes as per-variant sample-index sets keyed by ``var_idx``.

    Adapted from ``_encode_genotype_sets_by_var_idx`` in v4/compute_vp_counts.py. The
    upstream version selects between ``mt.adj`` and a freshly computed
    ``get_adj_expr(...)`` via a ``use_precomputed_adj`` flag; this one takes the
    ``adj_expr`` parameter that the rest of this script already threads through, since
    both loaders here annotate ``adj`` onto the MatrixTable themselves. The upstream
    physical-phase sidecar and high-AB correction are dropped (see the section header).

    Genotypes are encoded as:

        - missing (None) = hom-ref (space saving)
        - 0 = missing data (no GT call, or failed adj for ``adj_gt``)
        - 1 = het
        - 2 = hom-var

    Each sample falls in exactly one of these 7 disjoint per-variant categories (the
    implicit "adj-PASS-0/0" majority is never stored):

        cat | GT       | adj  | raw_gt | adj_gt | stored in
        ----+----------+------+--------+--------+--------------------------------
         1  | no entry |  -   |   NA   |   NA   | all_samples (only)
         2  | 0/0      | FAIL |   NA   |    0   | raw_hr_adj_missing (only)
         3  | 0/1      | FAIL |    1   |    0   | all_samples, raw_het
         4  | 0/1      | PASS |    1   |    1   | all_samples, raw_het, adj_het
         5  | 1/1      | FAIL |    2   |    0   | all_samples, raw_hv
         6  | 1/1      | PASS |    2   |    2   | all_samples, raw_hv, adj_hv
         7  | 0/0      | PASS |   NA   |   NA   | (implicit majority)

    ``all_samples`` stores cats 1, 3-6 directly (positive form), disjoint from
    ``raw_hr_adj_missing`` (cat 2). cat 7 is reconstructed by :func:`_count_from_sets`
    as ``n_samples - (cats 1-6)``.

    Keeping cat 1 (no entry) inside ``all_samples`` is what stops uncallable samples
    leaking into the hom-ref/hom-ref cell: they land in ``D'_v`` and so fall outside
    the hom-ref set the count kernel derives by exclusion.

    :param mt: MatrixTable with variant data. Row key must be ``(locus, alleles)``,
        with a ``GT`` entry field.
    :param adj_expr: Boolean expression on ``mt`` indicating high-quality genotypes.
    :return: Table keyed by ``v_idx``, with a ``samples`` global, holding per-variant
        sample-index sets plus the ``locus`` / ``alleles`` they came from.
    """
    gt_count_expr = (
        hl.case(missing_false=True)
        .when(~hl.is_missing(mt.GT) & ~mt.GT.is_non_ref(), hl.missing(hl.tint32))
        .when(mt.GT.is_het(), 1)
        .when(mt.GT.is_hom_var(), 2)
        .default(0)
    )
    adj_gt_count_expr = hl.if_else(
        adj_expr, gt_count_expr, 0, missing_false=True
    )

    mt = mt.select_entries(raw_gt=gt_count_expr, adj_gt=adj_gt_count_expr)
    ht = mt.localize_entries("_entries", "samples")

    gt = hl.enumerate(ht._entries)
    # cats 1-6, used only to derive n_with_data; never stored.
    not_adj_hom_ref = gt.filter(
        lambda x: hl.is_missing(x[1])
        | hl.is_defined(x[1].raw_gt)
        | hl.is_defined(x[1].adj_gt)
    )
    # cats 1, 3-6: cat 1 (no entry) plus cats 3-6 (raw_gt defined only for het /
    # hom-var). Excludes cat 2, which is stored separately in raw_hr_adj_missing.
    not_adj_hom_ref_no_F = gt.filter(
        lambda x: hl.is_missing(x[1]) | hl.is_defined(x[1].raw_gt)
    )
    raw_hr_adj_missing = gt.filter(
        lambda x: hl.is_defined(x[1])
        & hl.is_missing(x[1].raw_gt)
        & (x[1].adj_gt == 0)
    )

    ht = ht.select(
        all_samples=hl.set(not_adj_hom_ref_no_F.map(lambda x: x[0])),
        n_with_data=hl.int32(not_adj_hom_ref.length()),
        raw_hr_adj_missing=hl.set(raw_hr_adj_missing.map(lambda x: x[0])),
        n_raw_hr_adj_missing=hl.int32(raw_hr_adj_missing.length()),
        raw_het=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].raw_gt == 1).map(lambda x: x[0])
        ),
        raw_hv=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].raw_gt == 2).map(lambda x: x[0])
        ),
        adj_het=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].adj_gt == 1).map(lambda x: x[0])
        ),
        adj_hv=hl.set(
            not_adj_hom_ref.filter(lambda x: x[1].adj_gt == 2).map(lambda x: x[0])
        ),
    )

    # Index in place rather than joining a separately-computed var_idx Table, so the
    # MatrixTable (and therefore the VDS densify behind it) is only ever walked once.
    # add_index over these localized rows gives exactly what mt.rows().add_index()
    # would: same row set, same order. locus/alleles stay so the (locus, alleles) ->
    # v_idx lookup can be projected off the checkpoint; _project_count_fields drops
    # them again before the per-pair join, so they cost nothing there.
    ht = ht.add_index("v_idx")
    return ht.key_by("v_idx")


def _count_from_sets(
    v1_het: hl.expr.SetExpression,
    v1_hv: hl.expr.SetExpression,
    v1_all: hl.expr.SetExpression,
    v1_n: hl.expr.Int32Expression,
    v1_F: hl.expr.SetExpression,
    v1_n_F: hl.expr.Int32Expression,
    v2_het: hl.expr.SetExpression,
    v2_hv: hl.expr.SetExpression,
    v2_all: hl.expr.SetExpression,
    v2_n: hl.expr.Int32Expression,
    v2_F: hl.expr.SetExpression,
    v2_n_F: hl.expr.Int32Expression,
    n_samples: hl.expr.Int32Expression,
    *,
    include_raw_hr_adj_missing: bool,
) -> hl.expr.ArrayExpression:
    """
    Compute the 9-element genotype count array from per-variant sample sets.

    Borrowed verbatim from ``_count_from_sets`` in v4/compute_vp_counts.py.

    Count array layout: ``[AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]``,
    matching GENOTYPE_CLASSES, where A/a = v1 ref/alt and B/b = v2 ref/alt.

    Per-variant storage (all sets positive form -- sample indices held directly):

      - ``v_all`` = ``A_v`` = cats 1, 3-6, disjoint from ``v_F``.
        ``v_n`` = ``|cats 1-6|`` = n_with_data.
      - ``v_F`` = ``raw_hr_adj_missing`` = cat 2. ``v_n_F`` = |cat 2|.

    The "real" not-adj-hom-ref set is ``D'_v = A_v u F_v = cats 1-6``; intersections
    over D' decompose into the A/F cross terms.

    Two modes:

      - ``False`` (ADJ cells): hom-ref = adj-PASS-0/0 = ``N \\ D'``. Pass
        ``v_het = adj_het``, ``v_hv = adj_hv``.
      - ``True`` (RAW cells): hom-ref = raw-0/0 = ``(adj-PASS-0/0) u F``. Pass
        ``v_het = raw_het``, ``v_hv = raw_hv``. F adds the adj-fail hom-ref samples
        to the hom-ref pool.

    :param v1_het, v1_hv: Sample sets for v1 het / hom-var.
    :param v1_all, v1_n: v1 ``all_samples`` (= A_v) set + ``n_with_data``.
    :param v1_F, v1_n_F: v1 ``raw_hr_adj_missing`` (= cat 2) set + size.
    :param v2_het, v2_hv, v2_all, v2_n, v2_F, v2_n_F: same for v2.
    :param n_samples: Total number of samples in the cohort.
    :param include_raw_hr_adj_missing: ``True`` for raw cells, ``False`` for adj.
    :return: 9-element genotype count array.
    """
    # D'_v = A_v u F_v (disjoint per variant). Decompose |D'_v1 n D'_v2| into the four
    # cross terms; reuse the A-F and F-F pieces in the raw-cells branch below.
    a1_isect_a2 = v1_all.intersection(v2_all).length()
    a1_isect_f2 = v1_all.intersection(v2_F).length()
    f1_isect_a2 = v1_F.intersection(v2_all).length()
    f1_isect_f2 = v1_F.intersection(v2_F).length()
    d1_isect_d2 = a1_isect_a2 + a1_isect_f2 + f1_isect_a2 + f1_isect_f2

    # |H_v1 n H_v2| where H = N \ D'  (adj-PASS-0/0 at both).
    #   = N - |D'_v1 u D'_v2| = N - |D'_v1| - |D'_v2| + |D'_v1 n D'_v2|
    h1_isect_h2 = n_samples - v1_n - v2_n + d1_isect_d2

    # Carrier-carrier cells.
    het_het = v1_het.intersection(v2_het).length()
    het_hv = v1_het.intersection(v2_hv).length()
    hv_het = v1_hv.intersection(v2_het).length()
    hv_hv = v1_hv.intersection(v2_hv).length()

    # Edge cells (adj component): |carrier_v n H_other|
    #   = |carrier_v| - |carrier_v n D'_other|
    # where |carrier_v n D'_other| = |carrier_v n A_other| + |carrier_v n F_other|.
    v1_het_in_d2 = (
        v1_het.intersection(v2_all).length() + v1_het.intersection(v2_F).length()
    )
    v1_hv_in_d2 = (
        v1_hv.intersection(v2_all).length() + v1_hv.intersection(v2_F).length()
    )
    v2_het_in_d1 = (
        v2_het.intersection(v1_all).length() + v2_het.intersection(v1_F).length()
    )
    v2_hv_in_d1 = (
        v2_hv.intersection(v1_all).length() + v2_hv.intersection(v1_F).length()
    )
    v1_het_in_h2 = v1_het.length() - v1_het_in_d2
    v1_hv_in_h2 = v1_hv.length() - v1_hv_in_d2
    v2_het_in_h1 = v2_het.length() - v2_het_in_d1
    v2_hv_in_h1 = v2_hv.length() - v2_hv_in_d1

    if include_raw_hr_adj_missing:
        # RAW cells: hom-ref pool extended by F (cat 2) at each variant.
        # AABB_raw = |(H1 u F1) n (H2 u F2)| (H_v and F_v disjoint per variant, so the
        # union sums) = |H1nH2| + |H1nF2| + |F1nH2| + |F1nF2|. Reuse the cross terms:
        #   |D'_v1 n F_v2| = |A_v1 n F_v2| + |F_v1 n F_v2|
        #   |F_v1 n D'_v2| = |F_v1 n A_v2| + |F_v1 n F_v2|
        h1_isect_f2 = v2_n_F - (a1_isect_f2 + f1_isect_f2)
        f1_isect_h2 = v1_n_F - (f1_isect_a2 + f1_isect_f2)
        hom_ref_both = h1_isect_h2 + h1_isect_f2 + f1_isect_h2 + f1_isect_f2

        # Edge cells (raw): also extend hom-ref-at-other by F_other.
        v1_het_in_f2 = v1_het.intersection(v2_F).length()
        v1_hv_in_f2 = v1_hv.intersection(v2_F).length()
        v2_het_in_f1 = v2_het.intersection(v1_F).length()
        v2_hv_in_f1 = v2_hv.intersection(v1_F).length()

        v1_het_homref_v2 = v1_het_in_h2 + v1_het_in_f2
        v1_hv_homref_v2 = v1_hv_in_h2 + v1_hv_in_f2
        v2_het_homref_v1 = v2_het_in_h1 + v2_het_in_f1
        v2_hv_homref_v1 = v2_hv_in_h1 + v2_hv_in_f1
    else:
        # ADJ cells: hom-ref = H only.
        hom_ref_both = h1_isect_h2
        v1_het_homref_v2 = v1_het_in_h2
        v1_hv_homref_v2 = v1_hv_in_h2
        v2_het_homref_v1 = v2_het_in_h1
        v2_hv_homref_v1 = v2_hv_in_h1

    return hl.array([
        hom_ref_both,       # AABB
        v2_het_homref_v1,   # AABb
        v2_hv_homref_v1,    # AAbb
        v1_het_homref_v2,   # AaBB
        het_het,            # AaBb
        het_hv,             # Aabb
        v1_hv_homref_v2,    # aaBB
        hv_het,             # aaBb
        hv_hv,              # aabb
    ])


_COUNT_FROM_SETS_FIELDS = (
    "raw_het", "raw_hv",
    "all_samples", "n_with_data",
    "raw_hr_adj_missing", "n_raw_hr_adj_missing",
    "adj_het", "adj_hv",
)
"""Fields on the encoded Table that the counting path actually reads."""


def _project_count_fields(encoded_ht: hl.Table) -> hl.Table:
    """
    Restrict an encoded Table to the fields the count path reads.

    Borrowed from v4/compute_vp_counts.py (minus its ``phased_het`` field, which this
    script does not encode). Not a no-op here: it drops the ``locus`` / ``alleles`` the
    encoder retains for the var_idx lookup, so they never ride through the per-pair
    join.

    :param encoded_ht: Encoded genotype Table keyed by ``v_idx``.
    :return: ``encoded_ht`` with only the count fields.
    """
    return encoded_ht.select(*_COUNT_FROM_SETS_FIELDS)


def _drop_pairs_missing_v_idx(vp_ht: hl.Table, caller: str) -> hl.Table:
    """
    Drop pairs whose ``v1_idx`` or ``v2_idx`` is missing, logging the count.

    Borrowed verbatim from v4/compute_vp_counts.py.

    Pairs land here when the pair list's ``(locus, alleles)`` doesn't appear in
    ``var_idx_ht`` (the index over the dense MatrixTable's rows) -- most often because
    the pair list and the callset split multi-allelics differently, or because the
    variant simply isn't in the callset. Such pairs are dropped here and therefore end
    up with *missing* counts on the output table rather than zeros: the row is
    preserved (`add_genotype_matrix` annotates onto the input list), but a missing
    count means "this variant was never found", which is a different statement from a
    zero count meaning "found, and nobody carried it".

    :param vp_ht: Pair Table annotated with ``v1_idx`` / ``v2_idx``.
    :param caller: Caller name for the log message.
    :return: ``vp_ht`` filtered to rows with both v_idx fields defined.
    """
    n_missing = vp_ht.aggregate(
        hl.agg.count_where(
            hl.is_missing(vp_ht.v1_idx) | hl.is_missing(vp_ht.v2_idx)
        )
    )
    if n_missing > 0:
        logger.warning(
            "%s: %d pairs have v_idx missing for v1 and/or v2 (variant not present "
            "in the callset); dropping them here, so they will carry missing (not "
            "zero) genotype counts on the output.",
            caller, n_missing,
        )
        vp_ht = vp_ht.filter(
            hl.is_defined(vp_ht.v1_idx) & hl.is_defined(vp_ht.v2_idx)
        )
    return vp_ht


def count_pairs_via_index(
    vp_ht: hl.Table,
    var_idx_ht: hl.Table,
    encoded_gt_ht: hl.Table,
) -> hl.Table:
    """
    Count every pair by indexed lookup of the two endpoints' encoded sample sets.

    Borrowed from v4/compute_vp_counts.py: ``count_all_pairs_via_index`` plus the
    full-cohort half of ``_count_pairs_via_index``, merged (the upstream split exists
    to share the inner function with the per-population and light/heavy paths, neither
    of which is ported here).

    There is no explicit shuffle and no per-pair sample array: each pair row carries
    only two int64 indices, and Hail's planner picks the join against the (small,
    one-row-per-variant) encoded table.

    Pair rows are keyed by ``(v1_idx, v2_idx)`` alone for the duration of the count,
    and ``locus1``/``alleles1``/``locus2``/``alleles2`` are reconstructed from the
    encoded Table at the very end. This is "design C" from the upstream benchmark of
    pair-table key layouts (5 genes, 1.67M pairs): carrying both the indices *and* the
    locus/alleles through -- which is what this function used to do -- pays for the
    int64 keys without getting the slim rows, and was the only layout that couldn't
    clear the light count step on a disk-constrained cluster. Upstream is adopting C
    for the same reason. The end remap costs one sort of narrow rows (two int64 keys
    and two 9-element count arrays) instead of dragging four locus/allele fields
    through every join. At gene scale the difference is unmeasurable; it is there for
    the large pair lists this script also gets pointed at.

    :param vp_ht: Variant pair Table keyed by locus1/alleles1/locus2/alleles2.
    :param var_idx_ht: ``(locus, alleles) -> v_idx`` lookup, projected off the
        encoded Table. (Upstream names this field ``var_idx`` because it comes from a
        separate ``_create_var_idx_ht`` pass; here it is the encoder's own ``v_idx``.)
    :param encoded_gt_ht: Encoded genotype Table keyed by ``v_idx``.
    :return: Pair Table with ``gt_counts_raw`` / ``gt_counts_adj`` 9-element arrays.
    """
    n_samples = hl.int32(encoded_gt_ht.index_globals().samples.length())

    # v_idx -> (locus, alleles), for rebuilding the pair key after counting. Taken
    # before _project_count_fields drops those columns; it is a two-field projection
    # off an already-materialized checkpoint, not a recompute.
    variant_ht = encoded_gt_ht.select("locus", "alleles")

    encoded_gt_ht = _project_count_fields(encoded_gt_ht)

    vp_ht = vp_ht.annotate(
        v1_idx=var_idx_ht[vp_ht.locus1, vp_ht.alleles1].v_idx,
        v2_idx=var_idx_ht[vp_ht.locus2, vp_ht.alleles2].v_idx,
    )
    vp_ht = _drop_pairs_missing_v_idx(vp_ht, "count_pairs_via_index")

    # Drop to the two indices only. The locus/alleles fields are the previous key, so
    # they can only be dropped by keying away from them first.
    vp_ht = vp_ht.key_by("v1_idx", "v2_idx").select()

    v1 = encoded_gt_ht[vp_ht.v1_idx]
    v2 = encoded_gt_ht[vp_ht.v2_idx]
    vp_ht = vp_ht.select(
        gt_counts_raw=_count_from_sets(
            v1.raw_het, v1.raw_hv, v1.all_samples, v1.n_with_data,
            v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
            v2.raw_het, v2.raw_hv, v2.all_samples, v2.n_with_data,
            v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
            n_samples,
            include_raw_hr_adj_missing=True,
        ),
        gt_counts_adj=_count_from_sets(
            v1.adj_het, v1.adj_hv, v1.all_samples, v1.n_with_data,
            v1.raw_hr_adj_missing, v1.n_raw_hr_adj_missing,
            v2.adj_het, v2.adj_hv, v2.all_samples, v2.n_with_data,
            v2.raw_hr_adj_missing, v2.n_raw_hr_adj_missing,
            n_samples,
            include_raw_hr_adj_missing=False,
        ),
    )

    # Rebuild the pair key the caller expects. Orientation survives the round trip:
    # v1_idx was looked up from locus1/alleles1, so it maps back to them. Listing only
    # the non-key count fields in select() is deliberate -- Hail's check_keys rejects
    # a select() that names the table's own key fields.
    v1_variant = variant_ht[vp_ht.v1_idx]
    v2_variant = variant_ht[vp_ht.v2_idx]
    vp_ht = vp_ht.key_by(
        locus1=v1_variant.locus,
        alleles1=v1_variant.alleles,
        locus2=v2_variant.locus,
        alleles2=v2_variant.alleles,
    ).select("gt_counts_raw", "gt_counts_adj")

    return vp_ht.cache()


def create_variant_pair_genotype_counts_ht(
    mt: hl.MatrixTable,
    vp_ht: hl.Table,
    adj_expr: hl.expr.BooleanExpression,
    n_variants: Optional[int] = None,
    tmp_dir: str = DEFAULT_TMP_DIR,
    run_tag: Optional[str] = None,
    resume: bool = True,
) -> hl.Table:
    """
    Create a variant pair genotype counts Table from a MatrixTable and a pair list.

    Encodes each variant once into sample-index sets, then counts every pair off those
    sets (see the section header above).

    :param mt: MatrixTable with variant data. Row key must be ``(locus, alleles)``.
    :param vp_ht: Table of variant pairs, keyed by locus1/alleles1/locus2/alleles2 and
        carrying no other fields -- callers should strip extra columns off first.
    :param adj_expr: Boolean expression on ``mt`` indicating high-quality genotypes.
    :param n_variants: Number of variants being encoded, if the caller already knows it
        (`add_genotype_matrix` counts them anyway). Used only to size the encoded
        Table's partitioning; when omitted, whatever partitioning falls out of the
        dense MatrixTable is left alone.
    :param tmp_dir: Base temporary directory, used for checkpointing (see
        `_checkpoint`).
    :param run_tag: Stable identifier for this input/run, for resumable checkpointing
        (see `_checkpoint`); None disables resuming.
    :param resume: See `_checkpoint`.
    :return: Pair Table with ``gt_counts_raw`` / ``gt_counts_adj`` 9-element arrays.
    """
    encoded_gt_ht = _encode_genotype_sets_by_var_idx(mt, adj_expr)

    if n_variants:
        target_partitions = min(
            MAX_ENCODED_PARTITIONS,
            max(1, -(-n_variants // ENCODED_ROWS_PER_PARTITION)),
        )
        if target_partitions > encoded_gt_ht.n_partitions():
            logger.info(
                "Repartitioning the encoded genotype Table from %d to %d partitions "
                "so the per-pair join can read it in parallel.",
                encoded_gt_ht.n_partitions(), target_partitions,
            )
            encoded_gt_ht = encoded_gt_ht.repartition(target_partitions, shuffle=True)

    # The densify + encode is by far the most expensive stage, and the one worth never
    # redoing on a resumed run. This checkpoint is also what forces the repartition
    # above, hence try_new_shuffle.
    encoded_gt_ht = _checkpoint(
        encoded_gt_ht, tmp_dir, run_tag, "encoded_gt_sets",
        resume=resume, try_new_shuffle=True,
    )

    # Projected straight off the checkpoint rather than computed from `mt`, so this
    # costs a small read of two key fields instead of a second densify.
    var_idx_ht = encoded_gt_ht.key_by("locus", "alleles").select("v_idx")

    return count_pairs_via_index(vp_ht, var_idx_ht, encoded_gt_ht)


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

    logger.warning(
        "The gnomAD v4 high-AB het -> hom-alt correction is NOT applied to these "
        "counts. For variants with adj AF > 1% and allele balance > 0.9, hom-var "
        "counts will run low and het counts high relative to gnomAD v4 release "
        "frequencies. Rare variants are unaffected (the correction is AF-gated). "
        "See the 'Known gnomAD v4 differences' section at the top of this file."
    )

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

    # Only keep the pair key for internal processing so the helpers above (which add
    # their own fields, e.g. 'v1_idx') don't clash with any extra columns the user's
    # input table already has. Note: vp_ht's *actual* key may be something else
    # entirely (e.g. locus/alleles/gene/gene_id) -- re-key by the pair fields first,
    # then select() with no args to drop everything else.
    pair_key_ht = vp_ht.key_by("locus1", "alleles1", "locus2", "alleles2").select()

    # One count action, used for the log line below and for the canonical-order
    # check. The repartition/skew-salting that used to live here went away with the
    # shuffle it existed to balance -- pairs now carry only two integer indices
    # through an indexed join, so there is no large per-pair payload left to spread.
    #
    # The canonical-order check only warns: counts are computed for whatever
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

    logger.info("Counting genotypes for %d variant pairs...", n_pairs)
    counts_ht = create_variant_pair_genotype_counts_ht(
        mt, pair_key_ht, mt.adj, n_variants=n_variants,
        tmp_dir=tmp_dir, run_tag=run_tag, resume=resume,
    )

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
        help=(
            "Genome build of the input variant pair list. Default is grch37. "
            "grch38 counts against gnomAD v4; see the 'Known gnomAD v4 differences' "
            "section at the top of this file for where those counts intentionally "
            "diverge from v4 release frequencies (high-AB het correction, sex "
            "ploidy)."
        ),
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

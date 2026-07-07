"""Trio phasing and trio-vs-gnomAD phase comparison for gnomAD v4.

Consolidates the v2 trio/PBT analysis (``phase_by_transmission.py`` plus the
``create_pbt_*`` steps and the gnomAD-comparison export) into one script, on
the v4 VDS-based pipeline. It reuses, rather than duplicates, the v4 counts
machinery (``create_vp_matrix``) and EM phasing (``phase_gnomad``).

Trio side (high-quality samples, incl. unreleasable):

* ``--create-pbt-trio-matrix`` — densify the gnomAD v4 VDS over the finalized
  trio samples, compute ``adj``, build a trio MatrixTable, and phase each trio
  by transmission (``phase_trio_matrix_by_transmission``).
* ``--explode-pbt`` — explode the trio matrix into a per-sample MatrixTable.
* ``--phase-multi-families`` — consensus phase for samples in multiple trios.
* ``--derive-trio-vps`` — variant pairs co-carried by trio probands.
* ``--call-trio-chet`` — per-pair cis/trans (``n_same_hap`` / ``n_chet``)
  counts from the proband transmission phase (the trio "truth").

Comparison side (gnomAD = release samples, PBT members removed):

* ``--gnomad-counts-no-pbt`` — gnomAD genotype counts on the trio VPs over
  release samples minus all PBT trio members.
* ``--phase-gnomad-counts`` — EM-phase those counts.
* ``--export-comparison`` — join trio truth vs gnomAD ``p_chet`` → HT + TSV.

The v4 VDS is read split, so no ``split_multi_hts`` pass is needed. The
finalized pedigree is exomes-derived in gnomad_qc v4.
"""

import argparse
import logging
import os
import tempfile
import timeit

import hail as hl
from gnomad.utils.annotations import get_adj_expr
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.sample_qc import pedigree, trios

from gnomad_chets.v4.create_vp_matrix import (
    compute_counts_by_pop,
    count_all_pairs_via_index,
    create_variant_pair_filter_ht,
    create_variant_pair_ht,
    encode_genotypes,
    filter_pairs_by_an_pct,
)
from gnomad_chets.v4.phase_gnomad import get_phased_gnomad_ht
from gnomad_chets.v4.resources import (
    DATA_TYPE_CHOICES,
    DEFAULT_DATA_TYPE,
    DEFAULT_TMP_DIR,
    GLOBAL_POP,
    TEST_INTERVALS,
    _get_output_postfix,
    get_pops,
    get_sample_pop_ht,
    get_trio_phasing_resources,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("trio_phasing")
logger.setLevel(logging.INFO)


def _complete_trio_samples(ped: hl.Pedigree) -> list:
    """Return the sorted unique sample IDs of all complete trios in ``ped``."""
    return sorted(
        {
            s
            for trio in ped.complete_trios()
            for s in (trio.s, trio.pat_id, trio.mat_id)
            if s is not None
        }
    )


def _samples_ht(samples: list) -> hl.Table:
    """Build a Table keyed by ``s`` from a Python list of sample IDs."""
    return hl.Table.parallelize(
        [{"s": s} for s in samples],
        schema=hl.tstruct(s=hl.tstr),
        key="s",
    )


def phase_trio_matrix(mt: hl.MatrixTable, ped: hl.Pedigree) -> hl.MatrixTable:
    """Build a PBT-phased trio MatrixTable from a dense sample MatrixTable.

    Annotates per-entry ``adj``, keeps only the ``GT`` / ``adj`` entry fields,
    drops monomorphic rows, builds the trio MatrixTable for complete trios and
    phases each trio by transmission.

    :param mt: Dense per-sample MatrixTable with ``GT``, ``GQ``, ``DP``, ``AD``
        entries, keyed by ``s``.
    :param ped: Pedigree to build trios from. ``hl.trio_matrix`` filters it to
        the samples present in ``mt`` and keeps only complete trios.
    :return: Trio MatrixTable with ``PBT_GT`` added to each member entry.
    """
    mt = mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD))
    mt = mt.select_entries("GT", "adj")
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    tm = hl.trio_matrix(mt, ped, complete_trios=True)
    return hl.experimental.phase_trio_matrix_by_transmission(tm)


def explode_pbt_trio_matrix(tm: hl.MatrixTable) -> hl.MatrixTable:
    """Explode a PBT-phased trio MatrixTable into a per-sample MatrixTable.

    The trio-level ``trio_adj`` (all three members adj) is computed before the
    explode and carried onto each per-sample entry.

    :param tm: PBT-phased trio MatrixTable (from :func:`phase_trio_matrix`).
    :return: Per-sample MatrixTable with ``GT`` / ``adj`` / ``PBT_GT`` /
        ``trio_adj`` entries and a ``source_trio`` column struct.
    """
    tm = tm.annotate_entries(
        trio_adj=tm.proband_entry.adj & tm.father_entry.adj & tm.mother_entry.adj
    )
    pmt = hl.experimental.explode_trio_matrix(tm, keep_trio_entries=True)
    return pmt.transmute_entries(trio_adj=pmt.source_trio_entry.trio_adj)


def phase_multi_offspring_families(pmt: hl.MatrixTable) -> hl.MatrixTable:
    """Compute consensus PBT phase for samples in multiple consistent trios.

    Keeps samples that appear in more than one trio where every trio shares
    the same parents (the handful of samples with conflicting parents are
    dropped), then collapses each sample's per-trio phased calls into a
    consensus genotype. The consensus prefers phased calls, then the call with
    the most votes.

    :param pmt: Exploded per-sample PBT MatrixTable (from
        :func:`explode_pbt_trio_matrix`).
    :return: MatrixTable grouped to one column per sample with
        ``consensus_gt``, ``phase_concordance`` and ``discordant_gts`` entries.
    """
    nt_samples = pmt.cols()
    nt_samples = nt_samples.group_by("s").aggregate(
        trios=hl.agg.collect(nt_samples.source_trio)
    )
    # Drop samples in >1 trio whose parents differ across trios.
    nt_samples = nt_samples.filter(
        (hl.len(nt_samples.trios) > 1)
        & nt_samples.trios[1:].any(
            lambda x: (x.mother.s != nt_samples.trios[0].mother.s)
            | (x.father.s != nt_samples.trios[0].father.s)
        ),
        keep=False,
    )
    pmt = pmt.filter_cols(hl.is_defined(nt_samples[pmt.col_key]))

    # Group each sample's per-trio entries, keeping all phased calls in an
    # array; the consensus prefers phased calls then the most-voted call.
    pmt = pmt.group_cols_by("s").aggregate(
        PBT_GTs=hl.agg.filter(hl.is_defined(pmt.PBT_GT), hl.agg.collect(pmt.PBT_GT))
    )
    gt_counter = hl.sorted(
        hl.array(pmt.PBT_GTs.group_by(lambda x: x).map_values(lambda x: hl.len(x))),
        key=lambda x: x[0].phased * 100 + x[1],
        reverse=True,
    )
    phased_gt_counts = gt_counter.filter(lambda x: x[0].phased).map(lambda x: x[1])
    return pmt.annotate_entries(
        consensus_gt=gt_counter.map(lambda x: x[0]).find(lambda x: True),
        phase_concordance=phased_gt_counts.find(lambda x: True)
        / hl.sum(phased_gt_counts),
        discordant_gts=hl.len(
            hl.set(
                pmt.PBT_GTs.map(
                    lambda x: hl.if_else(x.phased, hl.call(x[0], x[1]), x)
                )
            )
        )
        > 1,
    )


def _trio_chet_counts(shared):
    """(n_same_hap, n_chet) over a shared-proband array — raw and adj."""
    return hl.struct(
        raw=hl.struct(
            n_same_hap=shared.filter(lambda x: x.same_hap).length(),
            n_chet=shared.filter(lambda x: ~x.same_hap).length(),
        ),
        adj=hl.struct(
            n_same_hap=shared.filter(lambda x: x.adj & x.same_hap).length(),
            n_chet=shared.filter(lambda x: x.adj & ~x.same_hap).length(),
        ),
    )


def call_trio_chet(
    pmt: hl.MatrixTable, vp_ht: hl.Table, pop_ht: hl.Table = None, pops=None
) -> hl.Table:
    """Count cis/trans probands per variant pair from transmission phase.

    For each proband phased het/het at both variants, the pair is *same_hap*
    (cis) when the alt is on the same parental haplotype at both sites
    (``PBT_GT[0]`` equal), else *chet* (trans). Counts are reported for raw and
    for adj (all three trio members adj).

    :param pmt: Per-sample PBT MatrixTable restricted to deduplicated probands,
        with ``PBT_GT`` / ``trio_adj`` entries.
    :param vp_ht: Trio variant-pair list keyed by
        ``(locus1, alleles1, locus2, alleles2)``.
    :param pop_ht: When provided (``s`` → ``pop``; see
        :func:`resources.get_sample_pop_ht`), also stratify the cis/trans
        counts by proband genetic-ancestry group into ``raw_by_pop`` /
        ``adj_by_pop`` dicts keyed by pop (incl. ``GLOBAL_POP`` = all probands),
        so the comparison can be run per pop against the matching per-pop gnomAD
        EM. Probands with no pop label contribute only to ``GLOBAL_POP``.
    :param pops: Optional subset of groups to keep in ``raw_by_pop`` /
        ``adj_by_pop`` (``GLOBAL_POP`` always kept); ``None`` = all present.
    :return: ``vp_ht`` with ``raw`` / ``adj`` structs (``n_same_hap``,
        ``n_chet``) — plus ``raw_by_pop`` / ``adj_by_pop`` when ``pop_ht`` is
        given — filtered to pairs with at least one phased-het proband.
    """
    et = pmt.select_entries("PBT_GT", "trio_adj").entries()
    et = et.filter(et.PBT_GT.phased & et.PBT_GT.is_het())
    carrier_fields = dict(s=et.s, hap0=et.PBT_GT[0], adj=et.trio_adj)
    if pop_ht is not None:
        carrier_fields["pop"] = pop_ht[et.s].pop
    carriers = et.group_by(et.locus, et.alleles).aggregate(
        carriers=hl.agg.collect(hl.struct(**carrier_fields))
    )
    carriers = carriers.checkpoint(
        hl.utils.new_temp_file("call_trio_chet.carriers", "ht")
    )
    empty = hl.empty_array(carriers.carriers.dtype.element_type)

    vp = vp_ht.annotate(
        _c1=hl.or_else(carriers[vp_ht.locus1, vp_ht.alleles1].carriers, empty),
        _c2=hl.or_else(carriers[vp_ht.locus2, vp_ht.alleles2].carriers, empty),
    )
    vp = vp.annotate(_c2_map=hl.dict(vp._c2.map(lambda x: (x.s, x))))

    def _shared_struct(x):
        fields = dict(
            same_hap=x.hap0 == vp._c2_map[x.s].hap0,
            adj=x.adj & vp._c2_map[x.s].adj,
        )
        if pop_ht is not None:
            fields["pop"] = x.pop
        return hl.struct(**fields)

    vp = vp.annotate(
        _shared=vp._c1.filter(lambda x: vp._c2_map.contains(x.s)).map(_shared_struct)
    )
    counts = _trio_chet_counts(vp._shared)
    vp = vp.annotate(raw=counts.raw, adj=counts.adj)
    if pop_ht is not None:
        # Per-pop cis/trans: group shared probands by pop, count each, then add
        # GLOBAL_POP = all probands (incl. any with a missing pop label).
        # --pops restricts to a requested subset of groups.
        shared_for_pop = vp._shared.filter(lambda x: hl.is_defined(x.pop))
        if pops is not None:
            allowed = hl.literal(set(pops))
            shared_for_pop = shared_for_pop.filter(lambda x: allowed.contains(x.pop))
        by_pop = hl.group_by(lambda x: x.pop, shared_for_pop).map_values(
            _trio_chet_counts
        )
        vp = vp.annotate(
            raw_by_pop=hl.dict(
                hl.array(by_pop.map_values(lambda c: c.raw)).append(
                    (GLOBAL_POP, vp.raw)
                )
            ),
            adj_by_pop=hl.dict(
                hl.array(by_pop.map_values(lambda c: c.adj)).append(
                    (GLOBAL_POP, vp.adj)
                )
            ),
        )
    vp = vp.drop("_c1", "_c2", "_c2_map", "_shared")
    return vp.filter(vp.raw.n_same_hap + vp.raw.n_chet > 0)


def _chet_call(counts: hl.expr.StructExpression) -> hl.expr.BooleanExpression:
    """cis (False) / trans (True) / missing from ``n_same_hap`` & ``n_chet``."""
    return (
        hl.case()
        .when((counts.n_same_hap > 0) & (counts.n_chet == 0), False)
        .when((counts.n_same_hap == 0) & (counts.n_chet > 0), True)
        .or_missing()
    )


def subtract_pbt_from_gnomad_counts(
    vp_ht: hl.Table, gnomad_all: hl.Table, pbt_counts: hl.Table
) -> hl.Table:
    """Subtract the PBT∩release contribution from precomputed gnomAD counts.

    ``gnomad_minus_pbt[cell] = gnomad_all[cell] - pbt[cell]`` element-wise on the
    9-element ``gt_counts_{raw,adj}`` arrays. Exact for the 8 carrier cells
    (per-sample counts, additive over the disjoint release∖PBT and PBT∩release
    partition); the hom-ref/hom-ref (AABB) cell is additive too as long as
    ``pbt_counts`` is produced by the same counting function as ``gnomad_all``
    (same inclusion-exclusion form over the same encoded sets). A ``max(·, 0)``
    clamp guards against any off-by-callability negative.

    :param vp_ht: Pairs to emit, keyed by ``(locus1, alleles1, locus2, alleles2)``.
    :param gnomad_all: Precomputed release counts with ``gt_counts_raw/adj``.
    :param pbt_counts: Counts over PBT∩release for the same pairs.
    :return: ``vp_ht`` with PBT-subtracted ``gt_counts_raw`` / ``gt_counts_adj``.
    """
    zeros = hl.range(9).map(lambda _: 0)

    def _sub(a, b):
        b = hl.or_else(b, zeros)
        return hl.zip(a, b).map(lambda x: hl.max(x[0] - x[1], 0))

    g = gnomad_all[vp_ht.key]
    p = pbt_counts[vp_ht.key]
    fields = dict(
        gt_counts_raw=_sub(g.gt_counts_raw, p.gt_counts_raw),
        gt_counts_adj=_sub(g.gt_counts_adj, p.gt_counts_adj),
    )
    # Per-pop subtraction (same element-wise logic per pop), when both sides
    # carry the --stratify-by-pop breakdown.
    if "gt_counts_by_pop" in gnomad_all.row and "gt_counts_by_pop" in pbt_counts.row:
        zeros_struct = hl.struct(raw=zeros, adj=zeros)
        fields["gt_counts_by_pop"] = hl.dict(
            hl.array(g.gt_counts_by_pop).map(
                lambda kv: (
                    kv[0],
                    hl.bind(
                        lambda pv: hl.struct(
                            raw=_sub(kv[1].raw, pv.raw),
                            adj=_sub(kv[1].adj, pv.adj),
                        ),
                        hl.or_else(p.gt_counts_by_pop.get(kv[0]), zeros_struct),
                    ),
                )
            )
        )
    return vp_ht.select(**fields)


def build_trio_comparison(trio_ht: hl.Table, gnomad_ht: hl.Table) -> hl.Table:
    """Join trio truth with gnomAD-minus-PBT EM phase for each pair.

    :param trio_ht: Trio phase-count Table (from :func:`call_trio_chet`).
    :param gnomad_ht: EM-phased gnomAD-minus-PBT Table with ``em`` /
        ``em_plus_one`` / ``gt_counts_*`` fields.
    :return: ``trio_ht`` annotated with the gnomAD EM phase, ``distance`` and a
        ``trio_chet`` cis/trans call (``raw`` / ``adj``).
    """
    # trio_ht is keyed by vp_ht_idx (from create_variant_pair_ht); the gnomAD
    # counts are keyed by the variant-pair locus tuple. Re-key to join.
    trio_ht = trio_ht.key_by("locus1", "alleles1", "locus2", "alleles2")
    gnomad_fields = ["em", "em_plus_one", "gt_counts_raw", "gt_counts_adj"]
    # Carry the per-pop gnomAD EM / counts when present (--stratify-by-pop).
    gnomad_fields += [
        f for f in ("em_by_pop", "em_plus_one_by_pop", "gt_counts_by_pop")
        if f in gnomad_ht.row
    ]
    g = gnomad_ht.select(*gnomad_fields)
    ht = trio_ht.annotate(**g[trio_ht.key])
    ht = ht.annotate(
        distance=ht.locus2.position - ht.locus1.position,
        trio_chet=hl.struct(raw=_chet_call(ht.raw), adj=_chet_call(ht.adj)),
    )
    # Per-pop trio cis/trans call, keyed by the same pop labels as the trio
    # counts (incl. GLOBAL_POP). Paired with em_by_pop in the output so the
    # accuracy analysis can compare per pop.
    if "raw_by_pop" in ht.row:
        ht = ht.annotate(
            trio_chet_by_pop=hl.dict(
                hl.array(ht.raw_by_pop).map(
                    lambda kv: (
                        kv[0],
                        hl.struct(
                            raw=_chet_call(kv[1]),
                            adj=_chet_call(ht.adj_by_pop[kv[0]]),
                        ),
                    )
                )
            )
        )
    return ht


def main(args):
    """Trio phasing and trio-vs-gnomAD phase comparison."""
    start = timeit.default_timer()
    tmp_dir = args.tmp_dir
    overwrite = args.overwrite
    output_postfix = args.output_postfix
    data_type = args.data_type
    min_an_pct = args.min_an_pct
    trio_set = args.trio_set
    test_chrom = args.test_chrom
    if test_chrom and not test_chrom.startswith("chr"):
        test_chrom = f"chr{test_chrom}"
    test = args.test or bool(args.gene) or bool(test_chrom)
    test_intervals = (
        {args.gene: TEST_INTERVALS[args.gene]} if args.gene else TEST_INTERVALS
    )

    hl.init(
        log=os.path.join(tempfile.gettempdir(), "trio_phasing.log"),
        tmp_dir=tmp_dir,
    )

    logger.info(
        f"""
        Running script with the following parameters:

            Data type: {data_type}
            Test: {test}
            Gene: {args.gene or 'all test intervals'}
            Test chrom: {test_chrom or 'n/a'}
            Output postfix: {output_postfix}
            Trio set: {trio_set}
            Min AN pct: {min_an_pct}
            Overwrite: {overwrite}
            Tmp dir: {tmp_dir}
        """
    )

    if data_type != "exomes":
        logger.warning(
            "The gnomad_qc v4 finalized pedigree is exomes-derived; running "
            "with data_type=%s pairs the %s VDS with exome trio sample IDs.",
            data_type,
            data_type,
        )

    resources = get_trio_phasing_resources(
        data_type=data_type,
        test=test,
        tmp_dir=tmp_dir if test else None,
        output_postfix=output_postfix,
        overwrite=overwrite,
        trio_set=trio_set,
    )

    # The pedigree resource for the chosen trio set: "pedigree" = all trios
    # (multiple per family), "trios" = one trio per family.
    ped_resource = pedigree(finalized=True) if trio_set == "pedigree" else trios()

    get_vds_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )

    # Per-pop stratification (shared across the trio-truth and gnomAD-count
    # steps): --pops restricts to a subset of groups and implies
    # --stratify-by-pop; the meta pop HT is looked up lazily where needed.
    requested_pops = None
    if args.pops:
        requested_pops = [p.strip() for p in args.pops.split(",") if p.strip()]
        valid = set(get_pops(data_type))
        unknown = [p for p in requested_pops if p not in valid]
        if unknown:
            raise ValueError(
                f"--pops has unknown group(s) {unknown}; valid groups for "
                f"{data_type}: {sorted(valid)}"
            )
        requested_pops = [p for p in requested_pops if p != GLOBAL_POP]
    stratify_by_pop = args.stratify_by_pop or requested_pops is not None
    pop_ht = get_sample_pop_ht(data_type) if stratify_by_pop else None

    if test_chrom:
        filter_intervals = [
            hl.parse_locus_interval(test_chrom, reference_genome="GRCh38")
        ]
    elif test:
        filter_intervals = list(test_intervals.values())
    else:
        filter_intervals = None
    count_output_dir = (
        f"{tmp_dir}/trio_genotype_count_intermediates"
        f"{_get_output_postfix(output_postfix, test)}.{trio_set}"
    )

    if args.create_pbt_trio_matrix:
        logger.info("Creating PBT-phased trio MatrixTable...")
        res = resources.create_pbt_trio_matrix
        res.check_resource_existence()

        ped = ped_resource.pedigree()
        trio_samples = _complete_trio_samples(ped)
        logger.info(
            "Pedigree has %d complete trios (%d unique samples).",
            len(ped.complete_trios()),
            len(trio_samples),
        )

        vds = get_vds_func(
            high_quality_only=True,
            split=True,
            filter_intervals=filter_intervals,
            entries_to_keep=["GT", "GQ", "DP", "AD"],
            split_reference_blocks=False,
        )
        vds = hl.vds.filter_samples(vds, _samples_ht(trio_samples), keep=True)
        mt = hl.vds.to_dense_mt(vds)

        tm = phase_trio_matrix(mt, ped)
        tm = tm.checkpoint(res.pbt_trio_matrix.path, overwrite=overwrite)
        logger.info(
            "The PBT trio MatrixTable has been written to %s. "
            "Number of trios: %d, number of variants: %d",
            res.pbt_trio_matrix.path,
            tm.count_cols(),
            tm.count_rows(),
        )

    if args.explode_pbt:
        logger.info("Exploding the PBT trio MatrixTable into a per-sample MT...")
        res = resources.explode_pbt
        res.check_resource_existence()

        tm = res.pbt_trio_matrix.mt()
        pmt = explode_pbt_trio_matrix(tm)
        pmt = pmt.checkpoint(res.pbt_mt.path, overwrite=overwrite)
        logger.info(
            "The exploded PBT MatrixTable has been written to %s. "
            "Number of sample-trio columns: %d",
            res.pbt_mt.path,
            pmt.count_cols(),
        )

    if args.phase_multi_families:
        logger.info("Computing consensus phase for multi-offspring families...")
        res = resources.phase_multi_families
        res.check_resource_existence()

        pmt = res.pbt_mt.mt()
        pmt = phase_multi_offspring_families(pmt)
        pmt = pmt.checkpoint(res.pbt_multi_families_mt.path, overwrite=overwrite)
        logger.info(
            "The multi-offspring-family consensus MatrixTable has been written "
            "to %s. Number of samples: %d",
            res.pbt_multi_families_mt.path,
            pmt.count_cols(),
        )

    if args.derive_trio_vps:
        logger.info("Deriving trio variant-pair list from probands...")
        res = resources.derive_trio_vps
        if args.variant_filter_path:
            logger.info(
                "Using override variant filter HT: %s", args.variant_filter_path
            )
            filter_ht = hl.read_table(args.variant_filter_path)
        else:
            res.check_resource_existence()
            filter_ht = res.variant_filter_ht.ht()

        pmt = res.pbt_mt.mt()
        # Probands are the trio-id columns of the chosen trio set.
        pmt = pmt.filter_cols(pmt.s == pmt.source_trio.id)
        ht = create_variant_pair_ht(
            pmt, filter_ht, drop_oe_only_pairs=args.include_in_trans_oe_candidates
        )
        ht = ht.checkpoint(res.trio_vp_list_ht.path, overwrite=overwrite)
        logger.info(
            "The trio variant-pair list has been written to %s. "
            "Number of pairs: %d",
            res.trio_vp_list_ht.path,
            ht.count(),
        )

    if args.call_trio_chet:
        logger.info("Calling trio cis/trans (chet) per pair...")
        res = resources.call_trio_chet
        res.check_resource_existence()

        pmt = res.pbt_mt.mt()
        # Probands are the trio-id columns of the chosen trio set.
        pmt = pmt.filter_cols(pmt.s == pmt.source_trio.id)
        vp_ht = filter_pairs_by_an_pct(res.trio_vp_list_ht.ht(), min_an_pct)
        ht = call_trio_chet(pmt, vp_ht, pop_ht=pop_ht, pops=requested_pops)
        ht = ht.checkpoint(res.trio_phase_counts_ht.path, overwrite=overwrite)
        logger.info(
            "The trio phase-count Table has been written to %s. "
            "Number of pairs with a phased-het proband: %d",
            res.trio_phase_counts_ht.path,
            ht.count(),
        )

    if args.gnomad_counts_no_pbt:
        res = resources.gnomad_counts_no_pbt
        res.check_resource_existence()
        vp_ht = filter_pairs_by_an_pct(res.trio_vp_list_ht.ht(), min_an_pct)
        trio_samples = _complete_trio_samples(ped_resource.pedigree())

        if args.gnomad_counts_path:
            # Reuse precomputed release counts (join) and subtract the
            # PBT∩release contribution, computed by REUSING the already-densified
            # PBT MT — no new densify. The PBT MT carries GT + a precomputed adj
            # for every trio member, so the subtraction is exact (same encode +
            # count path, AABB anchor n_samples = |PBT∩release| picked up
            # automatically).
            logger.info(
                "Reusing precomputed gnomAD counts (%s); subtracting the "
                "PBT∩release contribution from the densified PBT MT...",
                args.gnomad_counts_path,
            )
            if not args.gnomad_encoded_path:
                raise ValueError(
                    "--gnomad-encoded-path (the encoded_gt_sets HT for the "
                    "precomputed counts) is required with --gnomad-counts-path: "
                    "its `samples` global is the exact release cohort to subtract."
                )
            gnomad_all = hl.read_table(args.gnomad_counts_path)
            vp_ht = vp_ht.key_by("locus1", "alleles1", "locus2", "alleles2")
            covered = vp_ht.semi_join(gnomad_all).checkpoint(
                hl.utils.new_temp_file("trio_covered_vps", "ht")
            )
            logger.info(
                "Trio pairs with a precomputed gnomAD count: %d", covered.count()
            )

            # Exact gnomad_all cohort = the encoded-sets `samples` global; the
            # PBT∩release samples are the trio members in that cohort.
            release_set = {
                r.s
                for r in hl.eval(
                    hl.read_table(args.gnomad_encoded_path).index_globals().samples
                )
            }
            pbt_in_release = sorted(set(trio_samples) & release_set)
            logger.info(
                "PBT∩release samples to subtract: %d of %d trio members",
                len(pbt_in_release),
                len(trio_samples),
            )

            # Restrict the exploded PBT MT to those samples (dedup multi-trio
            # columns to one per sample) and covered-pair variants, then encode
            # with the precomputed adj and count.
            release_lit = hl.literal(set(pbt_in_release))
            pbt_mt = res.pbt_mt.mt()
            pbt_mt = pbt_mt.filter_cols(release_lit.contains(pbt_mt.s))
            pbt_mt = pbt_mt.group_cols_by(pbt_mt.s).aggregate(
                GT=hl.agg.take(pbt_mt.GT, 1)[0],
                adj=hl.agg.take(pbt_mt.adj, 1)[0],
            )
            pbt_mt = pbt_mt.semi_join_rows(create_variant_pair_filter_ht(covered))
            pbt_dir = f"{count_output_dir}/pbt"
            pbt_mt = pbt_mt.annotate_globals(min_an_pct=min_an_pct).checkpoint(
                f"{pbt_dir}/pbt.mt", overwrite=overwrite
            )
            encode_genotypes(
                pbt_mt,
                covered,
                output_dir=pbt_dir,
                min_an_pct=min_an_pct,
                use_precomputed_adj=True,
            )
            # count_all_pairs_via_index rebuilds locus1/… via select, so pass
            # the pairs with those as non-key fields (covered stays locus-keyed
            # for the join/subtract below). With --stratify-by-pop the
            # subtraction below needs --gnomad-counts-path to also carry a
            # per-pop breakdown (same pop labels).
            pbt_var_idx = hl.read_table(f"{pbt_dir}/var_idx.ht")
            pbt_encoded = hl.read_table(f"{pbt_dir}/encoded_gt_sets_by_var_idx.ht")
            if pop_ht is not None:
                pbt_counts = compute_counts_by_pop(
                    covered.key_by(), pbt_var_idx, pbt_encoded, pop_ht, requested_pops,
                )
            else:
                pbt_counts = count_all_pairs_via_index(
                    covered.key_by(), pbt_var_idx, pbt_encoded,
                )
            ht = subtract_pbt_from_gnomad_counts(covered, gnomad_all, pbt_counts)
        else:
            # Count from scratch over release-minus-PBT (no precomputed counts).
            logger.info(
                "Computing gnomAD genotype counts (release minus PBT members) "
                "from scratch..."
            )
            pbt_members = _samples_ht(trio_samples)
            filter_ht = create_variant_pair_filter_ht(vp_ht)
            vds = get_vds_func(
                release_only=True,
                split=True,
                filter_intervals=filter_intervals,
                filter_variant_ht=filter_ht,
                entries_to_keep=["GT", "GQ", "DP", "AD"],
                split_reference_blocks=False,
            )
            vds = hl.vds.filter_samples(vds, pbt_members, keep=False)
            mt = hl.vds.to_dense_mt(vds).annotate_globals(min_an_pct=min_an_pct)
            mt = mt.checkpoint(f"{count_output_dir}/dense.mt", overwrite=overwrite)
            encode_genotypes(
                mt, vp_ht, output_dir=count_output_dir, min_an_pct=min_an_pct,
            )
            scratch_var_idx = hl.read_table(f"{count_output_dir}/var_idx.ht")
            scratch_encoded = hl.read_table(
                f"{count_output_dir}/encoded_gt_sets_by_var_idx.ht"
            )
            if pop_ht is not None:
                ht = compute_counts_by_pop(
                    vp_ht, scratch_var_idx, scratch_encoded, pop_ht, requested_pops,
                )
            else:
                ht = count_all_pairs_via_index(
                    vp_ht, scratch_var_idx, scratch_encoded,
                )

        ht = ht.checkpoint(res.gnomad_no_pbt_counts_ht.path, overwrite=overwrite)
        logger.info(
            "The gnomAD-minus-PBT genotype counts have been written to %s. "
            "Number of pairs: %d",
            res.gnomad_no_pbt_counts_ht.path,
            ht.count(),
        )

    if args.phase_gnomad_counts:
        logger.info("EM-phasing the gnomAD-minus-PBT counts...")
        res = resources.phase_gnomad_counts
        res.check_resource_existence()

        ht = res.gnomad_no_pbt_counts_ht.ht()
        ht = ht.annotate(**get_phased_gnomad_ht(ht))
        ht = ht.checkpoint(res.gnomad_no_pbt_phased_ht.path, overwrite=overwrite)
        logger.info(
            "The EM-phased gnomAD-minus-PBT Table has been written to %s.",
            res.gnomad_no_pbt_phased_ht.path,
        )

    if args.export_comparison:
        logger.info("Joining trio truth vs gnomAD phase and exporting...")
        res = resources.export_comparison
        res.check_resource_existence()

        ht = build_trio_comparison(
            res.trio_phase_counts_ht.ht(), res.gnomad_no_pbt_phased_ht.ht()
        )
        ht = ht.filter(ht.locus1.in_autosome())
        ht = ht.checkpoint(res.trio_comparison_ht.path, overwrite=overwrite)
        tsv_path = res.trio_comparison_ht.path[: -len(".ht")] + ".tsv"
        ht.flatten().export(tsv_path)
        logger.info(
            "The trio comparison has been written to %s (and %s). "
            "Number of pairs: %d",
            res.trio_comparison_ht.path,
            tsv_path,
            ht.count(),
        )

    stop = timeit.default_timer()
    logger.info(f"Time taken to run the script is {stop - start} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tmp-dir",
        default=DEFAULT_TMP_DIR,
        help="Temporary directory for intermediate files.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to test intervals (all genes in TEST_INTERVALS) for testing.",
    )
    parser.add_argument(
        "--gene",
        choices=list(TEST_INTERVALS),
        help=(
            "Run on a single gene; uses that gene's interval from TEST_INTERVALS and "
            "implies --test."
        ),
    )
    parser.add_argument(
        "--test-chrom",
        help=(
            "Scope the run to a single chromosome (e.g. '19' or 'chr19'). Implies "
            "--test (outputs go to the tmp dir)."
        ),
    )
    parser.add_argument(
        "--output-postfix",
        help=(
            'Postfix to append to output file names (e.g., "sgca_trio" for files like '
            "exomes.pbt.sgca_trio.mt)."
        ),
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Whether to overwrite existing files."
    )
    parser.add_argument(
        "--data-type",
        default=DEFAULT_DATA_TYPE,
        choices=DATA_TYPE_CHOICES,
        help=(
            f'Data type to use. Must be one of {", ".join(DATA_TYPE_CHOICES)}. Default '
            f"is {DEFAULT_DATA_TYPE}. The finalized pedigree is exomes-derived."
        ),
    )
    parser.add_argument(
        "--min-an-pct",
        type=int,
        default=0,
        help=(
            "AN_percent floor applied to the trio variant pairs at count time "
            "(must match across --gnomad-counts-no-pbt and --call-trio-chet)."
        ),
    )
    parser.add_argument(
        "--trio-set",
        default="pedigree",
        choices=("pedigree", "trios"),
        help=(
            "Which gnomad_qc v4 pedigree resource to use: 'pedigree' (all trios, "
            "multiple per family; n=15061, matches v2) or 'trios' (one trio per "
            "family; n=12731). Folded into every trio-phasing output path so the "
            "two sets never overwrite each other. Default is 'pedigree'."
        ),
    )
    parser.add_argument(
        "--create-pbt-trio-matrix",
        action="store_true",
        help="Build the PBT-phased trio MatrixTable from the gnomAD VDS + pedigree.",
    )
    parser.add_argument(
        "--explode-pbt",
        action="store_true",
        help="Explode the PBT trio MatrixTable into a per-sample MatrixTable.",
    )
    parser.add_argument(
        "--phase-multi-families",
        action="store_true",
        help="Compute consensus phase from PBT in families with multiple offspring.",
    )
    parser.add_argument(
        "--derive-trio-vps",
        action="store_true",
        help="Derive the trio variant-pair list from proband genotypes.",
    )
    parser.add_argument(
        "--variant-filter-path",
        help=(
            "Explicit path to a variant filter HT for --derive-trio-vps, overriding "
            "the postfix-derived resource path. Use to reuse the production "
            "genome-wide gs://gnomad/v4.1/variant_cooccurrence/exomes.variant_filter.ht."
        ),
    )
    parser.add_argument(
        "--include-in-trans-oe-candidates",
        action="store_true",
        help=(
            "In --derive-trio-vps, drop pairs where both sides are OE-candidate-only "
            "(drop_oe_only_pairs=True), matching the chr19 variant-pair-list run. "
            "Requires the variant filter HT to carry a 'source' field."
        ),
    )
    parser.add_argument(
        "--call-trio-chet",
        action="store_true",
        help="Count cis/trans probands per pair from transmission phase (trio truth).",
    )
    parser.add_argument(
        "--stratify-by-pop",
        action="store_true",
        help=(
            "Stratify the trio-vs-gnomAD comparison by genetic-ancestry group. "
            "In --call-trio-chet, adds raw_by_pop/adj_by_pop keyed by proband "
            "pop (meta.population_inference.pop). In --gnomad-counts-no-pbt, "
            "emits per-pop gnomAD counts (requires --gnomad-counts-path to also "
            "carry a gt_counts_by_pop breakdown when reusing precomputed counts). "
            "--phase-gnomad-counts and --export-comparison then pick up the "
            "per-pop EM / trio_chet automatically. Pass consistently across the "
            "steps. Use --pops to restrict to a subset of groups."
        ),
    )
    parser.add_argument(
        "--pops",
        help=(
            "Comma-separated genetic-ancestry groups to stratify by (e.g. "
            "'nfe,afr,eas'). Implies --stratify-by-pop and restricts the per-pop "
            "breakdown (trio counts + gnomAD counts) to these groups; 'all' is "
            "always included. Requested groups absent from the cohort are "
            "skipped. Default (with --stratify-by-pop, no --pops): all groups."
        ),
    )
    parser.add_argument(
        "--gnomad-counts-no-pbt",
        action="store_true",
        help="Compute gnomAD genotype counts on the trio VPs, release minus PBT members.",
    )
    parser.add_argument(
        "--gnomad-counts-path",
        help=(
            "Path to a precomputed gnomAD release counts HT (gt_counts_raw/adj). "
            "When set, --gnomad-counts-no-pbt reuses it by join and subtracts only "
            "the PBT∩release contribution (much cheaper than recounting all of "
            "release); the comparison is restricted to the covered pairs. Omit to "
            "count from scratch over release-minus-PBT."
        ),
    )
    parser.add_argument(
        "--gnomad-encoded-path",
        help=(
            "Path to the encoded_gt_sets HT that produced --gnomad-counts-path "
            "(required with it). Its `samples` global is the exact release cohort; "
            "the PBT∩release subset is subtracted by reusing the already-densified "
            "PBT MT (no new densify)."
        ),
    )
    parser.add_argument(
        "--phase-gnomad-counts",
        action="store_true",
        help="EM-phase the gnomAD-minus-PBT genotype counts.",
    )
    parser.add_argument(
        "--export-comparison",
        action="store_true",
        help="Join trio truth vs gnomAD EM phase and export an HT + TSV.",
    )

    args = parser.parse_args()
    main(args)

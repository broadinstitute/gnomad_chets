import hail as hl
from resources import phased_vp_count_ht_path
from gnomad.utils.liftover import get_liftover_genome
from gnomad.utils.slack import try_slack
import gnomad.resources.grch37.gnomad as gnomad
from phasing import get_em_expr, flatten_gt_counts
import argparse
from math import ceil
import logging

logger = logging.getLogger("compute_phase")
logger.setLevel(logging.INFO)


def liftover_expr(
        locus: hl.expr.LocusExpression,
        alleles: hl.expr.ArrayExpression,
        destination_ref: hl.ReferenceGenome
) -> hl.expr.StructExpression:
    lifted_over_locus = hl.liftover(locus, destination_ref, include_strand=True)
    lifted_over_alleles = alleles.map(
        lambda a: hl.if_else(lifted_over_locus.is_negative_strand, hl.reverse_complement(a), a)
    )
    return hl.struct(
        locus=lifted_over_locus.result,
        alleles=lifted_over_alleles
    )


def load_cmg(cmg_csv: str) -> hl.Table:
    cmg_ht = hl.import_table(cmg_csv, impute=True, delimiter=",", quote='"')

    cmg_ht = cmg_ht.transmute(
        locus1_b38=hl.locus("chr" + hl.str(cmg_ht.chrom_1), cmg_ht.pos_1, reference_genome='GRCh38'),
        alleles1_b38=[cmg_ht.ref_1, cmg_ht.alt_1],
        locus2_b38=hl.locus("chr" + hl.str(cmg_ht.chrom_2), cmg_ht.pos_2, reference_genome='GRCh38'),
        alleles2_b38=[cmg_ht.ref_2, cmg_ht.alt_2]
    )

    liftover_references = get_liftover_genome(cmg_ht.rename({'locus1_b38': 'locus'}))
    lifted_over_variants = hl.sorted(
        hl.array([
            liftover_expr(cmg_ht.locus1_b38, cmg_ht.alleles1_b38, liftover_references[1]),
            liftover_expr(cmg_ht.locus2_b38, cmg_ht.alleles2_b38, liftover_references[1])
        ]),
        lambda x: x.locus
    )

    cmg_ht = cmg_ht.key_by(
        locus1=lifted_over_variants[0].locus,
        alleles1=lifted_over_variants[0].alleles,
        locus2=lifted_over_variants[1].locus,
        alleles2=lifted_over_variants[1].alleles
    )

    return cmg_ht.annotate(
        bad_liftover=(
                hl.is_missing(cmg_ht.locus1) |
                hl.is_missing(cmg_ht.locus2) |
                (cmg_ht.locus1.sequence_context() != cmg_ht.alleles1[0][0]) |
                (cmg_ht.locus2.sequence_context() != cmg_ht.alleles2[0][0])
        )
    )


def get_sorted_variants_expr(
        locus1: hl.expr.LocusExpression,
        alleles1: hl.expr.ArrayExpression,
        locus2: hl.expr.LocusExpression,
        alleles2: hl.expr.ArrayExpression
) -> hl.expr.StructExpression:
    if locus1.dtype.reference_genome.name == 'GRCh37':
        variants = [
            hl.struct(
                locus=locus1,
                alleles=alleles1
            ),
            hl.struct(
                locus=locus2,
                alleles=alleles2
            )
        ]
    else:
        logger.warning("Variants are not on GRCh37; they will be lifted over.")
        _, destination_ref = get_liftover_genome(hl.struct(locus=locus1))
        variants = [
            liftover_expr(
                locus=locus1,
                alleles=alleles1,
                destination_ref=destination_ref
            ),
            liftover_expr(
                locus=locus2,
                alleles=alleles2,
                destination_ref=destination_ref
            )
        ]

    sorted_variants = hl.sorted(variants, key=lambda x: x.locus)
    return hl.struct(
        locus1=sorted_variants[0].locus,
        alleles1=sorted_variants[0].alleles,
        locus2=sorted_variants[1].locus,
        alleles2=sorted_variants[1].alleles
    )


def read_variants_ht(path: str) -> hl.Table:
    variants_ht = hl.read_table(path)

    # Make sure that types match
    assert (
            isinstance(variants_ht.key[0], hl.expr.LocusExpression) &
            (variants_ht.key[1].dtype == hl.tarray(hl.tstr)) &
            isinstance(variants_ht.key[2], hl.expr.LocusExpression) &
            (variants_ht.key[3].dtype == hl.tarray(hl.tstr))
    )

    variants_ht = variants_ht.key_by(
        **get_sorted_variants_expr(
            variants_ht.key[0],
            variants_ht.key[1],
            variants_ht.key[2],
            variants_ht.key[3]
        )
    ).persist()

    return variants_ht


def variants_ht_from_text(variants: str) -> hl.Table:
    v1, v2 = variants.split(",")
    chr1, pos1, ref1, alt1 = v1.strip().split(":")
    chr2, pos2, ref2, alt2 = v2.strip().split(":")

    ref = "GRCh38" if chr1.startswith("chr") else "GRCh37"

    variants_ht = hl.Table.parallelize(
        rows=[
            get_sorted_variants_expr(
                hl.locus(chr1, int(pos1), reference_genome=ref),
                hl.array([ref1, alt1]),
                hl.locus(chr2, int(pos2), reference_genome=ref),
                hl.array([ref2, alt2])
            )
        ],
        key=['locus1', 'alleles1', 'locus2', 'alleles2']
    ).persist()

    logger.info("Found the following variant pair to phase:")
    variants_ht.show()

    return variants_ht


# TODO: How to handle one variant absent from gnomAD?
def annotate_unphased_pairs(unphased_ht: hl.Table, n_variant_pairs: int):
    # unphased_ht = vp_ht.filter(hl.is_missing(vp_ht.all_phase))
    # unphased_ht = unphased_ht.key_by()

    # Explode variant pairs
    unphased_ht = unphased_ht.annotate(
        las=[
            hl.tuple([unphased_ht.locus1, unphased_ht.alleles1]),
            hl.tuple([unphased_ht.locus2, unphased_ht.alleles2])
        ]
    ).explode('las', name='la')

    unphased_ht = unphased_ht.key_by(
        locus=unphased_ht.la[0],
        alleles=unphased_ht.la[1]
    ).persist()  # .checkpoint('gs://gnomad-tmp/vp_ht_unphased.ht')

    # Annotate single variants with gnomAD freq
    gnomad_ht = gnomad.public_release('exomes').ht()
    gnomad_ht = gnomad_ht.semi_join(unphased_ht).repartition(
        ceil(n_variant_pairs / 10000),
        shuffle=True
    )

    missing_freq = hl.struct(
        AC=0,
        AF=0,
        AN=125748 * 2,  # set to no missing for now
        homozygote_count=0
    )

    logger.info(f"{gnomad_ht.count()}/{unphased_ht.count()} single variants from the unphased pairs found in gnomAD.")

    gnomad_freq = gnomad_ht[unphased_ht.key].freq
    unphased_ht = unphased_ht.annotate(
        adj_freq=hl.or_else(
            gnomad_freq[0],
            missing_freq
        ),
        raw_freq=hl.or_else(
            gnomad_freq[1],
            missing_freq
        ),
        # pop_max_freq=hl.or_else(
        #     gnomad_exomes.popmax[0],
        #     missing_freq.annotate(
        #         pop=hl.null(hl.tstr)
        #     )
        # )
    )
    unphased_ht = unphased_ht.persist()
    # unphased_ht = unphased_ht.checkpoint('gs://gnomad-tmp/unphased_ann.ht', overwrite=True)

    loci_expr = hl.sorted(
        hl.agg.collect(
            hl.tuple([
                unphased_ht.locus,
                hl.struct(
                    adj_freq=unphased_ht.adj_freq,
                    raw_freq=unphased_ht.raw_freq,
                    # pop_max_freq=unphased_ht.pop_max_freq
                )
            ])
        ),
        lambda x: x[0]  # sort by locus
    ).map(
        lambda x: x[1]  # get rid of locus
    )

    vp_freq_expr = hl.struct(
        v1=loci_expr[0],
        v2=loci_expr[1]
    )

    # [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]
    def get_gt_counts(freq: str):
        return hl.array([
            hl.min(vp_freq_expr.v1[freq].AN, vp_freq_expr.v2[freq].AN),  # AABB
            vp_freq_expr.v2[freq].AC - (2 * vp_freq_expr.v2[freq].homozygote_count),  # AABb
            vp_freq_expr.v2[freq].homozygote_count,  # AAbb
            vp_freq_expr.v1[freq].AC - (2 * vp_freq_expr.v1[freq].homozygote_count),  # AaBB
            0,  # AaBb
            0,  # Aabb
            vp_freq_expr.v1[freq].homozygote_count,  # aaBB
            0,  # aaBb
            0  # aabb
        ])

    gt_counts_raw_expr = get_gt_counts('raw_freq')
    gt_counts_adj_expr = get_gt_counts('adj_freq')

    # gt_counts_pop_max_expr = get_gt_counts('pop_max_freq')
    return unphased_ht.group_by(
        unphased_ht.locus1,
        unphased_ht.alleles1,
        unphased_ht.locus2,
        unphased_ht.alleles2
    ).aggregate(
        pop='all',  # TODO Add option for multiple pops?
        phase_info=hl.struct(
            gt_counts=hl.struct(
                raw=gt_counts_raw_expr,
                adj=gt_counts_adj_expr
            ),
            em=hl.struct(
                raw=get_em_expr(gt_counts_raw_expr),
                adj=get_em_expr(gt_counts_raw_expr)
            )
        )
        # pop_max_gt_counts_adj=gt_counts_raw_expr,
        # pop_max_em_p_chet_adj=get_em_expr(gt_counts_raw_expr).p_chet,
    )  # .key_by()


def flatten_phased_ht(phased_ht: hl.Table) -> hl.Table:
    phased_ht = phased_ht.key_by()

    # phase_ht = phase_ht.key_by()

    def flatten_phase_dict(
            expr: hl.expr.StructExpression
    ) -> hl.expr.StructExpression:
        return hl.struct(
            raw=flatten_gt_counts(expr.gt_counts.raw),
            adj=flatten_gt_counts(expr.gt_counts.adj),
            em_p_chet_raw=expr.em.raw.p_chet,
            em_p_chet_adj=expr.em.adj.p_chet,
            # em1_p_chet_raw=expr.em_plus_one.raw.p_chet,
            # em1_p_chet_adj=expr.em_plus_one.adj.p_chet
        )

    return phased_ht.transmute(
        chrom=phased_ht.locus1.contig,
        pos1=phased_ht.locus1.position,
        ref1=phased_ht.alleles1[0],
        alt1=phased_ht.alleles1[1],
        pos2=phased_ht.locus2.position,
        ref2=phased_ht.alleles2[0],
        alt2=phased_ht.alleles2[1],
        **{
            k: v for k, v in flatten_phase_dict(phased_ht.phase_info).items()
        }
    ).flatten()


def explode_phase_info(ht: hl.Table, remove_all_ref: bool = True) -> hl.Table:
    ht = ht.transmute(phase_info=hl.array(ht.phase_info))
    ht = ht.explode('phase_info')
    ht = ht.transmute(
        pop=ht.phase_info[0],
        phase_info=ht.phase_info[1]
    )

    if remove_all_ref:
        ht = ht.filter(hl.sum(ht.phase_info.gt_counts.raw[1:]) > 0)

    return ht


def compute_phase(variants_ht: hl.Table) -> hl.Table:
    n_variant_pairs = variants_ht.count()
    logger.info(f"Looking up phase for {n_variant_pairs} variant pair(s).")

    # Join with gnomad phased variants
    vp_ht = hl.read_table(phased_vp_count_ht_path('exomes'))
    phased_ht = vp_ht.semi_join(variants_ht)
    n_phased = phased_ht.count()
    phased_ht = explode_phase_info(phased_ht)  # explodes phase_info by pop
    phased_ht = phased_ht.transmute(
        phase_info=phased_ht.phase_info.select('gt_counts', 'em')
    ).repartition(ceil(n_variant_pairs / 10000), shuffle=True)
    phased_ht = phased_ht.persist()  # .checkpoint("gs://gnomad-tmp/vp_ht.ht")

    # If not all pairs had at least one carrier of both, then compute phase estimate from single variants
    logger.info(f"{n_phased}/{n_variant_pairs} variant pair(s) found with carriers of both in gnomAD.")

    if n_phased < n_variant_pairs:
        unphased_ht = variants_ht.anti_join(vp_ht)
        unphased_ht = annotate_unphased_pairs(unphased_ht, n_variant_pairs)
        phased_ht = phased_ht.union(
            unphased_ht,
            unify=True
        )

    return phased_ht


def main(args):
    # Load data
    variants_ht = read_variants_ht(args.ht) if args.ht else variants_ht_from_text(args.variants)

    # Add phase
    phased_ht = compute_phase(variants_ht)

    # Write results
    if args.out.endswith(".ht"):
        phased_ht.write(args.out, overwrite=args.overwrite)
    else:
        phased_ht = flatten_phased_ht(phased_ht)
        phased_ht.export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--ht', help='HT containing variants. Needs to be keyed by locus1, alleles1, locus2, alleles2.')
    data_grp.add_argument('--variants', help='Variants to phase in format chr1:pos1:ref1:alt1,chr2:pos2:ref2:alt2. Note that chromosome needs to start with "chr" for GRCh38 variants')
    parser.add_argument('--out', help="Output file path. Output file format depends on extension (.ht, .csv or .csv.gz)")
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

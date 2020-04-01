import hail as hl
from resources import *
from gnomad.utils.liftover import lift_data, get_liftover_genome
import gnomad.resources.grch37.gnomad as gnomad
from typing import  List
from phasing import *
import argparse

CMG_CSV = "gs://gnomad/projects/compound_hets/Feb_2020_CMG_Compound_het_list.csv"


def liftover_expr(
        locus: hl.expr.LocusExpression,
        alleles: hl.expr.ArrayExpression,
        destination_ref: hl.ReferenceGenome
) -> hl.expr.StructExpression:

    lifted_over_locus=hl.liftover(locus, destination_ref, include_strand=True)
    lifted_over_alleles=alleles.map(
        lambda a: hl.if_else(lifted_over_locus.is_negative_strand,  hl.reverse_complement(a), a)
    )
    return hl.struct(
        locus=lifted_over_locus.result,
        alleles=lifted_over_alleles
    )


def load_cmg(cmg_csv: str) -> hl.Table:
    cmg_ht = hl.import_table(cmg_csv, impute=True,  delimiter=",", quote='"')

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

    cmg_ht =  cmg_ht.key_by(
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

def main(args):
    # Load CMG data
    cmg_ht = load_cmg(CMG_CSV)
    # print(f"{cmg_ht.aggregate(hl.agg.count_where(~cmg_ht.bad_liftover))} sites lifted over successfully.")
    vp_ht = hl.read_table(phased_vp_count_ht_path('exomes'))
    vp_ht = cmg_ht.join(vp_ht, how="left")
    vp_ht = vp_ht.transmute(
        all_phase=vp_ht.phase_info['all'],
        pop_max=hl.sorted(
            hl.array(vp_ht.phase_info),
            key=lambda x: -x[1].em.adj.p_chet
        )[0]
    ).repartition(2, shuffle=True).checkpoint("gs://gnomad-tmp/vp_ht.ht")

    unphased_ht = vp_ht.filter(hl.is_missing(vp_ht.all_phase))
    unphased_ht = unphased_ht.key_by()
    unphased_ht = unphased_ht.annotate(
        las=[
            hl.tuple([unphased_ht.locus1, unphased_ht.alleles1]),
            hl.tuple([unphased_ht.locus2, unphased_ht.alleles2])
        ]
    ).explode('las', name='la')

    unphased_ht = unphased_ht.key_by(locus=unphased_ht.la[0], alleles=unphased_ht.la[1]).checkpoint('gs://gnomad-tmp/vp_ht_unphased.ht')
    loci=unphased_ht.locus.collect(_localize=False)
    gnomad_exomes = gnomad.public_release('exomes').ht()
    gnomad_exomes =  gnomad_exomes.filter(
        loci.contains(gnomad_exomes.locus)
    ).repartition(2, shuffle=True)[unphased_ht.key]

    missing_freq=hl.struct(
        AC=0,
        AF=0,
        AN=125748*2, # set to no missing for now
        homozygote_count=0
    )

    unphased_ht = unphased_ht.annotate(
        adj_freq=hl.or_else(
            gnomad_exomes.freq[0],
            missing_freq
        ),
        raw_freq=hl.or_else(
            gnomad_exomes.freq[0],
            missing_freq
        ),
        # pop_max_freq=hl.or_else(
        #     gnomad_exomes.popmax[0],
        #     missing_freq.annotate(
        #         pop=hl.null(hl.tstr)
        #     )
        # )
    )
    unphased_ht = unphased_ht.checkpoint('gs://gnomad-tmp/unphased_ann.ht', overwrite=True)

    unphased_ht = hl.read_table('gs://gnomad-tmp/unphased_ann.ht')
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
        lambda x: x[0] # sort by locus
    ).map(
        lambda x: x[1] # get rid of locus
    )

    vp_freq_expr = hl.struct(
        v1=loci_expr[0],
        v2=loci_expr[1]
    )

    # [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]
    def get_gt_counts(freq: str):
        return hl.array([
            hl.min(vp_freq_expr.v1[freq].AN, vp_freq_expr.v2[freq].AN), # AABB
            vp_freq_expr.v2[freq].AC -  (2*vp_freq_expr.v2[freq].homozygote_count), # AABb
            vp_freq_expr.v2[freq].homozygote_count, # AAbb
        vp_freq_expr.v1[freq].AC -  (2*vp_freq_expr.v1[freq].homozygote_count), #AaBB
            0, # AaBb
            0, # Aabb
            vp_freq_expr.v1[freq].homozygote_count, # aaBB
            0, # aaBb
            0 # aabb
        ])

    gt_counts_raw_expr = get_gt_counts('raw_freq')
    gt_counts_adj_expr = get_gt_counts('adj_freq')
    # gt_counts_pop_max_expr = get_gt_counts('pop_max_freq')

    # [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]
    def flatten_gt_counts(gt_counts: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
        return hl.struct(
            ref_ref=gt_counts[0],
            ref_het=gt_counts[1],
            ref_hom=gt_counts[2],
            het_ref=gt_counts[3],
            het_het=gt_counts[4],
            het_hom=gt_counts[5],
            hom_ref=gt_counts[6],
            hom_het=gt_counts[7],
            hom_hom=gt_counts[8]
        )

    unphased_ht = unphased_ht.group_by(
        unphased_ht.locus1,
        unphased_ht.alleles1,
        unphased_ht.locus2,
        unphased_ht.alleles2
    ).aggregate(
        raw=flatten_gt_counts(gt_counts_raw_expr),
        em_p_chet_raw=get_em_expr(gt_counts_raw_expr).p_chet,
        adj=flatten_gt_counts(gt_counts_adj_expr),
        em_p_chet_adj=get_em_expr(gt_counts_adj_expr).p_chet,
        # pop_max_gt_counts_adj=gt_counts_raw_expr,
        # pop_max_em_p_chet_adj=get_em_expr(gt_counts_raw_expr).p_chet,
    ).key_by()

    phase_ht = vp_ht.filter(hl.is_defined(vp_ht.all_phase))
    phase_ht = phase_ht.key_by()

    def flatten_phase_dict(
            expr: hl.expr.StructExpression
    ) -> hl.expr.StructExpression:
        return hl.struct(
            raw=flatten_gt_counts(expr.gt_counts.raw),
            adj=flatten_gt_counts(expr.gt_counts.adj),
            em_p_chet_raw=expr.em.raw.p_chet,
            em_p_chet_adj=expr.em.adj.p_chet,
            em1_p_chet_raw=expr.em_plus_one.raw.p_chet,
            em1_p_chet_adj=expr.em_plus_one.adj.p_chet
        )
    phase_ht.describe()

    phase_ht = phase_ht.transmute(
        pop_max=phase_ht.pop_max,
        **{
            k: v for k,v in flatten_phase_dict(phase_ht.all_phase).items()
        },
        **{
            f"pmax_{k}": v for k, v in flatten_phase_dict(phase_ht.pop_max[1]).items()
        }
    )
    phase_ht = phase_ht.union(
        unphased_ht,
        unify=True
    )

    phase_ht = phase_ht.transmute(
        chrom=phase_ht.locus1.contig,
        pos1=phase_ht.locus1.position,
        ref1=phase_ht.alleles1[0],
        alt1=phase_ht.alleles1[1],
        pos2=phase_ht.locus2.position,
        ref2=phase_ht.alleles2[0],
        alt2=phase_ht.alleles2[1],
        pop_max=phase_ht.pop_max[0]
    )

    phase_ht.flatten().export('gs://gnomad/projects/compound_hets/Feb_2020_CMG_Compound_het_list_phased.tsv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--file', help='File containing variants to phase.',
                        action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--pbt', help='If --pbt is specified, then only sites present in PBT samples are used and counts exclude PBT samples.',
                        action='store_true')
    parser.add_argument('--least_consequence', help=f'Least consequence for the input (just to get the path right). (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--no_em', help=f'Do not compute EM phase', action='store_true')
    parser.add_argument('--no_lr', help=f'Do not compute likelihood-ratio phase', action='store_true')
    parser.add_argument('--no_shr', help=f'Do not compute single het ratio phase', action='store_true')
    parser.add_argument('--max_freq', help=f'Maximum global adj AF for the input (just to get the path right). (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


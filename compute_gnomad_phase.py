from resources import *
from phasing import *
import argparse
from phasing import get_phased_gnomad_ht

from gnomad.resources.grch37.gnomad import EXOME_POPS

POPS = ["all"] + [pop.lower() for pop in EXOME_POPS]


def main(args):
    hl.init(log="/tmp/phasing.hail.log")

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]

    if args.create_phased_vp_summary:
        ht = hl.read_table(vp_count_ht_path(*path_args))
        ht = get_phased_gnomad_ht(
            ht,
            not args.no_em,
            not args.no_lr,
            not args.no_shr
        )

        ht.write(phased_vp_count_ht_path(*path_args), overwrite=args.overwrite)

    if args.create_phased_vp_summary_release:
        ht = hl.read_table(phased_vp_count_ht_path(*path_args))
        ht = ht.annotate(
            phase_info={
                pop: hl.struct(
                    gt_counts=ht.phase_info[pop].gt_counts.adj,
                    em=ht.phase_info[pop].em.adj.annotate(
                        same_haplotype=ht.phase_info[pop].em.adj.p_chet < args.same_haplotype_em_cutoff,
                        different_haplotype=ht.phase_info[pop].em.adj.p_chet > args.different_haplotypes_em_cutoff,
                    )
                )
                for pop in POPS
            },
        )

        ht = ht.annotate_globals(
            max_freq=args.max_freq,
            least_consequence=args.least_consequence,
            same_haplotype_em_probability_cutoff=args.same_haplotype_em_cutoff,
            different_haplotypes_em_probability_cutoff=args.different_haplotypes_em_cutoff,
            global_annotation_descriptions=hl.struct(
                max_freq="Maximum global variant allele frequency (adj filtered) used as the filter for inclusion in the Hail Table.",
                least_consequence="Consequence used to determine variant inclusion in the Hail Table. The table includes all variants for which the VEP worst_consequence is at least as bad as the least_consequence. The order of consequences for this determination is is taken from gnomad_hail.constants.",
                same_haplotype_em_probability_cutoff="Expectation-Maximization probability cutoff used for the same_haplotype annotation. Variant pairs with an EM probability (em.p_chet) less than this value are likely found on the same haplotype in most individuals in gnomAD.",
                different_haplotypes_em_probability_cutoff="Expectation-Maximization probability cutoff used for the different_haplotype annotation. Variant pairs with an EM probability (em.p_chet) greater than this value are likely found on different haplotypes in most individuals in gnomAD.",
            ),
            row_annotation_descriptions=hl.struct(
                locus1="Locus of the first variant in the variant pair. Contains contig and position information.",
                alleles1="Alleles of the first variant in the variant pair.",
                locus2="Locus of the second variant in the variant pair. Contains contig and position information.",
                alleles2="Alleles of the second variant in the variant pair.",
                phase_info=hl.struct(
                    description="The phase_info annotation is a dictionary of phase information broken down by population. The keys are the populations represented in the gnomAD v2.1.1 exomes (afr, amr, asj, eas, fin, nfe, sas, oth) as well as a key for phasing information for `all` populations combined.",
                    gt_counts="An array of the number of gnomAD v2.1.1 exomes that have the following combinations of genotypes for the two variants: [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]. A/a is the variant defined by locus1 and alleles1 and B/b is the variant defined by locus2 and alleles2. A/B indicates reference and a/b indicates alternate.",
                    em=hl.struct(
                        hap_counts="Array of estimated haplotype counts for the given population.",
                        p_chet="Expectation-Maximization probability that the pair of variants occur on different haplotypes. This is based on their co-occurrence pattern in gnomAD.",
                        same_haplotype="Based on their co-occurrence pattern in gnomAD, these variants are likely found on the same haplotype in most individuals in gnomAD.",
                        different_haplotype="Based on their co-occurrence pattern in gnomAD, these variants are likely found on different haplotypes in most individuals in gnomAD.",
                    ),
                ),
            )
        )
        ht.write(phased_vp_count_ht_path(*path_args, release=True), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('--create_phased_vp_summary', help='Adds gnomAD statistical phasing inforrmation onto the summarised VP table, with counts in release samples only.',
                        action='store_true')
    data_grp.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
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
    parser.add_argument('--chrom', help='Only run on given chromosome')
    parser.add_argument('--create_phased_vp_summary_release', help='Creates the release version of the gnomAD phased summarised VP table created with --create_phased_vp_summary.',
                        action='store_true')
    parser.add_argument('--same_haplotype_em_cutoff', help='EM probability cutoff for same haplotypes in the co-occurrence release HT.', default=0.164, type=float)
    parser.add_argument('--different_haplotypes_em_cutoff', help='EM probability cutoff for different haplotypes in the co-occurrence release HT.', default=0.505, type=float)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    main(args)

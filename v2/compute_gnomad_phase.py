from resources import *
from phasing import *
import argparse
from phasing import get_phased_gnomad_ht
import sys
from google.cloud import storage


from gnomad.resources.grch37.gnomad import EXOME_GEN_ANC_GROUPS

POPS = ["all"] + [pop.lower() for pop in EXOME_GEN_ANC_GROUPS]

def gcs_directory_exists(gcs_path: str) -> bool:
    """
    Checks if a "directory" on GCS contains at least one object.
    
    This is the most robust way to check for directory existence on GCS.
    """
    print(f"Checking for objects within GCS directory: {gcs_path}")
    try:
        # Create a client to interact with the GCS API.
        # Assumes you are authenticated (e.g., via `gcloud auth application-default login`).
        storage_client = storage.Client()

        # Parse the bucket and directory path from the full GCS URI.
        # e.g., 'gs://my-bucket/data/reports/' -> 'my-bucket', 'data/reports/'
        if not gcs_path.startswith("gs://"):
            raise ValueError("GCS path must start with 'gs://'")
        
        bucket_name, directory_prefix = gcs_path.replace("gs://", "").split("/", 1)
        
        # Ensure the prefix ends with a slash to treat it as a directory.
        if not directory_prefix.endswith('/'):
            directory_prefix += '/'

        # Use list_blobs with a prefix and ask for only one result for efficiency.
        blobs = storage_client.list_blobs(bucket_name, prefix=directory_prefix, max_results=1)
        
        # next(blobs, None) will return the first item if it exists, or None otherwise.
        # We just need to know if it's not None.
        return next(blobs, None) is not None
    except Exception as e:
        print(f"An error occurred while checking GCS: {e}", file=sys.stderr)
        return False


def main(args):
    hl.init(log="/tmp/phasing.hail.log")

    data_type = 'exomes' if args.exomes else 'genomes'
    path_args = [data_type, args.pbt, args.least_consequence, args.max_freq, args.chrom]
    
    
    if args.create_phased_vp_summary:
        if args.outfile:
            if gcs_directory_exists(f"{args.outfile}"):
                    print(f"File {args.outfile} already exists, no need to run this.")
                    sys.exit(0)
        if args.version=='v2':
            if not args.testing:
                ht = hl.read_table(vp_count_ht_path(*path_args))
            else:
                ht=hl.read_table(args.infile)
        ht = get_phased_gnomad_ht(
            ht,
            not args.no_em,
            not args.no_lr,
            not args.no_shr
        )
        if not args.outfile:
            ht.write(phased_vp_count_ht_path(*path_args), overwrite=args.overwrite)
        else:
            ht.write(args.outfile, overwrite=args.overwrite)

    if args.create_phased_vp_summary_release:
        if args.outfile:
            if gcs_directory_exists(f"{args.outfile}"):
                    print(f"File {args.outfile} already exists, no need to run this.")
                    sys.exit(0)       
        if args.version=='v2':
            if not args.testing:
                ht = hl.read_table(vp_count_ht_path(*path_args))
            else:
                ht=hl.read_table(args.infile)
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
                least_consequence="Consequence used to determine variant inclusion in the Hail Table. The table includes all variants for which the VEP worst_consequence is at least as bad as the least_consequence. The order of consequences for this determination is is taken from gnomad.utils.vep.CSQ_ORDER.",
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
                        hap_counts="Array of haplotype counts estimated with Expectation-Maximization. Ordered [AB, aB, Ab, ab], where A/a is the variant defined by locus1 and alleles1 and B/b is the variant defined by locus2 and alleles2. A/B indicates reference and a/b indicates alternate.",
                        p_chet="Expectation-Maximization probability that the pair of variants occur on different haplotypes. This is based on their co-occurrence pattern in gnomAD.",
                        same_haplotype="Based on their co-occurrence pattern in gnomAD, these variants are likely found on the same haplotype in most individuals in gnomAD.",
                        different_haplotype="Based on their co-occurrence pattern in gnomAD, these variants are likely found on different haplotypes in most individuals in gnomAD.",
                    ),
                ),
            )
        )
        
        if not args.outfile:
            ht.write(phased_vp_count_ht_path(*path_args, release=True), overwrite=args.overwrite)
        else:
            ht.write(args.outfile, overwrite=args.overwrite)


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
    parser.add_argument('--infile', help='File to be read for create_phased_vp_summary')
    parser.add_argument('--outfile', help='phased outfile, if want to specify')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--version', help='gnomAD version (default: v2)', default='v2', choices=['v2', 'v4'])
    parser.add_argument('--testing', help='For gnomad team only', action='store_false')

    args = parser.parse_args()
    main(args)
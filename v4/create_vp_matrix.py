import argparse
from gnomad.utils.file_utils import file_exists
import logging
import sys


def create_variant_pair_ht(vds):
   return 0 

#check if should overwrite existing file, if so see if file exists. If it does, exit program.
def check_overwrite(outfile, overwrite,logger):
    if not overwrite:
        if file_exists(outfile):
            logger.info(f"{outfile} already exists, exiting program.")
            sys.exit(0)
        else:
            logger.info(f"{outfile} does not exist, running program.")
    else:
        logger.info(f"overwrite is set to True, running program.")

def main(args):
    tmp_dir=args.tmp_dir
    infile_vds=args.infile_vds
    overwrite=args.overwrite
    name=args.name
    data_type = 'exomes' if args.exomes else 'genomes'
  
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    #if you are not overwriting file, check if it exists

   
    #generate just the list of possible variant pairs
    if args.create_vp_list:
        outfile=f"{tmp_dir}/{data_type}_{name}_list.ht"
        check_overwrite(outfile, overwrite,logger)
     
        import hail as hl
        
        #read in vds file
        vds=hl.vds.read_vds(infile_vds)
        vds=create_variant_pair_ht(vds)
                

#The order this should be run in is first create_vp_list (or vp_list_by_chrom), then create_full_vp, then create_vp_summary.

if __name__ == "__main__":    
    
    # Argument parsing for exomes or genomes, testing, and tmp-dir.
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument(
        '--exomes', 
        action='store_true',
        help='Run on exomes. One and only one of --exomes or --genomes is required.')
    data_grp.add_argument('--genomes',
        action='store_true',
        help='Run on genomes. One and only one of --exomes or --genomes is required.')
    parser.add_argument(
        '--testing',
        action='store_true',
        help='if you are testing the pipeline, for developers only')
    parser.add_argument(
        '--tmp-dir',
        required=True,
        help='Temporary directory for intermediate files.'
    )
    parser.add_argument(
        '--create_vp_list',
        action='store_true',
        help='first create just the list of possible variant pairs.'
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Whether to overwrite existing files.')
    parser.add_argument(
        '--name',
        required=True,
        help='unique name to be used in file naming.'
    )
    parser.add_argument(
        '--infile_vds',
        required=False,
        help='Path to input VDS file.'
    )

    args = parser.parse_args()
    main(args)

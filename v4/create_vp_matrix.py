import hail as hl
import argparse
from gnomad.utils.file_utils import file_exists
import logging
import sys
from gnomad_qc.v4.resources.annotations import get_vep, get_freq
from typing import Iterable #for create_full_vp_vds_efficient
from gnomad.utils.vep import CSQ_ORDER, filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.variant_qc import final_filter

def _get_ordered_vp_struct(v1: hl.expr.StructExpression, v2: hl.expr.StructExpression):
    return hl.if_else(
        v1.locus.position < v2.locus.position,
        hl.struct(v1=v1, v2=v2),
        hl.if_else(
            v1.locus.position == v2.locus.position,  # If positions are equal, sort on alt allele
            hl.if_else(
                v1.alleles[1] < v2.alleles[1],
                hl.struct(v1=v1, v2=v2),
                hl.struct(v1=v2, v2=v1)

            ),
            hl.struct(v1=v2, v2=v1)
        )

    )


def create_variant_pair_ht(vds, least_csq, max_freq, genotype_field='LGT'):
    """
    Create a Hail Table of unique ordered variant pairs per sample per gene.

    Parameters
    - vds: object with attribute `variant_data` (a MatrixTable)
    - least_csq: lowest-severity consequence to keep (must be in CSQ_ORDER)
    - max_freq: maximum global AF to keep (inclusive)
    - genotype_field: entry field that indicates genotype (default 'LGT'; use 'GT' if appropriate)

    Returns:
    - hail Table keyed by locus2, alleles2, locus1, alleles1 with one distinct row per unique pair.
    """
    mt = vds.variant_data
    
    #subset to variants that pass QC
    filt_resources=final_filter("exomes")
    filt_resources=filt_resources.ht()
    mt = mt.filter_rows(hl.is_defined(filt_resources[mt.locus, mt.alleles]))


    if least_csq not in CSQ_ORDER:
        raise ValueError(f"least_csq '{least_csq}' not in CSQ_ORDER")

    AF_table=get_freq(data_type="exomes").ht()

    vep_ht = get_vep(data_type="exomes").ht()
    vep_ht = vep_ht.key_by('locus', 'alleles')
    allowed_csqs = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index("3_prime_UTR_variant") + 1]))

    # Filter transcripts using gnomAD helper (first arg must be transcript_consequences expr)
    vep_ht = vep_ht.annotate(
            csqs = filter_vep_transcript_csqs_expr(
                vep_ht.vep.transcript_consequences,
                protein_coding=True,
                ensembl_only=True,
                additional_filtering_criteria=[
                    lambda tc: tc.consequence_terms.any(lambda c: allowed_csqs.contains(c))
                ]
            )
        )

    
    mt_anno = mt.select_rows(
        vep=vep_ht[mt.row_key].csqs.gene_id,
        af=hl.float32(AF_table[mt.row_key].freq[0].AF)
    )
    mt_anno = mt_anno.filter_rows(hl.is_defined(mt_anno.vep) & (hl.len(mt_anno.vep) > 0) & (mt_anno.af > 0) & (mt_anno.af <= max_freq))
    mt_anno = mt_anno.explode_rows(mt_anno.vep)
    mt_anno = mt_anno.rename({'vep': 'gene_id'})
    mt = mt_anno.filter_rows(hl.is_defined(mt_anno.gene_id))

    et = mt.select_cols().select_rows("gene_id").entries()
    et = et.group_by("gene_id", *mt.col_key)._set_buffer_size(5).aggregate(
        vgt=(hl.agg.collect(hl.struct(locus=et.locus, alleles=et.alleles)))
    )

    et = et.annotate(
        vgt=hl.range(0, hl.len(et.vgt))
                      .flatmap(lambda i1: hl.range(i1+1, hl.len(et.vgt))
                               .map(lambda i2: _get_ordered_vp_struct(et.vgt[i1], et.vgt[i2])))
    )

    et = et.explode(et.vgt)
    et = et.key_by(locus2=et.vgt.v2.locus, alleles2=et.vgt.v2.alleles, locus1=et.vgt.v1.locus, alleles1=et.vgt.v1.alleles)
    et = et.select().distinct()


    return et

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
    least_csq=args.least_csq
    max_freq=args.max_freq
  
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
        
        variant_pair_ht=create_variant_pair_ht(vds, least_csq, max_freq, genotype_field='LGT')
        variant_pair_ht.write(outfile, overwrite=overwrite)
        


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
        '--tmp_dir',
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
    parser.add_argument(
        '--least_csq',
        required=False,
        default='3_prime_UTR_variant',
        help='Lowest-severity consequence to keep (must be in CSQ_ORDER). Default is 3_prime_UTR_variant.'
    )
    parser.add_argument(
        '--max_freq',
        type=float,
        required=False,
        default=0.05,
        help='Maximum global AF to keep (inclusive). Default is 0.05.'
    )

    args = parser.parse_args()
    main(args)

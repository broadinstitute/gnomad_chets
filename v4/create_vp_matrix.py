import hail as hl
import argparse
from gnomad.utils.file_utils import file_exists
import logging
import sys
from gnomad_qc.v4.resources.annotations import get_vep, get_freq
from typing import Iterable #for create_full_vp_vds_efficient
from gnomad.utils.vep import CSQ_ORDER, filter_vep_transcript_csqs_expr
from gnomad_qc.v4.resources.variant_qc import final_filter
import timeit


# def create_full_vp(
#         mt: hl.MatrixTable,
#         vp_list_ht: hl.Table,
#         data_type: str,
#         tmp_dir: str
# ):
#     # TODO: This implementation was causing memory challenges.

#     vp_list_ht = vp_list_ht.key_by('locus2', 'alleles2')
#     vp_list_ht = vp_list_ht.select(locus1=vp_list_ht.locus1, alleles1=vp_list_ht.alleles1)
#     vp_mt = mt.annotate_rows(v1=vp_list_ht.index(mt.row_key, all_matches=True))
#     vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)
#     vp_mt = vp_mt.rename({x: f'{x}2' for x in vp_mt.entry})

#     vp_mt = vp_mt.explode_rows(vp_mt.v1)
#     vp_mt = vp_mt.transmute_rows(**vp_mt.v1)
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp0.mt', overwrite=True)

#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp1.mt', overwrite=True)

#     mt_joined = mt[vp_mt.row_key, vp_mt.col_key]
#     vp_mt = vp_mt.annotate_entries(**{f'{x}1': mt_joined[x] for x in mt.entry})
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp2.mt', overwrite=True)
#     vp_mt = vp_mt.repartition(10000, shuffle=True)
#     vp_mt = vp_mt.checkpoint(f'{tmp_dir}/{data_type}_vp_mt_tmp3.mt', overwrite=True)
#     vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})
#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')
    
#     return vp_mt


# # An efficient implementation of creating a VariantDataset of variant pairs. 
# def create_full_vp_vds_efficient(
#     vds_in: hl.vds.VariantDataset,
#     vp_list_ht: hl.Table,                # expected fields: locus1, alleles1, locus2, alleles2
#     tmp_dir: str,
#     entry_fields: Iterable[str] = None,  # None => keep all entry fields; otherwise list like ['GT','DP','GQ']
#     n_partitions: int = 2000,
#     flatten_pairs: bool = False          # if True: create GT1,GT2,... flattened fields (slower)
# ) -> hl.vds.VariantDataset:
#     """
#     Create a VDS of variant pairs efficiently from an input VDS.

#     Returns the VariantDataset object and writes it to out_path_vds.
#     """

#     # Work on variant_data (MatrixTable)
#     mt = vds_in.variant_data

#     # Choose which entry fields to carry. By default, use all.
#     if entry_fields is None:
#         entry_fields = list(mt.entry)

#     # 0) Minimal preselection on the original mt:
#     #    keep only row keys, row fields used by vp_list lookup (locus/alleles), and minimal entry fields.
#     #    This avoids carrying large row/col annotations through the pipeline.
#     mt_min = mt.select_rows().select_cols().select_entries(**{f: mt[f] for f in entry_fields})

#     # 1) Prepare vp_list keyed by locus2/alleles2
#     vp_list_ht = vp_list_ht.key_by('locus2', 'alleles2')
#     vp_list_ht = vp_list_ht.select(locus1=vp_list_ht.locus1, alleles1=vp_list_ht.alleles1)

#     # 2) Attach to mt_min an array of v1 partners (for each v2)
#     #    This requires indexing by mt_min.row_key (which is (locus, alleles))
#     vp_mt = mt_min.annotate_rows(v1 = vp_list_ht.index(mt_min.row_key, all_matches=True))

#     # 3) Keep only v2 rows that actually appear in pairs
#     vp_mt = vp_mt.filter_rows(hl.len(vp_mt.v1) > 0)

#     # 4) Turn existing entry fields into a single nested struct v2 to avoid many columns
#     #    This is much more compact in Spark and Hail than duplicating many top-level fields.
#     vp_mt = vp_mt.annotate_entries(v2 = hl.struct(**{f: vp_mt[f] for f in entry_fields}))
#     # drop the original entry fields, keep only v2
#     vp_mt = vp_mt.select_entries('v2')

#     # 5) Explode v1 list so each row is one (v1, v2) pair (v2 = current row)
#     vp_mt = vp_mt.explode_rows(vp_mt.v1)
#     vp_mt = vp_mt.transmute_rows(**vp_mt.v1)  # brings locus1/alleles1 to top-level row fields

#     # checkpoint an intermediate MT to break lineage (helps memory)
#     cp0 = f'{tmp_dir}/vp_pairs_tmp0.mt'
#     vp_mt = vp_mt.checkpoint(cp0, overwrite=True)

#     # 6) Re-key by locus1/alleles1 so we can slice original mt_min to fetch v1 entries cheaply
#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1')

#     # checkpoint again to persist the rekeyed state
#     cp1 = f'{tmp_dir}/vp_pairs_tmp1.mt'
#     vp_mt = vp_mt.checkpoint(cp1, overwrite=True)

#     # 7) Create a tiny mt that contains entry struct 'v1' for the original variants to slice from
#     mt_for_v1 = mt_min.annotate_entries(v1 = hl.struct(**{f: mt_min[f] for f in entry_fields}))
#     mt_for_v1 = mt_for_v1.select_entries('v1')

#     # 8) Slice mt_for_v1 by the vp_mt keys to get v1 entries aligned to vp_mt rows/cols.
#     #    This is a partition-local slice and is much cheaper than full join/repartition.
#     mt_joined = mt_for_v1[vp_mt.row_key, vp_mt.col_key]

#     # 9) Annotate vp_mt entries with v1 struct
#     vp_mt = vp_mt.annotate_entries(v1 = mt_joined.v1)

#     # At this point entry has two nested structs: v1 and v2

#     # checkpoint a compact paired MT
#     cp2 = f'{tmp_dir}/vp_pairs_tmp2.mt'
#     vp_mt = vp_mt.checkpoint(cp2, overwrite=True)

#     # 10) Repartition moderately (n_partitions) to make subsequent writing / VDS packaging efficient.
#     #     Avoid excessively large shuffles; n_partitions should match cluster size.
#     vp_mt = vp_mt.repartition(n_partitions, shuffle=True)

#     cp3 = f'{tmp_dir}/vp_pairs_tmp3.mt'
#     vp_mt = vp_mt.checkpoint(cp3, overwrite=True)

#     # 11) Rename locus/alleles -> locus2/alleles2 (current row fields represent the original v2)
#     vp_mt = vp_mt.rename({'locus': 'locus2', 'alleles': 'alleles2'})

#     # 12) Final canonical row-key ordering
#     vp_mt = vp_mt.key_rows_by('locus1', 'alleles1', 'locus2', 'alleles2')

#     # 13) Optionally flatten the nested structs into top-level fields like GT1, GT2, DP1, DP2, ...
#     if flatten_pairs:
#         # Flatten v1 and v2 structs into top-level entry fields with suffixes
#         # This increases width of the entry schema but may be necessary for downstream code.
#         def flatten_struct_to_entries(mt_in, struct_name, suffix):
#             # mt_in: MatrixTable with entry struct named struct_name (e.g., 'v1' or 'v2')
#             # returns a MatrixTable with additional top-level entry fields like GT1 = mt_in.entry.v1.GT
#             fields = list(entry_fields)
#             # annotate_entries with each flattened field
#             mt_out = mt_in.annotate_entries(**{f'{f}{suffix}': getattr(mt_in.entry[struct_name], f) for f in fields})
#             # we do not drop the nested struct because downstream code might expect it; you can drop it if desired:
#             mt_out = mt_out.select_entries(*[f'{f}{suffix}' for f in fields])
#             return mt_out

#         vp_mt = flatten_struct_to_entries(vp_mt, 'v1', '1')
#         vp_mt = flatten_struct_to_entries(vp_mt, 'v2', '2')

#         # optionally keep nested structs as well or drop them (we dropped them above)

#         # checkpoint after flatten
#         cp_flat = f'{tmp_dir}/vp_pairs_flattened.mt'
#         vp_mt = vp_mt.checkpoint(cp_flat, overwrite=True)

#     # 14) Build final VariantDataset with original sample_data (unchanged)
#     vds_out = hl.vds.VariantDataset(variant_data=vp_mt, sample_data=vds_in.sample_data)
#     return vds_out


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

def create_variant_pair_ht(vds, least_csq, max_freq,tmp_dir,name, genotype_field='LGT'):
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
    #logging.info(f"When first reading in the VDS, count number of variants: %s", mt.count())
    
    
    #subset to variants that pass QC
    filt_resources=final_filter("exomes").ht()
    AF_table=get_freq(data_type="exomes").ht()
    vep_ht = get_vep(data_type="exomes").ht().key_by('locus', 'alleles')
    
    allowed_csqs = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_csq) + 1]))
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
    

    #logging.info("After filtering to final variants, count number of variants: %s",mt.count())

    if least_csq not in CSQ_ORDER:
        raise ValueError(f"least_csq '{least_csq}' not in CSQ_ORDER")

    mt = mt.annotate_rows(
        passes_qc=hl.is_defined(filt_resources[mt.locus, mt.alleles]),
        gene_id=vep_ht[mt.row_key].csqs.gene_id,
        af=hl.float32(AF_table[mt.row_key].freq[0].AF)
    )
    
    mt = mt.filter_rows(
                        mt.passes_qc &
                        hl.is_defined(mt.gene_id) & 
                        (hl.len(mt.gene_id) > 0) & 
                        (mt.af > 0) & 
                        (mt.af <= max_freq)
    )

    #mt = mt.filter_rows(hl.is_defined(filt_resources[mt.locus, mt.alleles]))
    # mt_anno = mt.select_rows(
    #     vep=vep_ht[mt.row_key].csqs.gene_id,
    #     af=hl.float32(AF_table[mt.row_key].freq[0].AF)
    # )    
    
    mt=mt.select_rows("gene_id", "af")
    mt = mt.explode_rows(mt.gene_id)
    #logging.info("After exploding on gene_id, count number of variants %s", mt.count())
    

    filtered_vds=hl.vds.VariantDataset(vds.reference_data, mt)
    logging.info("Now turning to dense MatrixTable")
    mt=hl.vds.to_dense_mt(filtered_vds)
    logging.info("Dense MatrixTable created")
    
    mt = mt.checkpoint(f'{tmp_dir}/filtered_{name}_mt.mt', overwrite=True) #avoids recomputation

    et = mt.select_cols().select_rows("gene_id").entries()
    #per gene_id x sample, collect list of variants
    # Step 1: Collect unique variants per gene/sample
    et = et.group_by("gene_id", *mt.col_key)._set_buffer_size(20).aggregate(
        variants=hl.array(hl.agg.collect_as_set(hl.struct(locus=et.locus, alleles=et.alleles)))
    )
    

    # Step 2: Filter samples with <2 variants
    et = et.filter(hl.len(et.variants) >= 2)

    # Step 3: Create pairs
    et = et.annotate(
        pairs=hl.range(0, hl.len(et.variants))
            .flatmap(lambda i1: hl.range(i1+1, hl.len(et.variants))
                    .map(lambda i2: _get_ordered_vp_struct(et.variants[i1], et.variants[i2])))
    )

    # Step 4: Rename for consistency if needed
    et = et.transmute(vgt=et.pairs)

    #logging.info("After generating variant pairs, count number of variant pairs: %s", et.count())
    
    et = et.explode(et.vgt)
    
    et = et.transmute(
        locus1=et.vgt.v1.locus,
        alleles1=et.vgt.v1.alleles,
        locus2=et.vgt.v2.locus,
        alleles2=et.vgt.v2.alleles
    )
    
    #logging.info("After exploding variant pairs, count number of variant pairs: %s", et.count())
    
    et = et.key_by('locus2', 'alleles2', 'locus1', 'alleles1')
    et = et.select().distinct()
    
    #logging.info("After distinct(), count number of unique variant pairs: %s", et.count())


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
    start = timeit.default_timer()

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
        
        #variant_pair_ht=create_variant_pair_ht_ultra_fast(vds, least_csq, max_freq,tmp_dir)
        variant_pair_ht=create_variant_pair_ht(vds, least_csq, max_freq,tmp_dir,name, genotype_field='LGT')
        variant_pair_ht.write(outfile, overwrite=overwrite)
        
    if args.create_full_vp:
        outfile=f"{tmp_dir}/{data_type}_{name}_full.vds"
        check_overwrite(outfile, overwrite,logger)
        
        import hail as hl
        
        #read in vds file
        vds=hl.vds.read_vds(infile_vds)
        
        #read in variant pair ht
        vp_list_ht=hl.read_table(f"{tmp_dir}/{data_type}_{name}_list.ht")
        
        #add variant pairs to vds
        full_vp_vds=create_full_vp_vds_efficient(vds, vp_list_ht, tmp_dir)
        
        hl.vds.write_vds(full_vp_vds, outfile, overwrite=True)

    stop = timeit.default_timer()
    logger.info(f"Time taken to run the script is {stop - start} seconds.")

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
    parser.add_argument(
        '--create_full_vp',
        action='store_true',
        help='then create the full variant pair VDS.'
    )

    args = parser.parse_args()
    main(args)

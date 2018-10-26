#!/bin/bash

#Drop table if exists already
bq rm -f -t gnomad.exome_variant_pairs
#
##Create empty table
bq mk --force --schema gene:string,pop:string,prob_same_haplotype:float,chrom1:string,pos1:integer,ref1:string,alt1:string,chrom2:string,pos2:integer,ref2:string,alt2:string,g0000:integer,g0100:integer,g1100:integer,g0001:integer,g0101:integer,g1101:integer,g0011:integer,g0111:integer,g1111:integer,h00,h01,h10,h11 -t gnomad.exome_variant_pairs

##Load genotypes data
bq load --source_format=CSV -F "\t" --null_marker=NA --skip_leading_rows 1 gnomad.exome_variant_pairs gs://gnomad/projects/compound_hets/browser/gnomad_exomes_adj_all_vp.part0.tsv.bgz.tsv.bgz/part-*.bgz



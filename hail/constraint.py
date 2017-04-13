from utils import *


hc = hail.HailContext(log="/hail.log")

context_vds_path = "gs://gnomad-resources/Homo_sapiens_assembly19.fasta.vds"
mega_annotations_path = 'gs://hail-common/annotation_v0.2.vds'

def import_fasta():
    vds = hc.import_fasta("gs://gnomad-resources/Homo_sapiens_assembly19.fasta",
                          filter_Ns=True, flanking_context=3, create_snv_alleles=True, create_deletion_size=2,
                          create_insertion_size=2, line_limit = 5000)
    vds.write(context_vds_path, overwrite=True)





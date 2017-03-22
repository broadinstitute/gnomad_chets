
from utils import *

# Inputs

vds_path = "gs://gnomad/gnom.ad.vds"
meta_path = "gs://gnomad/gnomad.final.all_meta.txt"

# Outputs
odonovan_prefix = "gs://gnomad-genomes/subsets/odonovan_quads/odonovan_quads"

# Actions
odonovan = True


hc = hail.HailContext(log="/hail_log/subset.log")


if odonovan:
    vds = hc.read(vds_path)
    vds = vds.annotate_samples_table(meta_path, 'Sample', root='sa.meta', config=hail.TextTableConfig(impute=True))
    vds = vds.filter_samples_expr('sa.meta.Title == "PX10 O\'Donovan Quads G4L")')
    vds







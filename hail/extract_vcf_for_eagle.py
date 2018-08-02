from hail import *
from utils import *


hc = HailContext()
for data_type in ['exomes','genomes']:
    vds = get_gnomad_data(hc, data_type, hardcalls="adj", split=True, release_annotations=CURRENT_RELEASE)
    vds = vds.filter_intervals(Interval.parse("21:1-200000000"))
    vds = vds.filter_variants_expr('va.info.AC > 1')
    vds = vds.filter_samples_expr('isMissing(sa.fam.famID)')
    vds.export_vcf("gs://gnomad/projects/compound_hets/eagle/vcf/{}.chr21.no_singletons.no_fams.split.vcf.bgz".format(data_type))

from utils import *

gnomad_path = "gs://gnomad/gnom.ad.vds"

hc = HailContext(log='/variantqc.log')

(
    hc.read(gnomad_path)
    .
    .sample_qc()
    .export_samples(output='gs://gnomad/sampleqc.')

 )






truth_dir = 'gs://gnomad-public/truth-sets'
omni_path = "%s/1000G_omni2.5.b37.vds" % truth_dir
mills_path = "%s/Mills_and_1000G_gold_standard.indels.b37.vds" % truth_dir
hapmap_path = "%s/hapmap_3.3.b37.vds" % truth_dir
dbsnp_vcf = "%s/vcf/All_20160601.vcf.bgz" % truth_dir
NA12878_path = "%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vds" % truth_dir
NA12878_high_conf_regions_path = "%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed" % truth_dir


lcr_path = "gs://gnomad-public/intervals/LCR.interval_list"
decoy_path = "gs://gnomad-public/intervals/mm-2-merged.bed.gz"
syndip_path = "gs://gnomad/truth-sets/hybrid.m37m.vds"
syndip_high_conf_regions_path = "gs://gnomad/truth-sets/hybrid.m37m.bed"
NA12878_high_conf_exome_regions_path = "gs://exac2/intervals/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
exomes_high_conf_regions_path = "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list"

exome_calling_intervals = 'gs://gnomad-exomes/intervals/exome_calling_regions.v1.interval_list'
evaluation_intervals = 'gs://gnomad-exomes/intervals/exome_evaluation_regions.v1.intervals'
high_coverage_intervals = 'gs://gnomad-exomes/intervals/high_coverage.auto.interval_list'

additional_vcf_header = "gs://gnomad/gnomad.extra_header_fields.vcf"

vep_config = "/vep/vep-gcloud.properties"
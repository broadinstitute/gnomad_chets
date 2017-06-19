#Evaluation data
truth_dir = 'gs://gnomad-public/truth-sets'
omni_path = "%s/1000G_omni2.5.b37.vds" % truth_dir
mills_path = "%s/Mills_and_1000G_gold_standard.indels.b37.vds" % truth_dir
hapmap_path = "%s/hapmap_3.3.b37.vds" % truth_dir
dbsnp_vcf = "%s/vcf/All_20160601.vcf.bgz" % truth_dir
kgp_high_conf_snvs = "%s/vcf/1000G_phase1.snps.high_confidence.b37.vds" % truth_dir
NA12878_path = "%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vds" % truth_dir
NA12878_high_conf_regions_path = "%s/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed" % truth_dir
NA12878_high_conf_exome_regions_path = "gs://exac2/intervals/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
syndip_path = "%s/hybrid.m37m.vds" % truth_dir
syndip_high_conf_regions_path = "%s/hybrid.m37m.bed" % truth_dir
clinvar_variants = "gs://gnomad-resources/annotations/clinvar_alleles.single.b37.tsv.gz"
clinvar_vds = "gs://gnomad-resources/annotations/clinvar_alleles.single.b37.vds"

#Exome/genome duplicate samples
exomes_to_combined_IDs = "gs://gnomad-resources/exomes_to_combined.IDs.txt"
exomes_qc_pass_samples = "gs://gnomad-resources/exomes_qc_pass_samples.txt.gz" #Samples that didn't fail QC-metric filters (contains relateds and non-releasable samples)
genomes_to_combined_IDs = "gs://gnomad-resources/genomes_to_combined.IDs.txt"
genomes_qc_pass_samples = "gs://gnomad-resources/genomes_qc_pass_samples.txt.gz" #Samples that didn't fail QC-metric filters (contains relateds and non-releasable samples)

#Usefult intervals
lcr_path = "gs://gnomad-public/intervals/LCR.interval_list"
decoy_path = "gs://gnomad-public/intervals/mm-2-merged.bed.gz"
purcell5k_path = "gs://gnomad-public/intervals/purcell5k.interval_list"

#Exome intervals
exomes_high_conf_regions_path = "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list"
exome_calling_intervals = 'gs://gnomad-exomes/intervals/exome_calling_regions.v1.interval_list'
exome_calling_noheader_intervals = 'gs://gnomad-exomes/intervals/exome_calling_regions.v1.nohead.interval_list'
evaluation_intervals = 'gs://gnomad-exomes/intervals/exome_evaluation_regions.v1.intervals'
high_coverage_intervals = 'gs://gnomad-exomes/intervals/high_coverage.auto.interval_list'

additional_vcf_header = "gs://gnomad-resources/gnomad.extra_header_fields.vcf"

vep_config = "/vep/vep-gcloud.properties"

# Full VDSs
full_exac_v1_vds = 'gs://gnomad-exomes-raw/exacv1/exac.all.vds'
full_exome_vds = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
full_exome_hardcalls_vds = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw_hardcalls.vds'
full_exome_hardcalls_split_vds = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.raw_hardcalls.split.vds'
full_genome_vds = 'gs://gnomad-genomes-raw/full/gnomad.genomes.all.vds'
full_genome_hardcalls_vds = "gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.vds"
full_genome_hardcalls_split_vds = "gs://gnomad-genomes-raw/hardcalls/gnomad.raw_hardcalls.split.vds"
full_genomes_vep_split_vds = "gs://gnomad-genomes-raw/gnomad.genomes.vep.split.vds"
full_exomes_vep_split_vds = "gs://gnomad-exomes-raw/gnomad.exomes.vep.split.vds"

#Filtering VDSs
vqsr_vds_path = 'gs://gnomad-exomes/variantqc/gnomad.exomes.vqsr.unsplit.vds'
genomes_rf_path = "gs://gnomad-genomes/variantqc/RF/gnomad.sites.RF.newStats24.vds"
exomes_rf_path = "gs://gnomad-exomes/variantqc/gnomad.exomes.rf.vds"

# Release Sites VDSs
final_exac_sites_vds = 'gs://gnomad-exomes-raw/exacv1/exac.sites.vds'
final_exome_autosomes = 'gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.autosomes.vds'
final_genome_autosomes = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.autosomes.vds'
final_exome_vds = 'gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.vds'
final_exome_split_vds = 'gs://gnomad-exomes/sites/vds/gnomad.exomes.r2.0.1.sites.split.vds'
final_genome_vds = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.vds'
final_genome_split_vds = 'gs://gnomad-genomes/sites/vds/gnomad.genomes.r2.0.1.sites.split.vds'

#Meta data
genomes_meta = "gs://gnomad-genomes-raw/gnomad.final.all_meta.txt.bgz"
genomes_fam = "gs://gnomad-genomes-raw/gnomad.final.goodTrios.fam"
exomes_meta = 'gs://gnomad-exomes-raw/super_meta_april_01_2017.txt.bgz'
exomes_fam = "gs://gnomad-exomes/variantqc/gnomad_exomes.qctrios.fam"

#PCA
gnomad_pca = "gs://gnomad-genomes/sampleqc/gnomad.pca.vds"

#Annotations
methylation_kt = "gs://gnomad-resources/methylation.kt"
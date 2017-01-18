from pyhail import *

#Load data in VDS format
hc = HailContext(log='/tmp/hail.log')
( hc.import_vcf('gs://gnomad/genomes.*.vcf.bgz')
  .write('gs://gnomad/genomes.vds')
)

vds = hc.read('gs://gnomad/genomes.vds')

#Remove LCR and SegDup intervals
(vds.filter_variants_intervals('gs://gnomad/lcr.intervals', keep=False)
 .filter_variants_intervals('gs://gnomad/decoy.intervals', keep=False)
 )

#Compute sample QC stats
vds.sample_qc()

#Impute sex
vds.impute_sex()

#Export sample QC metrics
vds.export_samples(
    output='gs://gnomad/genomes.sampleqc.txt',
    condition='s.id,'
              'sa.qc.*,'
              'sa.imputesex.*'
)

#PCA
#Sample well behaved sites in exomes
exomes = hc.read('gs://gnomad/exomes.vds')
(
    exomes.filter_variants_intervals('gs://gnomad/exome_evaluation.intervals')
    .annotate_samples_table('gs://gnomad/exomes.meta.txt',
                            sample_expr='Sample',
                            root='sa.meta')
    .filter_samples('sa.meta.keep')
    .annotate_variants_expr("va.callstats = gs.callStats(g => v) ")
    .filter_variants('v.isBiallelic && '
                     'va.callstats.AF[1] > 0.001 && '
                     'gs.fraction(g => g.isCalled) > 0.99')
    .export_plink('gs://gnomad/exomes.good_sites')
)
##LD-pruning using PLINK (soon in Hail!)

#Load genomes and filter samples / variants
genomes = (
    hc.read('gs://gnomad/genomes.vds')
    .annotate_samples_table('gs://gnomad/genomes.meta.txt',
                            sample_expr='Sample', root='sa.meta')
    .filter_samples('sa.meta.keep')
    .filter_variants_intervals('gs://gnomad/pca_sites.txt')
)
exomes = (
    hc.read('gs://gnomad/exomes.vds')
    .annotate_samples_table('gs://gnomad/exomes.meta.txt',
                            sample_expr='Sample', root='sa.meta')
    .filter_samples('sa.meta.keep')
    .filter_variants_intervals('gs://gnomad/pca_sites.txt')
)

#Run PCA jointly
all = (
    genomes.join(exomes)
    .filter_multi
    .pca(scores='sa.pca', components=10)
)
all.export_samples(
    output='gs://gnomad/gnomad.pca.txt',
    condition='s.id, sa.pca.*'
)

##Variant QC
#Compute Mendel errors and transmission
genomes.mendel_errors('gs://gnomad/gnomad.mendel')
genomes.tdt(fam='gs://gnomad/genomes.fam', root='va.tdt')

#Create allele-specific metrics
genomes.annotate_alleles_expr(
    'va.stats.gq = gs.filter(g => g.isCalledNonRef).map(g => g.gq).stats()',
    'va.stats.dp = gs.filter(g => g.isCalledNonRef).map(g => g.dp).stats()',
    'va.stats.nrq = gs.filter(g => g.isCalledNonRef).map(g => g.dosage[0]).stats()',
    'va.stats.ab = gs.filter(g => g.isHet).map(g => g.ad[1]/g.dp).stats()'
)

features = ['va.variantType',
             'va.info.MQ',
             'va.info.MQRankSum',
             'va.info.SOR',
             'va.info.InbreedingCoeff',
             'va.info.ReadPosRankSum',
             'va.stats.raw.nrq_median',
             'va.stats.raw.ab_median',
             'va.stats.raw.dp_median',
             'va.stats.raw.gq_median']

(
    genomes.split_multi
    .annotate_variants_vds('gs://gnomad/hapmap.vds', code='va.hapmap = isDefined(vds)')
    .annotate_variants_vds('gs://gnomad/omni.vds', code='va.omni = isDefined(vds)')
    .annotate_variants_vds('gs://gnomad/mills.vds', code='va.mills = isDefined(vds)')
    .annotate_variants_expr('va.transmitted_singleton = va.tdt.nTransmitted == 1 && '
                            'va.calldata.raw.AC[va.aIndex]==2')
    .annotate_variants_expr('va.TP = va.omni || va.mills || va.transmitted_singleton, '
                               'va.FP = va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30')
    .annotate_variants_expr('va.label = if(!isMissing(va.FP) && va.FP) "FP" '
                            'else if(va.TP) "TP" '
                            'else NA: String, '
                            'va.train = v.contig != "20" && (va.TP || va.FP)')
    .random_forests(training='va.train', label='va.label', root='va.RF',
                    features=features, num_trees=500, max_depth=5)
    .write('gs://gnomad/genomes.RF.vds')
)


### Other examples for slides
vds.filter_variants('v.nAltAlleles > 1', keep=True)
vds.annotate_variants_expr('va.allelesTxt = v.altAlleles.map(a => a.alt).mkString(",")')

vds.annotate_variants_expr('va.hweCase = gs.filter(g => sa.isCase).hardyWeinberg()')

vds.annotate_gloabl_expr('global.chr_indels = variants.filter(v => v.altAllele.isIndel)'
                         '                            .map(v => v.contig).counter()')




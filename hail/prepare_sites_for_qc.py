
from variantqc import *

hc = HailContext()

get_exome_data = False
get_genome_data = False
join_exomes_genome = False
export_to_plink = False
ld_prune = False
pca = True
export_pca = True

bucket = 'gs://gnomad/'
root = '%s/sampleqc' % bucket

cpg_vds_path = 'gs://gnomad-public/cpg.vds'

joined_vds_path = 'gs://gnomad/hardcalls/gnomad.all.common_SNPs.vds'
plink_output_path = '%s/gnomad' % root

common_criteria = 'v.isBiallelic && !va.cpg && v.altAllele.isSNP && v.isAutosomal && va.calldata.all_samples_raw.AF[1] > 0.001 && va.callrate > 0.99'


def main():
    full_exome_vds_path = 'gs://gnomad-exomes-raw/full/gnomad.exomes.all.vds'
    common_sites_exome_vds_path = 'gs://gnomad-exomes-raw/hardcalls/gnomad.exomes.common_SNPs.vds'

    if get_exome_data:
        cpg_vds = hc.read(cpg_vds_path)
        full_vds = hc.read(full_exome_vds_path)
        write_common_variants(cpg_vds, full_vds, common_sites_exome_vds_path)

    full_genome_vds_path = 'gs://gnomad/gnom.ad.vds'
    common_sites_genome_vds_path = 'gs://gnomad/hardcalls/gnomad.genomes.common_SNPs.vds'
    genome_sample_list = 'gs://gnomad/gnomad.hardfiltered_samples.txt'

    if get_genome_data:
        cpg_vds = hc.read(cpg_vds_path)
        genome_vds = hc.read(full_genome_vds_path).filter_samples_list(genome_sample_list, keep=False)
        write_common_variants(cpg_vds, genome_vds, common_sites_genome_vds_path)

    if join_exomes_genome:
        genome_sample_rename = '/tmp/genome_ids.txt'
        exome_sample_rename = '/tmp/exome_ids.txt'

        exome_vds = hc.read(common_sites_exome_vds_path)
        samples = exome_vds.query_samples('samples.map(s => s.id).collect()')[0]
        with open(exome_sample_rename, 'w') as g:
            for sample in samples:
                g.write('%s\texome_%s\n' % (sample, sample.replace(' ', '_')))

        genome_vds = hc.read(common_sites_genome_vds_path)
        samples = genome_vds.query_samples('samples.map(s => s.id).collect()')[0]
        with open(genome_sample_rename, 'w') as g:
            for sample in samples:
                g.write('%s\tgenome_%s\n' % (sample, sample.replace(' ', '_')))

        (exome_vds.rename_samples('file://%s' % exome_sample_rename).annotate_samples_expr('sa = {}')
         .join(genome_vds.rename_samples('file://%s' % genome_sample_rename).annotate_samples_expr('sa = {}'))
         .annotate_variants_expr(['va.calldata.combined = gs.callStats(g => v)', 'va.callrate = gs.fraction(g => g.isCalled)'])
         .filter_variants_expr('va.callrate > 0.99 && va.calldata.combined.AF[1] > 0.001')
         .write(joined_vds_path))

    if export_to_plink:
        vds = hc.read(joined_vds_path)
        vds.filter_multi().export_plink(plink_output_path)

    if ld_prune:
        # We used plink for now (2/2/17), but this is for the future
        vds = hc.read(joined_vds_path)
        vds.ld_prune(0.1, num_cores=80).export_variants("%s/ldpruned.variants" % root, "v")

    # perl -p -e 's/ /_/g' /humgen/atgu1/fs03/kristen/exac/relatedness/removals/working_cp4_samples_to_keep.txt | awk '{print "exome_"$0}' > working_cp4_samples_to_keep_translated.txt
    # perl -p -e 's/ /_/g' /humgen/atgu1/fs03/kristen/exac/relatedness/removals/working_genome_cp4_samples_to_keep.txt | awk '{print "genome_"$0}' > working_genome_cp4_samples_to_keep_translated.txt
    if pca:
        vds = (hc.read(joined_vds_path)
               .annotate_samples_list('%s/working_cp4_samples_to_keep_translated.txt' % root, root='sa.exome_keep')
               .annotate_samples_list('%s/working_genome_cp4_samples_to_keep_translated.txt' % root, root='sa.genome_keep')
               .filter_samples_expr('sa.exome_keep || sa.genome_keep')
               .filter_multi())
        vds = (vds.filter_variants_list('%s/gnomad_0.1_0.001.variants.txt' % root))
        vds = (vds.filter_variants_expr('va.callrate > 0.999'))
        print vds.count()

        (vds.pca(scores='sa.pca', loadings='va.pca_loadings', eigenvalues='global.pca_evals', k=20)
         .write('%s/gnomad.pca.vds' % root, overwrite=True))

    if export_pca:
        pca_vds = hc.read('%s/gnomad.pca.vds' % root)
        pca_vds.export_samples('%s/gnomad.pca.txt.bgz' % root, 'Sample = s, sa.pca.*')


def write_common_variants(cpg_vds, full_vds, output_path):
    (full_vds.annotate_variants_expr(['va.calldata.all_samples_raw = gs.callStats(g => v)', 'va.callrate = gs.fraction(g => g.isCalled)'])
     .annotate_variants_vds(cpg_vds, code='va.cpg = !isMissing(vds)')
     .filter_variants_expr(common_criteria)
     .hardcalls()
     .write(output_path))

if __name__ == '__main__':
    main()
from utils import *

def get_kt_sparse_vsm(vds, keys):
    kt = vds.filter_genotypes('!g.isHomRef').genotypes_table()
    kt = kt.annotate('stuff = {v:v, va:va, s:s, sa:sa, gt:g.gt}')  # Dump sa ?
    kt = kt.aggregate_by_key(keys, 'index(stuff.collect(), v)')


def annotate_methylation(vds):
    mkt = vds.hc.read_table(methylation_kt_path).select(['locus', 'MEAN'])
    vds = vds.annotate_variants_table(mkt, root='va.methylation.value')
    vds = vds.annotate_variants_expr(
        ['va.methylated_cpg = v.altAllele().isTransition() && va.methylation.value >= 0.25'])
    #logger.debug('Number of methylated CpGs: {}'.format(
    #    str(vds.query_variants(['variants.map(v => va.methylated_cpg).counter()'])[0])))
    return vds


def flatten_counts(kt, gc_ann="genotype_counts", hc_ann="haplotype_counts", gt_anns=None, out_prefix=""):
    gc_cols = ['AABB', 'AaBB', 'aaBB', 'AABb', 'AaBb', 'aaBb', 'AAbb', 'Aabb', 'aabb']
    hc_cols = ['AB', 'Ab', 'aB', 'ab']

    gt_expr = {}
    for a in gt_anns:
        gt_expr[a + '_gt'] = a + '.gt'
        gt_expr[a + '_dp'] = a + '.dp'
        gt_expr[a + '_gq'] = a + '.gq'
        gt_expr[a + '_ad0'] = a + '.ad[0]'
        gt_expr[a + '_ad1'] = a + '.ad[1]'
        gt_expr[a + '_pl0'] = a + '.pl[0]'
        gt_expr[a + '_pl1'] = a + '.pl[1]'
        gt_expr[a + '_pl2'] = a + '.pl[2]'

    return (
        kt.annotate(['{0}{1} = {2}[{3}]'.format(out_prefix, gt, gc_ann, gc_cols.index(gt)) for gt in gc_cols] +
                    ['{0}{1} = {2}[{3}]'.format(out_prefix, hp, hc_ann, hc_cols.index(hp)) for hp in hc_cols] +
                    ['{0}{1} = {2}'.format(out_prefix, ann, expr) for ann,expr in gt_expr.iteritems()]),
        [out_prefix + gt for gt in gc_cols] + [out_prefix + hp for hp in hc_cols] + [out_prefix + gt for gt in gt_expr.keys()]
    )


def annotate_gene_impact(vds):
    vds = filter_vep_to_canonical_transcripts(vds)
    vds = process_consequences(vds)

    impact = {x: "high" for x in CSQ_CODING_HIGH_IMPACT  }
    impact.update({x: "medium" for x in CSQ_CODING_MEDIUM_IMPACT})
    impact.update({x: "low" for x in CSQ_CODING_LOW_IMPACT})

    vds = vds.annotate_global('global.impact', impact, TDict(TString(), TString()))

    vds = vds.annotate_variants_expr([
        'va.impact = if("-LC" ~ va.vep.worst_csq_suffix ) "medium" else if (global.impact.contains(va.vep.worst_csq)) global.impact[va.vep.worst_csq] else NA:String',
        'va.gene = if("-HC" ~ va.vep.worst_csq_suffix) '
        '   va.vep.transcript_consequences.find(x => x.lof == "HC").gene_symbol '
        'else'
        '   va.vep.transcript_consequences.find(x => x.consequence_terms.toSet.contains(va.vep.worst_csq)).gene_symbol',
        'va.alleleType = if(v.altAllele.isSNP) "SNP" else "indel"'
    ])

    #logger.debug("Variant impact counts: {0}, Variants with gene annotations: {1}".format(*vds.query_variants(
    #    ['variants.map(v => va.impact).counter()', 'variants.map(v => isDefined(va.gene)).counter()'])))
    vds = vds.filter_variants_expr('isDefined(va.impact) && isDefined(va.gene)')
    return vds

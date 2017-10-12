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
    if gt_anns is not None:
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


def annotate_gene_impact(vds, vep_root='va.vep'):
    vds = filter_vep_to_canonical_transcripts(vds, vep_root)
    vds = process_consequences(vds, vep_root)

    impact = {x: "high" for x in CSQ_CODING_HIGH_IMPACT  }
    impact.update({x: "medium" for x in CSQ_CODING_MEDIUM_IMPACT})
    impact.update({x: "low" for x in CSQ_CODING_LOW_IMPACT})

    vds = vds.annotate_global('global.impact', impact, TDict(TString(), TString()))

    vds = vds.annotate_variants_expr([
        'va.impact = if("-LC" ~ {0}.worst_csq_suffix ) "medium" else if (global.impact.contains({0}.worst_csq)) global.impact[{0}.worst_csq] else NA:String'.format(vep_root),
        'va.gene = if("-HC" ~ {0}.worst_csq_suffix) '
        '   {0}.transcript_consequences.find(x => x.lof == "HC").gene_symbol '
        'else'
        '   {0}.transcript_consequences.find(x => x.consequence_terms.toSet.contains({0}.worst_csq)).gene_symbol'.format(vep_root),
        'va.alleleType = if(v.altAllele.isSNP) "SNP" else "indel"'
    ])

    #logger.debug("Variant impact counts: {0}, Variants with gene annotations: {1}".format(*vds.query_variants(
    #    ['variants.map(v => va.impact).counter()', 'variants.map(v => isDefined(va.gene)).counter()'])))
    vds = vds.filter_variants_expr('isDefined(va.impact) && isDefined(va.gene)')
    return vds


def get_variants_phase_from_ref(kt, vds, va_agg, sa_agg=None, num_partitions=None, va_to_keep=None):
    """

    Runs EM and aggregate results based on a KeyTable with 2 variant columns

    :param KeyTable kt: Input KeyTable needs to be keyed by exactly to columns of type TVariant
    :param VariantDataset vds: VDS to compte EM on
    :param list of str va_agg: VA aggregation fields for EM (usually gene)
    :param list of str sa_agg: SA aggregation fields for EM (usually None or population)
    :param num_partitions int: Number of partitions for EM aggregation. If None then the same number of partitions than the input KeyTable is used.
    :param dict of str:str va_to_keep: Variant annotations to keep ({name: expr})
    :return: KeyTable with
    """

    all_variants = {v for vx in kt.query(['{0}.collectAsSet()'.format(k) for k in kt.key]) for v in vx}
    vds = vds.filter_variants_list(all_variants)

    agg_keys = {x: x for x in va_agg}
    if sa_agg is not None:
        vds = vds.annotate_samples_expr('sa = {{{0}}}'.format(",".join([expr for expr in sa_agg])))
        agg_keys.update({x: x for x in sa_agg})
    else:
        vds = vds.annotate_samples_expr('sa = {}')

    if num_partitions is None:
        num_partitions = kt.num_partitions()

    em_kt = vds.phase_em(va_agg, num_partitions,
                                      sa_keys=sa_agg,
                                      variant_pairs= kt)

    return aggregate_em_ref(em_kt,
                            agg_keys,
                            va_to_keep)


def aggregate_em_ref(reference_kt, keys , va_fields=None):
    """

    Aggregates reference after EM.

    :param KeyTable reference_kt: Reference KeyTable to aggregate
    :param dict of str:str keys: Aggregation key (usually v1, v2 and pop) names and expr
    :param dict of str:str va_fields: Dict with keys being the target column name and values the variant annotation to get. Note that they do not need to be preceded by 'va'
    :return: Aggregated KeyTable
    :rtype: KeyTable
    """

    if va_fields is None:
        va_fields = {}

    return reference_kt.aggregate_by_key(['{} = {}'.format(k, v) for k, v in keys.iteritems()],
                                         [
                                             'haplotype_counts = haplotype_counts.takeBy(x => isMissing(x).toInt,1)[0]',
                                             'genotype_counts = genotype_counts.takeBy(x => isMissing(x).toInt,1)[0]',
                                             'prob_same_haplotype = prob_same_haplotype.takeBy(x => isMissing(x).toInt,1)[0]'
                                         ] +
                                         ['{0}{1} = va{1}.takeBy(x => isMissing(x).toInt,1)[0].{2}'.format(name, num,
                                                                                                        field[3:] if field.startswith("va.") else field) for num
                                          in [1, 2] for name, field in va_fields.iteritems()])

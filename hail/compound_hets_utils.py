from gnomad_hail import *

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


def flatten_genotypes(kt, gt_columns, out_prefix=""):
    """

    Flattens genotype counts from the columns specified.
    Example flatten_genotypes(kt, ['g'], out_prefix = "dad_") will produce columns
    dad_g_gt, dad_g_dq, dad_g_gq, dad_g_ad0, dad_g_ad0, dad_g_ad1, dad_g_pl0, dad_g_pl1, dad_g_pl2

    Currently only supports bi-allelic genotypes

    :param KeyTable kt: input KT
    :param list of str gt_columns: Columns of type TGenotype to flatten
    :return: Keytable with flatten columns, name of columns added
    :rtype: (KeyTable, list of str)
    """
    gt_expr = {}
    for c in gt_columns:
        gt_expr['{0}{1}_gt'.format(out_prefix, c)] = '{}.gt'.format(c)
        gt_expr['{0}{1}_dp'.format(out_prefix, c)] = '{}.dp'.format(c)
        gt_expr['{0}{1}_gq'.format(out_prefix, c)] = '{}.gq'.format(c)
        gt_expr['{0}{1}_ad0'.format(out_prefix, c)] ='{}.ad[0]'.format(c)
        gt_expr['{0}{1}_ad1'.format(out_prefix, c)] = '{}.ad[1]'.format(c)
        gt_expr['{0}{1}_pl0'.format(out_prefix, c)] = '{}.pl[0]'.format(c)
        gt_expr['{0}{1}_pl1'.format(out_prefix, c)] = '{}.pl[1]'.format(c)
        gt_expr['{0}{1}_pl2'.format(out_prefix, c)] = '{}.pl[2]'.format(c)

    return kt.annotate(['{0} = {1}'.format(ann, expr) for ann, expr in gt_expr.iteritems()]), gt_expr.keys()


def flatten_variant(kt, variant_col, output_prefix="", output_suffix=""):
    """

    Flattens a column of type TVariant into chrom, pos, ref, alt.
    Currently only support bi-allelic variants

    :param KeyTable kt: Input KT
    :param str variant_col: Variant column
    :param str output_prefix: Output columns prefix
    :param str output_suffix: Output columns suffix
    :return: KeyTable with added flatten columns, list of columns added
    :rtype: (KeyTable, list of str)
    """

    expr = {
        '{}chrom{}'.format(output_prefix, output_suffix): '{}.contig'.format(variant_col),
        '{}pos{}'.format(output_prefix, output_suffix): '{}.start'.format(variant_col),
        '{}ref{}'.format(output_prefix, output_suffix): '{}.ref'.format(variant_col),
        '{}alt{}'.format(output_prefix, output_suffix): '{}.alt'.format(variant_col)
    }

    return kt.annotate(['{} = {}'.format(k,v) for k,v in expr.iteritems()]), expr.keys()



def flatten_haplotype_counts(kt, gc_col="genotype_counts", hc_col="haplotype_counts", out_prefix="", numeric_naming=False):
    """

    Flattens genotype pairs and haplotype counts.

    :param KeyTable kt: Input KT
    :param str gc_col: Column containing genotype pairs Array
    :param str hc_col: Columns containing haplotype counts Array
    :param str out_prefix: output column prefix
    :param bool numeric_naming: If set, then alleles are coded `0/1` instead of `AaBb`
    :return: Keytable with flatten counts and list of columns added to KT
    :rtype: (KeyTable, list of str)
    """
    if not gc_col:
        gc_cols = []
    elif numeric_naming:
        gc_cols = ['0000', '0100', '1100', '0001', '0101', '1101', '0011', '0111', '1111']
    else:
        gc_cols = ['AABB', 'AaBB', 'aaBB', 'AABb', 'AaBb', 'aaBb', 'AAbb', 'Aabb', 'aabb']

    if not hc_col:
        hc_cols = []
    elif numeric_naming:
        hc_cols = ['00', '01', '10', '11']
    else:
        hc_cols = ['AB', 'Ab', 'aB', 'ab']

    return (
        kt.annotate(['{0}{1} = {2}[{3}]'.format(out_prefix, gt, gc_col, gc_cols.index(gt)) for gt in gc_cols] +
                    ['{0}{1} = {2}[{3}]'.format(out_prefix, hp, hc_col, hc_cols.index(hp)) for hp in hc_cols]),
        [out_prefix + gt for gt in gc_cols] + [out_prefix + hp for hp in hc_cols]
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


def conditional_column_swap(kt, swap_expr, columns, gt_counts_col=None, hc_counts_col=None):
    """

    Swaps values of columns of a KeyTable if a swap condition is met.

    :param KeyTable kt: Input KT
    :param str swap_expr: Expression returning a boolean that is evaluated to swap columns or not
    :param list of (str, str) columns: Ordered corresponding columns (e.g. [('v1','v2'),('va1','va2')])
    :param str gt_counts_col: Columns with genotype counts to reorder when swapping variants in pairs
    :param str hc_counts_col: Columns with haplotype counts to reorder when swapping variants in pairs
    :return: KT with swapped columns at rows where condition is met
    :rtype: KeyTable
    """

    kt = kt.annotate(['swap_columns = {}'.format(swap_expr)])

    annotations = ["{0} = if(swap_columns) {1} else {0}".format(col[0], col[1]) for col in columns]
    annotations.extend(["{1} = if(swap_columns) {0} else {1}".format(col[0], col[1]) for col in columns])

    if gt_counts_col is not None:
        annotations.extend(['{0} = if(swap_columns) [{0}[0], {0}[3], {0}[6], {0}[1], {0}[4], {0}[7], {0}[2], {0}[5], {0}[8]] else {0}'.format(gt_counts_col)])

    if hc_counts_col is not None:
        annotations.extend(['{0} = if(swap_columns) [{0}[0], {0}[2], {0}[1], {0}[3]] else {0}'.format(hc_counts_col)])

    kt = kt.annotate(annotations)

    return kt.drop(['swap_columns'])


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

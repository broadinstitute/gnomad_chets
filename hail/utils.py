__author__ = 'konrad'


POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']


def magic(hc):
    from pyhail.java import scala_object
    scala_object(hc.jvm.org.apache.spark.deploy, 'SparkHadoopUtil').get().conf().setLong("parquet.block.size", 1099511627776L)

def write_interval_files(file_path):
    CHROMS = map(str, range(1, 23))
    #CHROMS.extend(['X', 'Y'])

    # Prepare some intervals files
    with open('%s/chrX.txt' % file_path, 'w') as f:
        f.write('X:1-1000000000')
    with open('%s/intervals/chrY.txt' % file_path, 'w') as f:
        f.write('Y:1-1000000000')
    with open('%s/intervals/autosomes.txt' % file_path, 'w') as f:
        f.write('\n'.join(['%s:1-1000000000' % x for x in CHROMS]))


def popmax_text(pops, skip_other=True):
    if skip_other:
        if 'oth' in pops: pops.remove('oth')
        if 'OTH' in pops: pops.remove('OTH')
    af_pops = ','.join(['va.info.AF_%s' % pop for pop in pops])
    skip_other_text = '.filter(x => x != "oth" && x != "OTH")' if skip_other else ''

    get_af_max = 'va.AF_max = let af = [%s] and pops = global.pops%s in range(v.nAltAlleles).map(a => range(pops.size).sortBy(x => af[x][a],false)[0])' % (af_pops, skip_other_text)

    command = []
    for pop in pops:
        this_pop = '''if(pops[va.AF_max[a]] == "%(pop_lower)s" && va.info.AC_%(pop)s[a]>0)
        {pop: "%(pop)s", AF: va.info.AF_%(pop)s[a], AC: va.info.AC_%(pop)s[a], AN: va.info.AN_%(pop)s}''' % {'pop': pop, 'pop_lower': pop.lower()}
        command.append(this_pop)

    command.append('na')
    get_popmax = """va.popmax = let na = NA : Struct{ pop: String, AF: Double, AC: Int, AN: Int} and pops = global.pops%s in
        range(va.AF_max.size).map(a => %s)""" % (skip_other_text, '\n else'.join(command))

    extract_popmax = """va.info.POPMAX = va.popmax.map(x => x.pop), va.info.AC_POPMAX = va.popmax.map(x => x.AC),
        va.info.AN_POPMAX = va.popmax.map(x => x.AN), va.info.AF_POPMAX = va.popmax.map(x => x.AF)"""

    return get_af_max, get_popmax, extract_popmax


def get_hom_from_gc(destination, target):
    return '%s = range(v.nAltAlleles).map(i => let n = i + 2 in %s[(n * (n + 1) / 2).toInt - 1])' % (destination, target)


def unfurl_hom_text(pops, simple_hom=True):
    expressions = [get_hom_from_gc('va.info.Hom_%s' % pop, 'va.info.GC_%s' % pop) for pop in pops]
    if simple_hom: expressions.append('va.info.Hom = range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC[(n * (n + 1) / 2).toInt - 1])')
    return ',\n'.join(expressions)


def unfurl_callstats_text(criteria_pops, lower=True, gc=True, additional_annotations=None):
    expression = []
    for criterion, pop in criteria_pops:
        input_dict = {'pop': pop.lower() if lower else pop, 'pop_upper': pop, 'criterion': criterion}
        expression.append('va.calldata.%(pop_upper)s = gs.filter(g => %(criterion)s == "%(pop)s").callStats(g => v)' % input_dict)
    callstats_command = ',\n'.join(expression)

    expression = []
    metrics = ['AC', 'AN', 'AF']
    if gc: metrics.append('GC')
    for metric in metrics:
        end = '[1:]' if metric in ('AC', 'AF') else ''
        for _, pop in criteria_pops:
            input_dict = {'pop': pop.lower() if lower else pop, 'pop_upper': pop, 'metric': metric, 'end': end}
            expression.append('va.info.%(metric)s_%(pop_upper)s = va.calldata.%(pop_upper)s.%(metric)s%(end)s' % input_dict)
    right_shift_command = ',\n'.join(expression)

    return callstats_command, right_shift_command


def cut_allele_from_g_array(target, destination=None):
    if destination is None: destination = target
    return ('%s = let removed_alleles = range(1, v.nAltAlleles + 1).filter(i => !aIndices.toSet.contains(i)).toSet in\n'
            'range(%s.size).filter(i => !removed_alleles.contains(gtj(i)) && !removed_alleles.contains(gtk(i)))\n'
            '.map(i => %s[i])' % (destination, target, target))


def unfurl_filter_alleles_annotation(a_based=None, r_based=None, g_based=None, additional_annotations=None):

    annotations = []
    if r_based:
        for ann in r_based:
            annotations.append('%s = aIndices.map(i => %s[i])' % (ann, ann))

    if a_based:
        for ann in a_based:
            annotations.append('%s = aIndices[1:].map(i => %s[i - 1])' % (ann, ann))

    if g_based:
        for ann in g_based:
            annotations.append(cut_allele_from_g_array(ann))

    if additional_annotations:
        if type(additional_annotations) == str:
            annotations.append(additional_annotations)
        else:
            annotations.extend(additional_annotations)

    return ',\n'.join(annotations)


def konrad_special_text(destination, template, na_struct='hist', reference=True):
    """

    :param str destination: Variant annotation to write to. For instance, va.info.hist
    :param str template: Basic command with placeholder for allele. For instance: gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.gq).hist(0, 100, 20)
    :param str na_struct: One of 'hist', 'stats', or a custom format for the NA structure
    :param bool reference: Whether the reference allele should be included (i.e. Number=R)
    :return: Unfurled command
    """
    post_process = ''
    if na_struct == 'hist':
        na_struct = "NA: Struct { binEdges: Array[Double], binFrequencies: Array[Long], nLess: Long, nGreater: Long }"
        post_process = '.map(x => x.binFrequencies.map(y => str(y)).mkString("|"))'
    elif na_struct == 'stats':
        na_struct = 'NA: Struct { mean: Double, stdev: Double, min: Double, max: Double, nNotMissing: Long, sum: Double }'
    elif na_struct == 'counter':
        na_struct = 'NA: Array[ Struct { key: String, count: Long } ]'

    if reference:
        cut = 'v.nAlleles'
        start = 0
    else:
        cut = 'v.nAltAlleles'
        start = 1

    template = template.replace('%s', '%(allele)s')
    command = []
    for i in range(start, 7):
        allele_command = template % {'allele': i}
        if i > 1:
            command.append('if (v.nAltAlleles > %s) %s else na' % (i - 1, allele_command))
        else:
            command.append(allele_command)
    command = ',\n'.join(command)

    full_command_text = """%(destination)s = let na = %(na_struct)s in [
            %(command)s
        ][:%(cut)s]%(post_process)s""" % {'destination': destination, 'na_struct': na_struct, 'command': command, 'cut': cut, 'post_process': post_process}

    return full_command_text


import pyhail.dataset
import pyhail.context
from pyhail.java import jarray
import pyspark.sql
import json


class VariantDataset(pyhail.dataset.VariantDataset):
    def konrad_special(self, destination, template, na_struct='hist', reference=True):
        return self.annotate_variants_expr(konrad_special_text(destination, template, na_struct, reference))

    def unfurl_callstats(self, pops, lower=False, gc=True):
        callstats_command, right_shift_command = unfurl_callstats_text(pops, lower, gc)
        return (self.annotate_variants_expr(callstats_command)
                .annotate_variants_expr(right_shift_command))

    def unfurl_hom(self, pops, simple_hom=True):
        hom_command = unfurl_hom_text(pops, simple_hom)
        return self.annotate_variants_expr(hom_command)

    def popmax(self, pops, skip_other=True):
        get_af_max, get_popmax, extract_popmax = popmax_text(pops, skip_other)
        return (self.annotate_variants_expr(get_af_max)
                .annotate_variants_expr(get_popmax)
                .annotate_variants_expr(extract_popmax))

    def projectmax(self):
        return ( self.annotate_alleles_expr('let nNonRef = gs.filter(g => g.isCalledNonRef).map(g => sa.meta.project_description).counter() and'
                                   'nSamples = gs.filter(g => g.isCalled).map(g => sa.meta.project_description).counter() in'
                                   'va.projectmax = nNonRef.map(x => {key: x.key, count: x.count, nsamples: nSamples.find(y => x.key == y.key).count}).sortBy(x =>x.count / x.nsamples,false)[0:5]')
                 .annotate_variants_expr('va.info.PROJECTMAX = va.projectmax.map(a => a.map(x => x.key).mkString("|")), ' \
                            'va.info.PROJECTMAX_NSamples = va.projectmax.map(a => a.map(x => str(x.nsamples)).mkString("|")), ' \
                            'va.info.PROJECTMAX_NonRefSamples = va.projectmax.map(a => a.map(x => str(x.count)).mkString("|")), ' \
                            'va.info.PROJECTMAX_PropNonRefSamples = va.projectmax.map(a => a.map(x => str(x.count / x.nsamples)).mkString("|"))')
                 )

    def filter_to_adj(self):
        return self.filter_genotypes('g.gq >= 20 && g.dp >= 10 && (!g.isHet || (g.gtj > 0 || g.ad[g.gtk]/g.dp > 0.2))') ##Assumes gtj <= gtk

    def filter_star(self, a_based=None, r_based=None, g_based=None, additional_annotations=None):
        annotation = unfurl_filter_alleles_annotation(a_based=a_based, r_based=r_based, g_based=g_based, additional_annotations=additional_annotations)
        return self.filter_alleles('v.altAlleles[aIndex - 1].alt == "*"', annotation=annotation, keep=False)

    def head(self):
        return json.loads(pyspark.sql.DataFrame(self.jvds.variantsDF(self.hc.jsql_context), self.hc.sql_context).toJSON().first())

    def remove_filter_status(self, criteria):
        return self.annotate_variants_expr('')


# class HailContext(pyhail.context.HailContext):
#     def run_command(self, vds, pargs):
#         jargs = jarray(self.gateway, self.jvm.java.lang.String, pargs)
#         t = self.jvm.org.broadinstitute.hail.driver.ToplevelCommands.lookup(jargs)
#         cmd = t._1()
#         cmd_args = t._2()
#         result = cmd.run(self._jstate(vds.jvds if vds != None else None),
#                          cmd_args)
#         return VariantDataset(self, result.vds())

def annotate_non_split_from_split(hc, non_split_vds_path, split_vds, annotations, annotation_exp_out_path):

    variant_annotated_vds = (
        hc.read(non_split_vds_path,sites_only=True)
        .annotate_variants_expr('va.variant = str(v)')
        .split_multi()
    )

    ann_agg_codes = ["`%s` = g.map(g => %s).collect()" % (a,a) for a in annotations]
    print(",".join(ann_agg_codes))

    x = (
        split_vds
            .annotate_variants_vds(variant_annotated_vds, 'va.variant = vds.variant')
            .filter_variants_expr('isDefined(va.variant)')
     )

    print(x.count())

    agg = x.aggregate_by_key(key_code='variant = va.variant', agg_code = ",".join(ann_agg_codes))
    print(agg.nrows())
    agg.export(output=annotation_exp_out_path, types_file=annotation_exp_out_path + '.types')

    ann_codes = ['%s = table.`%s`' % (a,a) for a in annotations]
    return(
        hc.read(non_split_vds_path)
        .annotate_variants_table(annotation_exp_out_path,'Variant(variant)',code=",".join(ann_codes),
                                 config= pyhail.TextTableConfig(types=annotation_exp_out_path + '.types'))
    )



def get_variant_type_expr(code="va.variantType"):
    return(['''%s =
    let non_star = v.altAlleles.filter(a => a.alt != "*") in
        if (non_star.forall(a => a.isSNP))
            if (non_star.length > 1)
                "multi-snv"
            else
                "snv"
        else if (non_star.forall(a => a.isIndel))
            if (non_star.length > 1)
                "multi-indel"
            else
                "indel"
        else
            "mixed"''' % code])

def get_hists_expr(root="va.hists"):
    hists = [
        '%s.GQ_HIST = gs.map(g => g.gq).hist(0, 100, 20)',
        '%s.DP_HIST = gs.map(g => g.dp).hist(0, 100, 20)',
        '%s.AB_HIST = gs.map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)'
    ]
    return(
        [x % root for x in hists]
    )

def get_stats_expr(root="va.stats", medians=False):
    stats = ['%s.gq = gs.filter(g => g.isCalledNonRef).map(g => g.gq).stats()',
             '%s.dp = gs.filter(g => g.isCalledNonRef).map(g => g.dp).stats()',
             '%s.nrq = gs.filter(g => g.isCalledNonRef).map(g => g.dosage[0]).stats()',
             '%s.ab = gs.filter(g => g.isCalledNonRef).map(g => g.ad[1]/g.dp).stats()']

    stats_expr = [x % root for x in stats]

    if(medians):
        template = (
            '%(destination)s = let sorted_vals = gs.filter(g => g.isCalledNonRef && !isMissing(%(metric)s)).map(g => %(metric)s).collect().sort() in '
            'if (sorted_vals.size == 0) NA: Double else '
            'if (sorted_vals.size %% 2 == 1) sorted_vals[(sorted_vals.size/2).toInt] else '
            '(sorted_vals[(sorted_vals.size/2).toInt] + sorted_vals[(sorted_vals.size/2).toInt - 1])/2.0')
        medians = [('g.gq', '%s.gq_median' % root), ('g.dp', '%s.dp_median' % root),
                   ('g.dosage[0]', '%s.nrq_median' % root), ('g.ad[1]/g.dp', '%s.ab_median' % root)]
        stats_expr.extend(
            [template % {'metric': metric, 'destination': destination} for (metric, destination) in medians])

    return stats_expr

def get_add_filter_annotation(filtername, filterexpr):
    return ('va.filters = \n'
    'if (%s)\n'
        'if (va.filters.contains("PASS"))\n'
            '["%s"].toSet\n'
        'else\n'
            '[va.filters, ["%s"].toSet].toSet.flatten()\n'
    'else\n'
        'va.filters' % (filterexpr, filtername, filtername)
    )

def create_sites_vds_annotations(vds, pops, dbsnp_path=None, npartitions=1000, shuffle=True):
    sexes = ['Male', 'Female']
    cuts = pops
    cuts.extend(sexes)

    g_based_annotations = ['va.info.GC', 'va.info.GC_Adj']
    g_based_annotations.extend(['va.info.GC_%s' % x for x in cuts])
    a_based_annotations = ['va.info.AC', 'va.info.AC_Adj', 'va.info.AF', 'va.info.AF_Adj']
    a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                                'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
    a_based_annotations.extend(['va.info.AC_%s' % x for x in cuts])
    a_based_annotations.extend(['va.info.AF_%s' % x for x in cuts])
    r_based_annotations = ['va.info.GQ_HIST', 'va.info.DP_HIST', 'va.info.AB_HIST']

    criterion_pops = [('sa.meta.population', x) for x in pops]
    criterion_pops.extend([('sa.meta.sex', x) for x in sexes])

    star_annotations = ['va.info.STAR_AC = va.info.AC[aInddices[0]]']

    ##TODO: Add dbSNP annotation once --comment option is available
    # if(dbsnp is not None):
    #     vds = vds.annotate_variants_loci(dbsnp_path,
    #                                      locus_expr='Locus(_0,_1)',
    #                                      )

    return (vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
            .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
            .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
            .annotate_alleles_expr(get_hists_expr('va.info'))
            .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
            .filter_to_adj()
            .projectmax()
            .annotate_variants_expr('va.calldata.Adj = gs.callStats(g => v)')
            .unfurl_callstats(criterion_pops, lower=True)
            .filter_samples_all()
            .annotate_variants_expr('va.info.AC = va.calldata.raw.AC[1:], '
                                    'va.info.AN = va.calldata.raw.AN, '
                                    'va.info.AF = va.calldata.raw.AF[1:], '
                                    'va.info.GC = va.calldata.raw.GC, '
                                    'va.info.AC_Adj = va.calldata.Adj.AC[1:], '
                                    'va.info.AN_Adj = va.calldata.Adj.AN, '
                                    'va.info.AF_Adj = va.calldata.Adj.AF[1:], '
                                    'va.info.GC_Adj = va.calldata.Adj.GC')
            .unfurl_hom(cuts, simple_hom=True)
            .filter_star(a_based_annotations, r_based_annotations, g_based_annotations, additional_annotations=star_annotations)
            .popmax(pops)
            .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
            .repartition(npartitions, shuffle=shuffle)
            )

    #TODO:
    #AB_HIST before or after Adj? Was after -- now is before
    #removed the following since I think it's taken care of by .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False):
    #
    ##TODO: Needs fixing
    ##                  .annotate_variants_loci('%s/sites/exac2.RF.allsites.txt.bgz' % root, 'Locus(chrom, pos)',
     ##                                    code='va.info.RF_PROB = table.rfprob, va.info.RF_TYPE = table.type', impute=True)
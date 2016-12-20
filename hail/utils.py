__author__ = 'konrad'
import re

from subprocess import check_output

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
        return ( self.annotate_alleles_expr('va.projectmax = let nNonRef = gs.filter(g => g.isCalledNonRef).map(g => sa.meta.project_description).counter() and '
                                   'nSamples = gs.filter(g => g.isCalled).map(g => sa.meta.project_description).counter() in '
                                   'nNonRef.map(x => {key: x.key, count: x.count, nsamples: nSamples.find(y => x.key == y.key).count}).sortBy(x =>x.count / x.nsamples,false)[0:5]')
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

    def histograms(self, root='va.info'):
        return (
            self.annotate_alleles_expr(
                ['%s.GQ_HIST_ALT = gs.filter(g => g.isCalledNonRef).map(g => g.gq).hist(0, 100, 20)' % root,
                 '%s.DP_HIST_ALT = gs.filter(g => g.isCalledNonRef).map(g => g.dp).hist(0, 100, 20)' % root,
                 '%s.AB_HIST_ALT = gs.filter(g => g.isCalledNonRef).map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)' % root])
                .annotate_variants_expr(
                ['%s.GQ_HIST_ALL = gs.map(g => g.gq).hist(0, 100, 20)' % root,
                 '%s.DP_HIST_ALL = gs.map(g => g.dp).hist(0, 100, 20)' % root,
                 '%s.AB_HIST_ALL = gs.filter(g => g.isCalledNonRef).map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)' % root])
        )


class HailContext(pyhail.context.HailContext):
    def run_command(self, vds, pargs):
        jargs = jarray(self.gateway, self.jvm.java.lang.String, pargs)
        t = self.jvm.org.broadinstitute.hail.driver.ToplevelCommands.lookup(jargs)
        cmd = t._1()
        cmd_args = t._2()
        result = cmd.run(self._jstate(vds.jvds if vds != None else None),
                         cmd_args)
        return VariantDataset(self, result.vds())

def annotate_non_split_from_split(hc, non_split_vds_path, split_vds, annotations, annotation_exp_out_path):

    variant_annotated_vds = (
        hc.read(non_split_vds_path, sites_only=True)
        .annotate_variants_expr('va.variant = str(v)')
        .split_multi()
    )

    ann_agg_codes = ["`%s` = index(va.map(x => {val: %s, aIndex: str(va.aIndex)}).collect(), aIndex)" % (a,a) for a in annotations]
    agg = (
        split_vds
            .annotate_variants_vds(variant_annotated_vds, 'va.variant = vds.variant')
            .filter_variants_expr('isDefined(va.variant)')
            .variants_keytable()
            .aggregate_by_key(key_code='variant = va.variant', agg_code=",".join(ann_agg_codes))

     )

    agg.export(output=annotation_exp_out_path, types_file=annotation_exp_out_path + '.types')

    schema_command = ['gsutil'] if(annotation_exp_out_path.startswith('gs://')) else []
    schema_command.extend(['cat', annotation_exp_out_path + '.types'])
    schema = check_output(schema_command).decode('ascii')

    #Get types from schema
    #types are named with names surrounded by backticks, comma delimited and nested in a Dict as follows where
    # t1 and t2 are the types and n1 and n2 their respective names
    # `n1` = Dict[Struct{val: t1}], `n2` = Dict[Struct{val: t2}]
    types = list(map(lambda x: x[:-2],re.split(",`[^`]+`:Dict\[Struct\{val:",schema)))[1:]

    ann_codes = ['%s = table.`%s`' % (a,a) for a in annotations]
    sort_ann = ['%s = range(v.nAltAlleles).map(i => if(%s.contains(str(i+1))) %s[str(i+1)].val else NA: %s)' % (a, a, a, b)
                for (a,b) in zip(annotations,types)]

    return(
        hc.read(non_split_vds_path)
        .annotate_variants_table(annotation_exp_out_path,'Variant(variant)',code=",".join(ann_codes),
                                 config= pyhail.TextTableConfig(types=schema))
        .annotate_variants_expr(sort_ann)
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


def get_stats_expr(root="va.stats", medians=False):
    stats = ['%s.gq = gs.filter(g => g.isCalledNonRef).map(g => g.gq).stats()',
             '%s.dp = gs.filter(g => g.isCalledNonRef).map(g => g.dp).stats()',
             '%s.nrq = gs.filter(g => g.isCalledNonRef).map(g => g.dosage[0]).stats()',
             '%s.ab = gs.filter(g => g.isHet).map(g => g.ad[1]/g.dp).stats()']

    stats_expr = [x % root for x in stats]

    if medians:
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


def create_sites_vds_annotations(vds, pops, tmp_path="/tmp", dbsnp_path=None, npartitions=1000, shuffle=True):

    auto_intervals_path = '%s/autosomes.txt' % tmp_path
    with open(auto_intervals_path, 'w') as f:
        f.write('\n'.join(['%s:1-1000000000' % x for x in map(str, range(1, 23))]))

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
    a_based_annotations.extend(['va.info.GQ_HIST_ALT', 'va.info.DP_HIST_ALT', 'va.info.AB_HIST_ALT'])

    criterion_pops = [('sa.meta.population', x) for x in pops]
    criterion_pops.extend([('sa.meta.sex', x) for x in sexes])

    star_annotations = ['va.info.STAR_%s = let removed_allele = range(1, v.nAltAlleles + 1).find(i => !aIndices.toSet.contains(i)) \n' \
                        'in if(isDefined(removed_allele)) va.info.%s[removed_allele] else NA: Int' % (a, a) for a in ['AC', 'AC_Adj', 'Hom']]

    vds =  vds.filter_variants_intervals('file://' + auto_intervals_path)

    if(dbsnp_path is not None):
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code = 'va.rsid=table._2',
                                         config=pyhail.TextTableConfig(noheader=True,comment="#",types='_0: String, _1: Int')
                                         )

    return (vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
            .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
            .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
            .histograms('va.info')
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
            .filter_star(a_based_annotations, None, g_based_annotations,
                         additional_annotations=star_annotations)
            .popmax(pops)
            .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
            .repartition(npartitions, shuffle=shuffle)
            )


def create_sites_vds_annotationsX(vds, pops, tmp_path="/tmp", dbsnp_path=None, npartitions=1000, shuffle=True):
    x_intervals = '%s/chrX.txt' % tmp_path
    with open(x_intervals, 'w') as f:
        f.write('X:1-1000000000')

    sexes = ['Male', 'Female']
    cuts = pops
    cuts.extend(sexes)

    g_based_annotations = ['va.info.GC', 'va.info.GC_Adj']
    g_based_annotations.extend(['va.info.GC_%s' % x for x in sexes])
    g_based_annotations.extend(['va.info.GC_%s_%s' % (y, x) for x in sexes for y in pops])

    a_based_annotations = ['va.info.AC', 'va.info.AC_Adj', 'va.info.AF', 'va.info.AF_Adj']
    a_based_annotations.extend(['va.info.AC_Male', 'va.info.AC_Female', 'va.info.AF_Male', 'va.info.AF_Female'])
    a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                                'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
    a_based_annotations.extend(['va.info.AC_%s_%s' % (y, x) for x in sexes for y in pops])
    a_based_annotations.extend(['va.info.AF_%s_%s' % (y, x) for x in sexes for y in pops])

    a_based_annotations.extend(['va.info.GQ_HIST_ALT', 'va.info.DP_HIST_ALT', 'va.info.AB_HIST_ALT'])

    generate_callstats_expression = []
    for sex in sexes:
        for pop in pops:
            input_dict = {'sex': sex.lower(), 'sex_label': sex.capitalize(), 'pop': pop.lower(),
                          'pop_upper': pop.upper()}
            generate_callstats_expression.append(
                'va.calldata.%(pop_upper)s_%(sex_label)s = gs.filter(g => sa.meta.population == "%(pop)s" && sa.meta.sex == "%(sex)s").callStats(v)' % input_dict)
        input_dict = {'sex': sex.lower(), 'sex_label': sex.capitalize()}
        generate_callstats_expression.append(
            'va.calldata.%(sex_label)s = gs.filter(g => sa.meta.sex == "%(sex)s").callStats(v)' % input_dict)
    generate_callstats_expression = ',\n'.join(generate_callstats_expression)

    rearrange_callstats_expression = []
    for metric in ['AC', 'AN', 'AF', 'GC']:
        for sex in sexes:
            start = ''
            end = '[1:]' if metric in ('AC', 'AF') else ''
            if sex == 'male' and metric in ('AC', 'AN'):
                start = 'let div_factor = if (v.inXNonPar) 2 else 1 in ('
                end += '/div_factor)'
                end += '.map(x => x.toInt)' if metric == 'AC' else '.toInt'
            for pop in pops:
                input_dict = {'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric, 'pop': pop,
                              'pop_upper': pop.upper(), 'start': start, 'end': end}
                rearrange_callstats_expression.append(
                    'va.info.%(metric)s_%(pop_upper)s_%(sex_label)s = %(start)sva.calldata.%(pop_upper)s_%(sex_label)s.%(metric)s%(end)s' % input_dict)
            input_dict = {'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric, 'start': start, 'end': end}
            rearrange_callstats_expression.append(
                'va.info.%(metric)s_%(sex_label)s = %(start)sva.calldata.%(sex_label)s.%(metric)s%(end)s' % input_dict)

    rearrange_callstats_expression.extend(['va.info.GC = va.calldata.raw.GC',
                                           'va.info.AC = (va.calldata.raw.AC - (va.calldata.hemi_raw.AC/2)).map(x => x.toInt)[1:]',
                                           'va.info.AN = (va.calldata.raw.AN - (va.calldata.hemi_raw.AN/2)).toInt',
                                           'va.info.GC_Adj = va.calldata.Adj.GC',
                                           'va.info.AC_Adj = (va.calldata.Adj.AC - (va.calldata.Hemi_Adj.AC/2)).map(x => x.toInt)[1:]',
                                           'va.info.AN_Adj = (va.calldata.Adj.AN - (va.calldata.Hemi_Adj.AN/2)).toInt'])
    rearrange_callstats_expression = ',\n'.join(rearrange_callstats_expression)

    ac_an_expression = []
    for metric in ['AC', 'AN']:
        for pop in pops:
            input_dict = {'metric': metric, 'pop': pop, 'pop_upper': pop.upper()}
            ac_an_expression.append(
                'va.info.%(metric)s_%(pop_upper)s = va.info.%(metric)s_%(pop_upper)s_Male + va.info.%(metric)s_%(pop_upper)s_Female' % input_dict)
    ac_an_expression = ',\n'.join(ac_an_expression)

    af_expression = []
    for pop in pops:
        input_dict = {'pop': pop, 'pop_upper': pop.upper()}
        af_expression.append(
            'va.info.AF_%(pop_upper)s = va.info.AC_%(pop_upper)s / va.info.AN_%(pop_upper)s' % input_dict)
    af_expression = ',\n'.join(af_expression)

    hom_hemi_expression = []
    for sex in sexes:
        metric = 'Hom' if sex == 'female' else 'Hemi'
        for pop in pops:
            input_dict = {'pop': pop, 'pop_upper': pop.upper(), 'sex': sex, 'sex_label': sex.capitalize(),
                          'metric': metric}
            hom_hemi_expression.append(
                'va.info.%(metric)s_%(pop_upper)s = if (v.inXNonPar) range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(pop_upper)s_%(sex_label)s[(n * (n + 1) / 2).toInt - 1]) else NA: Array[Int]' % input_dict)
        input_dict = {'sex': sex, 'sex_label': sex.capitalize(), 'metric': metric}
        hom_hemi_expression.append(
            'va.info.%(metric)s = if (v.inXNonPar) range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(sex_label)s[(n * (n + 1) / 2).toInt - 1]) else NA: Array[Int]' % input_dict)
    hom_hemi_expression = ',\n'.join(hom_hemi_expression)

    star_annotations = [
        'va.info.STAR_%s = let removed_allele = range(1, v.nAltAlleles + 1).find(i => !aIndices.toSet.contains(i)) \n' \
        'in if(isDefined(removed_allele)) va.info.%s[removed_allele] else NA: Int' % (a, a) for a in
        ['AC', 'AC_Adj', 'Hom','Hemi']]

    if (dbsnp_path is not None):
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=pyhail.TextTableConfig(noheader=True, comment="#",
                                                                       types='_0: String, _1: Int')
                                         )

    return (vds.filter_variants_intervals('file://' + x_intervals)
            .filter_genotypes('sa.meta.sex == "male" && g.isHet && v.inXNonPar', keep=False)
            .annotate_variants_expr('va.calldata.raw = gs.callStats(v)')
            .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
            .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
            .histograms('va.info')
            .annotate_variants_expr
            .annotate_variants_expr('va.calldata.raw = gs.callStats(v), '
                                    'va.calldata.hemi_raw = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(v)')
            .filter_to_adj()
            .projectmax()
            .annotate_variants_expr('va.calldata.Adj = gs.callStats(v), '
                                    'va.calldata.Hemi_Adj = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(v)')
            .annotate_variants_expr(generate_callstats_expression)
            .filter_samples_all()
            .annotate_variants_expr(rearrange_callstats_expression)
            .annotate_variants_expr(
        'va.info.AF = va.info.AC.map(x => if (va.info.AN > 0) x.toDouble/va.info.AN else NA: Double), '
        'va.info.AF_Adj = va.info.AC_Adj.map(x => if (va.info.AN_Adj > 0) x.toDouble/va.info.AN_Adj else NA: Double)')  # Got here
            .annotate_variants_expr(hom_hemi_expression)
            .annotate_variants_expr(ac_an_expression)
            .annotate_variants_expr(af_expression)
            .filter_star(a_based_annotations, None, g_based_annotations,
                         additional_annotations=star_annotations)
            .popmax(pops)
            .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
            )


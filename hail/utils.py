__author__ = 'konrad'


POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']


def magic(hc):
    from pyhail.java import scala_object
    scala_object(hc.jvm.org.apache.spark.deploy, 'SparkHadoopUtil').get().conf().setLong("parquet.block.size", 1099511627776L)


def popmax_text(pops=POPS, skip_other=True):
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


def projectmax_text():
    get_nonref_samples_by_project = konrad_special_text('va.nonref_samples_by_project',
                                                        'gs.filter(g => (g.gtj == %s || g.gtk == %s)).map(g => sa.project_or_cohort).counter()',
                                                        na_struct='counter', reference=False)
    flatten_proejcts = 'va.projects = va.nonref_samples_by_project.flatMap(a => a.map(x => x.key)).toSet'
    get_samples_by_project_command = 'va.samples_by_project = gs' \
                                     '.filter(g => g.isCalled && va.projects.contains(sa.project_or_cohort))' \
                                     '.map(g => sa.project_or_cohort).counter()'
    get_project_max_command = 'va.projectmax = va.nonref_samples_by_project' \
                              '.map(a => a.map(x => {key: x.key, count: x.count, nsamples: va.samples_by_project.find(y => x.key == y.key).count})' \
                              '.sortBy(x =>x.count / x.nsamples,false)[0:5])'
    format_projectmax = 'va.info.PROJECTMAX = va.projectmax.map(a => a.map(x => x.key).mkString("|")), ' \
                        'va.info.PROJECTMAX_NSamples = va.projectmax.map(a => a.map(x => str(x.nsamples)).mkString("|")), ' \
                        'va.info.PROJECTMAX_NonRefSamples = va.projectmax.map(a => a.map(x => str(x.count)).mkString("|")), ' \
                        'va.info.PROJECTMAX_PropNonRefSamples = va.projectmax.map(a => a.map(x => str(x.count / x.nsamples)).mkString("|"))'
    return get_nonref_samples_by_project, flatten_proejcts, get_samples_by_project_command, get_project_max_command, format_projectmax


def get_hom_from_gc(destination, target):
    return '%s = range(v.nAltAlleles).map(i => let n = i + 2 in %s[(n * (n + 1) / 2).toInt - 1])' % (destination, target)


def unfurl_hom_text(pops, simple_hom=True):
    expressions = [get_hom_from_gc('va.info.Hom_%s' % pop, 'va.info.GC_%s' % pop) for pop in pops]
    if simple_hom: expressions.append('va.info.Hom = range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC[(n * (n + 1) / 2).toInt - 1])')
    return ',\n'.join(expressions)


def unfurl_callstats_text(criteria_pops, lower=True, gc=True):
    expression = []
    for criterion, pop in criteria_pops:
        input_dict = {'pop': pop.lower() if lower else pop, 'pop_upper': pop, 'criterion': criterion}
        expression.append('va.calldata.%(pop_upper)s = gs.filter(g => %(criterion)s == "%(pop)s").callStats(v)' % input_dict)
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


def unfurl_filter_alleles_annotation(a_based=None, r_based=None, g_based=None):

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

    full_command_text = """
        %(destination)s = let na = %(na_struct)s in [
            %(command)s
        ][:%(cut)s]%(post_process)s
    """ % {'destination': destination, 'na_struct': na_struct, 'command': command, 'cut': cut, 'post_process': post_process}

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
        get_nonref_samples_by_project, flatten_proejcts, get_samples_by_project_command, get_project_max_command, format_projectmax = projectmax_text()
        return (self.annotate_variants_expr(get_nonref_samples_by_project)
                .annotate_variants_expr(flatten_proejcts)
                .annotate_variants_expr(get_samples_by_project_command)
                .annotate_variants_expr(get_project_max_command)
                .annotate_variants_expr(format_projectmax))

    def filter_to_adj(self):
        return self.filter_genotypes('g.gq >= 20 && g.dp >= 10')

    def filter_star(self, a_based=None, r_based=None, g_based=None):
        annotation = unfurl_filter_alleles_annotation(a_based=a_based, r_based=r_based, g_based=g_based)
        return self.filter_alleles('v.altAlleles[aIndex - 1].alt == "*"', annotation=annotation, keep=False)

    def head(self):
        return json.loads(pyspark.sql.DataFrame(self.jvds.variantsDF(self.hc.jsql_context), self.hc.sql_context).toJSON().first())


class HailContext(pyhail.context.HailContext):
    def run_command(self, vds, pargs):
        jargs = jarray(self.gateway, self.jvm.java.lang.String, pargs)
        t = self.jvm.org.broadinstitute.hail.driver.ToplevelCommands.lookup(jargs)
        cmd = t._1()
        cmd_args = t._2()
        result = cmd.run(self._jstate(vds.jvds if vds != None else None),
                         cmd_args)
        return VariantDataset(self, result.vds())
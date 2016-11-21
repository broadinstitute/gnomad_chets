__author__ = 'konrad'


def magic(hc):
    from pyhail.java import scala_object
    scala_object(hc.jvm.org.apache.spark.deploy, 'SparkHadoopUtil').get().conf().setLong("parquet.block.size", 1099511627776L)


def popmax(vds, pops='AFR,AMR,ASJ,EAS,FIN,NFE,SAS', skip_other=True):

    pops = pops.split(',')
    af_pops = ','.join(['va.info.AF_%s' % pop for pop in pops])
    skip_other_text = '.filter(x => x != "oth")' if skip_other else ''

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

    return vds.annotate_variants_expr(get_af_max)\
        .annotate_variants_expr(get_popmax)\
        .annotate_variants_expr(extract_popmax)


def unfurl_pop_callstats(vds=None, pops='AFR,AMR,ASJ,EAS,FIN,NFE,SAS', lower=True):
    pops = pops.split(',')

    expression = []
    for pop in pops:
        input_dict = {'pop': pop.lower() if lower else pop, 'pop_upper': pop}
        expression.append('va.calldata.%(pop_upper)s = gs.filter(g => sa.meta.population == "%(pop)s").callStats(v)' % input_dict)
    callstats_command = ',\n'.join(expression)

    expression = []
    for metric in ['AC', 'AN', 'AF', 'GC']:
        end = '[1:]' if metric in ('AC', 'AF') else ''
        for pop in pops:
            input_dict = {'pop': pop.lower() if lower else pop, 'pop_upper': pop, 'metric': metric, 'end': end}
            expression.append('va.info.%(metric)s_%(pop_upper)s = va.calldata.%(pop_upper)s.%(metric)s%(end)s' % input_dict)
    right_shift_command = ',\n'.join(expression)

    if vds:
        vds.annotate_variants_expr(callstats_command)\
            .annotate_variants_expr(right_shift_command)
    else:
        return (callstats_command, right_shift_command)


def konrad_special(destination, template, vds=None, na_struct='hist', reference=True):
    """

    :param destination: Variant annotation to write to. For instance, va.info.hist
    :param template: Basic command with placeholder for allele. For instance: gs.filter(g => g.gtj == %s || g.gtk == %s).map(g => g.gq).hist(0, 100, 20)
    :param na_struct: One of 'hist', 'stats', or a custom format for the NA structure
    :param reference: Whether the reference allele should be included (i.e. Number=R)
    :return: Unfurled command
    """
    post_process = ''
    if na_struct == 'hist':
        na_struct = "NA: Struct { binEdges:Array[Double], binFrequencies: Array[Long], nSmaller: Long, nGreater: Long }"
        post_process = '.map(x => x.binFrequencies.map(y => str(y)).mkString("|"))'
    elif na_struct == 'stats':
        na_struct = 'NA: Struct { mean: Double, stdev: Double, min: Double, max: Double, nNotMissing: Long, sum: Double }'

    if reference:
        cut = 'v.nAlleles'
        start = 0
    else:
        cut = 'v.nAltAlleles'
        start = 1

    command = []
    for i in range(start, 7):
        if i > 1:
            command.append('if (v.nAltAlleles > %s) %s' % (i - 1, template % i))
        else:
            command.append(template % i)
    command = ','.join(command)

    full_command_text = """
        %(destination)s = let na = %(na_struct)s in [
            %(command)s
        ][%(cut)s]%(post_process)s
    """ % {'destination': destination, 'na_struct': na_struct, 'command': command, 'cut': cut, 'post_process': post_process}

    if vds:
        return vds.annotate_variants_expr(full_command_text)
    else:
        return full_command_text
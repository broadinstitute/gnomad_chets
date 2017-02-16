__author__ = 'konrad'
import re
import sys
import hail
from hail.java import jarray, raise_py4j_exception
import pyspark.sql
import json
import copy
import time
from py4j.protocol import Py4JJavaError
from subprocess import check_output
from pprint import pprint

from resources import *
from hail.type import *
from hail.representation import *

POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
POP_NAMES = {'AFR': "African/African American",
             'AMR': "Admixed American",
             'ASJ': "Ashkenazi Jewish",
             'EAS': "East Asian",
             'FIN': "Finnish",
             'NFE': "Non-Finnish European",
             'OTH': "Other (population not assigned)",
             'SAS': "South Asian"
             }

SEXES = {
    'Male': 'Male',
    'Female': 'Female'
}

ANNOTATION_DESC = {
    'AC': ('A', 'Allele count in %sgenotypes, for each ALT allele, in the same order as listed'),
    'AF': ('A', 'Allele Frequency among %sgenotypes, for each ALT allele, in the same order as listed'),
    'AN': ('1', 'Total number of alleles in %scalled genotypes'),
    'Hom': ('A', 'Count of homozygous %sindividuals'),
    'Hemi': ('A', 'Count of hemizygous %sindividuals'),
    'GC': ('G', 'Count of %sindividuals for each genotype')
}

ADJ_GQ = 20
ADJ_DP = 10
ADJ_AB = 0.2

FILTERS_DESC = {
    'InbreedingCoeff': 'InbreedingCoeff < -0.3',
    'LCR': 'In a low complexity region',
    'LowQual': 'Low quality',
    'PASS': 'All filters passed for at least one of the alleles at that site (see AS_FilterStatus for allele-specific filter status)',
    'RF': 'Failed random forests filters for all alleles (SNV cutoff %s, indels cutoff %s)',
    'SEGDUP': 'In a segmental duplication region',
    'AC0': 'Allele Count is zero for all alleles (i.e. no high-confidence genotype (GQ >= %(gq)s, DP >= %(dp)s, AB => %(ab)s for het calls) was found for each alternate allele)' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}
}

adj_criteria = 'g.gq >= %(gq)s && g.dp >= %(dp)s && (' \
               '!g.isHet || ' \
               '(g.gtj == 0 && g.ad[1]/g.dp >= %(ab)s) || ' \
               '(g.gtj > 0 && g.ad[0]/g.dp >= %(ab)s && g.ad[1]/g.dp >= %(ab)s)' \
               ')' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}


def get_info_va_attr():
    va_attr = {
        'AC_raw': [("Number", "A"), ("Description",
                                     "Allele counts before filtering low-confidence genotypes, for each ALT allele, in the same order as listed")],
        'AF_raw': [("Number", "A"), ("Description",
                                     "Allele frequency before filtering low-confidence genotypes, for each ALT allele, in the same order as listed")],
        'AN_raw': [("Description", "Total number of alleles before filtering low-confidence genotypes")],
        'GC_raw': [("Number", "G"), ("Description",
                                     "Raw count of individuals for each genotype before filtering low-confidence genotypes")],
        'Hom_raw': [("Number", "A"), ("Description", "Count of homozygous individuals in raw genotypes before filtering low-confidence genotypes")],
        'Hemi_raw': [("Number", "A"), ("Description", "Count of hemizygous individuals in raw genotypes before filtering low-confidence genotypes")],
        'BaseQRankSum': [("Description", "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities")],
        'ClippingRankSum': [("Description", 'Z-score from Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases')],
        'DB': [("Description", "dbSNP Membership")],
        'DP': [("Description", "Approximate read depth; some reads may have been filtered")],
        'FS': [("Description", "Phred-scaled p-value using Fisher's exact test to detect strand bias")],
        'HaplotypeScore': [("Description", "Consistency of the site with at most two segregating haplotypes")],
        'InbreedingCoeff': [("Description", 'Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation')],
        'MQ': [('Description', 'RMS Mapping Quality')],
        'MQRankSum': [("Description", "Z-score from Wilcoxon rank sum test of Alt vs. Ref read mapping qualities")],
        'QD': [("Description", "Variant Confidence/Quality by Depth")],
        'ReadPosRankSum': [("Description", "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias")],
        'VQSLOD': [("Description", "Log odds ratio of being a true variant versus being false under the trained VQSR gaussian mixture model (deprecated; see AS_RF)")],
        'VQSR_culprit': [("Description",
                    "The annotation which was the worst performing in the VQSR Gaussian mixture model (deprecated; see AS_RF)")],
        'VQSR_POSITIVE_TRAIN_SITE': [("Description", "This variant was used to build the positive training set of good variants for VQSR (deprecated; see AS_RF_POSITIVE_TRAIN)")],
        'VQSR_NEGATIVE_TRAIN_SITE': [("Description", "This variant was used to build the negative training set of bad variants for VQSR (deprecated; see AS_RF_NEGATIVE_TRAIN)")],
        'POPMAX': [("Number", "A"), ("Description", "Population with max AF")],
        'AC_POPMAX': [("Number", "A"), ("Description", "AC in the population with the max AF")],
        'AF_POPMAX': [("Number", "A"), ("Description", "Maximum Allele Frequency across populations (excluding OTH)")],
        'AN_POPMAX': [("Number", "A"), ("Description", "AN in the population with the max AF")],
        'STAR_AC': [("Description", "AC of deletions spanning this position")],
        'STAR_AC_raw': [("Description", "Allele counts of deletions spanning this position before filtering low-confidence genotypes")],
        'STAR_Hom': [("Description", "Count of individuals homozygous for a deletion spanning this position")],
        'STAR_Hemi': [("Description", "Count of individuals hemizygous for a deletion spanning this position")],
        'AS_RF': [("Number", "A"),("Description", "Random Forests probability for each allele")],
        'AS_FilterStatus': [("Number", "A"), ("Description", "Random Forests filter status for each allele")],
        'AS_RF_POSITIVE_TRAIN': [("Number", "."), ("Description", "Contains the indices of all alleles used as positive examples during training of random forests")],
        'AS_RF_NEGATIVE_TRAIN': [("Number", "."), ("Description", "Contains the indices of all alleles used as negative examples during training of random forests")],
        'SOR': [('Description', 'Symmetric Odds Ratio of 2x2 contingency table to detect strand bias')],
        'AB_HIST_ALT': [('Number', 'A'), ('Description', 'Histogram for Allele Balance in heterozygous individuals for each allele; 100*AD[i_alt]/sum(AD); Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5')],
        'GQ_HIST_ALT': [("Number", 'A'), ("Description", "Histogram for GQ for each allele; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5")],
        'DP_HIST_ALT': [("Number", 'A'), ("Description",
                                          "Histogram for DP for each allele; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5")],
        'AB_HIST_ALL': [('Number', '1'), ('Description',
                                          'Histogram for Allele Balance in heterozygous individuals; 100*AD[i_alt]/sum(AD); Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5')],
        'GQ_HIST_ALL': [("Number", '1'), ("Description",
                                          "Histogram for GQ; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5")],
        'DP_HIST_ALL': [("Number", '1'), ("Description",
                                          "Histogram for DP; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5")],
        'GQ_MEDIAN': [("Number", "A"), ('Description', 'Median GQ in carriers of each allele')],
        'DP_MEDIAN': [("Number", "A"), ('Description', 'Median DP in carriers of each allele')],
        'AB_MEDIAN': [("Number", "A"), ('Description', 'Median allele balance in heterozygote carriers of each allele')],
        'DREF_MEDIAN': [("Number", "A"), ('Description', 'Median dosage of homozygous reference in carriers of each allele')],
        'PROJECTMAX': [("Number", "A"), ('Description', "Projects with the highest AF for each allele (up to 5 projects per allele)")],
        'PROJECTMAX_NSamples': [("Number", "A"),
                                ('Description', "Number of samples in projects with the highest proportion of non-ref samples for each allele (up to 5 projects per allele)")],
        'PROJECTMAX_NonRefSamples': [("Number", "A"),
                                     ('Description', "Number of non-ref samples in projects with the highest proportion of non-ref samples for each allele (up to 5 projects per allele)")],
        'PROJECTMAX_PropNonRefSamples': [("Number", "A"),
                                         ('Description', "Proportion of non-ref samples in projects with the highest proportion of non-ref samples for each allele (up to 5 projects per allele)")]
    }

    for ann, ann_desc in ANNOTATION_DESC.items():
        va_attr[ann] = [("Number", ann_desc[0]), ("Description", ann_desc[1] % "")]
        for pop, pop_name in POP_NAMES.items():
            va_attr[ann + "_" + pop] = [("Number", ann_desc[0]), ("Description", ann_desc[1] % (pop_name + " "))]
            for sex, sex_name in SEXES.items():
                va_attr[ann + "_" + pop + "_" + sex] = [("Number", ann_desc[0]), ("Description", ann_desc[1] % (pop_name + " " + sex_name + " "))]
        for sex, sex_name in SEXES.items():
            va_attr[ann + "_" + sex] = [("Number", ann_desc[0]), ("Description", ann_desc[1] % (sex_name + " "))]

    return va_attr


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


def popmax_text(input_pops, skip_other=True):
    pops = copy.deepcopy(input_pops)
    if skip_other:
        if 'oth' in pops: pops.remove('oth')
        if 'OTH' in pops: pops.remove('OTH')
    af_pops = ','.join(['va.info.AF_%s' % pop for pop in pops])
    skip_other_text = '.filter(x => x != "oth" && x != "OTH")' if skip_other else ''

    get_af_max = 'va.AF_max = let af = [%s] and pops = global.pops%s in range(v.nAltAlleles).map(a => range(pops.size).sortBy(x => af[x][a],false)[0])' % (af_pops, skip_other_text)

    command = []
    for pop in pops:
        this_pop = '''if(pops[va.AF_max[a]] == "%(pop_lower)s" && va.info.AC_%(pop)s[a]>0)
        {pop: "%(pop)s", AF: va.info.AF_%(pop)s[a], AC: va.info.AC_%(pop)s[a], AN: va.info.AN_%(pop)s} ''' % {'pop': pop, 'pop_lower': pop.lower()}
        command.append(this_pop)

    command.append('na')
    get_popmax = """va.popmax = let na = NA : Struct{ pop: String, AF: Double, AC: Int, AN: Int} and pops = global.pops%s in
        range(va.AF_max.size).map(a => %s) """ % (skip_other_text, '\n else '.join(command))

    extract_popmax = """va.info.POPMAX = va.popmax.map(x => x.pop), va.info.AC_POPMAX = va.popmax.map(x => x.AC),
        va.info.AN_POPMAX = va.popmax.map(x => x.AN), va.info.AF_POPMAX = va.popmax.map(x => x.AF) """

    return get_af_max, get_popmax, extract_popmax


def get_hom_from_gc(destination, target):
    return '%s = range(v.nAltAlleles).map(i => let n = i + 2 in %s[(n * (n + 1) / 2).toInt - 1])' % (destination, target)


def unfurl_hom_text(pops, simple_hom=True, hom_adj=True):
    expressions = [get_hom_from_gc('va.info.Hom_%s' % pop, 'va.info.GC_%s' % pop) for pop in pops]
    if simple_hom: expressions.append('va.info.Hom_raw = range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_raw[(n * (n + 1) / 2).toInt - 1])')
    if hom_adj: expressions.append('va.info.Hom = range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC[(n * (n + 1) / 2).toInt - 1])')
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


def index_into_arrays(a_based_annotations):
    annotations = []
    if a_based_annotations:
        for ann in a_based_annotations:
            annotations.append('%s = %s[va.aIndex - 1]' % (ann, ann))

    return annotations


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


class VariantDataset(hail.dataset.VariantDataset):
    """
    Custom extensions to VDS for gnomAD analyses
    """

    def unfurl_callstats(self, pops, lower=False, gc=True):
        callstats_command, right_shift_command = unfurl_callstats_text(pops, lower, gc)
        return (self.annotate_variants_expr(callstats_command)
                .annotate_variants_expr(right_shift_command))

    def unfurl_hom(self, pops, simple_hom=True, hom_adj=True):
        hom_command = unfurl_hom_text(pops, simple_hom, hom_adj)
        return self.annotate_variants_expr(hom_command)

    def popmax(self, pops, skip_other=True):
        get_af_max, get_popmax, extract_popmax = popmax_text(pops, skip_other)
        return (self.annotate_variants_expr(get_af_max)
                .annotate_variants_expr(get_popmax)
                .annotate_variants_expr(extract_popmax))

    def projectmax(self):
        return ( self.annotate_alleles_expr('va.projectmax = let nNonRef = gs.filter(g => g.isCalledNonRef).map(g => if(isDefined(g)) sa.meta.project_description else NA: String).counter() and '
                                   'nSamples = gs.filter(g => g.isCalled).map(g => if(isDefined(g)) sa.meta.project_description else NA: String).counter() in '
                                   'nNonRef.keys.map(x => {key: x, count: nNonRef[x], nsamples: nSamples[nSamples.keys.find(y => x == y)]}).sortBy(x => x.count / x.nsamples,false)[0:5]')
                 .annotate_variants_expr('va.info.PROJECTMAX = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => x.key).mkString("|")), '
                            'va.info.PROJECTMAX_NSamples = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => str(x.nsamples)).mkString("|")), '
                            'va.info.PROJECTMAX_NonRefSamples = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => str(x.count)).mkString("|")), '
                            'va.info.PROJECTMAX_PropNonRefSamples = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => str(x.count / x.nsamples)).mkString("|"))')
                 )

    def filter_to_adj(self):
        return self.filter_genotypes(adj_criteria)

    def filter_star(self, a_based=None, r_based=None, g_based=None, additional_annotations=None):
        annotation = unfurl_filter_alleles_annotation(a_based=a_based, r_based=r_based, g_based=g_based, additional_annotations=additional_annotations)
        return self.filter_alleles('v.altAlleles[aIndex - 1].alt == "*"', annotation=annotation, keep=False)

    def head(self):
        return json.loads(self.variants_keytable().to_dataframe().toJSON().first())

    def set_vcf_filters(self, filters_dict, filters_to_keep=[]):

        site_filters = ",".join(['if(%s) "%s" else NA: String' % (filter_expr,name) for (name,filter_expr) in filters_dict.items()])
        site_filters = '[%s].filter(x => isDefined(x)).toSet' %site_filters

        if len(filters_to_keep) > 0:
            let_stmt = 'let prev_filters = va.filters.filter(x => ["%s"].toSet.contains(x)) ' % '","'.join(filters_to_keep)
        else:
            let_stmt = 'let prev_filters = [""][:0].toSet '

        return( self.annotate_variants_expr('va.filters = ' + let_stmt +
               'if(prev_filters.isEmpty) %s \n' +
               'else [prev_filters,%s].toSet.flatten' %(site_filters,site_filters))
                )

    def histograms(self, root='va.info', AB=True, asText=True, extra_gs_filter=''):

        allele_hists = ['%s.GQ_HIST_ALT = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).hist(0, 100, 20)' % (root, extra_gs_filter),
                        '%s.DP_HIST_ALT = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).hist(0, 100, 20)' % (root, extra_gs_filter)]
        variants_hists = ['%s.GQ_HIST_ALL = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).hist(0, 100, 20)' % (root, extra_gs_filter),
                          '%s.DP_HIST_ALL = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).hist(0, 100, 20)' % (root, extra_gs_filter)]

        if AB:
            allele_hists.append('%s.AB_HIST_ALT = gs.filter(g => g.isHet %s).map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)' % (root, extra_gs_filter))
            variants_hists.append('%s.AB_HIST_ALL = gs.filter(g => g.isHet %s).map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)' % (root, extra_gs_filter))

        if asText:
            allele_hists = ['%s.binFrequencies.map(y => str(y)).mkString("|")' % x for x in allele_hists]
            variants_hists = ['%s.binFrequencies.map(y => str(y)).mkString("|")' % x for x in variants_hists]

        return (
            self.annotate_alleles_expr(allele_hists, propagate_gq=True)
                .annotate_variants_expr(variants_hists)
        )


class HailContext(hail.context.HailContext):
    def _run_command(self, vds, pargs):
        jargs = jarray(self._jvm.java.lang.String, pargs)
        t = self._hail.driver.ToplevelCommands.lookup(jargs)
        cmd = t._1()
        cmd_args = t._2()
        jstate = self._jstate(vds._jvds if vds != None else None)
        try:
            result = cmd.run(jstate, cmd_args)
        except Py4JJavaError as e:
            raise_py4j_exception(e)
        return VariantDataset(self, result.vds())


def getAnnType(annotation, schema):
    ann_path = annotation.split(".")[1:]
    ann_type = schema
    for p in ann_path:
        try:
            ann_type = [x for x in ann_type.fields if x.name == p][0].typ
        except Exception, e:
            print schema
            print 'ERROR: %s missing from schema above' % p
            sys.exit(1)
    return ann_type


def annotate_non_split_from_split(hc, non_split_vds_path, split_vds, annotations):

    ann_list = annotations.keys()

    ann_types = map(lambda x: str(getAnnType(x,split_vds.variant_schema)), ann_list)

    variant_annotated_vds = (
        hc.read(non_split_vds_path, sites_only=True)
        .annotate_variants_expr('va.variant = v')
        .split_multi()
    )

    ann_agg_codes = ["`%s` = index(va.map(x => {val: %s, aIndex: va.aIndex}).collect(), aIndex)" % (a, a) for a in ann_list]
    agg = (
        split_vds
            .annotate_variants_vds(variant_annotated_vds, 'va.variant = vds.variant, va.aIndex = vds.aIndex')
            .filter_variants_expr('isDefined(va.variant)')
            .variants_keytable()
            .aggregate_by_key(key_condition='variant = va.variant', agg_condition=",".join(ann_agg_codes))

     )

    ann_codes = ['%s = let x = table.`%s` in' \
                 ' range(table.variant.nAltAlleles).map(i => if(x.contains(i+1)) x[i+1].val else NA: %s)' % (annotations[ann], ann, typ)
                 for (ann, typ) in zip(ann_list, ann_types)]

    return (
        hc.read(non_split_vds_path)
            .annotate_variants_keytable(agg, ",".join(ann_codes))
    )


def get_variant_type_expr(code="va.variantType"):
    return('''%s =
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
            "mixed"''' % code)


def get_stats_expr(root="va.stats", medians=False, samples_filter_expr=''):

    if samples_filter_expr:
        samples_filter_expr = "&& " + samples_filter_expr

    stats = ['%s.gq = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).stats()',
             '%s.dp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).stats()',
             '%s.nrq = gs.filter(g => g.isCalledNonRef %s).map(g => g.dosage[0]).stats()',
             '%s.ab = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).stats()']

    medians_expr = ['%s.gq_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).collect().median',
                    '%s.dp_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).collect().median',
                    '%s.nrq_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.dosage[0]).collect().median',
                    '%s.ab_median = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).collect().median']

    stats_expr = [x % (root,samples_filter_expr) for x in stats]

    if medians:
        stats_expr.extend([x % (root,samples_filter_expr) for x in medians_expr])

    return stats_expr


def post_process_vds(hc, vds_path, rf_vds, rf_root, rf_train, rf_label, rf_snv_cutoff, rf_indel_cutoff, vep_config):
    print("Postprocessing %s\n" % vds_path)

    as_filters = {
        'AC0': 'isMissing(va.info.AC[i]) || va.info.AC[i]<1'
    }

    filters = {
        'RF': 'isMissing(va.info.AS_FilterStatus) || '
              '(va.info.AS_FilterStatus.forall(x => !x.contains("PASS")) && va.info.AS_FilterStatus.exists(x => x == "RF"))',
        'AC0': '(va.info.AS_FilterStatus.forall(x =>!x.contains("PASS")) && va.info.AS_FilterStatus.exists(x => x == "AC0"))',
        'SEGDUP': 'va.decoy',
        'LCR': 'va.lcr'
    }

    rf_annotations = {
        'va.stats.qc_samples_raw.nrq_median': 'va.info.DREF_MEDIAN',
        'va.stats.qc_samples_raw.gq_median': 'va.info.GQ_MEDIAN',
        'va.stats.qc_samples_raw.dp_median': 'va.info.DP_MEDIAN',
        'va.stats.qc_samples_raw.ab_median': 'va.info.AB_MEDIAN'
    }

    vds = annotate_from_rf(hc, vds_path, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, annotations=rf_annotations, train=rf_train, label=rf_label)

    vds = add_as_filters(vds,as_filters)

    vds = set_vcf_filters(vds, rf_snv_cutoff, rf_indel_cutoff, filters = filters, filters_to_keep = ['InbreedingCoeff'])

    vds = vds.vep(config=vep_config, csq=True, root='va.info.CSQ', force=True)

    return set_va_attributes(vds)


def write_vcfs(vds, contig, out_internal_vcf_prefix, out_external_vcf_prefix, intervals_tmp='/tmp', append_to_header=None, drop_fields=None):

    if contig != '':
        print 'Writing VCFs for chr%s' % contig
        interval_path = '%s/%s.txt' % (intervals_tmp, str(contig))
        with open(interval_path, 'w') as f:
            f.write('%s:1-1000000000' % str(contig))

        vds = vds.filter_variants_intervals('file://' + interval_path)
    else:
        contig = 'autosomes'

    if drop_fields is not None:
        vds = vds.annotate_variants_expr('va.info = drop(va.info, %s)' % ",".join(drop_fields))

    vds = vds.annotate_variants_expr(['va.filters = if(va.filters.isEmpty) ["PASS"].toSet else va.filters',
                                      'va.info.AS_FilterStatus = '
                                      'va.info.AS_FilterStatus.map(x => if(x.isEmpty) "PASS" else x.toArray.mkString("|")'])

    vds.export_vcf(out_internal_vcf_prefix + ".%s.vcf.bgz" % contig, append_to_header=append_to_header)

    (
        vds.annotate_variants_expr(
            'va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
            .export_vcf(out_external_vcf_prefix + ".%s.vcf.bgz" % contig, append_to_header=append_to_header)
    )


def common_sites_vds_annotations(vds):
    return(
        vds.annotate_variants_expr(['va.info.VQSR_culprit = va.info.culprit',
                                    'va.info.VQSR_NEGATIVE_TRAIN_SITE = va.info.NEGATIVE_TRAIN_SITE ',
                                    'va.info.VQSR_POSITIVE_TRAIN_SITE = va.info.POSITIVE_TRAIN_SITE'])
        .annotate_variants_expr('va.info = drop(va.info, culprit,NEGATIVE_TRAIN_SITE,POSITIVE_TRAIN_SITE, DS, END, HaplotypeScore)')
    )


def create_sites_vds_annotations(vds, pops, tmp_path="/tmp", dbsnp_path=None):

    auto_intervals_path = '%s/autosomes.txt' % tmp_path
    with open(auto_intervals_path, 'w') as f:
        f.write('\n'.join(['%s:1-1000000000' % x for x in map(str, range(1, 23))]))

    sexes = ['Male', 'Female']
    cuts = copy.deepcopy(pops)
    cuts.extend(sexes)

    g_based_annotations = ['va.info.GC', 'va.info.GC_raw']
    g_based_annotations.extend(['va.info.GC_%s' % x for x in cuts])
    a_based_annotations = ['va.info.AC', 'va.info.AC_raw', 'va.info.AF', 'va.info.AF_raw', 'va.info.Hom', 'va.info.Hom_raw']
    a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                                'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
    a_based_annotations.extend(['va.info.AC_%s' % x for x in cuts])
    a_based_annotations.extend(['va.info.AF_%s' % x for x in cuts])
    a_based_annotations.extend(['va.info.Hom_%s' % x for x in cuts])
    a_based_annotations.extend(['va.info.GQ_HIST_ALT', 'va.info.DP_HIST_ALT', 'va.info.AB_HIST_ALT'])

    criterion_pops = [('sa.meta.population', x) for x in pops]
    criterion_pops.extend([('sa.meta.sex', x) for x in sexes])

    star_annotations = ['va.info.STAR_%s = let removed_allele = range(1, v.nAltAlleles + 1).find(i => !aIndices.toSet.contains(i)) \n'
                        'in if(isDefined(removed_allele)) va.info.%s[removed_allele - 1] else NA: Int' % (a, a) for a in ['AC', 'AC_raw', 'Hom']]

    vds = vds.filter_variants_intervals('file://' + auto_intervals_path)

    vds = common_sites_vds_annotations(vds)

    if dbsnp_path is not None:
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=hail.TextTableConfig(noheader=True, comment="#", types='_0: String, _1: Int')
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
            .annotate_variants_expr('va.info.AC_raw = va.calldata.raw.AC[1:], '
                                    'va.info.AN_raw = va.calldata.raw.AN, '
                                    'va.info.AF_raw = va.calldata.raw.AF[1:], '
                                    'va.info.GC_raw = va.calldata.raw.GC, '
                                    'va.info.AC = va.calldata.Adj.AC[1:], '
                                    'va.info.AN = va.calldata.Adj.AN, '
                                    'va.info.AF = va.calldata.Adj.AF[1:], '
                                    'va.info.GC = va.calldata.Adj.GC')
            .unfurl_hom(cuts)
            .persist()
            .filter_star(a_based=a_based_annotations, g_based=g_based_annotations,
                         additional_annotations=star_annotations)
            .popmax(pops)
            .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)'))


def create_sites_vds_annotations_X(vds, pops, tmp_path="/tmp", dbsnp_path=None):
    x_intervals = '%s/chrX.txt' % tmp_path
    with open(x_intervals, 'w') as f:
        f.write('X:1-1000000000')

    sexes = ['Male', 'Female']

    g_based_annotations = ['va.info.GC', 'va.info.GC_raw']
    g_based_annotations.extend(['va.info.GC_%s' % x for x in sexes])
    g_based_annotations.extend(['va.info.GC_%s_%s' % (y, x) for x in sexes for y in pops])

    a_based_annotations = ['va.info.AC', 'va.info.AC_raw', 'va.info.AF', 'va.info.AF_raw', 'va.info.Hom',
                           'va.info.Hom_raw', 'va.info.Hemi', 'va.info.Hemi_raw']
    a_based_annotations.extend(['va.info.AC_Male', 'va.info.AC_Female', 'va.info.AF_Male', 'va.info.AF_Female'])
    a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                                'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
    a_based_annotations.extend(['va.info.AC_%s_%s' % (y, x) for x in sexes for y in pops])
    a_based_annotations.extend(['va.info.AF_%s_%s' % (y, x) for x in sexes for y in pops])
    a_based_annotations.extend(['va.info.AC_%s' % x for x in pops])
    a_based_annotations.extend(['va.info.AF_%s' % x for x in pops])
    a_based_annotations.extend(['va.info.Hom_%s' % x for x in pops])
    a_based_annotations.extend(['va.info.Hemi_%s' % x for x in pops])
    a_based_annotations.extend(['va.info.GQ_HIST_ALT', 'va.info.DP_HIST_ALT', 'va.info.AB_HIST_ALT'])

    generate_callstats_expression = []
    for sex in sexes:
        for pop in pops:
            input_dict = {'sex': sex.lower(), 'sex_label': sex.capitalize(), 'pop': pop.lower(),
                          'pop_upper': pop.upper()}
            generate_callstats_expression.append(
                'va.calldata.%(pop_upper)s_%(sex_label)s = gs.filter(g => sa.meta.population == "%(pop)s" && sa.meta.sex == "%(sex)s").callStats(g => v)' % input_dict)
        input_dict = {'sex': sex.lower(), 'sex_label': sex.capitalize()}
        generate_callstats_expression.append(
            'va.calldata.%(sex_label)s = gs.filter(g => sa.meta.sex == "%(sex)s").callStats(g => v)' % input_dict)
    generate_callstats_expression = ',\n'.join(generate_callstats_expression)

    rearrange_callstats_expression = []
    for metric in ['AC', 'AN', 'AF', 'GC']:
        for sex in ['male', 'female']:
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

    rearrange_callstats_expression.extend(['va.info.GC_raw = va.calldata.raw.GC',
                                           'va.info.AC_raw = if (v.inXNonPar) (va.calldata.raw.AC - (va.calldata.hemi_raw.AC/2)).map(x => x.toInt)[1:] else va.calldata.raw.AC[1:]',
                                           'va.info.AN_raw = if (v.inXNonPar) (va.calldata.raw.AN - (va.calldata.hemi_raw.AN/2)).toInt else va.calldata.raw.AN',
                                           'va.info.GC = va.calldata.Adj.GC',
                                           'va.info.AC = if (v.inXNonPar) (va.calldata.Adj.AC - (va.calldata.Hemi_Adj.AC/2)).map(x => x.toInt)[1:] else va.calldata.Adj.AC[1:]',
                                           'va.info.AN = if (v.inXNonPar) (va.calldata.Adj.AN - (va.calldata.Hemi_Adj.AN/2)).toInt else va.calldata.Adj.AN'])
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
            'va.info.AF_%(pop_upper)s = va.info.AC_%(pop_upper)s/va.info.AN_%(pop_upper)s' % input_dict)
    af_expression = ',\n'.join(af_expression)

    hom_hemi_expression = []
    for pop in pops:
        input_dict = {'pop': pop, 'pop_upper': pop.upper()}
        #Hom
        hom_hemi_expression.append(
            'va.info.Hom_%(pop_upper)s =  if (!v.inXNonPar) '
            '   let GC = va.info.GC_%(pop_upper)s_Male + va.info.GC_%(pop_upper)s_Female in range(v.nAltAlleles).map(i => let n = i + 2 in GC[(n * (n + 1) / 2).toInt - 1])'
            'else range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(pop_upper)s_Female[(n * (n + 1) / 2).toInt - 1])' % input_dict)
        #Hemi
        hom_hemi_expression.append('va.info.Hemi_%(pop_upper)s = if (v.inXNonPar) range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(pop_upper)s_Male[(n * (n + 1) / 2).toInt - 1]) else NA: Array[Int]' % input_dict)
    hom_hemi_expression.append(
        'va.info.Hom = let GC = if (v.inXNonPar) va.calldata.Adj.GC - va.calldata.Hemi_Adj.GC else va.calldata.Adj.GC'
        '   in range(v.nAltAlleles).map(i => let n = i + 2 in GC[(n * (n + 1) / 2).toInt - 1])')
    hom_hemi_expression.append(
        'va.info.Hemi = range(v.nAltAlleles).map(i => let n = i + 2 in va.calldata.Hemi_Adj.GC[(n * (n + 1) / 2).toInt - 1])')
    hom_hemi_expression.append(
        'va.info.Hom_raw = let GC = if (v.inXNonPar) va.calldata.raw.GC - va.calldata.hemi_raw.GC else va.calldata.raw.GC'
        '   in range(v.nAltAlleles).map(i => let n = i + 2 in GC[(n * (n + 1) / 2).toInt - 1])')
    hom_hemi_expression.append(
        'va.info.Hemi_raw = range(v.nAltAlleles).map(i => let n = i + 2 in va.calldata.hemi_raw.GC[(n * (n + 1) / 2).toInt - 1])')
    hom_hemi_expression = ',\n'.join(hom_hemi_expression)

    star_annotations = [
        'va.info.STAR_%s = let removed_allele = range(1, v.nAltAlleles + 1).find(i => !aIndices.toSet.contains(i)) \n'
        'in if(isDefined(removed_allele)) va.info.%s[removed_allele - 1] else NA: Int' % (a, a) for a in
        ['AC', 'AC_raw', 'Hom', 'Hemi']]

    vds = vds.filter_variants_intervals('file://' + x_intervals)

    vds = common_sites_vds_annotations(vds)

    if dbsnp_path is not None:
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=hail.TextTableConfig(noheader=True, comment="#",
                                                                       types='_0: String, _1: Int')
                                         )
    return (vds.filter_genotypes('sa.meta.sex == "male" && g.isHet && v.inXNonPar', keep=False)
            .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
            .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
            .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
            .histograms('va.info')
            .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v), '
                                    'va.calldata.hemi_raw = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(g => v)')
            .filter_to_adj()
            .projectmax()
            .annotate_variants_expr('va.calldata.Adj = gs.callStats(g => v), '
                                    'va.calldata.Hemi_Adj = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(g => v)')
            .annotate_variants_expr(generate_callstats_expression)
            .filter_samples_all()
            .annotate_variants_expr(rearrange_callstats_expression)
            .annotate_variants_expr(
        'va.info.AF_raw = va.info.AC_raw.map(x => if (va.info.AN_raw > 0) x.toDouble/va.info.AN_raw else NA: Double), '
        'va.info.AF = va.info.AC.map(x => if (va.info.AN > 0) x.toDouble/va.info.AN else NA: Double)')  # Got here
            .annotate_variants_expr(hom_hemi_expression)
            .annotate_variants_expr(ac_an_expression)
            .annotate_variants_expr(af_expression)
            .persist()
            .filter_star(a_based=a_based_annotations, g_based=g_based_annotations,
                         additional_annotations=star_annotations)
            .popmax(pops)
            .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
            )


def create_sites_vds_annotations_Y(vds, pops, tmp_path="/tmp", dbsnp_path=None):
    y_intervals = '%s/chrY.txt' % tmp_path
    with open(y_intervals, 'w') as f:
        f.write('Y:1-1000000000')

    criterion_pops = [('sa.meta.population', x) for x in pops]

    a_based_annotations = ['va.info.AC', 'va.info.AC_raw', 'va.info.AF', 'va.info.AF_raw']
    a_based_annotations.extend(['va.info.PROJECTMAX', 'va.info.PROJECTMAX_NSamples',
                                'va.info.PROJECTMAX_NonRefSamples', 'va.info.PROJECTMAX_PropNonRefSamples'])
    a_based_annotations.extend(['va.info.AC_%s' % x for x in pops])
    a_based_annotations.extend(['va.info.AF_%s' % x for x in pops])
    a_based_annotations.extend(['va.info.GQ_HIST_ALT', 'va.info.DP_HIST_ALT'])

    # Dividing AC and AN by 2 on the y chromosome
    correct_ac_an_command = ['va.info.AC = va.info.AC.map(x => (x/2).toInt), '
                             'va.info.AN = (va.info.AN/2).toInt']
    for pop in pops + ['raw']:
        correct_ac_an_command.append('va.info.AC_%(pop)s = va.info.AC_%(pop)s.map(x => (x/2).toInt), '
                                     'va.info.AN_%(pop)s = (va.info.AN_%(pop)s/2).toInt' % {'pop': pop})

    correct_ac_an_command = ',\n'.join(correct_ac_an_command)

    star_annotations = [
        'va.info.STAR_%s = let removed_allele = range(1, v.nAltAlleles + 1).find(i => !aIndices.toSet.contains(i)) \n'
        'in if(isDefined(removed_allele)) va.info.%s[removed_allele - 1] else NA: Int' % (a, a) for a in
        ['AC', 'AC_raw']]
    vds = vds.filter_variants_intervals('file://' + y_intervals)

    vds = common_sites_vds_annotations(vds)

    if dbsnp_path is not None:
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=hail.TextTableConfig(noheader=True, comment="#", types='_0: String, _1: Int')
                                         )

    return (vds.filter_variants_expr('v.inYNonPar')
                 .filter_samples_expr('sa.meta.sex == "male"')
                 .filter_genotypes('g.isHet', keep=False)
                 .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
                 .filter_alleles('va.calldata.raw.AC[aIndex] == 0', keep=False)  # change if default is no longer subset
                 .histograms('va.info', AB=False)
                 .annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
                 .filter_to_adj()
                 .projectmax()
                 .annotate_variants_expr('va.calldata.Adj = gs.callStats(g => v)')
                 .unfurl_callstats(criterion_pops, lower=True, gc=False)
                 .filter_samples_all()
                 .annotate_variants_expr('va.info.AC_raw = va.calldata.raw.AC[1:], '
                                         'va.info.AN_raw = va.calldata.raw.AN, '
                                         'va.info.AF_raw = va.calldata.raw.AF[1:], '
                                         'va.info.AC = va.calldata.Adj.AC[1:], '
                                         'va.info.AN = va.calldata.Adj.AN, '
                                         'va.info.AF = va.calldata.Adj.AF[1:]')
                 .annotate_variants_expr(correct_ac_an_command)
                 .persist()
                 .filter_star(a_based=a_based_annotations, additional_annotations=star_annotations)
                 .popmax(pops)
                 .annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')
                 )

def annotate_from_rf(hc, vds_path, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, annotations={}, train='va.train', label='va.label'):

    rf_ann_expr = (['va.info.AS_RF = if(isMissing(%s)) NA: Array[Double] '
                    '    else %s.map(x => if(isDefined(x)) x.probability["TP"] else NA: Double)' % (rf_root, rf_root),
                    'va.info.AS_FilterStatus = if(isMissing(%(root)s)) NA: Array[Set[String]]'
                    '    else range(v.nAltAlleles).map(i => '
                    '        if(isMissing(%(root)s[i])) ["RF"].toSet ' #Sets missing RF values to filtered...
                    '        else if(v.altAlleles[i].isSNP) '
                    '            if(%(root)s[i].probability["TP"] > %(snv).4f) ["PASS"].toSet else ["RF"].toSet '
                    '            else if(%(root)s[i].probability["TP"] > %(indel).4f) ["PASS"].toSet else ["RF"].toSet)' %
                    {'root': rf_root,
                     'snv': rf_snv_cutoff,
                     'indel': rf_indel_cutoff},
                    'va.info.AS_RF_POSITIVE_TRAIN = '
                    'range(v.nAltAlleles).filter(i => isDefined(%s) && isDefined(%s) && %s[i] && %s[i] == "TP")'
                    '.map(i => i+1)' %
                    (train, label, train, label),
                    'va.info.AS_RF_NEGATIVE_TRAIN = '
                    'range(v.nAltAlleles).filter(i => isDefined(%s) && isDefined(%s) && %s[i] && %s[i] == "FP")'
                    '.map(i => i+1)' %
                    (train, label, train, label)
                    ])

    annotations[train] = train
    annotations[label] = label
    annotations[rf_root] = rf_root

    vds = annotate_non_split_from_split(hc, non_split_vds_path=vds_path,
                                        split_vds=rf_vds,
                                        annotations=annotations)

    return vds.annotate_variants_expr(rf_ann_expr)


def ann_exists(vds, ann_path):
    ann_path = ann_path.split(".")[1:]
    ann_type = vds.variant_schema
    for p in ann_path:
        field = [x for x in ann_type.fields if x.name == p]
        if not field:
            return False
        else:
            ann_type = [x for x in ann_type.fields if x.name == p][0].typ
    return True


def add_as_filters(vds, filters, root='va.info.AS_FilterStatus'):
    """
    Filters should be a dict of name: filter_expr
    Where i in the filter_expr is the alternate allele index (a-based)
    """

    as_filters = ",".join(['if(%s) "%s" else NA: String' % (filter_expr, name) for (name, filter_expr) in
                           filters.items()])
    as_filters = '[%s].filter(x => isDefined(x)).toSet' % as_filters

    input_dict = {
        'root': root,
        'filters': as_filters,
    }
    if not ann_exists(vds, root):
        vds = vds.annotate_variants_expr('%(root)s = range(v.nAltAlleles)'
                                         '.map(i => %(filters)s)' % input_dict)
    else:
        vds = vds.annotate_variants_expr('%(root)s = range(v.nAltAlleles).map(i => '
                                         'if(isMissing(%(root)s[i])) %(filters)s '
                                         'else [%(root)s[i],%(filters)s].toSet.flatten'
                                          % input_dict)
    return vds


def set_vcf_filters(vds, rf_snv_cutoff, rf_indel_cutoff, filters = {}, filters_to_keep = []):

    if len(filters) > 0 or len(filters_to_keep) > 0:
        vds = vds.set_vcf_filters(filters, filters_to_keep)

    for filter in filters.keys():
        desc = ""
        if(FILTERS_DESC.has_key(filter)):
            desc = FILTERS_DESC[filter]
        else:
            print("WARN: description for filter %s not found in FILTERS_DESC" % filter)

        if(filter == "RF"):
            desc = desc % (rf_snv_cutoff, rf_indel_cutoff)

        vds = vds.set_va_attribute('va.filters',filter,desc)

    for filter in filters_to_keep:
        desc = ""
        if (FILTERS_DESC.has_key(filter)):
            desc = FILTERS_DESC[filter]
        else:
            print("WARN: description for filter %s not found in FILTERS_DESC" % filter)
        vds = vds.set_va_attribute('va.filters', filter, desc)

    vds = vds.set_va_attribute('va.filters','PASS',FILTERS_DESC['PASS'])

    return vds


def run_sanity_checks(vds, pops, verbose=True, sex_chrom=False, percent_missing_threshold=0.01):

    queries = []
    a_metrics = ['AC','Hom']
    one_metrics = ['STAR_AC','AN']

    #Filter counts
    queries.append('variants.map(v => va.filters.toArray.mkString(",")).counter()')

    #Check that raw is always larger than adj
    queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                    '.exists(i => va.info.%s[i] > va.info.%s_raw[i])).count()' % (x, x) for x in a_metrics])
    queries.extend(['variants.filter(v => va.info.%s > va.info.%s_raw).count()' % (x, x) for x in one_metrics])

    #Check that sum(pops) == total
    for metric in a_metrics:
        queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                        '.exists(i => %s != va.info.%s[i])).count()' % (" + ".join(["va.info.%s_%s[i]" % (metric, pop) for pop in pops]), metric)])

    for metric in one_metrics[1:]:
        queries.extend(['variants.filter(v => %s != va.info.%s).count()' % (
                        " + ".join(["va.info.%s_%s" % (metric, pop) for pop in pops]), metric)])

    #Check that male + female == total
    #Remove Hom for X
    if sex_chrom:
        a_metrics = a_metrics[:1]

    pop_strats = pops if sex_chrom else [None]
    for pop in pop_strats:
        pop_text = "" if pop is None else "_" + pop
        queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                        '.exists(i => va.info.%s%s_Male[i] + va.info.%s%s_Female[i] != va.info.%s%s[i])).count()' % (x, pop_text, x, pop_text, x, pop_text) for x in a_metrics])
        queries.extend(['variants.filter(v => va.info.%s%s_Male + va.info.%s%s_Female != va.info.%s%s).count()' % (x, pop_text, x, pop_text, x, pop_text) for x in one_metrics[1:]])

    end_counts = len(queries)

    missing_metrics = []
    va_info = [x for x in vds.variant_schema.fields if x.name == "info"][0].typ
    for field in va_info.fields:
        missing_metrics.append('va.info.%s' % field.name)
        queries.append('variants.map(v => va.info.%s).fraction(x => isMissing(x))' % field.name)

    stats = vds.query_variants(queries)

    #Print filters
    print("FILTERS CHECKS\n")
    print("Filter counts:\n")
    pprint(stats[0])

    #Check that all metrics sum as expected
    print("\n\nMETRICS COUNTS CHECK\n")
    nfail = 0
    for i in range(1,end_counts):
        if stats[i] != 0:
            print("FAILED METRICS CHECK for query: %s\n Expected: 0, Found: %s" % (queries[i], stats[i]))
            nfail += 1
        elif verbose:
            print("Success: %s\n" % queries[i])
    print("%s metrics count checks failed.\n" % nfail)

    #Check missing metrics
    print("MISSING METRICS CHECKS")
    nfail = 0
    missing_stats = stats[end_counts:]
    for i in range(len(missing_stats)):
        if missing_stats[i] > percent_missing_threshold:
            print("FAILED missing check for %s; %s%% missing." % (missing_metrics[i], 100*missing_stats[i]))
            nfail += 1
        elif verbose:
            print("SUCCESS missing check for %s; %s%% missing." % (missing_metrics[i], 100*missing_stats[i]))
    print("%s missing metrics checks failed.\n" % nfail)

    return vds

def set_va_attributes(vds):

    info_va_attr = get_info_va_attr()
    va_info = [x for x in vds.variant_schema.fields if x.name == "info"][0]

    for ann in va_info.typ.fields:
        if info_va_attr.has_key(ann.name):
            attributes = info_va_attr[ann.name]
            for att in attributes:
                vds = vds.set_va_attribute("va.info.%s" % ann.name, att[0], att[1])

        elif ann.name != "CSQ": print("WARN: No description found for va.info.%s\n" % ann.name)

    return vds


def send_message(channel=None, message=None):
    import getpass
    from slackclient import SlackClient
    # import os
    try:
        from slack_creds import slack_token
    except Exception, e:
        return 0

    # slack_token = os.environ["SLACK_API_TOKEN"]
    sc = SlackClient(slack_token)
    user = getpass.getuser()
    if user.startswith('konrad'): user = 'konradjk'
    users = [x['name'] for x in sc.api_call("users.list")['members']]
    if channel is None:
        channel = '#gnomad' if user not in users else '@' + user
    if message is None:
        message = "Hey %s! Your job is done :tada:" % user
    sc.api_call(
        "chat.postMessage",
        channel=channel,
        text=message,
        icon_emoji=':woohoo:'
    )



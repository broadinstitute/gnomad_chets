
import re
import sys
import hail
import pyspark.sql
import json
import copy
import time
import logging
from py4j.protocol import Py4JJavaError
from subprocess import check_output
from pprint import pprint, pformat
import gzip
import os

from resources import *
from hail.expr import *
from hail.representation import *
from slack_utils import *
from pyspark.sql.functions import bround
import subprocess

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)

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
    'PASS': 'All filters passed for at least one of the alleles at that site (see AS_FilterStatus for allele-specific filter status)',
    'RF': 'Failed random forests filters (SNV cutoff %s, indels cutoff %s)',
    'SEGDUP': 'In a segmental duplication region',
    'AC0': 'Allele Count is zero (i.e. no high-confidence genotype (GQ >= %(gq)s, DP >= %(dp)s, AB => %(ab)s for het calls))' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}
}

ADJ_CRITERIA = 'g.gq >= %(gq)s && g.dp >= %(dp)s && (' \
               '!g.isHet || ' \
               '(g.gtj == 0 && g.ad[g.gtk]/g.dp >= %(ab)s) || ' \
               '(g.gtj > 0 && g.ad[g.gtj]/g.dp >= %(ab)s && g.ad[g.gtk]/g.dp >= %(ab)s)' \
               ')' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_ORDER = ["transcript_ablation",
"splice_acceptor_variant",
"splice_donor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost",
"start_lost",  # new in v81
"initiator_codon_variant",  # deprecated
"transcript_amplification",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"protein_altering_variant",  # new in v79
"splice_region_variant",
"incomplete_terminal_codon_variant",
"stop_retained_variant",
"synonymous_variant",
"coding_sequence_variant",
"mature_miRNA_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"non_coding_transcript_exon_variant",
"non_coding_exon_variant",  # deprecated
"intron_variant",
"NMD_transcript_variant",
"non_coding_transcript_variant",
"nc_transcript_variant",  # deprecated
"upstream_gene_variant",
"downstream_gene_variant",
"TFBS_ablation",
"TFBS_amplification",
"TF_binding_site_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
"feature_elongation",
"regulatory_region_variant",
"feature_truncation",
"intergenic_variant",
""]


def get_info_va_attr():
    va_attr = {
        'AC_raw': {"Number": "A",
                   "Description": "Allele counts before filtering low-confidence genotypes, for each ALT allele, in the same order as listed"},
        "AF_raw": {"Number": "A", "Description":
            "Allele frequency before filtering low-confidence genotypes, for each ALT allele, in the same order as listed"},
        "AN_raw": {"Description": "Total number of alleles before filtering low-confidence genotypes"},
        "GC_raw": {"Number": "G", "Description":
            "Raw count of individuals for each genotype before filtering low-confidence genotypes"},
        "Hom_raw": {"Number": "A",
                    "Description": "Count of homozygous individuals in raw genotypes before filtering low-confidence genotypes"},
        "Hemi_raw": {"Number": "A",
                     "Description": "Count of hemizygous individuals in raw genotypes before filtering low-confidence genotypes"},
        "BaseQRankSum": {"Description": "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities"},
        "ClippingRankSum": {
            "Description": "Z-score from Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases"},
        "DB": {"Description": "dbSNP Membership"},
        "DP": {"Description": "Approximate read depth; some reads may have been filtered"},
        "FS": {"Description": "Phred-scaled p-value using Fisher's exact test to detect strand bias"},
        "HaplotypeScore": {
            "Description": "Consistency of the site with at most two segregating haplotypes"},
        "InbreedingCoeff": {
            "Description": "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"},
        "MQ": {
            "Description": "RMS Mapping Quality"},
        "MQRankSum": {
            "Description": "Z-score from Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"},
        "QD": {
            "Description": "Variant Confidence/Quality by Depth"},
        "ReadPosRankSum": {
            "Description": "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"},
        "VQSLOD": {
            "Description": "Log odds ratio of being a true variant versus being false under the trained VQSR gaussian mixture model deprecated; see AS_RF"},
        "VQSR_culprit": {
            "Description":
                "The annotation which was the worst performing in the VQSR Gaussian mixture model deprecated; see AS_RF"},
        "VQSR_POSITIVE_TRAIN_SITE": {
            "Description": "This variant was used to build the positive training set of good variants for VQSR deprecated; see AS_RF_POSITIVE_TRAIN"},
        "VQSR_NEGATIVE_TRAIN_SITE": {
            "Description": "This variant was used to build the negative training set of bad variants for VQSR deprecated; see AS_RF_NEGATIVE_TRAIN"},
        "POPMAX": {
            "Number": "A",
            "Description": "Population with max AF"},
        "AC_POPMAX": {
            "Number": "A",
            "Description": "AC in the population with the max AF"},
        "AF_POPMAX": {
            "Number": "A",
            "Description": "Maximum Allele Frequency across populations (excluding OTH)"},
        "AN_POPMAX": {
            "Number": "A",
            "Description": "AN in the population with the max AF"},
        "STAR_AC": {
            "Description": "AC of deletions spanning this position"},
        "STAR_AC_raw": {
            "Description": "Allele counts of deletions spanning this position before filtering low-confidence genotypes"},
        "STAR_Hom": {
            "Description": "Count of individuals homozygous for a deletion spanning this position"},
        "STAR_Hemi": {
            "Description": "Count of individuals hemizygous for a deletion spanning this position"},
        "AS_RF": {
            "Number": "A",
            "Description": "Random Forests probability for each allele"},
        "AS_FilterStatus": {
            "Number": "A",
            "Description": "Random Forests filter status for each allele"},
        "AS_RF_POSITIVE_TRAIN": {
            "Number": ".",
            "Description": "Contains the indices of all alleles used as positive examples during training of random forests"},
        "AS_RF_NEGATIVE_TRAIN": {
            "Number": ".",
            "Description": "Contains the indices of all alleles used as negative examples during training of random forests"},
        "SOR": {
            "Description": "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias"},
        "AB_HIST_ALT": {
            "Number": "A",
            "Description": "Histogram for Allele Balance in heterozygous individuals for each allele; 100*AD[i_alt]/sum(AD); Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5"},
        "GQ_HIST_ALT": {
            "Number": "A",
            "Description": "Histogram for GQ for each allele; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5"},
        "DP_HIST_ALT": {
            "Number": "A",
            "Description":
                "Histogram for DP for each allele; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5"},
        "AB_HIST_ALL": {
            "Number": "1",
            "Description":
                "Histogram for Allele Balance in heterozygous individuals; 100*AD[i_alt]/sumAD; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5"},
        "GQ_HIST_ALL": {
            "Number": "1",
            "Description":
                "Histogram for GQ; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5"},
        "DP_HIST_ALL": {
            "Number": "1",
            "Description":
                "Histogram for DP; Midpoints of histogram bins: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5"},
        "GQ_MEDIAN": {
            "Number": "A",
            "Description": "Median GQ in carriers of each allele"},
        "DP_MEDIAN": {
            "Number": "A",
            "Description": "Median DP in carriers of each allele"},
        "AB_MEDIAN": {
            "Number": "A",
            "Description": "Median allele balance in heterozygote carriers of each allele"},
        "DREF_MEDIAN": {
            "Number": "A",
            "Description": "Median dosage of homozygous reference in carriers of each allele"},
        "PROJECTMAX": {
            "Number": "A",
            "Description": "Projects with the highest AF for each allele (up to 5 projects per allele)"},
        "PROJECTMAX_NSamples": {
            "Number": "A",
            "Description": "Number of samples in projects with the highest proportion of non-ref samples for each allele (up to 5 projects per allele"},
        "PROJECTMAX_NonRefSamples": {"Number": "A",
                                     "Description": "Number of non-ref samples in projects with the highest proportion of non-ref samples for each allele (up to 5 projects per allele"},
        "PROJECTMAX_PropNonRefSamples": {"Number": "A",
                                         "Description": "Proportion of non-ref samples in projects with the highest proportion of non-ref samples for each allele (up to 5 projects per allele"}
    }

    for ann, ann_desc in ANNOTATION_DESC.items():
        va_attr[ann] = {"Number": ann_desc[0], "Description": ann_desc[1] % ""}
        for pop, pop_name in POP_NAMES.items():
            va_attr[ann + "_" + pop] = {"Number": ann_desc[0], "Description": ann_desc[1] % (pop_name + " ")}
            for sex, sex_name in SEXES.items():
                va_attr[ann + "_" + pop + "_" + sex] = {"Number": ann_desc[0], "Description": ann_desc[1] % (pop_name + " " + sex_name + " ")}
        for sex, sex_name in SEXES.items():
            va_attr[ann + "_" + sex] = {"Number": ann_desc[0], "Description": ann_desc[1] % (sex_name + " ")}

    return va_attr


def popmax_text(input_pops, skip_other=True):
    af_pops = ','.join(['va.info.AF_%s' % pop for pop in input_pops])
    skip_other_text = '.filter(x => x != "oth" && x != "OTH")' if skip_other else ''

    get_af_max = 'va.AF_max = let af = [%s] and pops = global.pops%s in range(v.nAltAlleles).map(a => range(pops.size).sortBy(x => af[x][a],false)[0])' % (af_pops, skip_other_text)

    command = []
    for pop in input_pops:
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


def unfurl_callstats_text(criteria_pops, lower=True, gc=True):
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


def index_into_arrays(a_based_annotations=None, r_based_annotations=None, vep_root=None):
    annotations = []
    if a_based_annotations:
        for ann in a_based_annotations:
            annotations.append('%s = %s[va.aIndex - 1]' % (ann, ann))
    if r_based_annotations:
        for ann in r_based_annotations:
            annotations.append('%s = [%s[0], %s[va.aIndex]]' % (ann, ann, ann))
    if vep_root:
        sub_fields = ['transcript_consequences', 'intergenic_consequences', 'motif_feature_consequences', 'regulatory_feature_consequences']
        annotations.extend(['%s.%s = %s.%s.filter(x => x.allele_num == va.aIndex)' % (vep_root, sub_field, vep_root, sub_field) for sub_field in sub_fields])

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


def unfurl_callstats(vds, pops, lower=False, gc=True):
    callstats_command, right_shift_command = unfurl_callstats_text(pops, lower, gc)
    return (vds.annotate_variants_expr(callstats_command)
            .annotate_variants_expr(right_shift_command))


def unfurl_hom(vds, pops, simple_hom=True, hom_adj=True):
        hom_command = unfurl_hom_text(pops, simple_hom, hom_adj)
        return vds.annotate_variants_expr(hom_command)


def popmax(vds, input_pops, skip_other=True):
    proc_pops = copy.deepcopy(input_pops)
    if skip_other:
        if 'oth' in proc_pops: proc_pops.remove('oth')
        if 'OTH' in proc_pops: proc_pops.remove('OTH')
    if len(proc_pops) < 2: return vds
    get_af_max, get_popmax, extract_popmax = popmax_text(proc_pops, skip_other)
    return (vds.annotate_variants_expr(get_af_max)
            .annotate_variants_expr(get_popmax)
            .annotate_variants_expr(extract_popmax))


def projectmax(vds):
    return (vds.annotate_alleles_expr(
        'va.projectmax = let nNonRef = gs.filter(g => g.isCalledNonRef).map(g => if(isDefined(g)) sa.meta.project_description else NA: String).counter() and '
        'nSamples = gs.filter(g => g.isCalled).map(g => if(isDefined(g)) sa.meta.project_description else NA: String).counter() in '
        'nNonRef.keys.map(x => {key: x, count: nNonRef[x], nsamples: nSamples[nSamples.keys.find(y => x == y)]}).sortBy(x => x.count / x.nsamples,false)[0:5]')
        .annotate_variants_expr(
        'va.info.PROJECTMAX = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => x.key).mkString("|")), '
        'va.info.PROJECTMAX_NSamples = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => str(x.nsamples)).mkString("|")), '
        'va.info.PROJECTMAX_NonRefSamples = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => str(x.count)).mkString("|")), '
        'va.info.PROJECTMAX_PropNonRefSamples = va.projectmax.map(a => if(a.isEmpty) NA:String else a.map(x => str(x.count / x.nsamples)).mkString("|"))')
    )


def filter_to_adj(vds):
    return vds.filter_genotypes(ADJ_CRITERIA)


def filter_star(vds, a_based=None, r_based=None, g_based=None, additional_annotations=None):
    annotation = unfurl_filter_alleles_annotation(a_based=a_based, r_based=r_based, g_based=g_based,
                                                  additional_annotations=additional_annotations)
    return vds.filter_alleles('v.altAlleles[aIndex - 1].alt == "*"', annotation=annotation, keep=False)


def histograms(vds, root='va.info', AB=True, asText=True, extra_gs_filter=''):
    allele_hists = ['%s.GQ_HIST_ALT = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).hist(0, 100, 20)' % (
    root, extra_gs_filter),
                    '%s.DP_HIST_ALT = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).hist(0, 100, 20)' % (
                    root, extra_gs_filter)]
    variants_hists = [
        '%s.GQ_HIST_ALL = gs.filter(g => g.isCalled %s).map(g => g.gq).hist(0, 100, 20)' % (root, extra_gs_filter),
        '%s.DP_HIST_ALL = gs.filter(g => g.isCalled %s).map(g => g.dp).hist(0, 100, 20)' % (root, extra_gs_filter)]

    if AB:
        allele_hists.append(
            '%s.AB_HIST_ALT = gs.filter(g => g.isHet %s).map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)' % (
            root, extra_gs_filter))
        variants_hists.append(
            '%s.AB_HIST_ALL = gs.filter(g => g.isHet %s).map(g => 100 - 100*g.ad[0]/g.dp).hist(0, 100, 20)' % (
            root, extra_gs_filter))

    if asText:
        allele_hists = ['%s.binFrequencies.map(y => str(y)).mkString("|")' % x for x in allele_hists]
        variants_hists = ['%s.binFrequencies.map(y => str(y)).mkString("|")' % x for x in variants_hists]

    return (
        vds.annotate_alleles_expr(allele_hists, propagate_gq=True)
        .annotate_variants_expr(variants_hists)
    )

def flatten_struct(struct, root ='', leaf_only = True):
    result = {}
    for f in struct.fields:
        path = '%s.%s' % (root, f.name)
        if isinstance(f.typ, hail.expr.TStruct):
            result.update(flatten_struct(f.typ, path))
            if not leaf_only:
                result[path] = f
        else:
            result[path] = f
    return result

def ann_exists(annotation, schema, root = 'va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    return anns.has_key(annotation)

def get_ann_field(annotation, schema, root = 'va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    if not annotation in anns:
        logger.error("%s missing from schema.", annotation)
        sys.exit(1)
    return anns[annotation]


def get_ann_type(annotation, schema, root = 'va'):
    return get_ann_field(annotation, schema, root).typ


def get_variant_type_expr(code="va.variantType"):
    return '''%s =
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
            "mixed"''' % code


def get_stats_expr(root="va.stats", medians=False, samples_filter_expr=''):

    if samples_filter_expr:
        samples_filter_expr = "&& " + samples_filter_expr

    stats = ['%s.gq = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).stats()',
             '%s.dp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).stats()',
             '%s.nrq = gs.filter(g => g.isCalledNonRef %s).map(g => -log10(g.dosage[0])).stats()',
             '%s.ab = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).stats()',
             '%s.best_ab = gs.filter(g => g.isHet %s).map(g => abs((g.ad[1]/g.dp) - 0.5)).min()']

    if medians:
        stats.extend(['%s.gq_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).collect().median',
                    '%s.dp_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).collect().median',
                    '%s.nrq_median = gs.filter(g => g.isCalledNonRef %s).map(g => -log10(g.dosage[0])).collect().median',
                    '%s.ab_median = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).collect().median'])

    stats_expr = [x % (root,samples_filter_expr) for x in stats]

    return stats_expr


def set_site_filters(vds, site_filters_dict, filters_to_keep=[], as_filters_root="va.info.AS_FilterStatus"):
    site_filters = ",".join(['orMissing(%s, "%s")' % (filter_expr,name) for (name,filter_expr) in site_filters_dict.items()])
    site_filters = '[%s].filter(x => isDefined(x)).toSet' % site_filters

    if len(filters_to_keep) > 0:
        prev_filters = 'va.filters.filter(x => ["%s"].toSet.contains(x))' % '","'.join(filters_to_keep)
    else:
        prev_filters = '[""][:0].toSet'

    input_dict = {
        'site_filters': site_filters,
        'prev_filters': prev_filters,
        'as_filters': as_filters_root
    }

    annotate_expr = ('va.filters = let prev_filters = %(prev_filters)s '
                     'and sites_filters = %(site_filters)s '
                     'and as_filters = %(as_filters)s.find(x => isDefined(x) && x.isEmpty())'
                     '.orElse(%(as_filters)s.find(x => isMissing(x))'
                     '.orElse(%(as_filters)s.toSet.flatten)) in '
                     'if(!prev_filters.isEmpty() || !sites_filters.isEmpty()) '
                     ' prev_filters.union(sites_filters).union(as_filters.orElse([""][:0].toSet)) '
                     'else as_filters' % input_dict)

    logger.debug(annotate_expr)

    return vds.annotate_variants_expr(annotate_expr)


def post_process_subset(subset_vds, release_vds_dict, as_filters_key, dot_annotations_dict = None):

    logger.info("Postprocessing %s", subset_vds)

    for release, release_dict in release_vds_dict.iteritems():
        annotations_to_ignore = ['DB', 'GQ_HIST_ALL', 'DP_HIST_ALL', 'AB_HIST_ALL', 'GQ_HIST_ALT', 'DP_HIST_ALT', 'AB_HIST_ALT', 'AF_.*', 'A[CN]_..._.*ale']
        if release == as_filters_key:
            annotations_to_ignore.extend([
                'BaseQRankSum', 'ClippingRankSum', 'DP', 'FS', 'InbreedingCoeff', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum',
                'SOR', 'VQSLOD', 'VQSR_culprit', 'VQSR_NEGATIVE_TRAIN_SITE', 'VQSR_POSITIVE_TRAIN_SITE'
            ])
        subset_vds = annotate_subset_with_release(subset_vds, release_dict, dot_annotations_dict=dot_annotations_dict, ignore= annotations_to_ignore)

    subset_vds = subset_vds.annotate_variants_expr("va.info.AS_FilterStatus = %sAS_FilterStatus" % release_vds_dict[as_filters_key]['out_root'])
    subset_vds = set_filters(subset_vds)

    as_filters_attributes = get_ann_field('va.info.AS_FilterStatus', release_vds_dict[as_filters_key]['vds'].variant_schema).attributes
    subset_vds = subset_vds.set_va_attributes("va.info.AS_FilterStatus", as_filters_attributes)

    return set_va_attributes(subset_vds, warn_if_not_found=False)


def post_process_vds(vds, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root,
                     rf_train='va.train', rf_label='va.label'):

    logger.info("Postprocessing %s", vds)

    rf_annotations = {
        'va.stats.qc_samples_raw.nrq_median': 'va.info.DREF_MEDIAN',
        'va.stats.qc_samples_raw.gq_median': 'va.info.GQ_MEDIAN',
        'va.stats.qc_samples_raw.dp_median': 'va.info.DP_MEDIAN',
        'va.stats.qc_samples_raw.ab_median': 'va.info.AB_MEDIAN'
    }

    vds = annotate_from_rf(vds, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, annotations=rf_annotations, train=rf_train, label=rf_label)
    vds = set_filters(vds, rf_snv_cutoff, rf_indel_cutoff)

    return set_va_attributes(vds)


def set_filters(vds, rf_snv_cutoff=None, rf_indel_cutoff=None):
    as_filters = {
        'AC0': 'isMissing(va.info.AC[i]) || va.info.AC[i]<1'
    }

    site_filters = {
        'InbreedingCoeff': 'isDefined(va.info.InbreedingCoeff) && va.info.InbreedingCoeff < -0.3',
        'SEGDUP': 'va.decoy',
        'LCR': 'va.lcr'
    }

    vds = add_as_filters(vds, as_filters)
    vds = set_filters_attributes(vds, rf_snv_cutoff, rf_indel_cutoff)
    return set_site_filters(vds, site_filters, filters_to_keep=['InbreedingCoeff'])


def write_vcfs(vds, contig, out_internal_vcf_prefix, out_external_vcf_prefix, rf_snv_cutoff, rf_indel_cutoff,
               as_filter_status_fields=['va.info.AS_FilterStatus'],
               append_to_header=None, drop_fields=None, nchunks=None):

    if contig != '':
        if not isinstance(contig, list):
            vds = vds.filter_variants_intervals(Interval.parse(str(contig)))
        else:
            vds = vds.filter_variants_intervals(IntervalTree.parse_all(contig))
    logger.info('Writing VCFs for chromosome: %s', contig if contig != '' else 'all')

    if drop_fields is not None:
        vds = vds.annotate_variants_expr('va.info = drop(va.info, %s)' % ",".join(drop_fields))

    parallel = False
    if nchunks is not None:
        parallel = True
        vds = vds.repartition(nchunks, shuffle=False)

    # AS_FilterStatus can be either:
    # Missing => no filtering was applied to this allele
    # {} => "PASS"
    # {RF|AC0} => this allele has a filter
    as_filter_status_attributes = flatten_struct(vds.variant_schema, root="va")
    as_filter_status_expression = ['%s = %s.map(x => orMissing(isDefined(x), if(x.isEmpty()) "PASS" else x.toArray.mkString("|")))' % (x, x) for x in as_filter_status_fields]
    vds = vds.annotate_variants_expr(as_filter_status_expression)
    for x in as_filter_status_fields:
        vds = vds.set_va_attributes(x, as_filter_status_attributes[x].attributes)
    vds = set_filters_attributes(vds, rf_snv_cutoff, rf_indel_cutoff)

    if out_internal_vcf_prefix:
        out_internal_vcf = "%s%s.vcf.bgz" % (out_internal_vcf_prefix, '' if contig == '' else '.%s' % contig)
        vds.export_vcf(out_internal_vcf, append_to_header=append_to_header, parallel=parallel)

    if out_external_vcf_prefix:
        out_external_vcf = "%s%s.vcf.bgz" % (out_external_vcf_prefix, '' if contig == '' else '.%s' % contig)
        (
            vds.annotate_variants_expr(
                'va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
            .export_vcf(out_external_vcf, append_to_header=append_to_header, parallel=parallel)
        )


def common_sites_vds_annotations(vds):
    return (
        vds.annotate_variants_expr(['va.info.VQSR_culprit = va.info.culprit',
                                    'va.info.VQSR_NEGATIVE_TRAIN_SITE = va.info.NEGATIVE_TRAIN_SITE ',
                                    'va.info.VQSR_POSITIVE_TRAIN_SITE = va.info.POSITIVE_TRAIN_SITE'])
        .annotate_variants_expr('va.info = drop(va.info, culprit,NEGATIVE_TRAIN_SITE,POSITIVE_TRAIN_SITE, DS, END, HaplotypeScore)')
    )


def create_sites_vds_annotations(vds, pops, dbsnp_path=None, drop_star=True, filter_alleles = True, generate_hists= True):

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

    vds = vds.filter_variants_intervals(IntervalTree.parse_all(['1-22']))

    vds = common_sites_vds_annotations(vds)

    if dbsnp_path is not None:
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=hail.TextTableConfig(noheader=True, comment="#", types='_0: String, _1: Int')
                                         )

    if filter_alleles:
        vds = (vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
               .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
               .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
        )
    if generate_hists:
        vds = histograms(vds, 'va.info')

    vds = vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
    vds = filter_to_adj(vds)
    vds = projectmax(vds)
    vds = vds.annotate_variants_expr('va.calldata.Adj = gs.callStats(g => v)')

    vds = unfurl_callstats(vds, criterion_pops, lower=True)

    vds = vds.drop_samples()
    vds = (vds.annotate_variants_expr('va.info.AC_raw = va.calldata.raw.AC[1:], '
                                      'va.info.AN_raw = va.calldata.raw.AN, '
                                      'va.info.AF_raw = va.calldata.raw.AF[1:], '
                                      'va.info.GC_raw = va.calldata.raw.GC, '
                                      'va.info.AC = va.calldata.Adj.AC[1:], '
                                      'va.info.AN = va.calldata.Adj.AN, '
                                      'va.info.AF = va.calldata.Adj.AF[1:], '
                                      'va.info.GC = va.calldata.Adj.GC')
    )
    vds = unfurl_hom(vds, cuts)

    if drop_star:
        vds = vds.persist()
        vds = filter_star(vds, a_based=a_based_annotations, g_based=g_based_annotations, additional_annotations=star_annotations)
    vds = popmax(vds, pops)
    return vds.annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')


def create_sites_vds_annotations_X(vds, pops, dbsnp_path=None, drop_star=True, filter_alleles = True, generate_hists= True):

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
        # Hom
        hom_hemi_expression.append(
            'va.info.Hom_%(pop_upper)s =  if (!v.inXNonPar) '
            '   let GC = va.info.GC_%(pop_upper)s_Male + va.info.GC_%(pop_upper)s_Female in range(v.nAltAlleles).map(i => let n = i + 2 in GC[(n * (n + 1) / 2).toInt - 1])'
            'else range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_%(pop_upper)s_Female[(n * (n + 1) / 2).toInt - 1])' % input_dict)
        # Hemi
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

    vds = vds.filter_variants_intervals(Interval.parse("X"))

    vds = common_sites_vds_annotations(vds)

    if dbsnp_path is not None:
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=hail.TextTableConfig(noheader=True, comment="#",
                                                                       types='_0: String, _1: Int')
                                         )
    vds = vds.filter_genotypes('sa.meta.sex == "male" && g.isHet && v.inXNonPar', keep=False)

    if filter_alleles:
        vds = (vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
               .filter_alleles('va.calldata.raw.AC[aIndex] == 0', subset=True, keep=False)
               .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
               )

    if generate_hists:
        vds = histograms(vds, 'va.info')

    vds = vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v), '
                                    'va.calldata.hemi_raw = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(g => v)')
    vds = filter_to_adj(vds)
    vds = projectmax(vds)

    vds = (vds.annotate_variants_expr('va.calldata.Adj = gs.callStats(g => v), '
                                      'va.calldata.Hemi_Adj = gs.filter(g => sa.meta.sex == "male" && v.inXNonPar).callStats(g => v)')
           .annotate_variants_expr(generate_callstats_expression))

    vds = vds.drop_samples()

    vds = (vds.annotate_variants_expr(rearrange_callstats_expression)
           .annotate_variants_expr(
        'va.info.AF_raw = va.info.AC_raw.map(x => if (va.info.AN_raw > 0) x.toDouble/va.info.AN_raw else NA: Double), '
        'va.info.AF = va.info.AC.map(x => if (va.info.AN > 0) x.toDouble/va.info.AN else NA: Double)')
           .annotate_variants_expr(hom_hemi_expression)
           .annotate_variants_expr(ac_an_expression)
           .annotate_variants_expr(af_expression)
    )
    if drop_star:
        vds = vds.persist()
        vds = filter_star(vds, a_based=a_based_annotations, g_based=g_based_annotations, additional_annotations=star_annotations)
    vds = popmax(vds,pops)

    return vds.annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')


def create_sites_vds_annotations_Y(vds, pops, dbsnp_path=None, drop_star=True, filter_alleles = True, generate_hists = True):

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

    vds = vds.filter_variants_intervals(Interval.parse("Y"))

    vds = common_sites_vds_annotations(vds)

    if dbsnp_path is not None:
        vds = vds.annotate_variants_loci(dbsnp_path,
                                         locus_expr='Locus(_0,_1)',
                                         code='va.rsid=table._2',
                                         config=hail.TextTableConfig(noheader=True, comment="#", types='_0: String, _1: Int')
                                         )

    vds = (vds.filter_variants_expr('v.inYNonPar')
           .filter_samples_expr('sa.meta.sex == "male"')
           .filter_genotypes('g.isHet', keep=False)
           )

    if filter_alleles:
        vds = (vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
               .filter_alleles('va.calldata.raw.AC[aIndex] == 0', keep=False)
               .filter_variants_expr('v.nAltAlleles == 1 && v.alt == "*"', keep=False)
               )

    if generate_hists:
        vds = histograms(vds, 'va.info', AB=False)

    vds = vds.annotate_variants_expr('va.calldata.raw = gs.callStats(g => v)')
    vds = filter_to_adj(vds)
    vds = projectmax(vds)
    vds = vds.annotate_variants_expr('va.calldata.Adj = gs.callStats(g => v)')
    vds = unfurl_callstats(vds, criterion_pops, lower=True, gc=False)

    vds = vds.drop_samples()
    vds = (vds.annotate_variants_expr('va.info.AC_raw = va.calldata.raw.AC[1:], '
                                      'va.info.AN_raw = va.calldata.raw.AN, '
                                      'va.info.AF_raw = va.calldata.raw.AF[1:], '
                                      'va.info.AC = va.calldata.Adj.AC[1:], '
                                      'va.info.AN = va.calldata.Adj.AN, '
                                      'va.info.AF = va.calldata.Adj.AF[1:]')
           .annotate_variants_expr(correct_ac_an_command)
    )
    if drop_star:
        vds = vds.persist()
        vds = filter_star(vds, a_based=a_based_annotations, additional_annotations=star_annotations)
    vds = popmax(vds, pops)

    return vds.annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')


def pre_calculate_metrics(vds, output_file):
    AF_BUCKETS = [0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    info_metrics = [
        ('BaseQRankSum', -20, 20),
        ('ClippingRankSum', -5, 5),
        ('FS', 0, 20),
        ('InbreedingCoeff', -1, 1),
        ('MQ', 0, 60),
        ('MQRankSum', -20, 20),
        ('QD', 0, 40),
        ('ReadPosRankSum', -20, 20),
        ('VQSLOD', -20, 20)
    ]
    query_command = ['variants.map(v => va.info.%s).hist(%s, %s, 40)' % (metric, start, end) for metric, start, end in info_metrics]
    logged_info_metrics = [
        ('DP', 8, 20)
    ]
    query_command.extend(['variants.map(v => log(va.info.%s)).hist(%s, %s, 40)' % (metric, start, end) for metric, start, end in logged_info_metrics])

    as_metrics = [
        ('AS_RF', 0, 1),
        ('DP_MEDIAN', 0, 200),
        ('GQ_MEDIAN', 0, 100),
        ('AB_MEDIAN', 0, 1)
    ]
    query_command.extend(['variants.flatMap(v => va.info.%s).hist(%s, %s, 40)' % (metric, start, end) for metric, start, end in as_metrics])

    logged_as_metrics = [
        ('DREF_MEDIAN', 0, 100)
    ]
    query_command.extend(['variants.flatMap(v => va.info.%s.map(x => -log(x))).hist(%s, %s, 40)' % (metric, start, end) for metric, start, end in logged_as_metrics])

    site_quality_criteria = [
        ('binned_singleton', 'va.info.AC.exists(x => x == 1)'),
        ('binned_doubleton', 'va.info.AC.exists(x => x == 2)'),
        ('binned_0.00005', 'va.info.AF.exists(x => x < 0.00005)')
    ]
    for i, x in enumerate(AF_BUCKETS[1:]):
        site_quality_criteria.append(('binned_%s' % x, 'va.info.AF.exists(x => x >= %s) && va.info.AF.exists(x => x < %s)' % (AF_BUCKETS[i], x)))

    query_command.extend(['variants.filter(v => %s).map(v => log(va.qual)).hist(4, 20, 40)' % criteria[1] for criteria in site_quality_criteria])

    all_metrics = zip(*(info_metrics + logged_info_metrics + as_metrics + logged_as_metrics + site_quality_criteria))[0]
    all_results = vds.query_variants(query_command)

    final_results = []
    for metric, results in zip(all_metrics, all_results):
        output_results = {'metric': metric, 'mids': results['binEdges'], 'hist': results['binFrequencies']}
        output_results['hist'][0] += results['nLess']
        output_results['hist'][-1] += results['nGreater']
        final_results.append(output_results)

    with open(output_file, 'w') as g:
        g.write(json.dumps(final_results))


def annotate_from_rf(vds, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, annotations={}, train='va.train', label='va.label'):

    #Strip va if present
    rf_root = rf_root[len('va.'):] if rf_root.startswith('va.') else rf_root
    train = train[len('va.'):] if train.startswith('va.') else train
    label = label[len('va.'):] if label.startswith('va.') else label

    rf_ann_expr = [
        'va.info.AS_RF = orMissing(vds.exists(x => isDefined(x) && isDefined(x.%s)), '
        'range(v.nAltAlleles).map(i => orMissing(isDefined(vds[i]), vds[i].%s.probability.get("TP"))))' % (rf_root, rf_root),
        'va.info.AS_FilterStatus = '
        '   if(vds.forall(x => isMissing(x) || isMissing(x.%(root)s ))) range(v.nAltAlleles).map(i => ["RF"].toSet)'
        '   else range(v.nAltAlleles).map(i => '
        '       if(isMissing(vds[i]) || isMissing(vds[i].%(root)s)) ["RF"].toSet'
        '       else'
        '           if(v.altAlleles[i].isSNP)'
        '               if(vds[i].%(root)s.probability["TP"] > %(snv).4f) [""][:0].toSet else ["RF"].toSet'
        '           else'
        '               if(vds[i].%(root)s.probability["TP"] > %(indel).4f) [""][:0].toSet else ["RF"].toSet'
        '       )' %
                    {'root': rf_root,
                     'snv': rf_snv_cutoff,
                     'indel': rf_indel_cutoff},
        'va.info.AS_RF_POSITIVE_TRAIN = let x = range(v.nAltAlleles).filter('
        'i => isDefined(vds[i]) && isDefined(vds[i].%s) && isDefined(vds[i].%s) && vds[i].%s && vds[i].%s == "TP")'
        '.map(i => i+1) in orMissing(!x.isEmpty, x)' % (train, label, train, label),
        'va.info.AS_RF_POSITIVE_TRAIN = let x = range(v.nAltAlleles).filter('
        'i => isDefined(vds[i]) && isDefined(vds[i].%s) && isDefined(vds[i].%s) && vds[i].%s && vds[i].%s == "FP")'
        '.map(i => i+1) in orMissing(!x.isEmpty, x)' % (train, label, train, label)
    ]

    for source, target in annotations.iteritems():
        # Strip va if present
        source = source[len('va.'):] if source.startswith('va.') else source
        rf_ann_expr.append('%s = orMissing(vds.exists(x => isDefined(x) && isDefined(x.%s)),'
                           ' range(v.nAltAlleles).map(i => orMissing(isDefined(vds[i]), '
                           ' vds[i].%s)))' % (target, source, source))

    return vds.annotate_alleles_vds(rf_vds, rf_ann_expr)


def add_as_filters(vds, filters, root='va.info.AS_FilterStatus'):
    """
    Filters should be a dict of name: filter_expr
    Where i in the filter_expr is the alternate allele index (a-based)
    """

    as_filters = ",".join(['orMissing(%s,"%s")' % (filter_expr, name) for (name, filter_expr) in
                           filters.items()])
    as_filters = '[%s].filter(x => isDefined(x)).toSet' % as_filters

    input_dict = {
        'root': root,
        'filters': as_filters,
    }
    if not ann_exists(root, vds.variant_schema):
        vds = vds.annotate_variants_expr('%(root)s = range(v.nAltAlleles)'
                                         '.map(i => %(filters)s)' % input_dict)
    else:
        vds = vds.annotate_variants_expr('%(root)s = range(v.nAltAlleles).map(i => '
                                         'let new_filters = %(filters)s in '
                                         'if(isMissing(%(root)s) || isMissing(%(root)s[i])) '
                                         '  orMissing(!new_filters.isEmpty, new_filters) '
                                         'else %(root)s[i].union(%(filters)s))'
                                          % input_dict)
    return vds


def set_filters_attributes(vds, rf_snv_cutoff, rf_indel_cutoff):
    filters_desc = copy.deepcopy(FILTERS_DESC)
    if rf_snv_cutoff is not None and rf_indel_cutoff is not None and "RF" in filters_desc:
        filters_desc["RF"] = filters_desc["RF"] % (rf_snv_cutoff, rf_indel_cutoff)

    return vds.set_va_attributes('va.filters', filters_desc)


def run_samples_sanity_checks(vds, reference_vds, n_samples = 10, verbose=True):
    logger.info("Running samples sanity checks on %d samples" % n_samples)

    comparison_metrics = ['nHomVar',
                          'nSNP',
                          'nTransition',
                          'nTransversion',
                          'nInsertion',
                          'nDeletion',
                          'nNonRef',
                          'nHet'
                          ]

    samples = vds.sample_ids[:n_samples]

    def get_samples_metrics(vds, samples):
        metrics = (vds.filter_samples_expr('["%s"].toSet.contains(s)' % '","'.join(samples))
                   .sample_qc()
                   .query_samples('samples.map(s => {sample: s, metrics: sa.qc }).collect()')
                   )
        return {x.sample: x.metrics for x in metrics}

    test_metrics = get_samples_metrics(vds, samples)
    ref_metrics = get_samples_metrics(reference_vds, samples)

    output = ''

    for s, m in test_metrics.iteritems():
        if s not in ref_metrics:
            output += "WARN: Sample %s not found in reference data.\n" % s
        else:
            rm = ref_metrics[s]
            for metric in comparison_metrics:
                if m[metric] == rm[metric]:
                    if verbose:
                        output += "SUCCESS: Sample %s %s matches (N = %d).\n" % (s, metric, m[metric])
                else:
                    output += "FAILURE: Sample %s, %s differs: Data: %s, Reference: %s.\n" % (
                    s, metric, m[metric], rm[metric])

    logger.info(output)
    return output


def run_sites_sanity_checks(vds, pops, verbose=True, contig='auto', percent_missing_threshold=0.01,
                            skip_star=False, split_lcr=False, split_star = False):
    logger.info("Running sites sanity checks on contig %s" % contig)

    output = ''

    #Grouped by filters
    ## By allele type
    pre_split_ann = get_variant_type_expr('va.final_variantType')
    pre_split_ann += ',va.hasStar = v.altAlleles.exists(a => a.isStar)'
    pre_split_ann += ',va.nAltAlleles = v.nAltAlleles'

    agg_key = 'type = if(v.altAllele.isSNP) "snv" else if(v.altAllele.isIndel) "indel" else "other"'
    select_cols = ['type',
                   'n',
                   bround('prop_filtered', 3).alias('All'),
                   bround('prop_hard_filtered', 3).alias('HF'),
                   bround('prop_AC0_filtered', 3).alias('AC0'),
                   bround('prop_RF_filtered', 3).alias('RF'),
                   bround('prop_hard_filtered_only', 3).alias('HF only'),
                   bround('prop_AC0_filtered_only', 3).alias('AC0 only'),
                   bround('prop_RF_filtered_only', 3).alias('RF only')]

    if split_lcr:
        agg_key += ', LCR = va.lcr'
        select_cols.append('LCR')

    df = (
        vds
            .annotate_variants_expr(pre_split_ann)
            .split_multi()
            .variants_keytable().aggregate_by_key(agg_key,
                                                  'n = va.count(), '
                                                                'prop_filtered = va.filter(x => !x.filters.isEmpty || !x.info.AS_FilterStatus[x.aIndex - 1].isEmpty).count(),'
                                                                'prop_hard_filtered = va.filter(x => x.filters.contains("LCR") || x.filters.contains("SEGDUP") || x.filters.contains("InbreedingCoeff")).count(),'
                                                                'prop_AC0_filtered = va.filter(x => x.info.AS_FilterStatus[x.aIndex - 1].contains("AC0")).count(),'
                                                                'prop_RF_filtered = va.filter(x => x.info.AS_FilterStatus[x.aIndex - 1].contains("RF")).count(),'
                                                                'prop_hard_filtered_only = va.filter(x => (x.filters.contains("LCR") || x.filters.contains("SEGDUP") || x.filters.contains("InbreedingCoeff")) && x.info.AS_FilterStatus[x.aIndex - 1].isEmpty).count(),'
                                                                'prop_AC0_filtered_only = va.filter(x => x.filters.forall(f => f == "AC0") &&  !x.info.AS_FilterStatus[x.aIndex - 1].isEmpty && x.info.AS_FilterStatus[x.aIndex - 1].forall(f => f == "AC0")).count(),'
                                                                'prop_RF_filtered_only = va.filter(x => x.filters.forall(f => f == "RF") &&  !x.info.AS_FilterStatus[x.aIndex - 1].isEmpty && x.info.AS_FilterStatus[x.aIndex - 1].forall(f => f == "RF")).count()')
            .annotate(['%s = %s / n' % (ann, ann) for ann in ['prop_filtered','prop_hard_filtered','prop_AC0_filtered','prop_RF_filtered','prop_hard_filtered_only','prop_AC0_filtered_only','prop_RF_filtered_only']])
            .to_dataframe()
    )

    df = df.select(select_cols)

    # At some point should print output too
    # pd = df.toPandas()
    # if return_string:
    #     output += "\nProportion of sites filtered by allele type:\n%s" % str(pd)

    print("\nProportion of sites filtered by allele type:\n")
    df.show()
    # By nAltAlleles

    agg_key = 'type = va.final_variantType, nAltAlleles = va.nAltAlleles'
    select_cols = ['type',
                   'nAltAlleles',
                   'n',
                   bround('prop_filtered', 3).alias('All'),
                   bround('prop_hard_filtered', 3).alias('HF'),
                   bround('prop_AC0_filtered', 3).alias('AC0'),
                   bround('prop_RF_filtered', 3).alias('RF'),
                   bround('prop_hard_filtered_only', 3).alias('HF only'),
                   bround('prop_AC0_filtered_only', 3).alias('AC0 only'),
                   bround('prop_RF_filtered_only', 3).alias('RF only')]

    order_by_cols = ["type","nAltAlleles"]

    if split_lcr:
        agg_key += ', LCR = va.lcr'
        select_cols.insert(0,'LCR')
        order_by_cols.append("LCR")

    if split_star:
        agg_key += ', hasStar = va.hasStar'
        select_cols.append('hasStar')
        order_by_cols.append("hasStar")

    df = (
        vds
            .annotate_variants_expr(pre_split_ann)
            .split_multi()
            .variants_keytable().aggregate_by_key(agg_key,
                                                  'n = va.count(), '
                                                                'prop_filtered = va.filter(x => !x.filters.isEmpty || !x.info.AS_FilterStatus[x.aIndex - 1].isEmpty).count(),'
                                                                'prop_hard_filtered = va.filter(x => x.filters.contains("LCR") || x.filters.contains("SEGDUP") || x.filters.contains("InbreedingCoeff")).count(),'
                                                                'prop_AC0_filtered = va.filter(x => x.info.AS_FilterStatus[x.aIndex - 1].contains("AC0")).count(),'
                                                                'prop_RF_filtered = va.filter(x => x.info.AS_FilterStatus[x.aIndex - 1].contains("RF")).count(),'
                                                                'prop_hard_filtered_only = va.filter(x => (x.filters.contains("LCR") || x.filters.contains("SEGDUP") || x.filters.contains("InbreedingCoeff")) && x.info.AS_FilterStatus[x.aIndex - 1].isEmpty).count(),'
                                                                'prop_AC0_filtered_only = va.filter(x => x.filters.forall(f => f == "AC0") &&  !x.info.AS_FilterStatus[x.aIndex - 1].isEmpty && x.info.AS_FilterStatus[x.aIndex - 1].forall(f => f == "AC0")).count(),'
                                                                'prop_RF_filtered_only = va.filter(x => x.filters.forall(f => f == "RF") &&  !x.info.AS_FilterStatus[x.aIndex - 1].isEmpty && x.info.AS_FilterStatus[x.aIndex - 1].forall(f => f == "RF")).count()')
            .annotate(['%s = %s / n' % (ann, ann) for ann in
                       ['prop_filtered', 'prop_hard_filtered', 'prop_AC0_filtered', 'prop_RF_filtered',
                        'prop_hard_filtered_only', 'prop_AC0_filtered_only', 'prop_RF_filtered_only']])
            .to_dataframe()
            .orderBy(order_by_cols)
    )

    df = df.select(select_cols)

    print("\nProportion of sites filtered by site type and number of alt alleles:\n")
    df.show(n = 200)


    df = (
        vds
            .annotate_variants_expr(pre_split_ann)
            .variants_keytable().aggregate_by_key(agg_key,
                                                  'n = va.count(), '
                                                                'prop_filtered = va.filter(x => !x.filters.isEmpty).count(),'
                                                                'prop_hard_filtered = va.filter(x => x.filters.contains("LCR") || x.filters.contains("SEGDUP") || x.filters.contains("InbreedingCoeff")).count(),'
                                                                'prop_AC0_filtered = va.filter(x => x.filters.contains("AC0")).count(),'
                                                                'prop_RF_filtered = va.filter(x => x.filters.contains("RF")).count(),'
                                                                'prop_hard_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => ["LCR","SEGDUP","InbreedingCoeff"].toSet.contains(f)) ).count(),'
                                                                'prop_AC0_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "AC0") ).count(),'
                                                                'prop_RF_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "RF") ).count()'
                                                                )
            .annotate(['%s = %s / n' % (ann, ann) for ann in
                       ['prop_filtered', 'prop_hard_filtered', 'prop_AC0_filtered', 'prop_RF_filtered',
                        'prop_hard_filtered_only', 'prop_AC0_filtered_only', 'prop_RF_filtered_only']])
            .to_dataframe()
            .orderBy(order_by_cols)
    )

    df = df.select(select_cols)

    print("\nProportion of sites filtered by variant type and number of alt alleles:\n")
    df.show(n = 200)

    #At some point should print output too
    # pd = df.toPandas()
    # if return_string:
    #     output += "\nProportion of sites filtered by variant type and number of alt alleles:\n%s" % str(pd)
    # print("\nProportion of sites filtered by variant type and number of alt alleles:\n%s" % str(pd))


    queries = ['variants.count()']
    a_metrics = ['AC', 'Hom']

    if contig == 'Y':
        a_metrics = a_metrics[:1]

    one_metrics = ['AN'] if skip_star else ['STAR_AC','AN']

    # Filter counts
    queries.extend(["let x = variants.map(v => !va.filters.isEmpty).counter() in orElse(x.get(true), 0L)/x.values().sum()",
                    'variants.map(v => va.filters.toArray.mkString(",")).counter()'])

    # Check number of samples
    queries.append('variants.map(v => va.info.AN).stats().max/2')
    queries.extend(['variants.map(v => va.info.AN_%s).stats().max/2' % pop for pop in pops])

    end_pop_counts = len(queries)

    # Check that raw is always larger than adj
    queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                    '.exists(i => va.info.%s[i] > va.info.%s_raw[i])).count()' % (x, x) for x in a_metrics])
    queries.extend(['variants.filter(v => va.info.%s > va.info.%s_raw).count()' % (x, x) for x in one_metrics])

    # Check that sum(pops) == total
    for metric in a_metrics:
        queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                        '.exists(i => %s != va.info.%s[i])).count()' % (" + ".join(["va.info.%s_%s[i]" % (metric, pop) for pop in pops]), metric)])

    queries.append('variants.filter(v => %s != va.info.AN).count()' % " + ".join(["va.info.AN_%s" % pop for pop in pops]))

    # Check that male + female == total
    # Remove Hom for X
    if contig == 'X':
        a_metrics = a_metrics[:1]

    if contig != 'Y':
        pop_strats = pops if contig == 'X' else [None]
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
        queries.append('let x = variants.map(v => isMissing(va.info.%s)).counter() in orElse(x.get(true), 0L)/x.values().sum()' % field.name)

    logger.debug('Queries: %s', '\n'.join(['%s: %s' % (i, x) for i, x in enumerate(queries)]))
    stats = vds.query_variants(queries)

    # Print filters

    # Double checking for no samples in VDS
    sample_count = vds.query_samples('samples.count()')
    if sample_count > 0:
        output += 'WARNING: %s samples found in VDS (should be 0)\n' % sample_count

    output += "Total number of sites:\n"
    output += pformat(int(stats[0])) + '\n'
    output += "FILTERS CHECKS\nTotal fraction sites filtered:\n"
    output += pformat(stats[1]) + '\n'
    output += "Filter counts:\n"
    output += pformat(stats[2]) + '\n'

    output += "\nPOPULATION COUNTS\n"
    output += "Total number of samples: %s\n" % stats[3]
    for i in range(4, end_pop_counts):
        output += '%s: %s\n' % (pops[i-4], stats[i])

    #Check that all metrics sum as expected
    output += "\nMETRICS COUNTS CHECK\n"
    nfail = 0
    for i in range(end_pop_counts, end_counts):
        if stats[i] != 0:
            output += "FAILED METRICS CHECK for query: %s\n Expected: 0, Found: %s\n" % (queries[i], stats[i])
            nfail += 1
        elif verbose:
            output += "Success: %s\n" % queries[i]
    output += "%s metrics count checks failed.\n" % nfail

    #Check missing metrics
    output += "\nMISSING METRICS CHECKS\n"
    nfail = 0
    missing_stats = stats[end_counts:]
    for i in range(len(missing_stats)):
        if missing_stats[i] > percent_missing_threshold:
            output += "FAILED missing check for %s; %s%% missing.\n" % (missing_metrics[i], 100*missing_stats[i])
            nfail += 1
        elif verbose:
            output += "SUCCESS missing check for %s; %s%% missing.\n" % (missing_metrics[i], 100*missing_stats[i])
    output += "%s missing metrics checks failed.\n" % nfail

    logger.info(output)
    return output


def set_va_attributes(vds, warn_if_not_found = True):

    info_va_attr = get_info_va_attr()
    va_info = [x for x in vds.variant_schema.fields if x.name == "info"][0]

    for ann in va_info.typ.fields:
        if ann.name in info_va_attr:
            vds = vds.set_va_attributes("va.info.%s" % ann.name, info_va_attr[ann.name])

        elif ann.name != "CSQ" and warn_if_not_found: logger.warn("No description found for va.info.%s", ann.name)

    return vds


def write_public_vds(vds, public_path, overwrite=False):
    vds = vds.annotate_samples_expr('sa = {}')
    vds = vds.annotate_variants_expr('va = select(va, rsid, qual, filters, pass, info)')
    vds = vds.annotate_variants_expr('va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
    vds.write(public_path, overwrite=overwrite)


def merge_schemas(vdses):

    vds_schemas = [vds.variant_schema for vds in vdses]

    for s in vds_schemas[1:]:
        if type(vds_schemas[0]) != type(s):
            logger.fatal("Cannot merge schemas as the root (va) is of different type: %s and %s", vds_schemas[0], s)
            sys.exit(1)

    if not isinstance(vds_schemas[0],hail.expr.TStruct):
        return(vdses)

    anns = [flatten_struct(s, root='va') for s in vds_schemas]

    all_anns = {}
    for i in reversed(range(len(vds_schemas))):
        common_keys = set(all_anns.keys()).intersection(anns[i].keys())
        for k in common_keys:
            if type(all_anns[k].typ) != type(anns[i][k].typ):
                logger.fatal(
                    "Cannot merge schemas as annotation %s type %s found in VDS %d is not the same as previously existing type %s"
                    % (k, anns[i][k].typ, i, all_anns[k].typ))
                sys.exit(1)
        all_anns.update(anns[i])

    for i in range(len(vdses)):
        vdses[i] = vdses[i].annotate_variants_expr(["%s = NA: %s" % (k, str(v.typ)) for k,v in
                                                    all_anns.iteritems() if k not in anns[i]])
        for ann, f in all_anns.iteritems():
            vdses[i] = vdses[i].set_va_attributes(ann, f.attributes)

    return vdses


def copy_schema_attributes(vds1, vds2):
    anns1 = flatten_struct(vds1.variant_schema, root='va')
    anns2 = flatten_struct(vds2.variant_schema, root='va')
    for ann in anns1.keys():
        if ann in anns2:
            vds1 = vds1.set_va_attributes(ann, anns2[ann].attributes)

    return vds1


def print_attributes(vds, path=None):
    anns = flatten_struct(vds.variant_schema, root='va')
    if path is not None:
        print "%s attributes: %s" % (path, anns[path].attributes)
    else:
        for ann, f in anns.iteritems():
            print "%s attributes: %s" % (ann, f.attributes)


def get_numbered_annotations(vds, root='va.info'):
    """
    Get all 1-, A-, G- numbered annotations from a VDS based on the Number va attributes.
    In addition returns arrays with no Number or Number=. va attribute separately
    :param vds: Input VDS
    :param root: Place to find annotations (defaults to va.info)
    :return: annotations, a_annotations, g_annotations, dot_annotations as list[Field]
    """
    a_annotations = []
    g_annotations = []
    dot_annotations = []
    annotations = []

    release_info = get_ann_type(root, vds.variant_schema)
    for field in release_info.fields:
        if isinstance(field.typ, TArray):
            if 'Number' in field.attributes:
                number = field.attributes['Number']
                if number == "A":
                    a_annotations.append(field)
                elif number == "G":
                    g_annotations.append(field)
                else:
                    dot_annotations.append(field)
        else:
            annotations.append(field)

    logger.info("Found the following fields:")
    logger.info("1-based annotations: " + ",".join([x.name for x in annotations]))
    logger.info("A-based annotations: " + ",".join([x.name for x in a_annotations]))
    logger.info("G-based annotations: " + ",".join([x.name for x in g_annotations]))
    logger.info("dot annotations: " + ",".join([x.name for x in dot_annotations]))

    return annotations, a_annotations, g_annotations, dot_annotations


def filter_annotations_regex(annotation_fields, ignore_list):
    def ann_in(name, list):
        # `list` is a list of regexes to ignore
        return any([x for x in list if re.search('^%s$' % x, name)])

    return [x for x in annotation_fields if not ann_in(x.name, ignore_list)]


def annotate_subset_with_release(subset_vds, release_dict, root="va.info", dot_annotations_dict = None, ignore = None, annotate_g_annotations = False):

    parsed_root = root.split(".")
    if parsed_root[0] != "va":
        logger.error("Found va annotation root not starting with va: %s", root)
    ann_root = ".".join(parsed_root[1:])

    annotations, a_annotations, g_annotations, dot_annotations = get_numbered_annotations(release_dict['vds'], root)

    if ignore is not None:
        annotations = filter_annotations_regex(annotations, ignore)
        a_annotations = filter_annotations_regex(a_annotations, ignore)
        g_annotations = filter_annotations_regex(g_annotations, ignore)
        dot_annotations = filter_annotations_regex(dot_annotations, ignore)

    annotation_expr = ['%s = vds.find(x => isDefined(x)).%s.%s' % (release_dict['out_root'] + ann.name, ann_root, ann.name) for ann in annotations]
    annotation_expr.extend(['%s = orMissing(vds.exists(x => isDefined(x)  && isDefined(x.%s.%s)), range(v.nAltAlleles)'
                            '.map(i => orMissing( isDefined(vds[i]), vds[i].%s.%s[aIndices[i]] )))'
                            % (release_dict['out_root'] + ann.name, ann_root, ann.name, ann_root, ann.name) for ann in a_annotations ])

    if annotate_g_annotations:
        annotation_expr.extend([
            '%s = orMissing(vds.exists(x => isDefined(x) && isDefined(x.%s.%s)), '
            'range(gtIndex(v.nAltAlleles,v.nAltAlleles)).map(i => let j = gtj(i) and k = gtk(i) and'
            'aj = if(j==0) 0 else aIndices[j-1]+1 and ak = if(k==0) 0 else aIndices[k-1]+1 in '
            'orMissing( isDefined(aj) && isDefined(ak),'
            'vds.find(x => isDefined(x)).%s.%s[ gtIndex(aj, ak)])))'
            % (release_dict['out_root'] + ann.name, ann_root, ann.name, ann_root, ann.name) for ann in g_annotations])

    if dot_annotations_dict is not None:
        for ann in dot_annotations:
            if ann in dot_annotations_dict:
                annotation_expr.append(dot_annotations_dict[ann.name] % (release_dict['out_root'] + ann.name))

    logger.debug("Annotating subset with the following expr:\n" + ",\n".join(annotation_expr))

    subset_vds = subset_vds.annotate_alleles_vds(release_dict['vds'], annotation_expr)

    #Set attributes for all annotations
    annotations.extend(a_annotations)
    if annotate_g_annotations:
        annotations.extend(g_annotations)

    if dot_annotations_dict is not None:
        for ann in dot_annotations:
            if ann in dot_annotations_dict:
                annotations.append(ann)

    for ann in annotations:
        attributes = {}
        for k,v in ann.attributes.iteritems():
            if k == "Description":
                v = "%s (source: %s)" % (v, v)
            attributes[k] = v

        subset_vds = subset_vds.set_va_attributes(release_dict['out_root'] + ann.name, attributes)

    return subset_vds


def pc_project(vds, pc_vds, pca_loadings_root = 'va.pca_loadings'):
    """
    Projects samples in `vds` on PCs computed in `pc_vds`
    :param vds: VDS containing the samples to project
    :param pc_vds: VDS containing the PC loadings for the variants
    :param pca_loadings_root: Annotation root for the loadings. Can be either an Array[Double] or a Struct{ PC1: Double, PC2: Double, ...}
    :return: VDS with
    """

    pca_loadings_type = get_ann_type(pca_loadings_root,pc_vds.variant_schema)


    pc_vds = pc_vds.annotate_variants_expr('va.pca.calldata = gs.callStats(g => v)')

    pcs_struct_to_array = ",".join(['vds.pca_loadings.PC%d' % x for x in range(1,21)])
    arr_to_struct_expr = ",".join(['PC%d: sa.pca[%d - 1]' % (x, x) for x in range(1, 21)])



    vds = (vds.filter_multi()
            .annotate_variants_vds(pc_vds, code = 'va.pca_loadings = [%s], va.pca_af = vds.pca.calldata.AF[1]' % pcs_struct_to_array)
           .filter_variants_expr('!isMissing(va.pca_loadings) && !isMissing(va.pca_af)')
     )

    n_variants = vds.query_variants(['variants.count()'])[0]

    return(vds
            .annotate_samples_expr('sa.pca = gs.filter(g => g.isCalled && va.pca_af > 0.0 && va.pca_af < 1.0).map(g => let p = va.pca_af in (g.gt - 2 * p) / sqrt(%d * 2 * p * (1 - p)) * va.pca_loadings).sum()' % n_variants)
            .annotate_samples_expr('sa.pca = {%s}' % arr_to_struct_expr)
    )


def read_list_data(input_file):
    if input_file.startswith('gs://'):
        hail.hadoop_copy(input_file, 'file:///' + input_file.split("/")[-1])
        f = gzip.open("/" + os.path.basename(input_file)) if input_file.endswith('gz') else open( "/" + os.path.basename(input_file))
    else:
        f = gzip.open(input_file) if input_file.endswith('gz') else open(input_file)
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def rename_samples(vds, input_file, filter_to_samples_in_file = False):
    names = {old: new for old,new in [x.split("\t") for x in read_list_data(input_file)]}
    logger.info("Found %d samples for renaming in input file.")
    logger.info("Renaming %d samples found in VDS" % len(set(names.keys()).intersection(set(vds.sample_ids)) ))

    if filter_to_samples_in_file:
        vds = vds.filter_samples_list(names.keys())
    return vds.rename_samples(names)



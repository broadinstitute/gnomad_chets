
import re
import sys
import json
import copy
import logging
import gzip
import os

from resources import *
from hail import *
from slack_utils import *
from pyspark.sql.functions import bround
from pprint import pprint, pformat

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)

GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
EXOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

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
CSQ_CODING_HIGH_IMPACT = ["transcript_ablation",
"splice_acceptor_variant",
"splice_donor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost"]

CSQ_CODING_MEDIUM_IMPACT = [
"start_lost",  # new in v81
"initiator_codon_variant",  # deprecated
"transcript_amplification",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"protein_altering_variant",  # new in v79
"splice_region_variant"
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
"stop_retained_variant",
"synonymous_variant",
"coding_sequence_variant"]

CSQ_NON_CODING = [
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
"intergenic_variant"
]

CSQ_ORDER = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING


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


def index_into_arrays(a_based_annotations=None, r_based_annotations=None, vep_root=None, drop_ref_ann = False):
    """

    Creates annotation expressions to get the correct values when splitting multi-allelics

    :param list of str a_based_annotations: A-based annotations
    :param list of str r_based_annotations: R-based annotations
    :param str vep_root: Root of the vep annotation
    :param bool drop_ref_ann: If set to True, then the reference value of R-based annotations is removed (effectively converting them in A-based annotations)
    :return: Annotation expressions
    :rtype: list of str
    """
    annotations = []
    if a_based_annotations:
        for ann in a_based_annotations:
            annotations.append('{0} = {0}[va.aIndex - 1]'.format(ann))
    if r_based_annotations:
        expr = '{0} = {0}[va.aIndex]' if drop_ref_ann else '{0} = [{0}[0], {0}[va.aIndex]]'
        for ann in r_based_annotations:
            annotations.append(expr.format(ann))
    if vep_root:
        sub_fields = ['transcript_consequences', 'intergenic_consequences', 'motif_feature_consequences', 'regulatory_feature_consequences']
        annotations.extend(['{0}.{1} = {0}.{1}.filter(x => x.allele_num == va.aIndex)'.format(vep_root, sub_field) for sub_field in sub_fields])

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
        if isinstance(additional_annotations, str):
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


def filter_to_adj(vds):
    return vds.filter_genotypes(ADJ_CRITERIA)


def filter_star(vds, a_based=None, r_based=None, g_based=None, additional_annotations=None):
    annotation = unfurl_filter_alleles_annotation(a_based=a_based, r_based=r_based, g_based=g_based,
                                                  additional_annotations=additional_annotations)
    return vds.filter_alleles('v.altAlleles[aIndex - 1].alt == "*"', annotation=annotation, keep=False)


def flatten_struct(struct, root='', leaf_only=True):
    result = {}
    for f in struct.fields:
        path = '%s.%s' % (root, f.name)
        if isinstance(f.typ, TStruct):
            result.update(flatten_struct(f.typ, path))
            if not leaf_only:
                result[path] = f
        else:
            result[path] = f
    return result


def ann_exists(annotation, schema, root='va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    return annotation in anns


def get_ann_field(annotation, schema, root='va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    if not annotation in anns:
        logger.error("%s missing from schema.", annotation)
        sys.exit(1)
    return anns[annotation]


def get_ann_type(annotation, schema, root='va'):
    return get_ann_field(annotation, schema, root).typ


def annotation_type_is_numeric(t):
    return (isinstance(t, TInt) or
            isinstance(t, TLong) or
            isinstance(t, TFloat) or
            isinstance(t, TDouble)
            )


def get_variant_type_expr(root="va.variantType"):
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
            "mixed"''' % root


def get_allele_stats_expr(root="va.stats", medians=False, samples_filter_expr=''):
    """

    Gets allele-specific stats expression: GQ, DP, NRQ, AB, Best AB, p(AB), NRDP, QUAL, combined p(AB)

    :param str root: annotations root
    :param bool medians: Calculate medians for GQ, DP, NRQ, AB and p(AB)
    :param str samples_filter_expr: Expression for filtering samples (e.g. "sa.keep")
    :return: List of expressions for `annotate_alleles_expr`
    :rtype: list of str
    """

    if samples_filter_expr:
        samples_filter_expr = "&& " + samples_filter_expr

    stats = ['%s.gq = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).stats()',
             '%s.dp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).stats()',
             '%s.nrq = gs.filter(g => g.isCalledNonRef %s).map(g => -log10(g.gp[0])).stats()',
             '%s.ab = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).stats()',
             '%s.best_ab = gs.filter(g => g.isHet %s).map(g => abs((g.ad[1]/g.dp) - 0.5)).min()',
             '%s.pab = gs.filter(g => g.isHet %s).map(g => g.pAB()).stats()',
             '%s.nrdp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).sum()',
             '%s.qual = -10*gs.filter(g => g.isCalledNonRef %s).map(g => if(g.pl[0] > 3000) -300 else log10(g.gp[0])).sum()',
             '%s.combined_pAB = let hetSamples = gs.filter(g => g.isHet %s).map(g => log(g.pAB())).collect() in orMissing(!hetSamples.isEmpty, -10*log10(pchisqtail(-2*hetSamples.sum(),2*hetSamples.length)))']

    if medians:
        stats.extend(['%s.gq_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.gq).collect().median',
                    '%s.dp_median = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).collect().median',
                    '%s.nrq_median = gs.filter(g => g.isCalledNonRef %s).map(g => -log10(g.gp[0])).collect().median',
                    '%s.ab_median = gs.filter(g => g.isHet %s).map(g => g.ad[1]/g.dp).collect().median',
                    '%s.pab_median = gs.filter(g => g.isHet %s).map(g => g.pAB()).collect().median'])

    stats_expr = [x % (root, samples_filter_expr) for x in stats]

    return stats_expr


def set_site_filters(vds, site_filters_dict, filters_to_keep=[], as_filters_root="va.info.AS_FilterStatus"):
    site_filters = ",".join(['orMissing(%s, "%s")' % (filter_expr, name) for (name, filter_expr) in site_filters_dict.items()])
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


def post_process_subset(subset_vds, release_vds_dict, as_filters_key, dot_annotations_dict=None):

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
        'va.info.AS_RF_NEGATIVE_TRAIN = let x = range(v.nAltAlleles).filter('
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


def run_samples_sanity_checks(vds, reference_vds, n_samples=10, verbose=True):
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


def set_va_attributes(vds, warn_if_not_found=True):

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
        if not isinstance(vds_schemas[0], type(s)):
            logger.fatal("Cannot merge schemas as the root (va) is of different type: %s and %s", vds_schemas[0], s)
            sys.exit(1)

    if not isinstance(vds_schemas[0], TStruct):
        return vdses

    anns = [flatten_struct(s, root='va') for s in vds_schemas]

    all_anns = {}
    for i in reversed(range(len(vds_schemas))):
        common_keys = set(all_anns.keys()).intersection(anns[i].keys())
        for k in common_keys:
            if not isinstance(all_anns[k].typ, type(anns[i][k].typ)):
                logger.fatal(
                    "Cannot merge schemas as annotation %s type %s found in VDS %d is not the same as previously existing type %s"
                    % (k, anns[i][k].typ, i, all_anns[k].typ))
                sys.exit(1)
        all_anns.update(anns[i])

    for i, vds in enumerate(vdses):
        vds = vds.annotate_variants_expr(["%s = NA: %s" % (k, str(v.typ)) for k, v in
                                          all_anns.iteritems() if k not in anns[i]])
        for ann, f in all_anns.iteritems():
            vds = vds.set_va_attributes(ann, f.attributes)

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
    def ann_in(name, lst):
        # `list` is a list of regexes to ignore
        return any([x for x in lst if re.search('^%s$' % x, name)])

    return [x for x in annotation_fields if not ann_in(x.name, ignore_list)]


def annotate_subset_with_release(subset_vds, release_dict, root="va.info", dot_annotations_dict = None, ignore = None, annotate_g_annotations=False):

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


def pc_project(vds, pc_vds, pca_loadings_root='va.pca_loadings'):
    """
    Projects samples in `vds` on PCs computed in `pc_vds`
    :param vds: VDS containing the samples to project
    :param pc_vds: VDS containing the PC loadings for the variants
    :param pca_loadings_root: Annotation root for the loadings. Can be either an Array[Double] or a Struct{ PC1: Double, PC2: Double, ...}
    :return: VDS with
    """

    pca_loadings_type = get_ann_type(pca_loadings_root, pc_vds.variant_schema)  # TODO: this isn't used?

    pc_vds = pc_vds.annotate_variants_expr('va.pca.calldata = gs.callStats(g => v)')

    pcs_struct_to_array = ",".join(['vds.pca_loadings.PC%d' % x for x in range(1, 21)])
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
        hadoop_copy(input_file, 'file:///' + input_file.split("/")[-1])
        f = gzip.open("/" + os.path.basename(input_file)) if input_file.endswith('gz') else open( "/" + os.path.basename(input_file))
    else:
        f = gzip.open(input_file) if input_file.endswith('gz') else open(input_file)
    output = []
    for line in f:
        output.append(line.strip())
    f.close()
    return output


def rename_samples(vds, input_file, filter_to_samples_in_file=False):
    names = {old: new for old, new in [x.split("\t") for x in read_list_data(input_file)]}
    logger.info("Found %d samples for renaming in input file %s." % (len(names.keys()), input_file))
    logger.info("Renaming %d samples found in VDS" % len(set(names.keys()).intersection(set(vds.sample_ids)) ))

    if filter_to_samples_in_file:
        vds = vds.filter_samples_list(names.keys())
    return vds.rename_samples(names)


def add_genomes_sa(vds):
    """
    Adds the genomes sample metadata to the VDS.

    :param VariantDataset vds: VDS to annotate
    :return: Annotated VDS.
    :rtype: VariantDataset
    """
    hc = vds.hc
    vds = vds.annotate_samples_table(KeyTable.import_fam(genomes_fam_path), root='sa.fam')
    vds = vds.annotate_samples_table(hc.import_table(genomes_meta_tsv_path, impute=True).key_by('Sample'), root='sa.meta')
    vds = vds.annotate_samples_table(
        hc.import_table(genomes_to_combined_IDs_tsv_path, impute=True, no_header=True).key_by('f0').select(['f0']),
        root='sa.in_exomes')
    vds = vds.annotate_samples_table(hc.import_table(genomes_qc_pass_samples_list_path, impute=True).key_by('sample'), root='sa.qc_pass')
    return vds


def add_exomes_sa(vds):
    """
    Adds the exomes sample metadata to the VDS.

    :param VariantDataset vds: VDS to annotate 
    :return: Annotated VDS.
    :rtype: VariantDataset
    """
    hc = vds.hc
    vds = vds.annotate_samples_table(KeyTable.import_fam(exomes_fam_path), root='sa.fam')
    vds = vds.annotate_samples_table(hc.import_table(exomes_meta_tsv_path, impute=True).key_by('sample'), root='sa.meta')
    vds = vds.annotate_samples_table(
        hc.import_table(exomes_to_combined_IDs_tsv_path, impute=True, no_header=True).key_by('f0').select(['f0']),
        root='sa.in_genomes')
    vds = vds.annotate_samples_table(hc.import_table(exomes_qc_pass_samples_list_path, impute=True).key_by('sample'), root='sa.qc_pass')
    return vds


def filter_low_conf_regions(vds, filter_lcr=True, filter_decoy=True, high_conf_regions=None):
    """
    Filters low-confidence regions

    :param VariantDataset vds: VDS to filter
    :param bool filter_lcr: Whether to filter LCR regions
    :param bool filter_decoy: Wheter to filter Segdup regions
    :param list of str high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return:
    """

    if filter_lcr:
        vds = vds.filter_variants_table(KeyTable.import_interval_list(lcr_intervals_path), keep=False)

    if filter_decoy:
        vds = vds.filter_variants_table(KeyTable.import_interval_list(decoy_intervals_path), keep=False)

    if high_conf_regions is not None:
        for region in high_conf_regions:
            vds = vds.filter_variants_table(KeyTable.import_interval_list(region), keep=True)

    return vds


def process_consequences(vds, vep_root='va.vep'):
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq and worst_csq_suffix (worst consequence across transcripts) into [vep_root]

    :param VariantDataset vds: Input VDS
    :param str vep_root: Root for vep annotation (probably va.vep)
    :return: VDS with better formatted consequences
    :rtype: VariantDataset
    """
    if vep_root + '.worst_csq' in flatten_struct(vds.variant_schema, root='va'):
        vds = (vds.annotate_variants_expr('%(vep)s.transcript_consequences = '
                                          ' %(vep)s.transcript_consequences.map('
                                          '     csq => drop(csq, most_severe_consequence)'
                                          ')' % {'vep': vep_root}))
    vds = (vds.annotate_global('global.csqs', CSQ_ORDER, TArray(TString()))
           .annotate_variants_expr(
        '%(vep)s.transcript_consequences = '
        '   %(vep)s.transcript_consequences.map(csq => '
        '   let worst_csq = global.csqs.find(c => csq.consequence_terms.toSet().contains(c)) in'
        # '   let worst_csq_suffix = if (csq.filter(x => x.lof == "HC").length > 0)'
        # '       worst_csq + "-HC" '
        # '   else '
        # '       if (csq.filter(x => x.lof == "LC").length > 0)'
        # '           worst_csq + "-LC" '
        # '       else '
        # '           if (csq.filter(x => x.polyphen_prediction == "probably_damaging").length > 0)'
        # '               worst_csq + "-probably_damaging"'
        # '           else'
        # '               if (csq.filter(x => x.polyphen_prediction == "possibly_damaging").length > 0)'
        # '                   worst_csq + "-possibly_damaging"'
        # '               else'
        # '                   worst_csq in'
        '   merge(csq, {most_severe_consequence: worst_csq'
        # ', most_severe_consequence_suffix: worst_csq_suffix'
        '})'
        ')' % {'vep': vep_root}
    ).annotate_variants_expr(
        '%(vep)s.worst_csq = global.csqs.find(c => %(vep)s.transcript_consequences.map(x => x.most_severe_consequence).toSet().contains(c)),'
        '%(vep)s.worst_csq_suffix = '
        'let csq = global.csqs.find(c => %(vep)s.transcript_consequences.map(x => x.most_severe_consequence).toSet().contains(c)) in '
        'if (%(vep)s.transcript_consequences.filter(x => x.lof == "HC").length > 0)'
        '   csq + "-HC" '
        'else '
        '   if (%(vep)s.transcript_consequences.filter(x => x.lof == "LC").length > 0)'
        '       csq + "-LC" '
        '   else '
        '       if (%(vep)s.transcript_consequences.filter(x => x.polyphen_prediction == "probably_damaging").length > 0)'
        '           csq + "-probably_damaging"'
        '       else'
        '           if (%(vep)s.transcript_consequences.filter(x => x.polyphen_prediction == "possibly_damaging").length > 0)'
        '               csq + "-possibly_damaging"'
        '           else'
        '               csq' % {'vep': vep_root}
    ))
    return vds


def filter_vep_to_canonical_transcripts(vds, vep_root='va.vep'):
    return vds.annotate_variants_expr(
        '%(vep)s.transcript_consequences = '
        '   %(vep)s.transcript_consequences.filter(csq => csq.canonical == 1)' % {'vep': vep_root})


def filter_to_pass(vds):
    """
    Does what it says

    :param VariantDataset vds: Input VDS
    :return: vds with only PASS
    :rtype: VariantDataset
    """
    return vds.annotate_variants_expr(index_into_arrays(['va.info.AS_FilterStatus'])).filter_variants_expr('va.filters.isEmpty && va.info.AS_FilterStatus.isEmpty')


def toSSQL(s):
    """
        Replaces `.` with `___`, since Spark ML doesn't support column names with `.`

    :param str s: The string in which the replacement should be done
    :return: string with `___`
    :rtype: str
    """
    return s.replace('.', '___')


def fromSSQL(s):
    """
        Replaces `___` with `.`, to go back from SSQL to hail annotations

    :param str s: The string in which the replacement should be done
    :return: string with `.`
    :rtype: str
    """
    return s.replace('___', '.')
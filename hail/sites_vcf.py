from utils import *
from hail import *
import time
import argparse


FILTERS_DESC = {
    'InbreedingCoeff': 'InbreedingCoeff < -0.3',
    'LCR': 'In a low complexity region',
    'PASS': 'All filters passed for at least one of the alleles at that site (see AS_FilterStatus for allele-specific filter status)',
    'RF': 'Failed random forests filters (SNV cutoff %s, indels cutoff %s)',
    'SEGDUP': 'In a segmental duplication region',
    'AC0': 'Allele Count is zero (i.e. no high-confidence genotype (GQ >= %(gq)s, DP >= %(dp)s, AB => %(ab)s for het calls))' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}
}

ANNOTATION_DESC = {
    'AC': ('A', 'Allele count in %sgenotypes, for each ALT allele, in the same order as listed'),
    'AF': ('A', 'Allele Frequency among %sgenotypes, for each ALT allele, in the same order as listed'),
    'AN': ('1', 'Total number of alleles in %scalled genotypes'),
    'Hom': ('A', 'Count of homozygous %sindividuals'),
    'Hemi': ('A', 'Count of hemizygous %sindividuals'),
    'GC': ('G', 'Count of %sindividuals for each genotype')
}


def preprocess_exomes_vds(vds, release=True):
    vds = (vds
           .annotate_samples_expr(['sa.meta.project_description = sa.meta.description'])  # Could be cleaner
           .annotate_variants_table(KeyTable.import_bed(decoy_intervals_path), root='va.decoy')
           .annotate_variants_table(KeyTable.import_interval_list(lcr_intervals_path), root='va.lcr')
    )
    return vds.filter_samples_expr('sa.meta.drop_status == "keep"') if release else vds


def preprocess_genomes_vds(vds, release=True):
    vds = (vds
           .annotate_samples_expr(['sa.meta.population = if(sa.meta.final_pop == "sas") "oth" else sa.meta.final_pop',
                                   'sa.meta.project_description = sa.meta.Title'])  # Could be cleaner
           .annotate_variants_table(KeyTable.import_bed(decoy_intervals_path), root='va.decoy')
           .annotate_variants_table(KeyTable.import_interval_list(lcr_intervals_path), root='va.lcr')
           .annotate_variants_expr('va.info = drop(va.info, MQ0, RAW_MQ)')
    )
    return vds.filter_samples_expr('sa.meta.keep') if release else vds


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


def set_va_attributes(vds, warn_if_not_found=True):

    info_va_attr = get_info_va_attr()
    va_info = [x for x in vds.variant_schema.fields if x.name == "info"][0]

    for ann in va_info.typ.fields:
        if ann.name in info_va_attr:
            vds = vds.set_va_attributes("va.info.%s" % ann.name, info_va_attr[ann.name])

        elif ann.name != "CSQ" and warn_if_not_found: logger.warn("No description found for va.info.%s", ann.name)

    return vds


def get_hom_from_gc(destination, target):
    return '%s = range(v.nAltAlleles).map(i => let n = i + 2 in %s[(n * (n + 1) / 2).toInt - 1])' % (destination, target)


def unfurl_callstats(vds, pops, lower=False, gc=True):
    callstats_command, right_shift_command = unfurl_callstats_text(pops, lower, gc)
    return (vds.annotate_variants_expr(callstats_command)
            .annotate_variants_expr(right_shift_command))


def unfurl_hom_text(pops, simple_hom=True, hom_adj=True):
    expressions = [get_hom_from_gc('va.info.Hom_%s' % pop, 'va.info.GC_%s' % pop) for pop in pops]
    if simple_hom: expressions.append('va.info.Hom_raw = range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC_raw[(n * (n + 1) / 2).toInt - 1])')
    if hom_adj: expressions.append('va.info.Hom = range(v.nAltAlleles).map(i => let n = i + 2 in va.info.GC[(n * (n + 1) / 2).toInt - 1])')
    return ',\n'.join(expressions)


def unfurl_hom(vds, pops, simple_hom=True, hom_adj=True):
    hom_command = unfurl_hom_text(pops, simple_hom, hom_adj)
    return vds.annotate_variants_expr(hom_command)


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


def histograms(vds, root='va.info', allele_balance=True, as_text=True, extra_gs_filter=''):
    allele_hists = [
        '{}.GQ_HIST_ALT = gs.filter(g => g.isCalledNonRef {}).map(g => g.gq).hist(0, 100, 20)'.format(root, extra_gs_filter),
        '{}.DP_HIST_ALT = gs.filter(g => g.isCalledNonRef {}).map(g => g.dp).hist(0, 100, 20)'.format(root, extra_gs_filter)]
    variants_hists = [
        '{}.GQ_HIST_ALL = gs.filter(g => g.isCalled {}).map(g => g.gq).hist(0, 100, 20)'.format(root, extra_gs_filter),
        '{}.DP_HIST_ALL = gs.filter(g => g.isCalled {}).map(g => g.dp).hist(0, 100, 20)'.format(root, extra_gs_filter)]

    if allele_balance:
        allele_hists.append(
            '{}.AB_HIST_ALT = gs.filter(g => g.isHet {}).map(g => 100*g.ad[1]/g.dp).hist(0, 100, 20)'.format(root, extra_gs_filter))
        variants_hists.append(
            '{}.AB_HIST_ALL = gs.filter(g => g.isHet {}).map(g => 100 - 100*g.ad[0]/g.dp).hist(0, 100, 20)'.format(root, extra_gs_filter))

    if as_text:
        allele_hists = ['{}.binFrequencies.map(y => str(y)).mkString("|")'.format(x) for x in allele_hists]
        variants_hists = ['{}.binFrequencies.map(y => str(y)).mkString("|")'.format(x) for x in variants_hists]

    return (
        vds.annotate_alleles_expr(allele_hists, propagate_gq=True)
        .annotate_variants_expr(variants_hists)
    )


def write_vcfs(vds, contig, out_internal_vcf_prefix, out_external_vcf_prefix, rf_snv_cutoff, rf_indel_cutoff,
               as_filter_status_fields=('va.info.AS_FilterStatus', ),
               append_to_header=None, drop_fields=None, nchunks=None):

    if contig != '':
        vds = vds.filter_intervals(Interval.parse(str(contig)))
    logger.info('Writing VCFs for chromosome: %s', contig if contig != '' else 'all')

    if drop_fields:
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


def create_sites_vds_annotations(vds, pops, dbsnp_kt=None, drop_star=True, filter_alleles=True, generate_hists=True):

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

    vds = vds.filter_intervals(Interval.parse('1-22'))

    vds = common_sites_vds_annotations(vds)

    if dbsnp_kt is not None:
        vds = vds.annotate_variants_table(dbsnp_kt, expr='va.rsid = table.f2')

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


def create_sites_vds_annotations_X(vds, pops, dbsnp_kt=None, drop_star=True, filter_alleles=True, generate_hists=True):

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

    # Abandon hope, all ye who enter here
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

    vds = vds.filter_intervals(Interval.parse('X'))

    vds = common_sites_vds_annotations(vds)

    if dbsnp_kt is not None:
        vds = vds.annotate_variants_table(dbsnp_kt, expr='va.rsid = table.f2')

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
    vds = popmax(vds, pops)

    return vds.annotate_variants_expr('va.info = drop(va.info, MLEAC, MLEAF)')


def create_sites_vds_annotations_Y(vds, pops, dbsnp_kt=None, drop_star=True, filter_alleles=True, generate_hists=True):

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

    vds = vds.filter_intervals(Interval.parse('Y'))

    vds = common_sites_vds_annotations(vds)

    if dbsnp_kt is not None:
        vds = vds.annotate_variants_table(dbsnp_kt, expr='va.rsid = table.f2')

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
        vds = histograms(vds, 'va.info', allele_balance=False)

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


def run_sites_sanity_checks(vds, pops, verbose=True, contig='auto', percent_missing_threshold=0.01,
                            skip_star=False, split_lcr=False, split_star=False):
    logger.info("Running sites sanity checks on contig %s" % contig)

    output = ''

    # Grouped by filters
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

    vds = vds.persist()

    split_vds = (
        vds
        .annotate_variants_expr(pre_split_ann)
    )
    rf_ann_expr = ['va.info.AS_RF_NEGATIVE_TRAIN = isDefined(va.info.AS_RF_NEGATIVE_TRAIN) && va.info.AS_RF_NEGATIVE_TRAIN.toSet.contains(va.aIndex)',
                   'va.info.AS_RF_POSITIVE_TRAIN = isDefined(va.info.AS_RF_POSITIVE_TRAIN) && va.info.AS_RF_POSITIVE_TRAIN.toSet.contains(va.aIndex)']
    split_vds = split_vds_and_annotations(split_vds, ['InbreedingCoeff'], 'va.info.AS_FilterStatus', rf_ann_expr, vep_root=None)
    split_vds = split_vds.persist()

    df = (
        split_vds
        .variants_table().aggregate_by_key(agg_key,
                                           'n = va.count(), '
                                           'prop_filtered = va.filter(x => !x.filters.isEmpty).count(),'
                                           'prop_hard_filtered = va.filter(x => x.filters.contains("InbreedingCoeff")).count(),'
                                           'prop_AC0_filtered = va.filter(x => x.filters.contains("AC0")).count(),'
                                           'prop_RF_filtered = va.filter(x => x.filters.contains("RF")).count(),'
                                           'prop_hard_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "InbreedingCoeff")).count(),'
                                           'prop_AC0_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "AC0")).count(),'
                                           'prop_RF_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "RF")).count()')
        .annotate(['{0} = {0} / n'.format(ann) for ann in
                   ['prop_filtered', 'prop_hard_filtered', 'prop_AC0_filtered', 'prop_RF_filtered',
                    'prop_hard_filtered_only', 'prop_AC0_filtered_only', 'prop_RF_filtered_only']])
        .to_dataframe()
    )

    df = df.select(select_cols)

    # At some point should print output too
    # pd = df.toPandas()
    # if return_string:
    #     output += "\nProportion of sites filtered by allele type:\n%s" % str(pd)

    print("\nProportion of alleles filtered by allele type:\n")
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

    order_by_cols = ["type", "nAltAlleles"]

    if split_lcr:
        agg_key += ', LCR = va.lcr'
        select_cols.insert(0, 'LCR')
        order_by_cols.append("LCR")

    if split_star:
        agg_key += ', hasStar = va.hasStar'
        select_cols.append('hasStar')
        order_by_cols.append("hasStar")

    df = (
        split_vds
        .variants_table().aggregate_by_key(agg_key,
                                           'n = va.count(), '
                                           'prop_filtered = va.filter(x => !x.filters.isEmpty).count(),'
                                           'prop_hard_filtered = va.filter(x => x.filters.contains("InbreedingCoeff")).count(),'
                                           'prop_AC0_filtered = va.filter(x => x.filters.contains("AC0")).count(),'
                                           'prop_RF_filtered = va.filter(x => x.filters.contains("RF")).count(),'
                                           'prop_hard_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "InbreedingCoeff")).count(),'
                                           'prop_AC0_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "AC0")).count(),'
                                           'prop_RF_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "RF")).count()')
        .annotate(['{0} = {0} / n'.format(ann) for ann in
                   ['prop_filtered', 'prop_hard_filtered', 'prop_AC0_filtered', 'prop_RF_filtered',
                    'prop_hard_filtered_only', 'prop_AC0_filtered_only', 'prop_RF_filtered_only']])
        .to_dataframe()
        .orderBy(order_by_cols)
    )

    df = df.select(select_cols)

    print("\nProportion of alleles filtered by site type and number of alt alleles:\n")
    df.show(n=200)

    df = (
        vds
        .annotate_variants_expr(pre_split_ann)
        .variants_table().aggregate_by_key(agg_key,
                                           'n = va.count(), '
                                           'prop_filtered = va.filter(x => !x.filters.isEmpty).count(),'
                                           'prop_hard_filtered = va.filter(x => x.filters.contains("InbreedingCoeff")).count(),'
                                           'prop_AC0_filtered = va.filter(x => x.filters.contains("AC0")).count(),'
                                           'prop_RF_filtered = va.filter(x => x.filters.contains("RF")).count(),'
                                           'prop_hard_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => ["InbreedingCoeff"].toSet.contains(f)) ).count(),'
                                           'prop_AC0_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "AC0") ).count(),'
                                           'prop_RF_filtered_only = va.filter(x => !x.filters.isEmpty && x.filters.forall(f => f == "RF") ).count()'
        )
        .annotate(['{0} = {0} / n'.format(ann) for ann in
                   ['prop_filtered', 'prop_hard_filtered', 'prop_AC0_filtered', 'prop_RF_filtered',
                    'prop_hard_filtered_only', 'prop_AC0_filtered_only', 'prop_RF_filtered_only']])
        .to_dataframe()
        .orderBy(order_by_cols)
    )

    df = df.select(select_cols)

    print("\nProportion of sites filtered by variant type and number of alt alleles:\n")
    df.show(n=200)

    #At some point should print output too
    # pd = df.toPandas()
    # if return_string:
    #     output += "\nProportion of sites filtered by variant type and number of alt alleles:\n%s" % str(pd)
    # print("\nProportion of sites filtered by variant type and number of alt alleles:\n%s" % str(pd))

    queries = ['variants.count()']
    a_metrics = ['AC', 'Hom']

    if contig == 'Y':
        a_metrics = a_metrics[:1]

    one_metrics = ['AN'] if skip_star else ['STAR_AC', 'AN']

    # Filter counts
    queries.extend(["let x = variants.map(v => !va.filters.isEmpty).counter() in orElse(x.get(true), 0L)/x.values().sum()",
                    'variants.map(v => va.filters.toArray.mkString(",")).counter()'])

    # Check number of samples
    queries.append('variants.map(v => va.info.AN).stats().max/2')
    queries.extend(['variants.map(v => va.info.AN_{}).stats().max/2'.format(pop) for pop in pops])

    end_pop_counts = len(queries)

    # Check that raw is always larger than adj
    queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                    '.exists(i => va.info.{0}[i] > va.info.{0}_raw[i])).count()'.format(x) for x in a_metrics])
    queries.extend(['variants.filter(v => va.info.{0} > va.info.{0}_raw).count()'.format(x) for x in one_metrics])

    # Check that sum(pops) == total
    for metric in a_metrics:
        queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                        '.exists(i => {} != va.info.{}[i])).count()'.format(" + ".join(["va.info.{}_{}[i]".format(metric, pop) for pop in pops]), metric)])

    queries.append('variants.filter(v => {} != va.info.AN).count()'.format(" + ".join(["va.info.AN_{}".format(pop) for pop in pops])))

    # Check that male + female == total
    # Remove Hom for X
    if contig == 'X':
        a_metrics = a_metrics[:1]

    if contig != 'Y':
        pop_strats = pops if contig == 'X' else [None]
        for pop in pop_strats:
            pop_text = "" if pop is None else "_" + pop
            queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                            '.exists(i => va.info.{x}{pop}_Male[i] + va.info.{x}{pop}_Female[i] != va.info.{x}{pop}[i])).count()'.format(x=x, pop=pop_text) for x in a_metrics])
            queries.extend(['variants.filter(v => va.info.{x}{pop}_Male + va.info.{x}{pop}_Female != va.info.{x}{pop}).count()'.format(x=x, pop=pop_text) for x in one_metrics[1:]])
    end_counts = len(queries)

    missing_metrics = []
    va_info = [x for x in vds.variant_schema.fields if x.name == "info"][0].typ
    for field in va_info.fields:
        missing_metrics.append('va.info.{}'.format(field.name))
        queries.append('let x = variants.map(v => isMissing(va.info.{})).counter() in orElse(x.get(true), 0L)/x.values().sum()'.format(field.name))

    logger.debug('Queries: %s', '\n'.join(['{}: {}'.format(i, x) for i, x in enumerate(queries)]))
    stats = vds.query_variants(queries)

    # Print filters

    # Double checking for no samples in VDS
    sample_count = vds.query_samples('samples.count()')
    if sample_count > 0:
        output += 'WARNING: {} samples found in VDS (should be 0)\n'.format(sample_count)

    output += "Total number of sites:\n{}\n".format(pformat(int(stats[0])))
    output += "FILTERS CHECKS\nTotal fraction sites filtered:\n{}\n".format(pformat(stats[1]))
    output += "Filter counts:{}\n".format(pformat(stats[2]))

    output += "\nPOPULATION COUNTS\nTotal number of samples: {}\n".format(stats[3])
    for i in range(4, end_pop_counts):
        output += '{}: {}\n'.format(pops[i - 4], stats[i])

    # Check that all metrics sum as expected
    output += "\nMETRICS COUNTS CHECK\n"
    nfail = 0
    for i in range(end_pop_counts, end_counts):
        if stats[i] != 0:
            output += "FAILED METRICS CHECK for query: {}\n Expected: 0, Found: {}\n".format(queries[i], stats[i])
            nfail += 1
        elif verbose:
            output += "Success: {}\n".format(queries[i])
    output += "{} metrics count checks failed.\n".format(nfail)

    # Check missing metrics
    output += "\nMISSING METRICS CHECKS\n"
    nfail = 0
    missing_stats = stats[end_counts:]
    for metric, stat in zip(missing_metrics, missing_stats):
        if stat > percent_missing_threshold:
            output += "FAILED missing check for {}; {}% missing.\n".format(metric, 100 * stat)
            nfail += 1
        elif verbose:
            output += "SUCCESS missing check for {}; {}% missing.\n".format(metric, 100 * stat)
    output += "{} missing metrics checks failed.\n".format(nfail)

    logger.info(output)
    return output


def post_process_vds(vds, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, rf_train='va.train', rf_label='va.label'):

    logger.info("Postprocessing...")
    rf_annotations = {
        'va.stats.qc_samples_raw.nrq_median': 'va.info.DREF_MEDIAN',
        'va.stats.qc_samples_raw.gq_median': 'va.info.GQ_MEDIAN',
        'va.stats.qc_samples_raw.dp_median': 'va.info.DP_MEDIAN',
        'va.stats.qc_samples_raw.ab_median': 'va.info.AB_MEDIAN'
    }

    vds = annotate_from_rf(vds, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, annotations=rf_annotations, train=rf_train, label=rf_label)
    vds = set_filters(vds, rf_snv_cutoff, rf_indel_cutoff)

    return set_va_attributes(vds)


def annotate_from_rf(vds, rf_vds, rf_snv_cutoff, rf_indel_cutoff, rf_root, annotations={}, train='va.train', label='va.label'):

    # Strip va if present
    rf_root = rf_root[len('va.'):] if rf_root.startswith('va.') else rf_root
    train = train[len('va.'):] if train.startswith('va.') else train
    label = label[len('va.'):] if label.startswith('va.') else label

    rf_ann_expr = [
        'va.info.AS_RF = orMissing(vds.exists(x => isDefined(x) && isDefined(x.{0})), '
        'range(v.nAltAlleles).map(i => orMissing(isDefined(vds[i]), vds[i].{0}.probability.get("TP"))))'.format(rf_root),
        'va.info.AS_FilterStatus = '
        '   if(vds.forall(x => isMissing(x) || isMissing(x.%(root)s ))) range(v.nAltAlleles).map(i => ["RF"].toSet)'
        '   else range(v.nAltAlleles).map(i => '
        '       if(isMissing(vds[i]) || isMissing(vds[i].%(root)s)) ["RF"].toSet'
        '       else'
        '           if(v.altAlleles[i].isSNP)'
        '               if(vds[i].%(root)s.probability["TP"] > %(snv).4f) [""][:0].toSet else ["RF"].toSet'
        '           else'
        '               if(vds[i].%(root)s.probability["TP"] > %(indel).4f) [""][:0].toSet else ["RF"].toSet'
        '       )' % {'root': rf_root,
                      'snv': rf_snv_cutoff,
                      'indel': rf_indel_cutoff},
        'va.info.AS_RF_POSITIVE_TRAIN = let x = range(v.nAltAlleles).filter('
        'i => isDefined(vds[i]) && isDefined(vds[i].{train}) && isDefined(vds[i].{label}) && vds[i].{train} && vds[i].{label} == "TP")'
        '.map(i => i+1) in orMissing(!x.isEmpty, x)'.format(train=train, label=label),
        'va.info.AS_RF_NEGATIVE_TRAIN = let x = range(v.nAltAlleles).filter('
        'i => isDefined(vds[i]) && isDefined(vds[i].{train}) && isDefined(vds[i].{label}) && vds[i].{train} && vds[i].{label} == "FP")'
        '.map(i => i+1) in orMissing(!x.isEmpty, x)'.format(train=train, label=label)
    ]

    for source, target in annotations.iteritems():
        # Strip va if present
        source = source[len('va.'):] if source.startswith('va.') else source
        rf_ann_expr.append('{target} = orMissing(vds.exists(x => isDefined(x) && isDefined(x.{source})),'
                           ' range(v.nAltAlleles).map(i => orMissing(isDefined(vds[i]), '
                           ' vds[i].{source})))'.format(target=target, source=source))

    return vds.annotate_alleles_vds(rf_vds, rf_ann_expr)


def set_filters_attributes(vds, rf_snv_cutoff, rf_indel_cutoff):
    filters_desc = copy.deepcopy(FILTERS_DESC)
    if rf_snv_cutoff is not None and rf_indel_cutoff is not None and "RF" in filters_desc:
        filters_desc["RF"] = filters_desc["RF"] % (rf_snv_cutoff, rf_indel_cutoff)

    return vds.set_va_attributes('va.filters', filters_desc)


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
    return set_site_filters(vds, site_filters, as_filters_root='va.info.AS_FilterStatus')


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


def write_public_vds(vds, public_path, overwrite=False):
    vds = vds.annotate_samples_expr('sa = {}')
    vds = vds.annotate_variants_expr('va = select(va, rsid, qual, filters, pass, info)')
    vds = vds.annotate_variants_expr('va.info = drop(va.info, PROJECTMAX, PROJECTMAX_NSamples, PROJECTMAX_NonRefSamples, PROJECTMAX_PropNonRefSamples)')
    vds.write(public_path, overwrite=overwrite)


def main(args):
    if args.debug: logger.setLevel(logging.DEBUG)
    hc = HailContext(log='/hail.sites_vcf.log')

    data_type = 'genomes' if args.genomes else 'exomes'
    pops = GENOME_POPS if args.genomes else EXOME_POPS
    RF_SNV_CUTOFF = None if args.genomes else 0.1
    RF_INDEL_CUTOFF = None if args.genomes else 0.1
    preprocess_vds = preprocess_genomes_vds if args.genomes else preprocess_exomes_vds
    rf_path = 'gs://gnomad-exomes/variantqc/170620_new/gnomad_exomes.rf.vds' if args.exomes else ''
    running = 'exomes' if args.exomes else 'genomes'

    if not (args.skip_preprocess_autosomes and args.skip_preprocess_X and args.skip_preprocess_Y):
        vds = get_gnomad_data(hc, data_type)
        vds = preprocess_vds(vds, True).annotate_samples_expr(
            'sa.meta = select(sa.meta, sex, population, project_description)'
        )
        if args.expr:
            vds = vds.filter_samples_expr(args.expr)
        logger.info('Found %s samples', vds.query_samples('samples.count()'))

        dbsnp_kt = (hc
                    .import_table(dbsnp_vcf_path, comment='#', no_header=True, types={'f0': TString(), 'f1': TInt()})
                    .annotate('locus = Locus(f0, f1)')
                    .key_by('locus')
        )
        if not args.skip_preprocess_autosomes:
            (
                create_sites_vds_annotations(vds, pops, dbsnp_kt=dbsnp_kt)
                .write(args.output + ".pre.autosomes.vds", overwrite=args.overwrite)
            )

        if not args.skip_preprocess_X:
            (
                create_sites_vds_annotations_X(vds, pops, dbsnp_kt=dbsnp_kt)
                .write(args.output + ".pre.X.vds", overwrite=args.overwrite)
            )

        if args.exomes and not args.skip_preprocess_Y:
            (
                create_sites_vds_annotations_Y(vds, pops, dbsnp_kt=dbsnp_kt)
                .write(args.output + ".pre.Y.vds", overwrite=args.overwrite)
            )

    if not args.skip_merge:
        vdses = [hc.read(args.output + ".pre.autosomes.vds"), hc.read(args.output + ".pre.X.vds")]
        if args.exomes: vdses.append(hc.read(args.output + ".pre.Y.vds"))
        vdses = merge_schemas(vdses)
        vds = vdses[0].union(vdses[1:])
        vds.write(args.output + '.pre.vds', overwrite=args.overwrite)

    if not args.skip_vep:
        (hc.read(args.output + ".pre.vds")
         .vep(config=vep_config, csq=True, root='va.info.CSQ')
         .write(args.output + ".pre.vep.vds", overwrite=args.overwrite)
         )

    if not args.skip_postprocess:
        vds = hc.read(args.output + ".pre.vep.vds")
        rf_vds = hc.read(rf_path)
        post_process_vds(vds, rf_vds, RF_SNV_CUTOFF, RF_INDEL_CUTOFF,
                         'va.rf').write(args.output + ".post.vds", overwrite=args.overwrite)

        vds = hc.read(args.output + ".post.vds")
        sanity_check = run_sites_sanity_checks(vds, pops)
        if args.slack_channel: send_snippet(args.slack_channel, sanity_check, 'sanity_%s.txt' % time.strftime("%Y-%m-%d_%H:%M"))

    if not args.skip_write:
        vds = hc.read(args.output + ".post.vds")
        if args.exomes:
            exome_intervals = KeyTable.import_interval_list(exome_calling_intervals_path)
            vds = vds.filter_variants_table(exome_intervals)
        write_vcfs(vds, '', args.output, False, RF_SNV_CUTOFF, RF_INDEL_CUTOFF, append_to_header=additional_vcf_header_path)
        write_public_vds(vds, args.output + ".vds", overwrite=args.overwrite)

    if not args.skip_pre_calculate_metrics:
        vds = hc.read(args.output + ".vds")
        fname = '{}_precalculated_metrics.txt'.format(running)
        pre_calculate_metrics(vds, fname)
        send_snippet('#gnomad_browser', open(fname).read())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--skip_preprocess_autosomes', help='Skip pre-processing autosomes (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_X', help='Skip pre-processing X (assuming already done)', action='store_true')
    parser.add_argument('--skip_preprocess_Y', help='Skip pre-processing Y (assuming already done)', action='store_true')
    parser.add_argument('--skip_postprocess', help='Skip merge and post-process (assuming already done)', action='store_true')
    parser.add_argument('--skip_merge', help='Skip merging data (assuming already done)', action='store_true')
    parser.add_argument('--skip_vep', help='Skip VEPping data (assuming already done)', action='store_true')
    parser.add_argument('--skip_write', help='Skip writing data (assuming already done)', action='store_true')
    parser.add_argument('--skip_pre_calculate_metrics', help='Skip pre-calculating metrics (assuming already done)', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    parser.add_argument('--expr', help='''Additional expression (e.g. "!sa.meta.remove_for_non_tcga)"''')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    main(args)
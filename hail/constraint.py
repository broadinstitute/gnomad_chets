#!/usr/bin/env python

import argparse
from utils import *
import statsmodels.formula.api as smf
import pandas as pd

try:
    hc = HailContext(log="/hail.log")
except Exception:
    pass

# Unprocessed files
# final_exome_vds = 'gs://gnomad-public/release-170228/gnomad.exomes.r2.0.1.sites.vds'
final_exome_vds = 'gs://gnomad-exomes/sites/170622_new/gnomad.exomes.sites.vds'
final_genome_vds = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.vds'

fasta_path = "gs://gnomad-resources/Homo_sapiens_assembly19.fasta"
gerp_annotations_path = 'gs://annotationdb/cadd/cadd.kt'  # Gerp is in here as gerpS
raw_context_vds_path = "gs://gnomad-resources/Homo_sapiens_assembly19.fasta.snps_only.vep.vds"
genome_coverage_kt_path = 'gs://gnomad-resources/genome_coverage.kt'
exome_coverage_kt_path = 'gs://gnomad-resources/exome_coverage.kt'

# Processed datasets
context_vds_path = 'gs://gnomad-resources/constraint/context_processed.vds'
genome_vds_path = 'gs://gnomad-resources/constraint/genome_processed.vds'
exome_vds_path = 'gs://gnomad-resources/constraint/exome_processed.vds'

CONTIG_GROUPS = ('1', '2', '3', '4', '5', '6', '7', '8-9', '10-11', '12-13', '14-16', '17-18', '19-20', '21', '22', 'X', 'Y')
# should have been: ('1', '2', '3', '4', '5', '6', '7', '8-9', '10-11', '12-13', '14-16', '17-19', '20-22', 'X', 'Y')
a_based_annotations = ['va.info.AC', 'va.info.AC_raw']

HIGH_COVERAGE_CUTOFF = 50
AF_CRITERIA = 'va.info.AN > 0 && va.info.AC/va.info.AN < 0.001'
GENOME_COVERAGE_CRITERIA = 'va.coverage.genome.mean >= 15 && va.coverage.genome.mean <= 60'


def remove_ttn(kt):
    """
    Remove TTN transcripts

    :param kt: input KeyTable
    :return: keytable without TTN transcripts
    :rtype KeyTable
    """
    ttn_transcripts = ["ENST00000342992", "ENST00000460472", "ENST00000589042", "ENST00000591111", "ENST00000342175", "ENST00000359218", "ENST00000414766", "ENST00000426232", "ENST00000446966", "ENST00000425332", "ENST00000448510", "ENST00000360870", "ENST00000436599", "ENST00000470257", "ENST00000412264"]
    # zcat gencode.v19.annotation.gtf.gz | awk '$3 == "transcript" {print}' | grep 'gene_name "TTN"' | awk '{print $12}' | perl -p -e 's/\..//g' | perl -p -e 's/;\n/, /g'
    return kt.filter('!({}.exists(x => x == transcript))'.format(str(ttn_transcripts).replace("'", '"')))


def import_fasta_and_vep(input_fasta_path, output_vds_path, overwrite=False):
    """
    Imports FASTA file with context and VEPs it. Only works with SNPs so far.
    WARNING: Very slow and annoying. Writes intermediate paths for now since VEPping the whole genome is hard.
    Some paths are hard-coded in here. Remove if needed down the line.

    :param str input_fasta_path: Input FASTA file
    :param str output_vds_path: Path to final output VDS
    :param bool overwrite: Whether to overwrite VDSes
    :return: None
    """
    input_vds_path = "gs://gnomad-resources/Homo_sapiens_assembly19.fasta.snps_only.vds"

    vds = hc.import_fasta(input_fasta_path,
                          filter_Ns=True, flanking_context=3, create_snv_alleles=True, create_deletion_size=0,
                          create_insertion_size=0, line_limit=2000)
    vds.split_multi().write(input_vds_path, overwrite=overwrite)

    for i in CONTIG_GROUPS:
        vds_name = "gs://gnomad-resources/chr_split/Homo_sapiens_assembly19.fasta.snps_only.split.%s.vds" % i
        print "Reading, splitting, and repartitioning contigs: %s..." % i
        hc.read(input_vds_path).filter_intervals(Interval.parse(i)).repartition(10000).write(vds_name,
                                                                                             overwrite=overwrite)
        print "Done! VEPping contigs: %s..." % i
        vds = hc.read(vds_name)
        vds = vds.vep(vep_config)
        vds.write(vds_name.replace('.vds', '.vep.vds'), overwrite=overwrite)

    print "Done! Unioning..."
    vds = hc.read(
        ["gs://gnomad-resources/chr_split/Homo_sapiens_assembly19.fasta.snps_only.split.%s.vep.vds" % i for i in
         CONTIG_GROUPS])
    vds.repartition(40000, shuffle=False).write(output_vds_path, overwrite=overwrite)


def pre_process_all_data():
    exome_coverage_kt = hc.read_table(exome_coverage_kt_path).select(['locus', 'mean', 'median'])
    genome_coverage_kt = hc.read_table(genome_coverage_kt_path).select(['locus', 'mean', 'median'])
    methylation_kt = hc.read_table(methylation_kt_path).select(['locus', 'MEAN'])

    gerp_kt = hc.read_table(gerp_annotations_path).select(['variant', 'GerpS'])
    context_vds = process_consequences(hc.read(raw_context_vds_path)
                                       .split_multi()
                                       .annotate_variants_expr(index_into_arrays(vep_root='va.vep'))
                                       .annotate_variants_table(gerp_kt, root='va.gerp')
                                       .annotate_variants_table(exome_coverage_kt, root='va.coverage.exome')
                                       .annotate_variants_table(genome_coverage_kt, root='va.coverage.genome')
                                       .annotate_variants_table(methylation_kt, root='va.methylation.value')
                                       .annotate_variants_expr('va.methylation.level = '
                                                               'if (v.altAllele().isTransition()) '
                                                               '    range(19, -1, -1).find(e => va.methylation.value > 20*e) '
                                                               'else NA: Int')
    )
    context_vds.write(context_vds_path, overwrite=args.overwrite)
    context_vds = hc.read(context_vds_path)
    genome_vds = (hc.read(final_genome_vds).split_multi()
                  .annotate_variants_expr(index_into_arrays(a_based_annotations))
                  .annotate_variants_vds(context_vds,
                                         expr='va.context = vds.context, va.gerp = vds.gerp, va.vep = vds.vep, '
                                              'va.coverage = vds.coverage, va.methylation = vds.methylation'))
    genome_vds.write(genome_vds_path, overwrite=args.overwrite)
    exome_vds = (hc.read(final_exome_vds).split_multi()
                 .annotate_variants_expr(index_into_arrays(a_based_annotations))
                 .annotate_variants_vds(context_vds,
                                        expr='va.context = vds.context, va.gerp = vds.gerp, va.vep = vds.vep, '
                                             'va.coverage = vds.coverage, va.methylation = vds.methylation'))
    exome_vds.write(exome_vds_path, overwrite=args.overwrite)


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases[::-1])


def rev_comp(ref, alt, context):
    if ref in ('G', 'T'):
        ref, alt, context = reverse_complement(ref), reverse_complement(alt), reverse_complement(context)
    return pd.Series({'new_ref': ref, 'new_alt': alt, 'new_context': context})


def variant_type(ref, alt, context):  # new_ref is only A and C
    if ref == 'C' and alt == 'T':
        return 'CpG' if context[2] == 'G' else 'transition'
    elif ref == 'A' and alt == 'G':
        return 'transition'
    elif ref == 'C':
        return 'CpG transversion' if context[2] == 'G' else 'transversion'
    return 'transversion'


def count_variants(vds, criteria=None, additional_groupings=None, trimer=False, explode=None, collapse_contexts=True, coverage=True, methylation=False):
    """
    Counts variants in VDS by context, ref, alt, and any other groupings provided

    :param VariantDataset vds: Input VDS
    :param str criteria: Any filtering criteria (e.g. non-coding, non-conserved), to be passed to filter_variants_expr
    :param additional_groupings: Whether to group further (e.g. by functional annotation)
    :type additional_groupings: str or list of str
    :param bool trimer: whether to use trimer context (default heptamer)
    :param str explode: criteria to explode by (most likely va.vep.transcript_consequences)
    :param bool or dict coverage: Whether to also aggregate coverage information (allowed to be dict for convenience)
    :return: keytable with counts as `variant_count`
    :rtype: KeyTable
    """
    if criteria is not None:
        vds = vds.filter_variants_expr(criteria)

    if trimer:
        vds = vds.annotate_variants_expr('va.context = va.context[2:5]')

    grouping = ['context = `va.context`', 'ref = `va.ref`',
                'alt = `va.alt`']  # va is flattened, so this is a little awkward, and now ref and alt are in va.
    if methylation: grouping.append('methylation_level = `va.methylation.level`')
    if additional_groupings is not None:
        if isinstance(additional_groupings, str):
            grouping.append(additional_groupings)
        else:
            grouping.extend(additional_groupings)

    kt = vds.variants_table().flatten()

    if collapse_contexts:
        kt = collapse_strand(kt)

    if explode:
        kt = kt.explode(explode).flatten()

    aggregation_functions = ['variant_count = v.count()']
    if coverage: aggregation_functions.append('sum_coverage = `va.coverage.exome.median`.sum()')

    return kt.aggregate_by_key(grouping, aggregation_functions)


def calculate_mutation_rate(possible_variants_vds, genome_vds, criteria=None, trimer=False, methylation=False):
    """
    Calculate mutation rate from all possible variants vds and observed variants vds
    Currently actually calculating more like "expected_proportion_variants"

    :param VariantDataset possible_variants_vds: synthetic VDS
    :param VariantDataset genome_vds: gnomAD WGS VDS
    :param bool trimer: whether to use trimer context (default heptamer)
    :return: keytable with mutation rates as `mutation_rate`
    :rtype: KeyTable
    """

    grouping = 'methylation_level = `va.methylation.level`' if methylation else None
    all_possible_kt = count_variants(possible_variants_vds, criteria=criteria, trimer=trimer, additional_groupings=grouping, coverage=False)
    observed_kt = count_variants(genome_vds, criteria=criteria, trimer=trimer, additional_groupings=grouping, coverage=False)

    kt = (all_possible_kt.rename({'variant_count': 'possible_variants'})
          .join(observed_kt, how='outer')
          .annotate('mutation_rate = variant_count/possible_variants'))

    return kt.filter('ref.length == 1 && alt.length == 1 && !("N" ~ context)')


def collapse_counts_by_exon(kt, methylation=False,
                            mutation_rate_weights=None, regression_weights=None, coverage_weights=None):
    """
    From context, ref, alt, transcript groupings, group by transcript and returns counts.
    Can optionally weight by mutation_rate in order to generate "expected counts"

    Note: coverage correction should only apply to expected counts (all possible VDS)

    :param KeyTable kt: key table to aggregate over
    :param KeyTable mutation_rate_weights: Mutation rate keytable (if provided, creates expected counts)
    :param dict regression_weights: dict of {'slope': float, 'intercept': float} to adjust aggregate_mutation_rate by
    :param dict coverage_weights: dict of coverage model weights
        e.g. (ExAC example): {'high_cutoff': 50, 'low_cutoff': 1, 'mid_beta': 0.217, 'mid_intercept': 0.089}
    :return: key table grouped by only transcript and exon
    :rtype: KeyTable
    """
    aggregation_expression = ['variant_count = variant_count.sum()']
    if mutation_rate_weights:
        keys = ['context', 'ref', 'alt']
        if methylation: keys.append('methylation_level')
        kt = (kt.key_by(keys)
              .join(mutation_rate_weights.select(keys + ['mutation_rate']), how='outer')
              .annotate(['aggregate_mutation_rate = mutation_rate * variant_count']))
        aggregation_expression.extend(['aggregate_mutation_rate = aggregate_mutation_rate.sum()',
                                       'sum_coverage = sum_coverage.sum()'])

    kt = kt.aggregate_by_key(['transcript = transcript', 'exon = exon', 'annotation = annotation'], aggregation_expression)

    if regression_weights is not None:
        kt = kt.annotate(['median_coverage = sum_coverage//variant_count',
                          'expected_variant_count = {slope} * aggregate_mutation_rate + {intercept}'.format(**regression_weights)])
        if coverage_weights is not None:
            kt = kt.annotate('expected_variant_count_adj = '
                             'if (median_coverage >= {high_cutoff}) '
                             '  expected_variant_count '
                             'else '
                             '  (expected_variant_count * (log(median_coverage) * {mid_beta} + {mid_intercept}))'.format(**coverage_weights))
    return kt


def get_coverage_weights(synonymous_kt_depth, high_cutoff=HIGH_COVERAGE_CUTOFF, low_cutoff=1):
    """
    Gets model for coverage from Keytable with observed and expected counts by coverage

    :param KeyTable synonymous_kt_depth: Keytable with median_coverage, observed, expected
    :return: dict with {'high_cutoff': 35, 'low_cutoff': 1, 'mid_beta': 0.217, 'mid_intercept': 0.089}
    :rtype dict
    """
    output_dict = {
        'high_cutoff': high_cutoff,
        'low_cutoff': low_cutoff
    }
    synonymous_kt_depth = (synonymous_kt_depth
                           .filter('median_coverage < {high_cutoff} && median_coverage > {low_cutoff}'.format(**output_dict))
                           .annotate(['oe = observed/expected', 'log_coverage = log(median_coverage)'])
    )
    synonymous_depth = synonymous_kt_depth.to_pandas()
    lm = smf.ols(formula='oe ~ log_coverage', data=synonymous_depth).fit()
    output_dict['mid_beta'] = lm.params['log_coverage']
    output_dict['mid_intercept'] = lm.params['Intercept']
    return output_dict


def build_synonymous_model(syn_kt):
    """
    Calibrates mutational model to synonymous variants in gnomAD exomes

    :param KeyTable syn_kt: Synonymous keytable
    :return:
    """
    syn_pd = syn_kt.to_pandas()
    lm = smf.ols(formula='variant_count ~ aggregate_mutation_rate', data=syn_pd).fit()
    return {
        'slope': lm.params['aggregate_mutation_rate'],
        'intercept': lm.params['Intercept']
    }


def build_synonymous_model_maps(syn_kt):
    """
    Calibrates mutational model to synonymous variants in gnomAD exomes

    :param KeyTable syn_kt: Synonymous keytable
    :return:
    """
    syn_pd = syn_kt.to_pandas()
    lm = smf.ols(formula='variant_count ~ expected_variant_count', data=syn_pd).fit()
    slope = lm.params['expected_variant_count']
    intercept = lm.params['Intercept']
    return slope, intercept


def get_observed_expected_kt(vds, all_possible_vds, mutation_kt, canonical=False,
                             criteria=None, coverage_cutoff=0, methylation=False,
                             coverage_weights=None, regression_weights=None):
    """
    Get a set of observed and expected counts based on some criteria (and optionally for only canonical transcripts)

    :param VariantDataset vds: VDS to generate observed counts
    :param VariantDataset all_possible_vds: VDS with all possible variants
    :param KeyTable mutation_kt: Mutation rate keytable
    :param bool canonical: Whether to only use canonical transcripts
    :param str criteria: Subset of data to use (e.g. only synonymous, to create model)
    :param int coverage_cutoff: Median coverage cutoff to apply (for creating first-pass model)
    :param bool methylation: Whether to add methylation status for mutational model
    :param dict coverage_weights: dict of coverage model weights to pass through to count_variants
        e.g. (ExAC example): {'high_cutoff': 50, 'low_cutoff': 1, 'mid_beta': 0.217, 'mid_intercept': 0.089}
    :param dict regression_weights: dict of {'slope': float, 'intercept': float} to adjust aggregate_mutation_rate by
    :return: KeyTable with observed (`variant_count`) and expected (`expected_variant_count`)
    :rtype KeyTable
    """
    if canonical:
        vds = filter_vep_to_canonical_transcripts(vds)
        all_possible_vds = filter_vep_to_canonical_transcripts(all_possible_vds)
    vds = process_consequences(vds)
    all_possible_vds = process_consequences(all_possible_vds)

    count_grouping = [
        'annotation = `va.vep.transcript_consequences.most_severe_consequence`',
        'transcript = `va.vep.transcript_consequences.transcript_id`',
        'exon = `va.vep.transcript_consequences.exon`']  # va gets flattened, so this a little awkward
    kt = count_variants(vds, additional_groupings=count_grouping,
                        explode='va.vep.transcript_consequences', trimer=True, methylation=methylation)
    all_possible_kt = count_variants(all_possible_vds,
                                     additional_groupings=count_grouping,
                                     explode='va.vep.transcript_consequences',
                                     trimer=True, methylation=methylation)
    if criteria:
        kt = kt.filter(criteria)
        all_possible_kt = all_possible_kt.filter(criteria)

    collapsed_kt = collapse_counts_by_exon(kt, methylation=methylation)
    collapsed_all_possible_kt = collapse_counts_by_exon(all_possible_kt,
                                                        methylation=methylation,
                                                        mutation_rate_weights=mutation_kt,
                                                        regression_weights=regression_weights,
                                                        coverage_weights=coverage_weights)
    # Calculating median_coverage only for "all possible" keytable means we get the exact mean for each exon
    collapsed_all_possible_kt = (collapsed_all_possible_kt
                                 .annotate('median_coverage = sum_coverage//variant_count')
                                 .rename({'variant_count': 'possible_variant_count'}))
    if coverage_weights is not None:
        coverage_cutoff = coverage_weights['low_cutoff']
    if coverage_cutoff:
        collapsed_all_possible_kt = collapsed_all_possible_kt.filter('median_coverage >= {}'.format(coverage_cutoff))
    return collapsed_kt.join(collapsed_all_possible_kt)


def get_proportion_observed(exome_vds, all_possible_vds, trimer=False, methylation=False):
    """
    Intermediate function to get proportion observed by context, ref, alt

    :param VariantDataset exome_vds: gnomAD exome VDS
    :param VariantDataset all_possible_vds: VDS with all possible variants
    :return: Key Table with context, ref, alt, proportion observed
    :rtype KeyTable
    """
    grouping = ['annotation = `va.vep.transcript_consequences.most_severe_consequence`']  # va gets flattened, so this a little awkward
    if methylation: grouping.append('methylation_level = `va.methylation.level`')
    exome_kt = count_variants(exome_vds,
                              additional_groupings=grouping,
                              explode='va.vep.transcript_consequences', trimer=trimer)
    all_possible_kt = count_variants(all_possible_vds,
                                     additional_groupings=grouping,
                                     explode='va.vep.transcript_consequences', trimer=trimer)

    full_kt = all_possible_kt.rename({'variant_count': 'possible_variants'}).join(exome_kt, how='outer')
    return full_kt.annotate('proportion_observed = variant_count/possible_variants')


def run_sanity_checks(vds, exome=True, csq_queries=False, return_data=False):
    """

    :param VariantDataset vds: Input VDS
    :param bool exome: Run and return exome queries
    :param bool csq_queries: Run and return consequence queries
    :return: whether VDS was split, queries, full sanity results, [exome sanity results], [csq results]
    """
    sanity_queries = ['variants.count()',
                      'variants.filter(v => isMissing(va.vep)).count()',
                      'variants.filter(v => isMissing(va.vep.transcript_consequences)).count()',
                      'variants.filter(v => isMissing(va.gerp)).count()',
                      'variants.filter(v => "CG" ~ va.context[2:5] && isMissing(va.methylation)).count()']

    additional_queries = ['variants.map(v => va.vep.worst_csq).counter()',
                          'variants.map(v => va.vep.worst_csq_suffix).counter()']

    # [x.replace('variants.', 'variants.filter(v => ).') for x in sanity_queries]
    # should probably annotate_variants_intervals and filter directly in query

    full_sanity_results = vds.query_variants(sanity_queries)

    if return_data:
        results = [vds.was_split(), sanity_queries, full_sanity_results]
    else:
        print 'Was split: %s' % vds.was_split()
        print "All data:\n%s" % zip(sanity_queries, full_sanity_results)

    if exome:
        exome_intervals = KeyTable.import_interval_list(exome_calling_intervals)
        exome_intervals_sanity_results = (vds.filter_variants_table(exome_intervals)
                                          .query_variants(sanity_queries))
        if return_data:
            results.append(exome_intervals_sanity_results)
        else:
            print "Exome Intervals:\n%s" % zip(sanity_queries, exome_intervals_sanity_results)

    if csq_queries:
        csq_query_results = vds.query_variants(additional_queries)
        if return_data:
            results.append(csq_query_results)
        else:
            print csq_query_results

    if return_data:
        return results


def collapse_strand(kt):
    """
    :param KeyTable kt: Keytable with context, ref, alt to collapse strands
    :return: Keytable with collapsed strands (puts ref and alt into va.ref and va.alt)
    :rtype KeyTable
    """

    def flip_text(root):
        return ('let base = %s in '
                'if (base == "A") "T" else '
                'if (base == "T") "A" else '
                'if (base == "C") "G" else '
                'if (base == "G") "C" else base' % root
        )

    return kt.annotate(['`va.ref` = if (v.ref == "T" || v.ref == "G") %s else v.ref' % flip_text('v.ref'),
                        '`va.alt` = if (v.ref == "T" || v.ref == "G") %s else v.alt' % flip_text('v.alt'),
                        '`va.context` = if (v.ref == "T" || v.ref == "G") (%s) + (%s) + (%s) '
                        'else `va.context`' % (
                            flip_text('`va.context`[2]'), flip_text('`va.context`[1]'), flip_text('`va.context`[0]'))])


def maps(vds, mutation_kt, additional_groupings=None, trimer=True, methylation=False):
    """

    :param VariantDataset vds: VDS
    :return:
    """

    if trimer: vds = vds.annotate_variants_expr('va.context = va.context[2:5]')

    mutation_kt = (mutation_kt.annotate('meth_int = (methylation_level*4 + 0.05).toInt')
                   .key_by(['context', 'ref', 'alt', 'meth_int']))

    aggregation = ['context = `va.context`', 'ref = `va.ref`', 'alt = `va.alt`']
    if methylation:
        aggregation.append('methylation_level = `va.methylation.level`')

    kt = (vds.repartition(1000, shuffle=False)
          .annotate_variants_expr('va.singleton = (va.info.AC == 1).toLong')
          .annotate_variants_expr('va = select(va, context, singleton, vep, methylation)')
          .annotate_variants_expr('va.vep = select(va.vep, transcript_consequences, worst_csq, worst_csq_suffix)')
          .variants_table()
          .flatten())
    kt = collapse_strand(kt).key_by([x.split('=')[-1].strip(' `') for x in aggregation])

    syn_kt = (kt.filter('`va.vep.worst_csq` == "synonymous_variant"')
              .aggregate_by_key(aggregation,
                                'proportion_singleton = `va.singleton`.sum()/v.count()'))

    syn_pd = syn_kt.join(mutation_kt).to_pandas()
    lm = smf.ols(formula='proportion_singleton ~ mutation_rate', data=syn_pd).fit()
    slope = lm.params['mutation_rate']
    intercept = lm.params['Intercept']

    grouping = ['worst_csq = `va.vep.worst_csq`']
    if additional_groupings is not None:
        if isinstance(additional_groupings, str):
            grouping.append(additional_groupings)
        else:
            grouping.extend(additional_groupings)

    return (syn_kt, kt.join(mutation_kt)
            .annotate('expected_ps = %s * mutation_rate + %s' % (slope, intercept))
            .aggregate_by_key(grouping,
                              ['num_singletons = `va.singleton`.sum()',
                               'num_variants = v.count()',
                               'total_expected_ps = expected_ps.sum()'])
            .annotate(['raw_ps = num_singletons/num_variants',
                       'maps = (num_singletons - total_expected_ps)/num_variants'])
            .annotate(['ps_sem = sqrt(raw_ps*(1-raw_ps)/num_variants)']))


def rebin_methylation(vds, bins=20):
    """
    Rebins va.methylation.level

    :param VariantDataset vds: VDS with original va.methylation.level
    :param int bins: Number of bins
    :return: VDS with new va.methylation.level number of bins
    :rtype: VariantDataset
    """
    return vds.annotate_variants_expr('va.methylation.level = '
                                      'if (v.altAllele().isTransition()) '
                                      '    range({}, -1, -1).find(e => va.methylation.value*{} >= e) '
                                      'else NA: Int'.format(str(bins - 1), bins))


def main(args):

    send_message(args.slack_channel, 'Started constraint process...')
    if args.methylation:
        mutation_rate_kt_path = 'gs://gnomad-resources/constraint/mutation_rate_methylation.kt'
        synonymous_kt_depth_path = 'gs://gnomad-resources/constraint/syn_depth_explore_methylation.kt'
        synonymous_kt_path = 'gs://gnomad-resources/constraint/syn_methylation.kt'
        full_kt_path = 'gs://gnomad-resources/constraint/constraint_methylation.kt'
    else:
        mutation_rate_kt_path = 'gs://gnomad-resources/constraint/mutation_rate.kt'
        synonymous_kt_depth_path = 'gs://gnomad-resources/constraint/syn_depth_explore.kt'
        synonymous_kt_path = 'gs://gnomad-resources/constraint/syn.kt'
        full_kt_path = 'gs://gnomad-resources/constraint/constraint.kt'

    if args.generate_fasta_vds:
        import_fasta_and_vep(fasta_path, context_vds_path, args.overwrite)

    if args.pre_process_data:
        pre_process_all_data()

    context_vds = hc.read(context_vds_path).filter_intervals(Interval.parse('1-22'))
    genome_vds = hc.read(genome_vds_path).filter_intervals(Interval.parse('1-22'))
    exome_vds = hc.read(exome_vds_path).filter_intervals(Interval.parse('1-22'))

    context_vds = rebin_methylation(context_vds)
    genome_vds = rebin_methylation(filter_to_pass(genome_vds))
    exome_vds = rebin_methylation(filter_to_pass(exome_vds))

    genome_vds = genome_vds.filter_variants_expr(GENOME_COVERAGE_CRITERIA)
    genome_vds = genome_vds.filter_variants_expr(AF_CRITERIA)
    exome_vds = exome_vds.filter_variants_expr(AF_CRITERIA)

    if args.run_sanity_checks:
        run_sanity_checks(context_vds)
        run_sanity_checks(genome_vds)
        run_sanity_checks(exome_vds, exome=False, csq_queries=True)
        proportion_observed = (
            get_proportion_observed(exome_vds, context_vds, trimer=True, methylation=args.methylation)
            .filter('"[ATCG]{3}" ~ context')
            .to_pandas().sort('proportion_observed', ascending=False)
        )
        print(proportion_observed)

    if args.calculate_mutation_rate:
        # TODO: PCR-free only
        mutation_kt = calculate_mutation_rate(context_vds,
                                              genome_vds,
                                              criteria='isDefined(va.gerp) && va.gerp < 0 && isMissing(va.vep.transcript_consequences)',
                                              methylation=args.methylation,
                                              trimer=True)
        mutation_kt.repartition(1).write(mutation_rate_kt_path, overwrite=args.overwrite)
        hc.read_table(mutation_rate_kt_path).export(mutation_rate_kt_path.replace('.kt', '.txt.bgz'))
        send_message(args.slack_channel, 'Mutation rate calculated!')

    mutation_kt = hc.read_table(mutation_rate_kt_path)

    if args.calibrate_raw_model:
        # First, get raw depth-uncorrected equation from only high coverage exons
        syn_kt = get_observed_expected_kt(exome_vds,
                                          context_vds, mutation_kt, canonical=True,
                                          criteria='annotation == "synonymous_variant"',
                                          coverage_cutoff=HIGH_COVERAGE_CUTOFF,
                                          methylation=args.methylation
        )
        syn_kt = remove_ttn(syn_kt)
        syn_kt.repartition(10).write(synonymous_kt_path, overwrite=args.overwrite)
        hc.read_table(synonymous_kt_path).export(synonymous_kt_path.replace('.kt', '.txt.bgz'))
        send_message(args.slack_channel, 'Raw model calibrated!')

    syn_kt = hc.read_table(synonymous_kt_path)
    syn_model = build_synonymous_model(syn_kt)
    print('\nRegression model weights: ')
    pprint(syn_model)

    if args.calibrate_coverage_model:
        # Then, get dependence of depth on O/E rate
        # Could try to save ~20 mins on 400 cores by combining this with raw model (need to apply regression_weights)
        syn_kt_by_coverage = get_observed_expected_kt(exome_vds,
                                                      context_vds, mutation_kt, canonical=True,
                                                      criteria='annotation == "synonymous_variant"',
                                                      methylation=args.methylation,
                                                      regression_weights=syn_model)
        (syn_kt_by_coverage.repartition(10).aggregate_by_key('median_coverage = median_coverage',
                                                             ['observed = variant_count.sum()',
                                                              'expected = expected_variant_count.sum()'])
         .write(synonymous_kt_depth_path, overwrite=args.overwrite))
        hc.read_table(synonymous_kt_depth_path).export(synonymous_kt_depth_path.replace('.kt', '.txt.bgz'))
        send_message(args.slack_channel, 'Coverage model calibrated!')

    synonymous_kt_depth = hc.read_table(synonymous_kt_depth_path)
    coverage_weights = get_coverage_weights(synonymous_kt_depth)
    print('\nCoverage model weights: ')
    pprint(coverage_weights)

    if args.build_full_model:
        full_kt = get_observed_expected_kt(exome_vds, context_vds, mutation_kt,
                                           methylation=args.methylation,
                                           regression_weights=syn_model,
                                           coverage_weights=coverage_weights)
        (full_kt.repartition(10)
         .aggregate_by_key(['transcript = transcript',
                            'annotation = annotation'],
                           ['variant_count = variant_count.sum()',
                            'sum_coverage = sum_coverage.sum()',
                            'possible_variant_count = possible_variant_count.sum()',
                            'aggregate_mutation_rate = aggregate_mutation_rate.sum()',
                            'expected_variant_count = expected_variant_count.sum()',
                            'expected_variant_count_adj = expected_variant_count_adj.sum()',
                            ])
         .write(full_kt_path, overwrite=args.overwrite))
        full_kt = hc.read_table(full_kt_path)
        full_kt.export(full_kt_path.replace('.kt', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--generate_fasta_vds', help='Generate FASTA VDS', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process all data (context, genome, exome)', action='store_true')
    parser.add_argument('--run_sanity_checks', help='Run sanity checks on all VDSes', action='store_true')
    parser.add_argument('--calculate_mutation_rate', help='Calculate mutation rate', action='store_true')
    parser.add_argument('--calibrate_raw_model', help='Re-calibrate model against synonymous variants', action='store_true')
    parser.add_argument('--calibrate_coverage_model', help='Calculate coverage model', action='store_true')
    parser.add_argument('--build_full_model', help='Build full model', action='store_true')
    parser.add_argument('--methylation', help='Use methylation model', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
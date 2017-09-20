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
final_exome_vds = 'gs://gnomad-exomes-hail01/sites/170622_new/gnomad.exomes.sites.vds'
final_genome_vds = 'gs://gnomad-public/release-170228/gnomad.genomes.r2.0.1.sites.vds'
exome_downample_vds_path = 'gs://gnomad-exomes-hail01/subsets/random_subsamples/gnomad.exomes.subsamples.sites.vds'
genome_downample_vds_path = 'gs://gnomad-genomes/subsets/random_subsamples/gnomad.genomes.subsamples.pcr_free.sites.vds'

fasta_path = "gs://gnomad-resources/Homo_sapiens_assembly19.fasta"
gerp_annotations_path = 'gs://annotationdb/cadd/cadd.kt'  # Gerp is in here as gerpS
raw_context_vds_path = "gs://gnomad-resources/Homo_sapiens_assembly19.fasta.snps_only.vep.vds"
genome_coverage_kt_path = 'gs://gnomad-resources/genome_coverage.kt'
exome_coverage_kt_path = 'gs://gnomad-resources/exome_coverage.kt'

# Processed datasets
context_vds_path = 'gs://gnomad-resources/constraint/context_processed.vds'
genome_vds_path = 'gs://gnomad-resources/constraint/genome_processed.vds'
exome_vds_path = 'gs://gnomad-resources/constraint/exome_processed.vds'

mutation_rate_kt_path = 'gs://gnomad-resources/constraint/mutation_rate.kt'
po_coverage_kt_path = 'gs://gnomad-resources/constraint/new/prop_observed_by_coverage.kt'
po_kt_path = 'gs://gnomad-resources/constraint/new/prop_observed.kt'

a_based_annotations = ['va.info.AC', 'va.info.AC_raw']

HIGH_COVERAGE_CUTOFF = 0.9
AF_CRITERIA = 'va.info.AN > 0 && va.info.AC/va.info.AN < 0.001'
GENOME_COVERAGE_CRITERIA = 'va.coverage.genome.mean >= 15 && va.coverage.genome.mean <= 60'
VARIANT_TYPES_FOR_MODEL = ('ACG', 'TCG', 'CCG', 'GCG', 'non-CpG')

EXOME_DOWNSAMPLINGS = ['10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '15000', '20000', '25000', '30000', '35000', '40000', '45000', '50000', '55000', '60000', '65000', '70000', '75000', '80000', '85000', '90000', '95000', '100000', '105000', '110000', '115000', '120000', '123136']
GENOME_DOWNSAMPLINGS = ['10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000']


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

    CONTIG_GROUPS = ('1', '2', '3', '4', '5', '6', '7', '8-9', '10-11', '12-13', '14-16', '17-18', '19-20', '21', '22', 'X', 'Y')
    # should have been: ('1', '2', '3', '4', '5', '6', '7', '8-9', '10-11', '12-13', '14-16', '17-19', '20-22', 'X', 'Y')

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
    # Upstream pre-processing of context VDS
    exome_coverage_kt = hc.read_table(exome_coverage_kt_path).select(['locus', '10'])
    genome_coverage_kt = hc.read_table(genome_coverage_kt_path).select(['locus', '10'])
    methylation_kt = hc.read_table(methylation_kt_path).select(['locus', 'MEAN'])
    gerp_kt = hc.read_table(gerp_annotations_path).select(['variant', 'GerpS'])
    context_vds = process_consequences(hc.read(raw_context_vds_path)
                                       .split_multi()
                                       .annotate_variants_expr(index_into_arrays(vep_root='va.vep'))
                                       .annotate_variants_table(gerp_kt, root='va.gerp')
                                       .annotate_variants_table(exome_coverage_kt, root='va.coverage.exome')
                                       .annotate_variants_table(genome_coverage_kt, root='va.coverage.genome')
                                       .annotate_variants_table(methylation_kt, root='va.methylation.value')
    )
    context_vds = rebin_methylation(context_vds)
    context_vds.write(context_vds_path, overwrite=args.overwrite)
    context_vds = hc.read(context_vds_path)

    # Genome processing
    genome_vds = (hc.read(final_genome_vds).split_multi()
                  .annotate_variants_expr(index_into_arrays(a_based_annotations))
                  .annotate_variants_vds(context_vds,
                                         expr='va.context = vds.context, va.gerp = vds.gerp, va.vep = vds.vep, '
                                              'va.coverage = vds.coverage, va.methylation = vds.methylation'))
    genome_ds_vds = hc.read(genome_downample_vds_path).split_multi()
    genome_ds_vds = genome_ds_vds.annotate_variants_expr(index_into_arrays(r_based_annotations=['va.calldata.raw.n{}.AC'.format(x) for x in GENOME_DOWNSAMPLINGS], drop_ref_ann=True))
    genome_vds = genome_vds.annotate_variants_vds(genome_ds_vds, expr='va.ds = vds.calldata.raw')
    genome_vds.write(genome_vds_path, overwrite=args.overwrite)

    # Exome processing
    exome_vds = (hc.read(final_exome_vds).split_multi()
                 .annotate_variants_expr(index_into_arrays(a_based_annotations))
                 .annotate_variants_vds(context_vds,
                                        expr='va.context = vds.context, va.gerp = vds.gerp, va.vep = vds.vep, '
                                             'va.coverage = vds.coverage, va.methylation = vds.methylation'))
    exome_ds_vds = hc.read(exome_downample_vds_path).split_multi()
    exome_ds_vds = exome_ds_vds.annotate_variants_expr(index_into_arrays(r_based_annotations=['va.calldata.raw.n{}.AC'.format(x) for x in EXOME_DOWNSAMPLINGS], drop_ref_ann=True))
    exome_vds = exome_vds.annotate_variants_vds(exome_ds_vds, expr='va.ds = vds.calldata.raw')
    exome_vds.write(exome_vds_path, overwrite=args.overwrite)


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
        exome_intervals = KeyTable.import_interval_list(exome_calling_intervals_path)
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
        return ('let base = {} in '
                'if (base == "A") "T" else '
                'if (base == "T") "A" else '
                'if (base == "C") "G" else '
                'if (base == "G") "C" else base'.format(root)
        )

    return kt.annotate(['`va.ref` = if (v.ref == "T" || v.ref == "G") {} else v.ref'.format(flip_text('v.ref')),
                        '`va.alt` = if (v.ref == "T" || v.ref == "G") {} else v.alt'.format(flip_text('v.alt')),
                        '`va.context` = if (v.ref == "T" || v.ref == "G") ({}) + ({}) + ({}) '
                        'else `va.context`'.format(
                            flip_text('`va.context`[2]'), flip_text('`va.context`[1]'), flip_text('`va.context`[0]'))])


def set_kt_cols_to_zero(kt, cols, variable_type=int):
    """
    Sets values to zero if missing

    :param KeyTable kt: Input Keytable
    :param list of str cols: List of columns to set to zero if missing
    :return: Keytable with columns set to zero if missing
    :rtype: KeyTable
    """
    else_string = 'L' if variable_type == int else '.0'
    return kt.annotate(['{0} = orElse({0}, 0{1})'.format(x, else_string) for x in cols])


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


def annotate_variant_types(kt):
    """
    Adds cpg, transition, and variant_type, variant_type_model columns

    :param KeyTable kt: input kt
    :return: Keytable with cpg, transition, variant_type, variant_type_model columns
    :rtype KeyTable
    """
    return (kt
            .annotate(['transition = (ref == "A" && alt == "G") || (ref == "G" && alt == "A") || (ref == "T" && alt == "C") || (ref == "C" && alt == "T")',
                       'cpg = (ref == "G" && alt == "A" && context[:2] == "CG") || (ref == "C" && alt == "T" && context[1:] == "CG")'])
            .annotate('variant_type = if (cpg) "CpG" else if (transition) "non-CpG transition" else "transversion"')
            # .annotate('variant_type_model = if (cpg) context else variant_type'))
            .annotate('variant_type_model = if (cpg) context else "non-CpG"'))


def filter_vep(vds, canonical=True, synonymous=True):
    if canonical: vds = filter_vep_to_canonical_transcripts(vds)
    vds = process_consequences(vds)
    if synonymous: vds = filter_vep_to_synonymous_variants(vds)

    return (vds.filter_variants_expr('!va.vep.transcript_consequences.isEmpty')
            .annotate_variants_expr('va.vep = select(va.vep, transcript_consequences)'))


def round_coverage(kt, keys):
    return (kt
            .annotate('coverage = (coverage*100).toInt()/100')
            .aggregate_by_key(['{0} = {0}'.format(x) for x in keys],
                              ['variant_count = variant_count.sum()',
                               'possible_variant_count = possible_variant_count.sum()'])
            .annotate('proportion_observed = variant_count/possible_variant_count'))


def get_high_coverage_kt(kt, keys):
    return (kt
            .filter('coverage >= {}'.format(HIGH_COVERAGE_CUTOFF))
            .aggregate_by_key(['{0} = {0}'.format(x) for x in keys],
                              'high_coverage_proportion_observed = variant_count.sum()/possible_variant_count.sum()')
            .select(keys + ['high_coverage_proportion_observed']))


def build_coverage_model(po_coverage_rounded_kt, keys):
    """

    :param KeyTable po_coverage_rounded_kt: kt
    :param KeyTable po_high_coverage_kt:
    :return:
    :rtype: dict
    """
    po_high_coverage_kt = get_high_coverage_kt(po_coverage_rounded_kt, keys)
    po_coverage_kt = po_coverage_rounded_kt.join(po_high_coverage_kt).annotate('scaled_proportion_observed = proportion_observed/high_coverage_proportion_observed')
    low_coverage_kt = po_coverage_kt.filter('coverage > 0 && coverage < {}'.format(HIGH_COVERAGE_CUTOFF))
    low_coverage_kt = low_coverage_kt.aggregate_by_key(['log_coverage = log(coverage)'],
                                                       ['log_mean_scaled_proportion_observed = log(scaled_proportion_observed.stats().mean)'])
    cov_slope, cov_intercept = calculate_coverage_model(low_coverage_kt)
    return {
        'coverage': cov_slope,
        'intercept': cov_intercept,
        'high_cutoff': HIGH_COVERAGE_CUTOFF
    }


def calculate_coverage_model(coverage_kt):
    """
    Calibrates mutational model to synonymous variants in gnomAD exomes

    :param KeyTable coverage_kt: Synonymous keytable
    :return:
    """
    coverage_pd = coverage_kt.to_pandas()
    lm = smf.ols(formula='log_mean_scaled_proportion_observed ~ log_coverage', data=coverage_pd).fit()
    slope = lm.params['log_coverage']
    intercept = lm.params['Intercept']
    return slope, intercept


def build_plateau_models(kt):
    """
    Calibrates mutational model to synonymous variants in gnomAD exomes

    :param KeyTable kt: Synonymous keytable
    :return:
    """
    output = {}
    for variant_type_model in VARIANT_TYPES_FOR_MODEL:
        high_coverage_pd = kt.filter('variant_type_model == "{}"'.format(variant_type_model)).to_pandas()
        lm = smf.ols(formula='high_coverage_proportion_observed ~ mutation_rate', data=high_coverage_pd).fit()
        output[variant_type_model] = dict(lm.params)
    return output


def count_variants_by_grouping(vds, grouping, criteria=None, trimer=False, singletons=False, downsample=False, explode=None, mutation_kt=None, regression_weights=None, coverage_weights=None, partitions=100):
    """
    Counts variants in VDS by provided groupings

    :param VariantDataset vds: Input VDS
    :param str criteria: Additional filtering criteria
    :param list of str grouping: List of variables to pass to key_exprs in aggregate_by_key
    :param bool trimer: whether to use trimer context (default heptamer)
    :param bool singletons: whether to split by singletons
    :param bool downsample: whether to use downsampled data
    :param str explode: criteria to explode by (most likely va.vep.transcript_consequences)
    :param KeyTable mutation_kt: Mutation rate KeyTable
    :param dict of dict regression_weights: Regression weights dict {'ACG': {'mutation_rate': slope, 'Intercept': intercept}, ...}
    :param dict coverage_weights: Coverage weights dict {'high_coverage': cutoff, 'coverage': slope, 'Intercept': intercept}
    :param int partitions: num_partitions to repartition by (if 0, no repartitioning)
    :return: keytable with counts as `variant_count`
    :rtype: KeyTable
    """
    if criteria is not None:
        vds = vds.filter_variants_expr(criteria)
    if trimer:
        vds = vds.annotate_variants_expr('va.context = va.context[2:5]')

    kt = collapse_strand(vds.variants_table().flatten())

    if explode:
        kt = kt.explode(explode).flatten()

    aggregation_functions = ['variant_count = v.count()']
    if singletons:
        aggregation_functions.append('singleton_count = v.filter(`va.info.AC_Raw` == 1).count()')
    if downsample:
        aggregation_functions.extend(['variant_count_n{0} = v.filter(v => `va.ds.n{0}.AC` > 0).count()'.format(x) for x in EXOME_DOWNSAMPLINGS])
        if singletons:
            aggregation_functions.extend(['singleton_count_n{0} = v.filter(v => `va.ds.n{0}.AC` == 1).count()'.format(x) for x in EXOME_DOWNSAMPLINGS])

    if mutation_kt:
        keys = ['context', 'ref', 'alt', 'methylation_level']
        kt = kt.annotate(['context = `va.context`', 'ref = `va.ref`', 'alt = `va.alt`', 'methylation_level = `va.methylation.level`'])
        kt = (kt.key_by(keys)
              .broadcast_left_join_distinct(mutation_kt.select(keys + ['mutation_rate'])))
        kt = annotate_variant_types(kt)
        kt = kt.annotate(
            'adjusted_mutation_rate = {} NA: Double'.format(' '.join(['if (variant_type_model == "{mut_type}") (mutation_rate * {mutation_rate} + {Intercept}) else'.format(mut_type=k, **v) for k, v in regression_weights.items()]))
        )
        aggregation_functions.extend(['mutation_rate = mutation_rate.sum()', 'adjusted_mutation_rate = adjusted_mutation_rate.sum()'])
    else:
        return kt.aggregate_by_key(grouping, aggregation_functions).repartition(partitions)

    kt = kt.aggregate_by_key(grouping, aggregation_functions).repartition(partitions)

    aggregation_functions[0] = 'variant_count = variant_count.sum()'
    if coverage_weights:
        kt = kt.annotate('expected_variant_count = '
                         'if (coverage >= {high_cutoff}) '
                         '  adjusted_mutation_rate '
                         'else '
                         '  (adjusted_mutation_rate * exp(log(coverage) * {coverage} + {intercept}))'.format(**coverage_weights))

        aggregation_functions.append('expected_variant_count = expected_variant_count.sum()')
    kt = kt.aggregate_by_key(['annotation = annotation', 'modifier = modifier', 'transcript = transcript'], aggregation_functions)
    return kt


def calculate_mutation_rate(possible_variants_vds, genome_vds, criteria=None, trimer=False):
    """
    Calculate mutation rate from all possible variants vds and observed variants vds
    Currently actually calculating more like "expected_proportion_variants"

    :param VariantDataset possible_variants_vds: synthetic VDS
    :param VariantDataset genome_vds: gnomAD WGS VDS
    :param bool trimer: whether to use trimer context (default heptamer)
    :return: keytable with mutation rates as `mutation_rate`
    :rtype: KeyTable
    """

    possible_variants_vds = (possible_variants_vds
                             .annotate_variants_vds(genome_vds
                                                    .filter_variants_expr('va.ds.n1000.AC > 1'),
                                                    expr='va.common_in_genome = !isMissing(vds)')
                             .filter_variants_expr('!va.common_in_genome'))
    genome_vds = genome_vds.filter_variants_expr('va.ds.n1000.AC == 1')
    # TODO: Switch to using common variants (remove singletons)

    grouping = ['context = `va.context`', 'ref = `va.ref`',
                'alt = `va.alt`', 'methylation_level = `va.methylation.level`']
    all_possible_kt = count_variants_by_grouping(possible_variants_vds, grouping, criteria=criteria, trimer=trimer)
    observed_kt = count_variants_by_grouping(genome_vds, grouping, criteria=criteria, trimer=trimer)

    kt = (all_possible_kt.rename({'variant_count': 'possible_variants'})
          .join(observed_kt, how='outer')
          .annotate('mutation_rate = variant_count/possible_variants'))

    return kt.filter('ref.length == 1 && alt.length == 1 && !("N" ~ context)')


def get_proportion_observed(vds, all_possible_vds, mutation_kt,
                            regression_weights, coverage_weights, mode='transcript',
                            canonical=False, synonymous=False, downsample=False):
    """
    Does what it says

    :param VariantDataset vds: VDS to generate observed counts
    :param VariantDataset all_possible_vds: VDS with all possible variants
    :param KeyTable mutation_kt: Mutation rate keytable
    :param bool downsample: Whether to only use downsampled data in addition to full dataset
    :return: KeyTable with coverage, mu, and proportion_observed
    :rtype KeyTable
    """
    vds = filter_vep(vds, canonical=canonical, synonymous=synonymous)
    all_possible_vds = filter_vep(all_possible_vds, canonical=canonical, synonymous=synonymous)

    vds = vds.annotate_variants_expr('va = select(va, ds, context, methylation, vep, coverage)')
    all_possible_vds = all_possible_vds.annotate_variants_expr('va = select(va, context, methylation, vep, coverage)')

    if mode == 'transcript':
        grouping = [
            'annotation = `va.vep.transcript_consequences.most_severe_consequence`',
            'modifier = if (`va.vep.transcript_consequences.most_severe_consequence` == "missense_variant") '
            ' `va.vep.transcript_consequences.polyphen_prediction` else '
            'if (`va.vep.transcript_consequences.lof` != "") `va.vep.transcript_consequences.lof` '
            'else NA: String',
            'transcript = `va.vep.transcript_consequences.transcript_id`',
            'coverage = `va.coverage.exome`']  # va is flattened, so this is a little awkward, and now ref and alt are in va.
        kt = count_variants_by_grouping(vds, grouping, explode='va.vep.transcript_consequences', trimer=True)
        all_possible_kt = count_variants_by_grouping(all_possible_vds, grouping, mutation_kt=mutation_kt,
                                                     regression_weights=regression_weights,
                                                     coverage_weights=coverage_weights,
                                                     explode='va.vep.transcript_consequences', trimer=True)
    elif mode == 'coverage':
        grouping = ['context = `va.context`', 'ref = `va.ref`',
                    'alt = `va.alt`', 'methylation_level = `va.methylation.level`',
                    'coverage = `va.coverage.exome`']
        kt = count_variants_by_grouping(vds, grouping=grouping,
                                        explode='va.vep.transcript_consequences', trimer=True, downsample=downsample)
        all_possible_kt = count_variants_by_grouping(all_possible_vds, grouping=grouping,
                                                     explode='va.vep.transcript_consequences', trimer=True)

    columns = ['variant_count']
    if downsample:
        columns.extend(['variant_count_n{}'.format(x) for x in EXOME_DOWNSAMPLINGS])
    final_kt = set_kt_cols_to_zero(
        kt.join(
            all_possible_kt.rename({'variant_count': 'possible_variant_count'}),
            how='right'), columns)
    if mode == 'coverage':
        keys = ['context', 'ref', 'alt', 'methylation_level']
        return final_kt.key_by(keys).join(mutation_kt.select(keys + ['mutation_rate']), how='outer')
    elif mode == 'transcript':
        return final_kt


def maps(vds, mutation_kt, additional_groupings=None, trimer=True):
    """

    :param VariantDataset vds: VDS
    :return:
    """

    if trimer: vds = vds.annotate_variants_expr('va.context = va.context[2:5]')

    aggregation = ['context = `va.context`', 'ref = `va.ref`', 'alt = `va.alt`', 'methylation_level = `va.methylation.level`']

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


def main(args):

    if args.generate_fasta_vds:
        import_fasta_and_vep(fasta_path, context_vds_path, args.overwrite)

    if args.pre_process_data:
        pre_process_all_data()

    context_vds = hc.read(context_vds_path).filter_intervals(Interval.parse('1-22'))
    exome_vds = hc.read(exome_vds_path).filter_intervals(Interval.parse('1-22'))

    exome_vds = filter_rf_variants(exome_vds.filter_variants_expr(AF_CRITERIA))

    if args.calculate_mutation_rate:
        # TODO: PCR-free only
        genome_vds = hc.read(genome_vds_path).filter_intervals(Interval.parse('1-22'))
        genome_vds = filter_rf_variants(genome_vds)

        mutation_kt = calculate_mutation_rate(context_vds.filter_variants_expr(GENOME_COVERAGE_CRITERIA),
                                              genome_vds.filter_variants_expr(GENOME_COVERAGE_CRITERIA),
                                              criteria='isDefined(va.gerp) && va.gerp < 0 && isMissing(va.vep.transcript_consequences)',
                                              trimer=True)
        mutation_kt.repartition(1).write(mutation_rate_kt_path, overwrite=args.overwrite)
        hc.read_table(mutation_rate_kt_path).export(mutation_rate_kt_path.replace('.kt', '.txt.bgz'))
        send_message(args.slack_channel, 'Mutation rate calculated!')

    mutation_kt = hc.read_table(mutation_rate_kt_path)

    if args.get_mu_coverage:
        po_coverage_kt = get_proportion_observed(exome_vds, context_vds, mutation_kt, mode='coverage', downsample=args.downsample)
        po_coverage_kt.write(po_coverage_kt_path, overwrite=True)
        hc.read_table(po_coverage_kt_path).export(po_coverage_kt_path.replace('.kt', '.txt.bgz'))

    po_coverage_kt = hc.read_table(po_coverage_kt_path)
    po_coverage_kt = annotate_variant_types(po_coverage_kt)

    keys = ['context', 'ref', 'alt', 'methylation_level', 'mutation_rate', 'cpg', 'transition', 'variant_type', 'variant_type_model']
    po_coverage_rounded_kt = round_coverage(po_coverage_kt, keys + ['coverage']).key_by(keys)

    plateau_models = build_plateau_models(get_high_coverage_kt(po_coverage_rounded_kt, keys))
    coverage_model = build_coverage_model(po_coverage_rounded_kt, keys)

    # Used to confirm model
    # get_proportion_observed_by_transcript(exome_vds, context_vds, mutation_kt, plateau_models, coverage_model, canonical=True, synonymous=True).write(po_kt_path)
    get_proportion_observed(exome_vds, context_vds, mutation_kt, plateau_models, coverage_model).write(po_kt_path, overwrite=args.overwrite)
    hc.read_table(po_kt_path).export(po_kt_path.replace('.kt', '.txt.bgz'))

    # Correlate mu from genomes with the plateau to impute plateau for unobserved contexts

    send_message('@konradjk', 'Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--generate_fasta_vds', help='Generate FASTA VDS', action='store_true')
    parser.add_argument('--downsample', help='Use downsampled data', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process all data (context, genome, exome)', action='store_true')
    parser.add_argument('--calculate_mutation_rate', help='Calculate mutation rate', action='store_true')
    parser.add_argument('--get_mu_coverage', help='Calculate proportion observed by mu by coverage', action='store_true')
    parser.add_argument('--build_coverage_model', help='Build coverage model', action='store_true')
    parser.add_argument('--calculate_mu_summary', help='Calculate proportion observed by mu', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
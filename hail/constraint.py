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

CONTIG_GROUPS = ('1', '2', '3', '4', '5', '6', '7', '8-9', '10-11', '12-13', '14-16', '17-18', '19-20', '21', '22', 'X', 'Y')
# should have been: ('1', '2', '3', '4', '5', '6', '7', '8-9', '10-11', '12-13', '14-16', '17-19', '20-22', 'X', 'Y')
a_based_annotations = ['va.info.AC', 'va.info.AC_raw']

HIGH_COVERAGE_CUTOFF = 0.9
AF_CRITERIA = 'va.info.AN > 0 && va.info.AC/va.info.AN < 0.001'
GENOME_COVERAGE_CRITERIA = 'va.coverage.genome.mean >= 15 && va.coverage.genome.mean <= 60'
VARIANT_TYPES_FOR_MODEL = ('ACG', 'TCG', 'CCG', 'GCG', 'non-CpG')

EXOME_DOWNSAMPLINGS = ['10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '15000', '20000', '25000', '30000', '35000', '40000', '45000', '50000', '55000', '60000', '65000', '70000', '75000', '80000', '85000', '90000', '95000', '100000', '105000', '110000', '115000', '120000', '123136']
GENOME_DOWNSAMPLINGS = ['10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000']
exome_downample_vds_path = 'gs://gnomad-exomes-hail01/subsets/random_subsamples/gnomad.exomes.subsamples.sites.vds'
genome_downample_vds_path = 'gs://gnomad-genomes/subsets/random_subsamples/gnomad.genomes.subsamples.pcr_free.sites.vds'


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


def count_variants(vds, criteria=None, additional_groupings=None, trimer=False, explode=None, coverage=True, singletons=False, downsample=False, partitions=100):
    """
    Counts variants in VDS by context, ref, alt, and any other groupings provided

    :param VariantDataset vds: Input VDS
    :param str criteria: Any filtering criteria (e.g. non-coding, non-conserved), to be passed to filter_variants_expr
    :param additional_groupings: Whether to group further (e.g. by functional annotation)
    :type additional_groupings: str or list of str
    :param bool trimer: whether to use trimer context (default heptamer)
    :param str explode: criteria to explode by (most likely va.vep.transcript_consequences)
    :param bool or dict coverage: Whether to also aggregate coverage information (allowed to be dict for convenience)
    :param bool singletons: Whether to also split out singletons
    :param bool downsample: Whether to use downsampled data in addition to full model
    :return: keytable with counts as `variant_count`
    :rtype: KeyTable
    """
    if criteria is not None:
        vds = vds.filter_variants_expr(criteria)

    if trimer:
        vds = vds.annotate_variants_expr('va.context = va.context[2:5]')

    grouping = ['context = `va.context`', 'ref = `va.ref`',
                'alt = `va.alt`', 'methylation_level = `va.methylation.level`']  # va is flattened, so this is a little awkward, and now ref and alt are in va.
    if additional_groupings is not None:
        if isinstance(additional_groupings, str):
            grouping.append(additional_groupings)
        else:
            grouping.extend(additional_groupings)

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
    if coverage: aggregation_functions.append('sum_coverage = `va.coverage.exome.median`.sum()')

    return kt.aggregate_by_key(grouping, aggregation_functions).repartition(partitions)


def count_variants_by_transcript(vds, criteria=None, trimer=False, explode=None, mutation_kt=None, regression_weights=None, coverage_weights=None, partitions=100):
    """
    Counts variants in VDS by context, ref, alt

    :param VariantDataset vds: Input VDS
    :param str criteria: Any filtering criteria (e.g. non-coding, non-conserved), to be passed to filter_variants_expr
    :param bool trimer: whether to use trimer context (default heptamer)
    :param str explode: criteria to explode by (most likely va.vep.transcript_consequences)
    :return: keytable with counts as `variant_count`
    :rtype: KeyTable
    """
    if criteria is not None:
        vds = vds.filter_variants_expr(criteria)

    if trimer:
        vds = vds.annotate_variants_expr('va.context = va.context[2:5]')

    grouping = [
        'annotation = `va.vep.transcript_consequences.most_severe_consequence`',
        'modifier = if (`va.vep.transcript_consequences.most_severe_consequence` == "missense_variant") '
        ' `va.vep.transcript_consequences.polyphen_prediction` else '
        'if (`va.vep.transcript_consequences.lof` != "") `va.vep.transcript_consequences.lof` '
        'else NA: String',
        'transcript = `va.vep.transcript_consequences.transcript_id`',
        'coverage = `va.coverage.exome`']  # va is flattened, so this is a little awkward, and now ref and alt are in va.

    kt = collapse_strand(vds.variants_table().flatten())

    if explode:
        kt = kt.explode(explode).flatten()

    aggregation_functions = ['variant_count = v.count()']

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

    all_possible_kt = count_variants(possible_variants_vds, criteria=criteria, trimer=trimer, coverage=False)
    observed_kt = count_variants(genome_vds, criteria=criteria, trimer=trimer, coverage=False)

    kt = (all_possible_kt.rename({'variant_count': 'possible_variants'})
          .join(observed_kt, how='outer')
          .annotate('mutation_rate = variant_count/possible_variants'))

    return kt.filter('ref.length == 1 && alt.length == 1 && !("N" ~ context)')


def collapse_counts_by_exon(kt, mutation_rate_weights=None, regression_weights=None, coverage_weights=None, split_singletons=False, fix_to_zero=False, downsample=False):
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
    if mutation_rate_weights:  # expected
        keys = ['context', 'ref', 'alt', 'methylation_level']
        kt = (kt.key_by(keys)
              .join(mutation_rate_weights.select(keys + ['mutation_rate']), how='outer')
              .annotate(['aggregate_mutation_rate = mutation_rate * variant_count']))
        aggregation_expression.extend(['aggregate_mutation_rate = aggregate_mutation_rate.sum()',
                                       'sum_coverage = sum_coverage.sum()'])
    else:  # observed
        if downsample:
            aggregation_expression.extend(['variant_count_n{0} = variant_count_n{0}.sum()'.format(x) for x in EXOME_DOWNSAMPLINGS])
    if split_singletons:
        aggregation_expression.append('singleton_count = singleton_count.sum()')
    kt = kt.aggregate_by_key(['transcript = transcript', 'exon = exon', 'annotation = annotation'], aggregation_expression).repartition(10)

    if regression_weights is not None:
        annotation = ['median_coverage = sum_coverage//variant_count']
        if downsample:
            annotation.extend(['expected_variant_count_n{n} = {slope} * aggregate_mutation_rate + {intercept}'.format(n=k, **v) for k, v in regression_weights.items() if k is not 'full'])
            regression_weights = regression_weights['full']
        annotation.append('expected_variant_count = {slope} * aggregate_mutation_rate + {intercept}'.format(**regression_weights))
        kt = kt.annotate(annotation)
        if fix_to_zero:
            kt = kt.annotate('expected_variant_count = if (expected_variant_count > 0.0) expected_variant_count else 0.0')
        if coverage_weights is not None:
            if downsample:
                kt = kt.annotate(['expected_variant_count_adj_n{n} = '
                                  'if (median_coverage >= {high_cutoff}) '
                                  '  expected_variant_count_n{n} '
                                  'else '
                                  '  (expected_variant_count_n{n} * (log(median_coverage) * {mid_beta} + {mid_intercept}))'.format(n=k, **v)
                for k, v in coverage_weights.items() if k is not 'full'])
                coverage_weights = coverage_weights['full']
            kt = kt.annotate('expected_variant_count_adj = '
                             'if (median_coverage >= {high_cutoff}) '
                             '  expected_variant_count '
                             'else '
                             '  (expected_variant_count * (log(median_coverage) * {mid_beta} + {mid_intercept}))'.format(**coverage_weights))
    return kt


def get_coverage_weights(synonymous_kt_depth, high_cutoff=HIGH_COVERAGE_CUTOFF, low_cutoff=1, downsample=False):
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
    if downsample:
        output = {'full': output_dict}
        for x in EXOME_DOWNSAMPLINGS:
            output_dict = {
                'high_cutoff': high_cutoff,
                'low_cutoff': low_cutoff
            }
            synonymous_kt_depth = (synonymous_kt_depth
                                   .filter('median_coverage < {high_cutoff} && median_coverage > {low_cutoff}'.format(**output_dict))
                                   .annotate(['oe = observed_n{0}/expected_n{0}'.format(x), 'log_coverage = log(median_coverage)'])
            )
            synonymous_depth = synonymous_kt_depth.to_pandas()
            lm = smf.ols(formula='oe ~ log_coverage', data=synonymous_depth).fit()
            output_dict['mid_beta'] = lm.params['log_coverage']
            output_dict['mid_intercept'] = lm.params['Intercept']
            output[x] = output_dict
        return output
    else:
        return output_dict


def build_synonymous_model(syn_kt, downsample=False):
    """
    Calibrates mutational model to synonymous variants in gnomAD exomes

    :param KeyTable syn_kt: Synonymous keytable
    :return:
    """
    syn_pd = syn_kt.to_pandas()
    lm = smf.ols(formula='variant_count ~ aggregate_mutation_rate', data=syn_pd).fit()
    output_dict = {
        'slope': lm.params['aggregate_mutation_rate'],
        'intercept': lm.params['Intercept']
    }
    if downsample:
        output = {'full': output_dict}
        for ds in EXOME_DOWNSAMPLINGS:
            lm = smf.ols(formula='variant_count_n{} ~ aggregate_mutation_rate'.format(ds), data=syn_pd).fit()
            output[ds] = {
                'slope': lm.params['aggregate_mutation_rate'],
                'intercept': lm.params['Intercept']
            }
        return output
    else:
        return output_dict


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
                             criteria=None, coverage_cutoff=0, split_singletons=False,
                             coverage_weights=None, regression_weights=None, fix_to_zero=False, downsample=False):
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
        'modifier = if (`va.vep.transcript_consequences.most_severe_consequence` == "missense_variant") '
        ' `va.vep.transcript_consequences.polyphen_prediction` else '
        'if (`va.vep.transcript_consequences.lof` != "") `va.vep.transcript_consequences.lof` '
        'else NA: String',
        'transcript = `va.vep.transcript_consequences.transcript_id`',
        'exon = `va.vep.transcript_consequences.exon`']  # va gets flattened, so this a little awkward

    columns = ['v', 'va.context', 'va.ref', 'va.alt', 'va.methylation.level', 'va.coverage.exome.median', 'va.coverage.genome.median',
               'va.vep.transcript_consequences.most_severe_consequence',
               'va.vep.transcript_consequences.polyphen_prediction',
               'va.vep.transcript_consequences.lof',
               'va.vep.transcript_consequences.transcript_id',
               'va.vep.transcript_consequences.exon']  # TODO: select columns here earlier
    if downsample:
        columns.extend(['va.ds.n{}.AC'.format(x) for x in EXOME_DOWNSAMPLINGS])
    columns.append('va.info.AC')
    # kt = kt.select(columns)  # Temporary for hail speedups
    kt = count_variants(vds, additional_groupings=count_grouping, singletons=split_singletons,
                        explode='va.vep.transcript_consequences', trimer=True, downsample=downsample)
    all_possible_kt = count_variants(all_possible_vds,
                                     additional_groupings=count_grouping,
                                     explode='va.vep.transcript_consequences',
                                     trimer=True)
    if criteria:
        kt = kt.filter(criteria)
        all_possible_kt = all_possible_kt.filter(criteria)

    collapsed_kt = collapse_counts_by_exon(kt, split_singletons=split_singletons, downsample=downsample)
    collapsed_all_possible_kt = collapse_counts_by_exon(all_possible_kt,
                                                        mutation_rate_weights=mutation_kt,
                                                        regression_weights=regression_weights,
                                                        coverage_weights=coverage_weights, fix_to_zero=fix_to_zero,
                                                        downsample=downsample)
    # Calculating median_coverage only for "all possible" keytable means we get the exact mean for each exon
    collapsed_all_possible_kt = (collapsed_all_possible_kt
                                 .annotate('median_coverage = sum_coverage//variant_count')
                                 .rename({'variant_count': 'possible_variant_count'}))
    if coverage_weights is not None:
        coverage_cutoff = coverage_weights['full']['low_cutoff'] if downsample else coverage_weights['low_cutoff']
    if coverage_cutoff:
        collapsed_all_possible_kt = collapsed_all_possible_kt.filter('median_coverage >= {}'.format(coverage_cutoff))
    # return collapsed_kt.join(collapsed_all_possible_kt)
    columns = ['variant_count']
    if downsample:
        columns.extend(['variant_count_n{}'.format(x) for x in EXOME_DOWNSAMPLINGS])
    final_kt = set_kt_cols_to_zero(collapsed_kt.join(collapsed_all_possible_kt, how='right'), columns)

    if coverage_cutoff:
        final_kt = final_kt.filter('median_coverage >= {}'.format(coverage_cutoff))
    return final_kt


def filter_vep(vds, canonical=True, synonymous=True):
    if canonical: vds = filter_vep_to_canonical_transcripts(vds)
    vds = process_consequences(vds)
    if synonymous: vds = filter_vep_to_synonymous_variants(vds)

    return (vds.filter_variants_expr('!va.vep.transcript_consequences.isEmpty')
            .annotate_variants_expr('va.vep = select(va.vep, transcript_consequences)'))


def get_proportion_observed_by_coverage(vds, all_possible_vds, mutation_kt,
                                        canonical=True, synonymous=True, downsample=True):
    """
    Does what it says

    :param VariantDataset vds: VDS to generate observed counts
    :param VariantDataset all_possible_vds: VDS with all possible variants
    :param KeyTable mutation_kt: Mutation rate keytable
    :param bool canonical: Whether to only use canonical transcripts
    :param bool synonymous: Whether to only use only synonymous variants
    :param bool downsample: Whether to only use downsampled data in addition to full dataset
    :return: KeyTable with coverage, mu, and proportion_observed
    :rtype KeyTable
    """
    vds = filter_vep(vds, canonical=canonical, synonymous=synonymous)
    all_possible_vds = filter_vep(all_possible_vds, canonical=canonical, synonymous=synonymous)

    vds = vds.annotate_variants_expr('va = select(va, ds, context, methylation, vep, coverage)')
    all_possible_vds = all_possible_vds.annotate_variants_expr('va = select(va, context, methylation, vep, coverage)')

    # count_grouping = ['coverage = `va.coverage.exome.median`']
    count_grouping = ['coverage = `va.coverage.exome`']

    # kt = kt.select(columns)  # Temporary for hail speedups
    kt = count_variants(vds, additional_groupings=count_grouping, coverage=False,
                        explode='va.vep.transcript_consequences', trimer=True, downsample=downsample)
    all_possible_kt = count_variants(all_possible_vds,
                                     additional_groupings=count_grouping, coverage=False,
                                     explode='va.vep.transcript_consequences',
                                     trimer=True)

    columns = ['variant_count']
    if downsample:
        columns.extend(['variant_count_n{}'.format(x) for x in EXOME_DOWNSAMPLINGS])
    final_kt = set_kt_cols_to_zero(
        kt.join(
            all_possible_kt.rename({'variant_count': 'possible_variant_count'}),
            how='right'), columns)
    keys = ['context', 'ref', 'alt', 'methylation_level']
    return final_kt.key_by(keys).join(mutation_kt.select(keys + ['mutation_rate']), how='outer')


def get_proportion_observed_by_transcript(vds, all_possible_vds, mutation_kt, regression_weights, coverage_weights,
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

    kt = count_variants_by_transcript(vds, explode='va.vep.transcript_consequences', trimer=True)
    all_possible_kt = count_variants_by_transcript(all_possible_vds, mutation_kt=mutation_kt,
                                                   regression_weights=regression_weights,
                                                   coverage_weights=coverage_weights,
                                                   explode='va.vep.transcript_consequences',
                                                   trimer=True)

    columns = ['variant_count']
    if downsample:
        columns.extend(['variant_count_n{}'.format(x) for x in EXOME_DOWNSAMPLINGS])
    return set_kt_cols_to_zero(
        kt.join(
            all_possible_kt.rename({'variant_count': 'possible_variant_count'}),
            how='right'), columns)


def set_kt_cols_to_zero_float(kt, cols):
    """
    Sets values to zero (float) if missing

    :param KeyTable kt: Input Keytable
    :param list of str cols: List of columns to set to zero if missing
    :return: Keytable with columns set to zero if missing
    :rtype: KeyTable
    """
    return kt.annotate(['{0} = orElse({0}, 0.0)'.format(x) for x in cols])


def set_kt_cols_to_zero(kt, cols):
    """
    Sets values to zero if missing

    :param KeyTable kt: Input Keytable
    :param list of str cols: List of columns to set to zero if missing
    :return: Keytable with columns set to zero if missing
    :rtype: KeyTable
    """
    return kt.annotate(['{0} = orElse({0}, 0L)'.format(x) for x in cols])


def get_proportion_observed(exome_vds, all_possible_vds, trimer=False, downsample=False):
    """
    Intermediate function to get proportion observed by context, ref, alt

    :param VariantDataset exome_vds: gnomAD exome VDS
    :param VariantDataset all_possible_vds: VDS with all possible variants
    :return: Key Table with context, ref, alt, proportion observed
    :rtype KeyTable
    """
    grouping = ['annotation = `va.vep.transcript_consequences.most_severe_consequence`']  # va gets flattened, so this a little awkward
    exome_kt = count_variants(exome_vds.repartition(2000, shuffle=False),
                              additional_groupings=grouping,
                              explode='va.vep.transcript_consequences', trimer=trimer, downsample=downsample)
    all_possible_kt = count_variants(all_possible_vds.repartition(2000, shuffle=False),
                                     additional_groupings=grouping,
                                     explode='va.vep.transcript_consequences', trimer=trimer)

    full_kt = all_possible_kt.rename({'variant_count': 'possible_variants'}).join(exome_kt.drop('sum_coverage'), how='outer')
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


def rebin_va(vds, root, bins=100):
    """
    Rebins value in root (assumes values between 0 and 1)

    :param VariantDataset vds: VDS
    :param str root: Path to va. annotation to rebin
    :param int bins: Number of bins
    :return: VDS with rebinned root value
    :rtype: VariantDataset
    """
    return vds.annotate_variants_expr('{0} = range({1}, -1, -1).find(e => {0}*{2} >= e)/{2}'.format(root, bins - 1, bins))


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


def main(args):

    if args.generate_fasta_vds:
        import_fasta_and_vep(fasta_path, context_vds_path, args.overwrite)

    if args.pre_process_data:
        pre_process_all_data()

    context_vds = hc.read(context_vds_path).filter_intervals(Interval.parse('1-22'))
    exome_vds = hc.read(exome_vds_path).filter_intervals(Interval.parse('1-22'))

    exome_coverage_kt = hc.read_table(exome_coverage_kt_path).select(['locus', '10'])
    context_vds = context_vds.annotate_variants_table(exome_coverage_kt, root='va.coverage.exome')
    exome_vds = exome_vds.annotate_variants_table(exome_coverage_kt, root='va.coverage.exome')
    # segdups = KeyTable.import_bed(decoy_intervals_path)
    # lcrs = KeyTable.import_interval_list(lcr_intervals_path)
    context_vds = rebin_methylation(context_vds)
    exome_vds = rebin_methylation(filter_rf_variants(exome_vds))

    exome_vds = exome_vds.filter_variants_expr(AF_CRITERIA)

    exome_ds_vds = hc.read(exome_downample_vds_path).filter_intervals(Interval.parse('1-22')).split_multi()
    exome_ds_vds = exome_ds_vds.annotate_variants_expr(index_into_arrays(r_based_annotations=['va.calldata.raw.n{}.AC'.format(x) for x in EXOME_DOWNSAMPLINGS], drop_ref_ann=True))
    exome_vds = exome_vds.annotate_variants_vds(exome_ds_vds, expr='va.ds = vds.calldata.raw')

    # TODO: need to blacklist LCR and SEGDUP from expected throughout
    if args.calculate_mutation_rate:
        # TODO: PCR-free only
        genome_vds = hc.read(genome_vds_path).filter_intervals(Interval.parse('1-22'))
        genome_vds = rebin_methylation(filter_rf_variants(genome_vds))

        genome_ds_vds = hc.read(genome_downample_vds_path).filter_intervals(Interval.parse('1-22')).split_multi()
        genome_ds_vds = genome_ds_vds.annotate_variants_expr(index_into_arrays(r_based_annotations=['va.calldata.raw.n{}.AC'.format(x) for x in GENOME_DOWNSAMPLINGS], drop_ref_ann=True))
        genome_vds = genome_vds.annotate_variants_vds(genome_ds_vds, expr='va.ds = vds.calldata.raw')
        mutation_kt = calculate_mutation_rate(context_vds.filter_variants_expr(GENOME_COVERAGE_CRITERIA),
                                              genome_vds.filter_variants_expr(GENOME_COVERAGE_CRITERIA),
                                              criteria='isDefined(va.gerp) && va.gerp < 0 && isMissing(va.vep.transcript_consequences)',
                                              trimer=True)
        mutation_kt.repartition(1).write(mutation_rate_kt_path, overwrite=args.overwrite)
        hc.read_table(mutation_rate_kt_path).export(mutation_rate_kt_path.replace('.kt', '.txt.bgz'))
        send_message(args.slack_channel, 'Mutation rate calculated!')

    mutation_kt = hc.read_table(mutation_rate_kt_path)

    if args.get_mu_coverage:
        po_coverage_kt = get_proportion_observed_by_coverage(exome_vds, context_vds, mutation_kt, downsample=args.downsample)
        po_coverage_kt.write(po_coverage_kt_path, overwrite=True)
        hc.read_table(po_coverage_kt_path).export(po_coverage_kt_path.replace('.kt', '.txt.bgz'))

    po_coverage_kt = hc.read_table(po_coverage_kt_path)
    # po_coverage_kt = hc.read_table('prop_observed_by_coverage.kt')
    po_coverage_kt = annotate_variant_types(po_coverage_kt)

    keys = ['context', 'ref', 'alt', 'methylation_level', 'mutation_rate', 'cpg', 'transition', 'variant_type', 'variant_type_model']
    po_coverage_rounded_kt = round_coverage(po_coverage_kt, keys + ['coverage']).key_by(keys)

    plateau_models = build_plateau_models(get_high_coverage_kt(po_coverage_rounded_kt, keys))
    coverage_model = build_coverage_model(po_coverage_rounded_kt, keys)

    # Used to confirm model
    # get_proportion_observed_by_transcript(exome_vds, context_vds, mutation_kt, plateau_models, coverage_model, canonical=True, synonymous=True).write(po_kt_path)
    get_proportion_observed_by_transcript(exome_vds, context_vds, mutation_kt, plateau_models, coverage_model).write(po_kt_path)
    hc.read_table(po_kt_path).export(po_kt_path.replace('.kt', '.txt.bgz'))

    # Correlate mu from genomes with the plateau to impute plateau for unobserved contexts

    send_message('@konradjk', 'Done!')


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
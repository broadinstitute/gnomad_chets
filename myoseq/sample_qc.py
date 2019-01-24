from gnomad_hail import *
from os.path import splitext, basename
import hail as hl
import numpy as np
import hdbscan


def assign_platform_pcs(platform_pc_table: hl.Table, pc_scores_ann: str = 'scores') -> hl.Table:

    # Read and format data for clustering
    data = platform_pc_table.to_pandas()
    callrate_data = np.matrix(data[pc_scores_ann].tolist())
    logger.info('Assigning platforms to {} exome samples in MT...'.format(len(callrate_data)))

    # Cluster data
    clusterer = hdbscan.HDBSCAN(min_cluster_size=100)
    cluster_labels = clusterer.fit_predict(callrate_data)
    n_clusters = len(set(cluster_labels)) - (-1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info('Found {} unique platforms during platform imputation...'.format(n_clusters))

    data['qc_platform'] = cluster_labels
    return hl.Table.from_pandas(data)


def compute_callrate_mt(mt: hl.MatrixTable,  intervals: hl.Table) -> hl.MatrixTable:
    callrate_mt = filter_to_autosomes(mt)
    callrate_mt = callrate_mt.annotate_rows(interval=intervals[callrate_mt.locus].target)
    callrate_mt = callrate_mt.filter_rows(hl.is_defined(callrate_mt.interval) & (hl.len(callrate_mt.alleles) == 2))
    callrate_mt = callrate_mt.select_entries(GT=hl.or_missing(hl.is_defined(callrate_mt.GT), hl.struct()))
    callrate_mt = callrate_mt.group_rows_by(callrate_mt.interval).aggregate(callrate=hl.agg.fraction(hl.is_defined(callrate_mt.GT)))
    return  callrate_mt

def run_platform_pca(callrate_mt: hl.MatrixTable) -> hl.Table:
    callrate_mt = callrate_mt.annotate_entries(callrate=hl.int(callrate_mt.callrate > 0.25))
    # Center until Hail's PCA does it for you
    callrate_mt = callrate_mt.annotate_rows(mean_callrate=hl.agg.mean(callrate_mt.callrate))
    callrate_mt = callrate_mt.annotate_entries(callrate=callrate_mt.callrate - callrate_mt.mean_callrate)
    eigenvalues, scores, _ = hl.pca(callrate_mt.callrate, compute_loadings=False)
    logger.info('Eigenvalues: {}'.format(eigenvalues))
    # [731282566.2824697, 78687228.90071851, 43837650.51729764, 33969298.61827205, 26308703.539534636, 21102437.512725923, 16949828.555817757, 12994894.187041137, 8372332.274295175, 8128326.814388647]
    return scores


def impute_sex(
        x_mt: hl.MatrixTable,
        y_mt: hl.MatrixTable,
        platform_ht: hl.Table,
        male_min_f_stat: float,
        female_max_f_stat: float,
        min_male_y_sites_called: int,
        max_y_call_rate_stdev: int
) -> hl.Table:
    x = hl.filter_intervals(x_mt, [hl.parse_locus_interval('X')])
    x_ht = hl.impute_sex(x.GT, aaf_threshold=0.05, female_threshold=female_max_f_stat, male_threshold=male_min_f_stat)
    y = hl.filter_intervals(y_mt, [hl.parse_locus_interval('Y')])
    sex_ht = y.annotate_cols(
        qc_platform=platform_ht[y.col_key].qc_platform,
        is_female=x_ht[y.col_key].is_female,
        y_call_rate=hl.agg.fraction(hl.is_defined(y.GT)),
        n_y_sites_called=hl.agg.count_where(hl.is_defined(y.GT)),
        **{ann: f'x_{ann}' for ann in x_ht.row_value if ann != 'is_female'}
    ).cols()

    mean_male_y_sites_called = sex_ht.aggregate(hl.agg.filter(~sex_ht.is_female, hl.agg.group_by(sex_ht.qc_platform, hl.agg.mean(sex_ht.n_y_sites_called))))
    y_call_rate_stats = sex_ht.aggregate(
        hl.agg.filter(
            hl.is_defined(sex_ht.is_female),
            hl.agg.group_by(hl.tuple([sex_ht.qc_platform, sex_ht.is_female]),
                            hl.agg.stats(sex_ht.y_call_rate)
                            )
        )
    )

    no_y_call_rate_platforms = set()
    for platform, sites_called in mean_male_y_sites_called.items():
        if sites_called < min_male_y_sites_called:
            logger.warn(f"Mean number of sites in males on Y chromosome for platform {platform} is < {min_male_y_sites_called} ({sites_called} sites found). Y call rate filter will NOT be applied for samples on platform {platform}.")
            no_y_call_rate_platforms.add(platform)

    sex_ht = sex_ht.annotate_globals(y_call_rate_stats=y_call_rate_stats,
                                     no_y_call_rate_platforms=no_y_call_rate_platforms if no_y_call_rate_platforms else hl.empty_set(hl.tstr))
    y_female_stats = sex_ht.y_call_rate_stats[(sex_ht.qc_platform, True)]
    y_male_stats = sex_ht.y_call_rate_stats[(sex_ht.qc_platform, False)]
    sex_ht = sex_ht.annotate(
        is_female=(
            hl.case()
                .when(sex_ht.no_y_call_rate_platforms.contains(sex_ht.qc_platform), sex_ht.is_female)
                .when(sex_ht.is_female & (sex_ht.y_call_rate < y_female_stats.mean + max_y_call_rate_stdev * y_female_stats.stdev), True)
                .when(~sex_ht.is_female & (sex_ht.y_call_rate > y_male_stats.mean - max_y_call_rate_stdev * y_male_stats.stdev), False)
                .or_missing()
        )
    )
    return(sex_ht)


def main(args):

    output_prefix = args.output_dir.rstrip("/") + "/" + splitext(basename(args.input_mt))[0]

    mt = hl.read_matrix_table(args.input_mt)

    if args.compute_qc_mt:
        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
        qc_mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                            (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                            (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99)).persist()
        pruned_ht = hl.ld_prune(qc_mt.GT, r2=0.1)
        qc_mt = qc_mt.filter_rows(hl.is_defined(pruned_ht[qc_mt.row_key]))
        qc_mt = qc_mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
        qc_mt.write(f'{output_prefix}.qc.mt', overwrite=args.overwrite)

    if args.compute_qc_metrics:
        qc_ht = hl.sample_qc(mt).cols()
        qc_ht.write(f'{output_prefix}.sample_qc.ht', overwrite=args.overwrite)

    if args.compute_callrate_mt:
        logger.info('Preparing data for platform PCA...')
        intervals = hl.import_locus_intervals(evaluation_intervals_path)
        callrate_mt = compute_callrate_mt(mt, intervals)
        callrate_mt.write(f'{output_prefix}.callrate.mt', args.overwrite)

    if args.run_platform_pca:
        logger.info("Running platform PCA")
        callrate_mt = hl.read_matrix_table(f'{output_prefix}.callrate.mt')
        scores_ht  = run_platform_pca(callrate_mt)
        scores_ht.write(f'{output_prefix}.platform_pca_scores.ht', overwrite=args.overwrite)

        scores_ht = hl.read_table(f'{output_prefix}.platform_pca_scores.ht')
        platform_ht = assign_platform_pcs(scores_ht)
        platform_ht.write(f'{output_prefix}.platform_pca_results.ht', overwrite=args.overwrite)

    if args.impute_sex:
        platform_ht = hl.read_table(f'{output_prefix}.platform_pca_results.ht')
        qc_mt = hl.read_matrix_table(f'{output_prefix}.qc.mt')
        sex_ht =  impute_sex(qc_mt, mt, platform_ht, args.male_threshold, args.female_threshold, args.min_male_y_sites_called, args.y_call_rate_stdev)
        sex_ht.write(f'{output_prefix}.sex.ht', overwrite=args.overwrite)

    if args.run_pc_relate:
        logger.info('Running PCA for PC-Relate...')
        qc_mt = hl.read_matrix_table(f'{output_prefix}.qc.mt')
        eig, scores, _ = hl.hwe_normalized_pca(qc_mt.GT, k=10, compute_loadings=False)
        scores.write(f'{output_prefix}.pruned.pca_scores.ht', args.overwrite)

        logger.info('Running PC-Relate...')
        scores = hl.read_table(f'{output_prefix}.pruned.pca_scores.ht')
        # NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes
        relatedness_ht = hl.pc_relate(qc_mt.GT, min_individual_maf=0.05, scores_expr=scores[qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=args.min_kinship, statistics='kin2')
        relatedness_ht.write(f'{output_prefix}.relatedness.ht', args.overwrite)

    if args.filter_related_samples:
        pass


    if args.filter_dups: # Note, this needs a prioritization scheme
        relatedness_ht = hl.read_table(f'{output_prefix}.relatedness.ht')
        logger.info("Getting duplicate samples")
        dups = get_duplicated_samples(relatedness_ht)
        dups_to_keep = [hl.struct(kept=list(d)[0], filtered=list(d)[1:]) for d in dups]  # Dict of kept -> filtered (random selection)
        dups_ht = hl.Table.parallelize(dups_to_keep, key='kept')
        dups_ht.write(f'{output_prefix}.dups.ht', overwrite=args.overwrite)

    if args.infer_families:
        relatedness_ht = hl.read_table(f'{output_prefix}.relatedness.ht')
        relatedness_ht = relatedness_ht.key_by( # Needed for familial inference at this point
            i=relatedness_ht.i.s,
            j=relatedness_ht.j.s
        )
        dups_ht = hl.read_table(f'{output_prefix}.dups.ht')
        dups_to_remove = dups_ht.aggregate(hl.agg.explode(lambda x: hl.agg.collect_as_set(x.s), dups_ht.filtered))
        logger.info(f"Found {len(dups)} duplicate sets; marking {len(dups_to_remove)} duplicates")
        x_ht = hl.read_table(f'{output_prefix}.sex.ht')
        sex = {row.s: row.is_female for row in x_ht.to_pandas().itertuples()}
        raw_ped = infer_families(relatedness_ht, sex, dups_to_remove)
        logger.info(f"Found {len(raw_ped.complete_trios())} complete trios.")
        raw_ped.write(f"{output_prefix}.ped")

    if args.run_pca:
        pass




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--input_mt', help='Input MT', default='gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.vcf.split.mt')
    parser.add_argument('--output_dir', help='Directory to output files to.', default='gs://gnomad/projects/compound_hets/myoseq/sample_qc')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    basic = parser.add_argument_group("Compute basic QC metrics and resources")
    basic.add_argument('--compute_qc_mt', help='Computes QC MT for PCA, sex, relationship inference, etc.', action='store_true')
    basic.add_argument('--compute_qc_metrics', help='Computes hail sample qc metrics', action='store_true')

    platform_pca = parser.add_argument_group("Platform PCA")
    platform_pca.add_argument('--compute_callrate_mt', help='Computes an interval by sample mt of callrate that will be used to compute platform PCs', action='store_true')
    platform_pca.add_argument('--run_platform_pca', help='Runs platform PCA (assumes callrate MT was computed)', action='store_true')

    sex_imputations = parser.add_argument_group("Sex imputation")
    sex_imputations.add_argument('--impute_sex', help='Imputes sex', action='store_true')
    sex_imputations.add_argument('--male_threshold', help='samples with X-chromosome F-stat > male_threshold will be annotated as males', default=0.7, type=float)
    sex_imputations.add_argument('--female_threshold', help='samples with X-chromosome F-stat < female_threshold will be annotated as females', default=0.5, type=float)
    sex_imputations.add_argument('--min_male_y_sites_called', help='Minimum number of sites called on Y for a male (based on X f-stat) on a platform to enable  Y call rate filtering.', default=500, type=int)
    sex_imputations.add_argument('--y_call_rate_stdev', help='Number of standard deviation from the mean of the platform Y call rate to remove a sample based on its Y call rate', default=4, type=int)

    relatedness = parser.add_argument_group("Relatedness")
    relatedness.add_argument('--run_pc_relate', help='Runs PC-relate on all samples', action='store_true')
    relatedness.add_argument('--min_kinship', help='Minimum kinship threshold for emitting a pair of samples in PC relate and filtering related individuals.', default=0.05, type=float)
    relatedness.add_argument('--filter_related_samples', help='Filter related samples (based on the pairs present from the --run_pc_relate, so using the --min_kinship value for that run)', action='store_true')
    relatedness.add_argument('--filter_dups', help='Filter duplicated samples', action='store_true')
    relatedness.add_argument('--infer_families', help='Extracts duplicate samples and infers families samples based on PC-relate results', action='store_true')

    pca = parser.add_argument_group("Population PCA)")
    pca.add_argument('--run_pca', help='Runs PCA on all samples', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
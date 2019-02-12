from gnomad_hail import *
from os.path import splitext, basename
import hail as hl
import numpy as np
import hdbscan
import pickle
from sklearn.ensemble import RandomForestClassifier

# TODO: Move to gnomad_hail common repo
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
    ht = hl.Table.from_pandas(data, key=[*platform_pc_table.key])
    ht = ht.annotate(qc_platform = hl.int32(ht.qc_platform))
    return ht


def compute_callrate_mt(mt: hl.MatrixTable,  intervals: hl.Table) -> hl.MatrixTable:
    callrate_mt = filter_to_autosomes(mt)
    callrate_mt = callrate_mt.annotate_rows(interval=intervals[callrate_mt.locus].target)
    callrate_mt = callrate_mt.filter_rows(hl.is_defined(callrate_mt.interval) & (hl.len(callrate_mt.alleles) == 2))
    callrate_mt = callrate_mt.select_entries(GT=hl.or_missing(hl.is_defined(callrate_mt.GT), hl.struct()))
    callrate_mt = callrate_mt.group_rows_by(callrate_mt.interval).aggregate(callrate=hl.agg.fraction(hl.is_defined(callrate_mt.GT)))
    return  callrate_mt


def run_platform_pca(callrate_mt: hl.MatrixTable) -> Tuple[List[float], hl.Table, hl.Table]:
    callrate_mt = callrate_mt.annotate_entries(callrate=hl.int(callrate_mt.callrate > 0.25))
    # Center until Hail's PCA does it for you
    callrate_mt = callrate_mt.annotate_rows(mean_callrate=hl.agg.mean(callrate_mt.callrate))
    callrate_mt = callrate_mt.annotate_entries(callrate=callrate_mt.callrate - callrate_mt.mean_callrate)
    eigenvalues, scores, loadings = hl.pca(callrate_mt.callrate, compute_loadings=True) # TODO:  Evaluate whether computing loadings is a good / worthy thing
    logger.info('Eigenvalues: {}'.format(eigenvalues))

    return eigenvalues, scores, loadings


def impute_sex(
        x_mt: hl.MatrixTable, # TODO: This feels somewhat unsatisfying to provide two MTs. Maybe just reapply the QC MT filters to both (minus callrate for Y)?
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
        **{f'x_{ann}': x_ht[y.col_key][ann] for ann in x_ht.row_value if ann != 'is_female'}
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
                                     no_y_call_rate_platforms=no_y_call_rate_platforms if no_y_call_rate_platforms else hl.empty_set(hl.tint32))
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
    sex_ht = sex_ht.annotate_globals(
        impute_sex_params=hl.struct(
            male_min_f_stat=male_min_f_stat,
            female_max_f_stat=female_max_f_stat,
            min_male_y_sites_called=min_male_y_sites_called,
            max_y_call_rate_stdev=max_y_call_rate_stdev
        )
    )
    return(sex_ht)


def get_related_samples_to_drop(
        relatedness_ht: hl.Table,
        min_filtering_kinship: float,
        rank_func: Callable[..., Tuple[hl.Table, Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression]]],
        rank_func_args: Optional[List]
) -> hl.Table:
    """
    Use the maximal independence function in Hail to intelligently prune clusters of related individuals, removing
    less desirable samples while maximizing the number of unrelated individuals kept in the sample set.

    :param relatedness_ht:
    :param min_filtering_kinship:
    :param rank_func:
    :param rank_func_args:
    :return:
    """
    # Define maximal independent set, using rank list
    related_pairs = relatedness_ht.filter(relatedness_ht.kin > min_filtering_kinship).key_by().select('i', 'j')
    n_related_samples = related_pairs.annotate(ij=[related_pairs.i, related_pairs.j]).explode('ij').key_by('ij').distinct().count()
    logger.info('{} samples with at least 2nd-degree relatedness found in callset'.format(n_related_samples))

    if rank_func_args is None:
        rank_func_args = []
    related_pairs, tie_breaker = rank_func(relatedness_ht, *rank_func_args)

    related_samples_to_drop_ranked = hl.maximal_independent_set(related_pairs.i, related_pairs.j,
                                                                keep=False, tie_breaker=tie_breaker)
    # FIXME: This is really quite awful, in particular there is no guarantee that rank_func will preserve all fields in i...
    related_samples_to_drop_ranked = related_samples_to_drop_ranked.key_by()
    related_samples_to_drop_ranked = related_samples_to_drop_ranked.select(**related_samples_to_drop_ranked.node)
    return related_samples_to_drop_ranked.key_by(*relatedness_ht.i)


def filter_dups(
        relatedness_ht: hl.Table,
        dup_ranking_func: Callable[..., Tuple[hl.Table, hl.expr.Expression]],
        dup_ranking_func_args: Optional[List]
):
    """
    Creates a HT with duplicated samples sets.
    Each row is indexed by the sample that is kept and also contains the set of duplicate samples that should be filtered.

    dup_ranking_func is a function to decide which duplicate to keep.
    It should take a table keyed by samples key and any number of additional arguments given through `dup_ranking_func_args`
    It should return a tuple containing:
        - the input table with any modifications (e.g. annotations) needed for ranking
        - A hail expression of a type that can be sorted giving the corresponding rank (where smaller is better)

    :param relatedness_ht: Input relatedness HT
    :param dup_ranking_func: Ranking function to pick amongst duplicates.
    :param dup_ranking_func_args: Optional additional arguments for `dup_ranking_func`
    :return:
    """
    logger.info("Getting duplicate samples")
    dups = get_duplicated_samples(relatedness_ht)
    logger.info(f"Found {len(dups)} duplicate sets.")
    dups_ht = hl.Table.parallelize([hl.struct(dup_set=i, dups=dups[i]) for i in range(0, len(dups))])
    dups_ht = dups_ht.explode(dups_ht.dups, name='_dup')
    dups_ht = dups_ht.key_by(**dups_ht._dup)
    if dup_ranking_func_args is None:
        dup_ranking_func_args  = []
    dups_ht, rank_expr = dup_ranking_func(dups_ht, *dup_ranking_func_args)
    dups_cols = hl.bind(
        lambda x: hl.struct(
            kept=x[0],
            filtered=x[1:]
        ),
        hl.sorted(hl.agg.collect(hl.tuple([dups_ht._dup, rank_expr])), key=lambda x: x[1]).map(lambda x: x[0])
    )
    dups_ht = dups_ht.group_by(dups_ht.dup_set).aggregate(
        **dups_cols
    )

    dups_ht = dups_ht.key_by(**{f'{x}_kept': expr for x, expr in dups_ht.kept.items()}).drop('kept')
    return dups_ht

def filter_rows_for_qc(mt: hl.MatrixTable, min_af: float = 0.001, min_callrate: float = 0.99):
    return mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                           (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > min_af) &
                           (hl.agg.fraction(hl.is_defined(mt.GT)) > min_callrate)).persist()


def get_qc_mt(mt: hl.MatrixTable, min_af: float = 0.001, min_callrate: float = 0.99, ld_r2: float = 0.1) -> hl.MatrixTable:
    qc_mt = filter_rows_for_qc(mt, min_af, min_callrate)
    pruned_ht = hl.ld_prune(qc_mt.GT, r2=ld_r2)
    qc_mt = qc_mt.filter_rows(hl.is_defined(pruned_ht[qc_mt.row_key]))
    qc_mt = qc_mt.annotate_globals(
        qc_mt_params=hl.struct(
            min_af=min_af,
            min_callrate=min_callrate,
            ld_r2=ld_r2
        )
    )
    return qc_mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))


def get_ped(relatedness_ht: hl.Table, dups_ht: hl.Table, sex_ht: hl.Table) -> hl.Pedigree:
    relatedness_ht = relatedness_ht.key_by(  # Needed for familial inference at this point -- should be generalized
        i=relatedness_ht.i.s,
        j=relatedness_ht.j.s
    )
    dups_to_remove = dups_ht.aggregate(hl.agg.explode(lambda x: hl.agg.collect_as_set(x.s), dups_ht.filtered))
    logger.info(f"Removing {len(dups_to_remove)} duplicates from family creation.")
    sex = {row.s: row.is_female for row in sex_ht.to_pandas().itertuples()}
    ped = infer_families(relatedness_ht, sex, dups_to_remove)
    logger.info(f"Found {len(ped.complete_trios())} complete trios.")
    return ped


def run_pca_with_relateds(mt: hl.MatrixTable, related_samples_to_drop: Optional[hl.Table], n_pcs: int = 10):
    if related_samples_to_drop:
        unrelated_mt = mt.filter_cols(hl.is_missing(related_samples_to_drop[mt.col_key]))

    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(unrelated_mt.GT, k=n_pcs, compute_loadings=True)
    pca_af_ht = unrelated_mt.annotate_rows(pca_af=hl.agg.mean(unrelated_mt.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)  # TODO: Evaluate if needed to write results at this point if relateds or not

    if not related_samples_to_drop:
        return pca_evals, pca_scores, pca_loadings
    else:
        pca_loadings = pca_loadings.persist()
        pca_scores = pca_scores.persist()
        related_mt = mt.filter_cols(hl.is_defined(related_samples_to_drop[mt.col_key]))
        related_scores = pc_project(related_mt, pca_loadings)
        pca_scores = pca_scores.union(related_scores)
        return pca_evals, pca_scores, pca_loadings


def assign_pops_from_pc(joint_samples_ht: hl.table, pca_scores_ht: hl.Table, known_col: str, min_assignment_prob: float) -> Tuple[hl.Table, RandomForestClassifier]: # TODO: Should we just integrate this in utils.assign_population_pcs ?
    pca_scores = pca_scores_ht[joint_samples_ht.key]
    pops_pd = joint_samples_ht.annotate(
        **{f'PC{i}': pca_scores.scores[i] for i in range(0, args.n_kgp_pcs)}
    ).to_pandas()

    pops_pd, pops_rf_model = assign_population_pcs(pops_pd, [f'PC{i}' for i in range(0, args.n_kgp_pcs)], known_col=known_col, min_prob=min_assignment_prob)

    pops_ht = hl.Table.from_pandas(pops_pd, key=list(joint_samples_ht.key))
    pops_ht.annotate_globals(
        assign_pops_from_pc_params=hl.struct(
            min_assignment_prob=min_assignment_prob
        )
    )
    return pops_ht, pops_rf_model


def compute_stratified_metrics_filter(ht: hl.Table, qc_metrics: List[str], strata: List[str] = None) -> hl.Table:
    """
    Compute median, MAD, and upper and lower thresholds for each metric used in pop- and platform-specific outlier filtering

    :param MatrixTable ht: HT containing relevant sample QC metric annotations
    :param list qc_metrics: list of metrics for which to compute the critical values for filtering outliers
    :param list of str strata: List of annotations used for stratification. These metrics should be discrete types!
    :return: Table grouped by pop and platform, with upper and lower threshold values computed for each sample QC metric
    :rtype: Table
    """

    def make_pop_filters_expr(ht: hl.Table, qc_metrics: List[str]) -> hl.expr.SetExpression:
        return hl.set(hl.filter(lambda x: hl.is_defined(x),
                                [hl.or_missing(ht[f'fail_{metric}'], metric) for metric in qc_metrics]))

    ht = ht.select(*strata, **ht.sample_qc.select(*qc_metrics)).key_by('s').persist()

    def get_metric_expr(ht, metric):
        return hl.bind(
            lambda x: x.annotate(
                upper=x.median + 4 * x.mad,
                lower=x.median - 4 * x.mad
            ),
            hl.bind(
                lambda elements, median: hl.struct(
                    median=median,
                    mad=1.4826 * hl.median(hl.abs(elements - median))
                ),
                *hl.bind(
                    lambda x: hl.tuple([x, hl.median(x)]),
                    hl.agg.collect(ht[metric])
                )
            )
        )

    agg_expr = hl.struct(**{metric: get_metric_expr(ht, metric) for metric in qc_metrics})
    if strata:
        ht = ht.annotate_globals(metrics_stats=ht.aggregate(hl.agg.group_by(hl.tuple([ht[x] for x in strata]), agg_expr)))
    else:
        ht = ht.annotate_globals(metrics_stats={(): ht.aggregate(agg_expr)})

    strata_exp = hl.tuple([ht[x] for x in strata]) if strata else hl.tuple([])

    fail_exprs = {
        f'fail_{metric}':
            (ht[metric] >= ht.metrics_stats[strata_exp][metric].upper) |
            (ht[metric] <= ht.metrics_stats[strata_exp][metric].lower)
        for metric in qc_metrics}
    ht = ht.transmute(**fail_exprs)
    pop_platform_filters = make_pop_filters_expr(ht, qc_metrics)
    return ht.annotate(pop_platform_filters=pop_platform_filters)


# MYOSEQ-specific methods

def rank_related_samples(
        related_pairs: hl.Table,
        meta_ht: hl.Table,
        qc_ht: hl.Table,
        fam_ht: hl.Table
) -> Tuple[hl.Table, Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression]]:

    # Load families and identify parents from cases as they will be thrown away anyways
    fam_ht = fam_ht.transmute(trio=[
        hl.struct(s=fam_ht.id, is_parent=False),
        hl.struct(s=fam_ht.pat_id, is_parent=True),
        hl.struct(s=fam_ht.mat_id, is_parent=True)
    ])
    fam_ht = fam_ht.explode(fam_ht.trio)
    fam_ht = fam_ht.key_by(s=fam_ht.trio.s)
    case_parents = fam_ht.filter(meta_ht[fam_ht.key].is_case & fam_ht.trio.is_parent)

    def annotate_related_pairs(related_pairs: hl.Table, index_col: str) -> hl.Table:
        related_pairs = related_pairs.key_by(**related_pairs[index_col])
        related_pairs = related_pairs.filter(hl.is_missing(case_parents[related_pairs.key]))
        return related_pairs.annotate(
            **{
                index_col: related_pairs[index_col].annotate(
                    case_rank=hl.or_else(hl.int(meta_ht[related_pairs.key].is_case), -1),
                    dp_mean=hl.or_else(qc_ht[related_pairs.key].sample_qc.dp_stats.mean, -1.0)
                )
            }
        ).key_by()

    related_pairs = annotate_related_pairs(related_pairs, "i")
    related_pairs = annotate_related_pairs(related_pairs, "j")

    def tie_breaker(l, r):
        return (
            hl.case()
            .when(l.case_rank != r.case_rank, r.case_rank - l.case_rank) # smaller is better
            .default(l.dp_mean - r.dp_mean) # larger is better
        )

    return related_pairs, tie_breaker


def rank_dup_samples(dups_ht: hl.Table, qc_ht: hl.Table) -> Tuple[hl.Table, hl.expr.Expression]:
    dups_ht = dups_ht.annotate(rank=-1*qc_ht[dups_ht.key].sample_qc.dp_stats.mean)
    return dups_ht, dups_ht.rank


def main(args):

    output_prefix = args.output_dir.rstrip("/") + "/" + splitext(basename(args.input_mt))[0]

    if args.compute_qc_mt:
        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
        qc_mt = get_qc_mt(hl.read_matrix_table(args.input_mt))
        qc_mt.write(f'{output_prefix}.qc.mt', overwrite=args.overwrite)

    if args.compute_qc_metrics:
        qc_ht = hl.sample_qc(hl.read_matrix_table(args.input_mt)).cols()
        qc_ht.write(f'{output_prefix}.sample_qc.ht', overwrite=args.overwrite)

    if args.compute_callrate_mt:
        logger.info('Preparing data for platform PCA...')
        intervals = hl.import_locus_intervals(evaluation_intervals_path)
        callrate_mt = compute_callrate_mt(hl.read_matrix_table(args.input_mt), intervals)
        callrate_mt.write(f'{output_prefix}.callrate.mt', args.overwrite)

    if args.run_platform_pca:
        logger.info("Running platform PCA")
        callrate_mt = hl.read_matrix_table(f'{output_prefix}.callrate.mt')
        eigenvalues, scores_ht, loadings_ht  = run_platform_pca(callrate_mt)
        scores_ht.write(f'{output_prefix}.platform_pca_scores.ht', overwrite=args.overwrite)
        loadings_ht.write(f'{output_prefix}.platform_pca_loadings.ht', overwrite=args.overwrite)

        scores_ht = hl.read_table(f'{output_prefix}.platform_pca_scores.ht')
        platform_ht = assign_platform_pcs(scores_ht)
        platform_ht.write(f'{output_prefix}.platform_pca_results.ht', overwrite=args.overwrite)

    if args.impute_sex:
        sex_ht = impute_sex(
            hl.read_matrix_table(f'{output_prefix}.qc.mt'),
            hl.read_matrix_table(args.input_mt),
            hl.read_table(f'{output_prefix}.platform_pca_results.ht'),
            args.male_threshold,
            args.female_threshold,
            args.min_male_y_sites_called,
            args.y_call_rate_stdev
        )
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
                                      block_size=4096, min_kinship=args.min_emission_kinship, statistics='all')
        relatedness_ht.write(f'{output_prefix}.relatedness.ht', args.overwrite)

    if args.filter_dups:
        dups_ht = filter_dups(
            hl.read_table(f'{output_prefix}.relatedness.ht'),
            rank_dup_samples,
            [hl.read_table(f'{output_prefix}.sample_qc.ht')]
        )
        dups_ht.write(f'{output_prefix}.dups.ht', overwrite=args.overwrite)

    if args.infer_families:
        ped = get_ped(
            hl.read_table(f'{output_prefix}.relatedness.ht'),
            hl.read_table(f'{output_prefix}.dups.ht'),
            hl.read_table(f'{output_prefix}.sex.ht')
        )
        ped.write(f"{output_prefix}.ped")

    if args.filter_related_samples:
        related_samples_to_drop = get_related_samples_to_drop(
            hl.read_table(f'{output_prefix}.relatedness.ht'),
            args.min_filtering_kinship,
            rank_related_samples,
            [
                hl.read_table(args.meta),
                hl.read_table(f'{output_prefix}.sample_qc.ht'),
                hl.import_fam(f"{output_prefix}.ped", delimiter="\t")
            ]
        )
        related_samples_to_drop.write(f'{output_prefix}.related_samples_to_drop.ht', overwrite=args.overwrite)

    if args.run_pca:
        qc_mt = filter_to_autosomes(hl.read_matrix_table(f'{output_prefix}.qc.mt'))
        related_samples_to_drop = hl.read_table(f'{output_prefix}.related_samples_to_drop.ht')
        pca_evals, pca_scores, pca_loadings = run_pca_with_relateds(qc_mt, related_samples_to_drop, args.n_pcs)
        pca_loadings.write(f'{output_prefix}.pca_loadings.ht', args.overwrite)
        pca_scores.write(f'{output_prefix}.pca_scores.ht', args.overwrite)

    if args.assign_pops_kgp:
        logger.info("Joining data with 1000 Genomes")
        qc_mt = filter_to_autosomes(hl.read_matrix_table(f'{output_prefix}.qc.mt')).select_rows().select_entries("GT")
        qc_mt = qc_mt.select_cols(known_pop=hl.null(hl.tstr), known_subpop=hl.null(hl.tstr))
        qc_mt = qc_mt.key_cols_by(_kgp=False,  *qc_mt.col_key)

        kgp_mt = filter_to_autosomes(hl.read_matrix_table(kgp_phase3_genotypes_mt_path())).select_rows()
        kgp_mt = kgp_mt.select_cols(
            known_pop=kgp_mt.super_pops.get(kgp_mt.population, "oth").lower(),
            known_subpop=kgp_mt.population.lower()
        )
        kgp_mt = kgp_mt.filter_rows(hl.is_defined(qc_mt.rows()[kgp_mt.row_key]))
        kgp_mt = filter_rows_for_qc(kgp_mt)
        kgp_mt = kgp_mt.key_cols_by(_kgp=True, *kgp_mt.col_key)

        joint_mt = qc_mt.union_cols(kgp_mt)
        joint_mt.write(f'{output_prefix}.union_kgp.qc.mt', overwrite=args.overwrite)

        logger.info("Computing PCA on data with 1000 Genomes")
        joint_mt = hl.read_matrix_table(f'{output_prefix}.union_kgp.qc.mt')
        related_samples_to_drop = hl.read_table(f'{output_prefix}.related_samples_to_drop.ht')
        related_samples_to_drop = related_samples_to_drop.key_by(_kgp=False,  *related_samples_to_drop.key)
        pca_evals, pca_scores, pca_loadings = run_pca_with_relateds(joint_mt, related_samples_to_drop, args.n_kgp_pcs)
        pca_loadings.write(f'{output_prefix}.union_kgp.pca_loadings.ht', args.overwrite)
        pca_scores.write(f'{output_prefix}.union_kgp.pca_scores.ht', args.overwrite)

        logger.info("Assigning populations based on 1000 Genomes labels")
        joint_cols = hl.read_matrix_table(f'{output_prefix}.union_kgp.qc.mt').cols()
        pca_scores = hl.read_table(f'{output_prefix}.union_kgp.pca_scores.ht')
        pops_ht, pops_rf_model = assign_pops_from_pc(joint_cols, pca_scores, 'known_pop', args.min_pop_prob)

        pops_ht.write(f'{output_prefix}.union_kgp.pops.ht', args.overwrite)

        with hl.hadoop_open(f'{output_prefix}.union_kgp.pops.rf_fit.pkl', 'wb') as out:
            pickle.dump(pops_rf_model, out)

    if args.assign_subpops_kgp:
        joint_mt = hl.read_matrix_table(f'{output_prefix}.union_kgp.qc.mt')
        pops_ht = hl.read_table(f'{output_prefix}.union_kgp.pops.ht')
        related_samples_to_drop = hl.read_table(f'{output_prefix}.related_samples_to_drop.ht')
        related_samples_to_drop = related_samples_to_drop.key_by(_kgp=False, *related_samples_to_drop.key)

        pops_for_subpop = pops_ht.aggregate(hl.agg.group_by(pops_ht.pop, hl.agg.count_where(~pops_ht._kgp)))
        pops_for_subpop = [pop for pop, n in pops_for_subpop.items() if n >= args.min_samples_for_subpop and pop is not None and pop != 'oth']
        logger.info(f"Assigning subpopulations based on 1000 Genomes for: {','.join(pops_for_subpop)}")

        for pop in pops_for_subpop:
            joint_mt = joint_mt.filter_cols(pops_ht[joint_mt.col_key].pop == pop)
            joint_samples_ht = joint_mt.cols().persist()
            joint_mt = filter_rows_for_qc(joint_mt).persist()
            pca_evals, pca_scores, pca_loadings = run_pca_with_relateds(joint_mt, related_samples_to_drop, args.n_kgp_pcs)
            pca_loadings.write(f'{output_prefix}.union_kgp.{pop}.pca_loadings.ht', args.overwrite)
            pca_scores.write(f'{output_prefix}.union_kgp.{pop}.pca_scores.ht', args.overwrite)

            subpops_ht, subpops_rf_model = assign_pops_from_pc(joint_samples_ht, pca_scores, 'known_subpop', args.min_pop_prob)
            subpops_ht.write(f'{output_prefix}.union_kgp.subpops.{pop}.ht', args.overwrite)

            with hl.hadoop_open(f'{output_prefix}.union_kgp.subpops.{pop}.rf_fit.pkl', 'wb') as out:
                pickle.dump(subpops_rf_model, out)

    if args.apply_stratified_filters:
        qc_ht = hl.read_table(f'{output_prefix}.sample_qc.ht')
        pops_ht = hl.read_table(f'{output_prefix}.union_kgp.pops.ht') # TODO: This should be customizable
        pops_ht = pops_ht.filter(~pops_ht._kgp)
        pops_ht = pops_ht.key_by(*[k for k in pops_ht.key if k != "_kgp"])
        platform_ht = hl.read_table(f'{output_prefix}.platform_pca_results.ht')
        qc_ht = qc_ht.annotate(
            qc_pop=pops_ht[qc_ht.key].pop,
            qc_platform=platform_ht[qc_ht.key].qc_platform
        )
        stratified_metrics_ht = compute_stratified_metrics_filter(qc_ht, args.filtering_qc_metrics.split(","), ['qc_pop', 'qc_platform'])
        stratified_metrics_ht.write(f'{output_prefix}.stratified_metrics_filters.ht', overwrite=args.overwrite)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--input_mt', help='Input MT', default='gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.mt')
    parser.add_argument('--meta', help='Location of the metadata HT', default="gs://gnomad/projects/compound_hets/myoseq/LGMD_VCF_Detail_Metrics_V2_01302019.ht")
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
    relatedness.add_argument('--run_pc_relate', help='Runs PC-relate on all samples. NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes', action='store_true')
    relatedness.add_argument('--min_emission_kinship', help='Minimum kinship threshold for emitting a pair of samples in PC relate and filtering related individuals.', default=0.05, type=float)
    relatedness.add_argument('--min_filtering_kinship', help='Minimum kinship threshold for filtering a pair of samples in PC relate and filtering related individuals. (Default = 0.08838835; 2nd degree relatives)', default= 0.08838835, type=float)
    relatedness.add_argument('--filter_dups', help='Filter duplicated samples', action='store_true')
    relatedness.add_argument('--infer_families', help='Extracts duplicate samples and infers families samples based on PC-relate results', action='store_true')
    relatedness.add_argument('--filter_related_samples', help='Filter related samples (based on the pairs present from the --run_pc_relate and using the --min_filtering_kinship value for that run)', action='store_true')

    pca = parser.add_argument_group("Population PCA")
    pca.add_argument('--run_pca', help='Runs PCA on all samples', action='store_true')
    pca.add_argument('--n_pcs', help='Number of PCs to compute (default: 20)', default=20, type=int)

    pop = parser.add_argument_group("Population assignment")
    pop.add_argument('--assign_pops_kgp', help='Runs a combined PCA with 1000 Genomes and assigns populations based on 1000 Genomes pops.', action='store_true')
    pop.add_argument('--n_kgp_pcs', help='Number of PCs to compute when joined with 1000 Genomes (default: 10)', default=10, type=int)
    pop.add_argument('--min_pop_prob', help='Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as "oth" (default: 0.6)', default=0.6, type=float) # TODO: Evaluate whether this is sensible. Also, should we consider the difference bewteen the two most likely pops instead?
    pop.add_argument('--assign_subpops_kgp', help='Runs a combined PCA with 1000 Genomes and assigns populations based on 1000 Genomes pops.', action='store_true')
    pop.add_argument('--min_samples_for_subpop', help='Minimum number of samples in a global population to run the subpopulation PCA / assignment (default: 500)', default=500, type=int)

    qc_metrics_filtering = parser.add_argument_group("Stratified (per population/platform) QC metrics filtering")
    qc_metrics_filtering.add_argument('--apply_stratified_filters', help="Compute per pop, per platform filtering and create a table with these annotations.", action='store_true')
    qc_metrics_filtering.add_argument('--filtering_qc_metrics', help="List of QC metrics for filtering.", default=",".join(['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']))

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
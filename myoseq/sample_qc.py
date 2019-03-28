from gnomad_hail import *
from os.path import splitext, basename
import hail as hl
from gnomad_hail.utils.sample_qc import *
import pickle


output_prefix = ""

def rank_related_samples(
        relatedness_ht: hl.Table,
        meta_ht: hl.Table,
        sample_qc_ht: hl.Table,
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
                    dp_mean=hl.or_else(sample_qc_ht[related_pairs.key].sample_qc.dp_stats.mean, -1.0)
                )
            }
        ).key_by()

    relatedness_ht = annotate_related_pairs(relatedness_ht, "i")
    relatedness_ht = annotate_related_pairs(relatedness_ht, "j")

    def tie_breaker(l, r):
        return (
            hl.case()
                .when(l.case_rank != r.case_rank, r.case_rank - l.case_rank)  # smaller  is better
                .default(l.dp_mean - r.dp_mean)  # larger is better
        )

    return relatedness_ht, tie_breaker


def assign_and_write_subpops(
        mt: hl.MatrixTable,
        related_samples_to_drop_ht: hl.Table,
        min_samples_for_subpop: int,
        n_pcs: int,
        min_pop_prob: float,
        overwrite: bool = False,
        pop_ann: str = 'pop',
        subpop_ann: str = 'subpop',
        files_prefix: str = '',
        include_in_pop_count: hl.expr.BooleanExpression = hl.expr.bool(True)
) -> None:
    logger.info("Assigning subpopulations{}".format(
        f'({files_prefix.rstrip("_")})' if files_prefix else ''
    ))

    pops_for_subpop = mt.aggregate_cols(hl.agg.group_by(mt[pop_ann], hl.agg.count_where(include_in_pop_count)))
    pops_for_subpop = [pop for pop, n in pops_for_subpop.items() if n >= min_samples_for_subpop and pop is not None and pop != 'oth']
    logger.info(f"Assigning subpopulations for: {','.join(pops_for_subpop)}")

    for pop in pops_for_subpop:
        logger.info(f"Running subpop pcs for {pop}.")
        mt = mt.filter_cols(mt[pop_ann] == pop)
        mt = filter_rows_for_qc(mt)
        pca_evals, subpop_pca_scores_ht, subpop_pca_loadings_ht = run_pca_with_relateds(
            mt,
            related_samples_to_drop_ht,
            n_pcs
        )
        subpop_pca_loadings_ht.write(path(f'{files_prefix}{pop}_pca_loadings.ht'), overwrite)
        subpop_pca_scores_ht.write(path(f'{files_prefix}{pop}_pca_scores.ht'), overwrite)

        subpop_pca_scores_ht = hl.read_table(path(f'{files_prefix}{pop}_pca_scores.ht'))
        subpop_pca_scores_ht = subpop_pca_scores_ht.annotate(
            **mt.cols()[subpop_pca_scores_ht.key].select(subpop_ann)
        )
        subpop_ht, subpop_rf_model = assign_population_pcs(
            subpop_pca_scores_ht,
            pc_cols=subpop_pca_scores_ht.scores[:n_pcs],
            known_col=subpop_ann,
            min_prob=min_pop_prob,
            output_col="subpop"
        )
        subpop_ht.write(path(f'{files_prefix}subpop_{pop}.ht'), overwrite)

        with hl.hadoop_open(path(f'{files_prefix}subpop_{pop}_rf_model.pkl'), 'wb') as out:
            pickle.dump(subpop_rf_model, out)


def path(file: str) -> str:
    return f'{output_prefix}.{file}'


def get_platform_specific_intervals(platform_pc_loadings_ht: hl.Table, threshold: float) -> List[hl.Interval]:
    """
    This takes the platform PC loadings and returns a list of intervals where the sum of the loadings above the given threshold.
    The experimental / untested idea behind this, is that those intervals may be problematic on some platforms.

    :param Table platform_pc_loadings_ht: Platform PCA loadings indexed by interval
    :param float threshold: Minimal threshold
    :param str intervals_path: Path to the intervals file to use (default: b37 exome calling intervals)
    :return: List of intervals with PC loadings above the given threshold
    :rtype: list of Interval
    """
    platform_specific_intervals = platform_pc_loadings_ht.filter(hl.sum(hl.abs(platform_pc_loadings_ht.loadings))>=threshold)
    return platform_specific_intervals.interval.collect()


def main(args):
    global output_prefix
    output_prefix = args.output_dir.rstrip("/") + "/" + splitext(basename(args.input_mt))[0]

    if args.compute_qc_mt:
        qc_mt = get_qc_mt(hl.read_matrix_table(args.input_mt))
        qc_mt = qc_mt.repartition(n_partitions=200)
        qc_mt.write(path('qc.mt'), overwrite=args.overwrite)

    if args.compute_qc_metrics:
        logger.info("Computing sample QC")
        mt = filter_to_autosomes(hl.read_matrix_table(args.input_mt))
        strats = {
            'bi_allelic': bi_allelic_expr(mt),
            'multi_allelic': ~bi_allelic_expr(mt)
        }
        for strat, filter_expr in strats.items():
            strat_sample_qc_ht =  hl.sample_qc(mt.filter_rows(filter_expr)).cols()
            strat_sample_qc_ht.write(path(f'{strat}_sample_qc.ht'), overwrite=args.overwrite)
        strat_hts = [hl.read_table(path(f'{strat}_sample_qc.ht')) for strat in strats]
        sample_qc_ht = strat_hts.pop()
        sample_qc_ht = sample_qc_ht.select(
            sample_qc=merge_sample_qc_expr([sample_qc_ht.sample_qc] + [strat_hts[i][sample_qc_ht.key].sample_qc for i in range(0,len(strat_hts))])
        )
        sample_qc_ht.write(path('sample_qc.ht'), overwrite=args.overwrite)

    if args.compute_callrate_mt:
        callrate_mt = compute_callrate_mt(hl.read_matrix_table(args.input_mt))
        callrate_mt.write(path('callrate.mt'), args.overwrite)

    if args.run_platform_pca:
        eigenvalues, scores_ht, loadings_ht = run_platform_pca(hl.read_matrix_table(path('callrate.mt')))
        scores_ht.write(path('platform_pca_scores.ht'), overwrite=args.overwrite)
        loadings_ht.write(path('platform_pca_loadings.ht'), overwrite=args.overwrite)

    if args.assign_platforms:
        platform_ht = assign_platform_from_pcs(
            hl.read_table(path('platform_pca_scores.ht')), 
            hdbscan_min_cluster_size=args.hdbscan_min_cluster_size, 
            hdbscan_min_samples=args.hdbscan_min_samples
        )
        platform_ht.write(f'{output_prefix}.platform_pca_results.ht', overwrite=args.overwrite)

    if args.impute_sex:
        sex_ht = infer_sex(
            hl.read_matrix_table(path('qc.mt')),
            hl.read_matrix_table(args.input_mt),
            hl.read_table(path('platform_pca_results.ht')),
            args.male_threshold,
            args.female_threshold,
            args.min_male_y_sites_called,
            args.max_y_female_call_rate,
            args.min_y_male_call_rate
        )
        sex_ht.write(path('sex.ht'), overwrite=args.overwrite)

    if args.run_pc_relate:
        logger.info('Running PCA for PC-Relate')
        qc_mt = hl.read_matrix_table(path('qc.mt')).unfilter_entries()
        eig, scores, _ = hl.hwe_normalized_pca(qc_mt.GT, k=10, compute_loadings=False)
        scores.write(path('pruned.pca_scores.ht'), args.overwrite)

        logger.info('Running PC-Relate')
        logger.warn("PC-relate requires SSDs and doesn't work with preemptible workers!")
        scores = hl.read_table(path('pruned.pca_scores.ht'))
        relatedness_ht = hl.pc_relate(qc_mt.GT, min_individual_maf=0.05, scores_expr=scores[qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=args.min_emission_kinship, statistics='all')
        relatedness_ht.write(path('relatedness.ht'), args.overwrite)

    if args.filter_dups:
        logger.info("Filtering duplicate samples")
        sample_qc_ht = hl.read_table(path('sample_qc.ht'))
        samples_rankings_ht = sample_qc_ht.select(rank = -1 * sample_qc_ht.sample_qc.dp_stats.mean)
        dups_ht = filter_duplicate_samples(
            hl.read_table(path('relatedness.ht')),
            samples_rankings_ht
        )
        dups_ht.write(path('duplicates.ht'), overwrite=args.overwrite)

    if args.infer_families:
        logger.info("Inferring families")
        duplicates_ht =hl.read_table(path('duplicates.ht'))
        dups_to_remove = duplicates_ht.aggregate(hl.agg.explode(lambda x: hl.agg.collect_as_set(x.s), duplicates_ht.filtered))
        ped = infer_families(
            hl.read_table(path('relatedness.ht')),
            hl.read_table(path('sex.ht')),
            dups_to_remove
        )
        ped.write(path('pedigree.ped'))

    if args.filter_related_samples:
        logger.info("Filtering related samples")
        related_pairs_ht, related_pairs_tie_breaker = rank_related_samples(
            hl.read_table(path('relatedness.ht')),
            hl.read_table(args.meta),
            hl.read_table(path('sample_qc.ht')),
            hl.import_fam(path('pedigree.ped'), delimiter="\t")
        )

        related_samples_to_drop_ht = hl.maximal_independent_set(related_pairs_ht.i, related_pairs_ht.j,
                                                             keep=False, tie_breaker=related_pairs_tie_breaker)
        related_samples_to_drop_ht = related_samples_to_drop_ht.key_by()
        related_samples_to_drop_ht = related_samples_to_drop_ht.select(**related_samples_to_drop_ht.node)
        related_samples_to_drop_ht = related_samples_to_drop_ht.key_by('s')
        related_samples_to_drop_ht.write(path('related_samples_to_drop.ht'), overwrite=args.overwrite)

    if args.run_pca:
        logger.info("Running population PCA")
        pca_evals, subpop_pca_scores_ht, subpop_pca_loadings_ht = run_pca_with_relateds(
            hl.read_matrix_table(path('qc.mt')), 
            hl.read_table(path('related_samples_to_drop.ht')), 
            args.n_pcs
        )
        subpop_pca_loadings_ht.write(path('pop_pca_loadings.ht'), args.overwrite)
        subpop_pca_scores_ht.write(path('pop_pca_scores.ht'), args.overwrite)

    if args.assign_pops:
        logger.info("Assigning global population labels")
        subpop_pca_scores_ht = hl.read_table(path("pop_pca_scores.ht"))
        gnomad_meta_ht = get_gnomad_meta('exomes').select("pop")[subpop_pca_scores_ht.key]
        subpop_pca_scores_ht = subpop_pca_scores_ht.annotate(
            known_pop=hl.or_missing(gnomad_meta_ht.pop != "oth", gnomad_meta_ht.pop)
        )
        pop_ht, pops_rf_model = assign_population_pcs(
            subpop_pca_scores_ht,
            pc_cols=subpop_pca_scores_ht.scores[:args.n_pcs],
            known_col='known_pop',
            min_prob=args.min_pop_prob
        )

        pop_ht.write(path('pop.ht'), args.overwrite)
        with hl.hadoop_open(path('pop_rf_model.pkl'), 'wb') as out:
            pickle.dump(pops_rf_model, out)

    if args.assign_subpops:
        qc_mt = hl.read_matrix_table(path('qc.mt'))
        pop_ht = hl.read_table(path('pop.ht'))
        meta_ht = hl.read_table(args.meta)[qc_mt.col_key]
        qc_mt = qc_mt.annotate_cols(
            pop=pop_ht[qc_mt.col_key].pop,
            is_case=meta_ht.is_case,
            country=meta_ht.country
        )

        platform_specific_intervals = get_platform_specific_intervals(hl.read_table(path('platform_pca_loadings.ht')), threshold=0.01)
        logger.info(f'Excluding {len(platform_specific_intervals)} platform-specific intervals for subpop PCA.')
        qc_mt = hl.filter_intervals(qc_mt, platform_specific_intervals, keep=False)

        assign_and_write_subpops(
            qc_mt ,
            hl.read_table(path('related_samples_to_drop.ht')),
            min_samples_for_subpop=args.min_samples_for_subpop,
            n_pcs=args.n_pcs,
            min_pop_prob=args.min_pop_prob,
            overwrite=args.overwrite,
            pop_ann='pop',
            subpop_ann='country',
            include_in_pop_count=qc_mt.is_case
        )

    if args.run_kgp_pca:
        logger.info("Joining data with 1000 Genomes")
        qc_mt = hl.read_matrix_table(path('qc.mt')).select_rows().select_entries("GT")
        qc_mt = qc_mt.select_cols(known_pop=hl.null(hl.tstr), known_subpop=hl.null(hl.tstr))
        qc_mt = qc_mt.key_cols_by(_kgp=False, *qc_mt.col_key)

        kgp_mt = hl.read_matrix_table(kgp_phase3_genotypes_mt_path()).select_rows()
        kgp_mt = kgp_mt.select_cols(
            known_pop=kgp_mt.super_pops.get(kgp_mt.population, "oth").lower(),
            known_subpop=kgp_mt.population.lower()
        )
        kgp_mt = kgp_mt.filter_rows(hl.is_defined(qc_mt.rows()[kgp_mt.row_key]))
        kgp_mt = filter_rows_for_qc(kgp_mt)
        kgp_mt = kgp_mt.key_cols_by(_kgp=True, *kgp_mt.col_key)

        union_kgp_qc_mt = qc_mt.union_cols(kgp_mt)
        union_kgp_qc_mt.write(path('union_kgp_qc.mt'), overwrite=args.overwrite)

        logger.info("Computing PCA on data with 1000 Genomes")
        union_kgp_qc_mt = hl.read_matrix_table(path('union_kgp_qc.mt'))
        related_samples_to_drop_ht = hl.read_table(path('related_samples_to_drop.ht'))
        related_samples_to_drop_ht = related_samples_to_drop_ht.key_by(_kgp=False, *related_samples_to_drop_ht.key)
        pca_evals, union_kgp_pca_scores_ht, union_kgp_pca_loadings_ht = run_pca_with_relateds(union_kgp_qc_mt, related_samples_to_drop_ht, args.n_kgp_pcs)
        union_kgp_pca_loadings_ht.write(path('union_kgp_pca_loadings.ht'), args.overwrite)
        union_kgp_pca_scores_ht.write(path('union_kgp_pca_scores.ht'), args.overwrite)

    if args.assign_pops_kgp:
        logger.info("Assigning populations based on 1000 Genomes labels")
        union_kgp_qc_mt = hl.read_matrix_table(path('union_kgp_qc.mt'))
        union_kgp_pca_scores_ht = hl.read_table(path('union_kgp_pca_scores.ht'))
        union_kgp_pca_scores_ht = union_kgp_pca_scores_ht.annotate(
            known_pop=union_kgp_qc_mt[union_kgp_pca_scores_ht.key].known_pop
        )
        union_kgp_pop_ht, union_kgp_pop_rf_model = assign_population_pcs(
            union_kgp_pca_scores_ht,
            pc_cols=union_kgp_pca_scores_ht.scores[:args.n_kgp_pcs],
            known_col='known_pop',
            min_prob=args.min_kgp_pop_prob
        )

        union_kgp_pop_ht.write(path('union_kgp_pop.ht'), args.overwrite)

        with hl.hadoop_open(path('union_kgp_pop_rf_model.pkl'), 'wb') as out:
            pickle.dump(union_kgp_pop_rf_model, out)

    if args.assign_subpops_kgp:
        union_kgp_qc_mt = hl.read_matrix_table(path('union_kgp_qc.mt'))
        meta_ht = hl.read_table(args.meta)
        union_kgp_pop_ht = hl.read_table(path('union_kgp_pop.ht'))
        union_kgp_qc_mt = union_kgp_qc_mt.annotate_cols(
            is_case=meta_ht[union_kgp_qc_mt.col_key].is_case,
            pop=union_kgp_pop_ht[union_kgp_qc_mt.col_key].pop
        )

        platform_specific_intervals = get_platform_specific_intervals(hl.read_table(path('platform_pca_loadings.ht')))
        logger.info(f'Excluding {len(platform_specific_intervals)} platform-specific intervals for subpop PCA.')
        union_kgp_qc_mt = hl.filter_intervals(union_kgp_qc_mt, platform_specific_intervals, keep=False)

        related_samples_to_drop_ht = hl.read_table(path('related_samples_to_drop.ht'))
        related_samples_to_drop_ht = related_samples_to_drop_ht.key_by(_kgp=False, *related_samples_to_drop_ht.key)

        assign_and_write_subpops(
            union_kgp_qc_mt,
            related_samples_to_drop_ht,
            min_samples_for_subpop=args.min_samples_for_subpop,
            n_pcs=args.n_kgp_pcs,
            min_pop_prob=args.min_kgp_pop_prob,
            overwrite=args.overwrite,
            pop_ann='pop',
            subpop_ann='known_subpop',
            include_in_pop_count=union_kgp_qc_mt.is_case,
            files_prefix='union_kgp_'
        )

    if args.apply_stratified_filters:
        logger.info("Computing stratified QC")
        for variant_class_prefix in ['', 'bi_allelic_', 'multi_allelic_']:
            sample_qc_ht = hl.read_table(path(f'{variant_class_prefix}sample_qc.ht'))
            pop_ht = hl.read_table(path('pops.ht'))
            platform_ht = hl.read_table(path('platform_pca_results.ht'))
            sample_qc_ht = sample_qc_ht.annotate(
                qc_pop=pop_ht[sample_qc_ht.key].pop,
                qc_platform=platform_ht[sample_qc_ht.key].qc_platform
            )
            stratified_metrics_ht = compute_stratified_metrics_filter(
                sample_qc_ht,
                args.filtering_qc_metrics.split(","),
                ['qc_pop', 'qc_platform']
            )
            stratified_metrics_ht.write(path(f'{variant_class_prefix}stratified_metrics_filters.ht'), overwrite=args.overwrite)

    if args.write_full_meta:
        logger.info("Writing metadata table")

        # List all tables to join with the base meta
        meta_annotation_hts = [
            hl.read_table(path('platform_pca_results.ht')).rename({'scores': 'platform_pc_scores'}),
            hl.read_table(path('sex.ht')),
            flatten_duplicate_samples_ht(hl.read_table(path('duplicates.ht'))),
            hl.read_table(path('related_samples_to_drop.ht')).select(related_filtered=True),
            hl.read_table(path('pca_scores.ht')).rename({'scores': 'pop_pc_scores'}),
            hl.read_table(path('pops.ht')).select('pop'),
            hl.read_table(path('nfe.pca_scores.ht')).rename({'scores': 'nfe_pc_scores'}),
            hl.read_table(path('subpops.nfe.ht')).select('subpop')
        ]

        # union_kgp_pops_ht = hl.read_table(path('union_kgp_pops.ht'))
        # union_kgp_pops_ht = union_kgp_pops_ht.filter(~union_kgp_pops_ht._kgp).key_by('s')
        # union_kgp_pops_ht = union_kgp_pops_ht.select(kgp_pop=union_kgp_pops_ht.pop)
        # meta_annotation_hts.append(union_kgp_pops_ht)
        #
        # union_kgp_pca_scores_ht = hl.read_table(path('union_kgp_pca_scores.ht')).rename({'scores': 'kgp_pop_pc_scores'})
        # union_kgp_pca_scores_ht = union_kgp_pca_scores_ht.filter(~union_kgp_pca_scores_ht._kgp).key_by('s')
        # meta_annotation_hts.append(union_kgp_pca_scores_ht)

        gnomad_meta_ht = get_gnomad_meta('exomes')
        gnomad_meta_ht = gnomad_meta_ht.select(gnomad_pop=gnomad_meta_ht.pop, gnomad_subpop=gnomad_meta_ht.subpop)
        meta_annotation_hts.append(gnomad_meta_ht)

        for variant_class_prefix in ['', 'bi_allelic_', 'multi_allelic_']:
            sample_qc_ht = hl.read_table(path(f'{variant_class_prefix}sample_qc.ht'))
            stratified_metrics_filters_ht = hl.read_table(path(f'{variant_class_prefix}stratified_metrics_filters.ht'))
            if variant_class_prefix:
                sample_qc_ht = sample_qc_ht.rename({'sample_qc': f'{variant_class_prefix}sample_qc'})
                stratified_metrics_filters_ht = stratified_metrics_filters_ht.rename(
                    {f: f'{variant_class_prefix}{f}' for f in list(stratified_metrics_filters_ht.globals) + list(stratified_metrics_filters_ht.row_value)}
                )
            meta_annotation_hts.extend([sample_qc_ht, stratified_metrics_filters_ht])

        meta_ht = hl.read_table(args.meta)
        meta_ht = meta_ht.annotate_globals(
            **{name: expr for ann_ht in meta_annotation_hts for name, expr in ann_ht.index_globals().items()}
        )

        meta_ht = meta_ht.annotate(
            **{name: expr for ann_ht in meta_annotation_hts for name, expr in ann_ht[meta_ht.key].items()}
        )

        filtering_col_prefix = '' if args.filtering_variant_class == 'all' else args.filtering_variant_class + "_"
        meta_ht = meta_ht.annotate_globals(filtering_variant_class=args.filtering_variant_class)
        meta_ht =  meta_ht.annotate(
            sample_filters=add_filters_expr(
                filters={
                    "ambiguous sex": hl.is_missing(meta_ht.is_female),
                    'call_rate': meta_ht.sample_qc.call_rate < args.min_call_rate,
                    'duplicate': hl.is_defined(meta_ht.dup_filtered) & meta_ht.dup_filtered,
                    'related': meta_ht.related_filtered
                },
                current_filters=meta_ht[f'{filtering_col_prefix}pop_platform_filters']
            )
        )

        meta_ht.write(path('full_meta.ht'), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--input_mt', help='Input MT', default='gs://gnomad/projects/compound_hets/myoseq/MacArthur_LGMD_Callset_Jan2019.mt')
    parser.add_argument('--meta', help='Location of the metadata HT', default="gs://gnomad/projects/compound_hets/myoseq/original_meta_files/myoseq_meta.ht")
    parser.add_argument('--output_dir', help='Directory to output files to.', default='gs://gnomad/projects/compound_hets/myoseq/sample_qc')
    parser.add_argument('--write_full_meta', help='Write metadata file with all information.', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    hard_filters = parser.add_argument_group("Hard filters")
    hard_filters.add_argument('--min_call_rate', help='Minimum call rate to keep sample (default: 0.95)', default=0.95, type=float)

    basic = parser.add_argument_group("Compute basic QC metrics and resources")
    basic.add_argument('--compute_qc_mt', help='Computes QC MT for PCA, sex, relationship inference, etc.', action='store_true')
    basic.add_argument('--compute_qc_metrics', help='Computes hail sample qc metrics', action='store_true')

    platform_pca = parser.add_argument_group("Platform PCA")
    platform_pca.add_argument('--compute_callrate_mt', help='Computes an interval by sample mt of callrate that will be used to compute platform PCs', action='store_true')
    platform_pca.add_argument('--run_platform_pca', help='Runs platform PCA (assumes callrate MT was computed)', action='store_true')
    platform_pca.add_argument('--assign_platforms', help='Assigns platforms based on callrate PCA results using HDBSCAN', action='store_true')
    platform_pca.add_argument('--hdbscan_min_samples', help='Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.', type=int, required=False)
    platform_pca.add_argument('--hdbscan_min_cluster_size', help='Minimum cluster size parameter for HDBSCAN.', type=int, default=50)

    sex_imputations = parser.add_argument_group("Sex imputation")
    sex_imputations.add_argument('--impute_sex', help='Imputes sex', action='store_true')
    sex_imputations.add_argument('--male_threshold', help='samples with X-chromosome F-stat > male_threshold will be annotated as males', default=0.7, type=float)
    sex_imputations.add_argument('--female_threshold', help='samples with X-chromosome F-stat < female_threshold will be annotated as females', default=0.5, type=float)
    sex_imputations.add_argument('--min_male_y_sites_called', help='Minimum number of sites called on Y for a male (based on X f-stat) on a platform to enable  Y call rate filtering.', default=500, type=int)
    sex_imputations.add_argument('--max_y_female_call_rate', help='Maximum normalized Y-chromosome callrate for females.', default=0.15, type=float)
    sex_imputations.add_argument('--min_y_male_call_rate', help='Minimum normalized Y-chromosome callrate for males.', default=0.8, type=float)

    relatedness = parser.add_argument_group("Relatedness")
    relatedness.add_argument('--run_pc_relate', help='Runs PC-relate on all samples. NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes', action='store_true')
    relatedness.add_argument('--min_emission_kinship', help='Minimum kinship threshold for emitting a pair of samples in PC relate and filtering related individuals.', default=0.05, type=float)
    relatedness.add_argument('--min_filtering_kinship', help='Minimum kinship threshold for filtering a pair of samples in PC relate and filtering related individuals. (Default = 0.08838835; 2nd degree relatives)', default=0.08838835, type=float)
    relatedness.add_argument('--filter_dups', help='Filter duplicated samples', action='store_true')
    relatedness.add_argument('--infer_families', help='Extracts duplicate samples and infers families samples based on PC-relate results', action='store_true')
    relatedness.add_argument('--filter_related_samples', help='Filter related samples (based on the pairs present from the --run_pc_relate and using the --min_filtering_kinship value for that run)', action='store_true')

    pca = parser.add_argument_group("Population PCA")
    pca.add_argument('--run_pca', help='Runs PCA on all samples', action='store_true')
    pca.add_argument('--n_pcs', help='Number of PCs to compute (default: 10)', default=10, type=int)

    pop = parser.add_argument_group("Population assignment")
    pop.add_argument('--assign_pops', help='Assign pops based on gnomAD samples known pops.', action='store_true')
    pop.add_argument('--assign_subpops', help='Assign pops based on gnomAD samples known pops.', action='store_true')
    pop.add_argument('--min_samples_for_subpop', help='Minimum number of samples in a global population to run the subpopulation PCA / assignment (default: 500)', default=500, type=int)
    pop.add_argument('--min_pop_prob', help='Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as "oth" (default: 0.9)', default=0.9,
                     type=float)  # TODO: Evaluate whether this is sensible. Also, should we consider the difference bewteen the two most likely pops instead?
    pop.add_argument('--run_kgp_pca', help='Runs a combined PCA with 1000 Genomes.', action='store_true')
    pop.add_argument('--n_kgp_pcs', help='Number of PCs to compute when joined with 1000 Genomes (default: 10)', default=10, type=int)
    pop.add_argument('--assign_pops_kgp', help='Assigns populations based on 1000 Genomes pca/pops.', action='store_true')
    pop.add_argument('--assign_pops_both', help='Assigns populations based on 1000 Genomes pca, using 1000 Genomes AND gnomAD pops.', action='store_true')
    pop.add_argument('--min_kgp_pop_prob', help='Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as "oth" (default: 0.6)', default=0.6,
                     type=float)  # TODO: Evaluate whether this is sensible. Also, should we consider the difference bewteen the two most likely pops instead?
    pop.add_argument('--assign_subpops_kgp', help='Runs a combined PCA with 1000 Genomes and assigns populations based on 1000 Genomes pops.', action='store_true')

    qc_metrics_filtering = parser.add_argument_group("Stratified (per population/platform) QC metrics filtering")
    qc_metrics_filtering.add_argument('--apply_stratified_filters', help="Compute per pop, per platform filtering and create a table with these annotations.", action='store_true')
    qc_metrics_filtering.add_argument('--filtering_qc_metrics', help="List of QC metrics for filtering.", default=",".join(['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']))
    qc_metrics_filtering.add_argument('--filtering_variant_class', help="Variant class to use for filtering; one of: 'all', 'bi_allelic', 'multi_allelic'", default='bi_allelic')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

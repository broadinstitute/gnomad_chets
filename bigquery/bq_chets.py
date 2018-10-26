import argparse
import sys
from bq_utils import *

def main(args):

    data_type = 'exomes' if args.exomes else 'genomes'
    use_cache = not args.no_query_cache

    client = bigquery.Client()
    dataset_ref = client.dataset(args.dataset)

    tables_prefix = f'broad-mpg-gnomad.{args.dataset}.{data_type}'
    test_prefix = '_test' if args.test else ''

    logger.info("Creating list of existing variant-pairs per transcript")
    sql = f"""
        SELECT DISTINCT v1.v as v1, v2.v as v2
        FROM 
        (   SELECT gt.s_id, tx.v, transcript_consequences_transcript_id as transcript_id FROM `{tables_prefix}_genotypes` as gt
            join `{tables_prefix}_transcripts` as tx on tx.v = gt.v
            join `{tables_prefix}_meta` as meta on meta.s_id = gt.s_id 
            where gt.is_adj and meta.release
        ) as v1
        join (
            SELECT gt.s_id, tx.v, transcript_consequences_transcript_id as transcript_id FROM `{tables_prefix}_genotypes` as gt
            join `{tables_prefix}_transcripts` as tx on tx.v = gt.v   
            join `{tables_prefix}_meta` as meta on meta.s_id = gt.s_id 
            where gt.is_adj and meta.release
        ) as v2
        on v1.transcript_id = v2.transcript_id and v1.s_id = v2.s_id
        where v1.v < v2.v
    """
    create_table(client, dataset_ref.table(f'{data_type}_existing_variant_pairs'), sql, args.overwrite, use_cache)

    if args.test:
        logger.info("Creating test existing variants table")
        sql = f"""
                SELECT * from `{tables_prefix}_existing_variant_pairs`
                LIMIT 100    
            """
        create_table(client, dataset_ref.table(f'{data_type}_test_existing_variant_pairs'), sql, overwrite=args.overwrite)

        logger.info("Creating test genotypes table")
        sql = f"""
                SELECT * FROM `{tables_prefix}_genotypes`
                WHERE v IN (
                    SELECT v1 FROM `{tables_prefix}_test_existing_variant_pairs` 
                    UNION ALL (
                        SELECT v2 FROM `{tables_prefix}_test_existing_variant_pairs`
                    )
                )    
            """
        create_table(client, dataset_ref.table(f'{data_type}_test_genotypes'), sql, overwrite=args.overwrite)

    # Works but produces 60TB table!!
    # for i in ['1', '2']:
    #     logger.info(f"Creating v{i} table")
    #     sql = f"""
    #         SELECT v1, v2, gt.pop, gt.s_id, is_het
    #         FROM
    #         (
    #             SELECT v, meta.pop, meta.s_id,
    #             CASE WHEN (is_adj is null or not is_adj) THEN null ELSE is_het END as is_het
    #             FROM `{tables_prefix}{test_prefix}_genotypes` as gt
    #             JOIN `{tables_prefix}_meta` as meta
    #             ON gt.s_id = meta.s_id
    #             WHERE meta.release
    #         ) as gt
    #         JOIN `{tables_prefix}{test_prefix}_existing_variant_pairs` as vp
    #         ON vp.v{i} = gt.v
    #     """
    #     create_table(client, dataset_ref.table(f'{data_type}{test_prefix}_v{i}_gt'), sql, overwrite=args.overwrite or args.test)
    #
    # logger.info("Creating variant-pairs table with genotypes")
    # sql = f"""
    #     SELECT
    #     g1.v1,
    #     g1.v2,
    #     pop,
    #     g1.s_id,
    #     g1.is_het as het1,
    #     g2.is_het as het2
    #     FROM {data_type}{test_prefix}_v1_gt as g1
    #     JOIN {data_type}{test_prefix}_v2_gt as g2
    #     ON g1.v1 = g2.v1 and g1.v2 = g2.v2 and g1.s_id = g2.s_id
    # """
    # create_table(client, dataset_ref.table(f'{data_type}{test_prefix}_vp_gt'), sql, overwrite=args.overwrite or args.test)

    logger.info("Grouping variant-pairs genotypes by v1, v2, pop")
    sql = f"""
        SELECT 
        v1,v2,pop,
        count(case when het1 and het2 then 1 end) as n11,
        count(case when het1 and het2 is not null and not het2 then 1 end) as n12,
        count(case when het1 is not null and not het1 and het2 then 1 end) as n21,
        count(case when het1 is not null and het2 is not null and not het1 and not het2 then 1 end) as n22,
        count(case when het1 is null and het2 then 1 end) as nx1,
        count(case when het1 is null and het2 is not null and not het2 then 1 end) as nx2,
        count(case when het1 and het2 is null then 1 end) as n1x,
        count(case when het1 is not null and not het1 and het2 is null then 1 end) as n2x
        FROM (
            SELECT
            g1.v1,
            g1.v2,
            pop,
            g1.s_id,
            g1.is_het as het1,
            g2.is_het as het2
            FROM (
                SELECT v1, v2, gt.pop, gt.s_id, is_het
                FROM
                (
                    SELECT v, meta.pop, meta.s_id,
                    CASE WHEN (is_adj is null or not is_adj) THEN null ELSE is_het END as is_het
                    FROM `{tables_prefix}{test_prefix}_genotypes` as gt
                    JOIN `{tables_prefix}_meta` as meta
                    ON gt.s_id = meta.s_id
                    WHERE meta.release
                ) as gt
                JOIN `{tables_prefix}{test_prefix}_existing_variant_pairs` as vp
                ON vp.v1 = gt.v
            ) as g1
            JOIN (
                SELECT v1, v2, gt.s_id, is_het
                FROM
                (
                    SELECT v, meta.s_id,
                    CASE WHEN (is_adj is null or not is_adj) THEN null ELSE is_het END as is_het
                    FROM `{tables_prefix}{test_prefix}_genotypes` as gt
                    JOIN `{tables_prefix}_meta` as meta
                    ON gt.s_id = meta.s_id
                    WHERE meta.release
                ) as gt
                JOIN `{tables_prefix}{test_prefix}_existing_variant_pairs` as vp
                ON vp.v2 = gt.v
            ) as g2
            ON g1.v1 = g2.v1 and g1.v2 = g2.v2 and g1.s_id = g2.s_id 
        )
        where het1 is not null or het2 is not null
        group by v1, v2, pop
    """
    create_table(client, dataset_ref.table(f'{data_type}{test_prefix}_grp_vp_gt'), sql, overwrite=args.overwrite or args.test)

    logger.info("Creating variants frequency view")
    sql = f"""
        SELECT
        idx as v,
        chrom,
        pos,
        ref,
        alt,
        freq.element.pop,
        freq.element.ac,
        freq.element.an,
        freq.element.af,
        freq.element.hom
        FROM `{tables_prefix}_variants` as var
        CROSS JOIN UNNEST(var.freq.list) as freq
    """
    create_table(client, dataset_ref.table(f'{data_type}_freq'), sql, overwrite=args.overwrite, view=True)

    logger.info("Creating variant-pairs table")
    sql = f"""
        SELECT v1, v2, vp.pop, 
        f1.an - f2.ac - f2.hom - n11 - n21 - nx1 - f2.hom - nx2 - n22 - n12 as n00,
        f2.ac - f2.hom - n11 - n21 - nx1 as n01,
        f2.hom - nx2 - n22 - n21 as n02,
        f1.ac - f1.hom - n11 - n12 - n1x as n10,
        n11, 
        n12, 
        f1.hom - n21 - n22 - n2x as n20,
        n21, 
        n22, 
        nx1, 
        nx2, 
        n1x, 
        n2x 
        FROM `{tables_prefix}{test_prefix}_grp_vp_gt` as vp
        join `{tables_prefix}_freq` as f1 on 
        v1 = f1.v and f1.pop = vp.pop
        join `{tables_prefix}_freq` as f2 on 
        v2 = f2.v and f2.pop = vp.pop
        order by v1,v2,pop
    """
    create_table(client, dataset_ref.table(f'{data_type}{test_prefix}_variant_pairs'), sql, overwrite=args.overwrite or args.test)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes (One of --exomes or --genomes is required)', action='store_true')
    parser.add_argument('--genomes', help='Run on genomes (One of --exomes or --genomes is required)', action='store_true')
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--test', help='Runs the script on a subset of 100 variant-pairs.', action='store_true')
    parser.add_argument('--overwrite', help='If set, will overwrite all existing tables.', action='store_true')
    parser.add_argument('--no_query_cache', help='If set, query cache is disabled.', action='store_true')

    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Exactly one of --exomes or --genomes is required.')

    main(args)

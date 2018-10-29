import argparse
from bq_utils import *
from bq_table_descriptions import *


def get_data_view(client: bigquery.Client, data_type: str, dataset: bigquery.DatasetReference) -> str:

    #Exclude data_type from meta, so we can use * in main  query
    meta_table = client.get_table(dataset.table(f'{data_type}_meta'))
    variants_table = client.get_table(dataset.table(f'{data_type}_variants'))
    genotypes_table = client.get_table(dataset.table(f'{data_type}_genotypes'))
    meta_cols = [f.name for f in meta_table.schema if f.name != 'data_type']
    genotypes_cols = [f.name for f in genotypes_table.schema]
    variants_cols = [f.name for f in variants_table.schema]
    first_cols = [f"'{data_type}' as data_type", "chrom", "pos", "ref", "alt"]

    return f"""
    
    SELECT {",".join(first_cols)},
           {",".join([f for f in genotypes_cols if f != 'v'])}, 
           {",".join([f for f in variants_cols if f not in first_cols])},
           {",".join([f for f in meta_cols if f not in genotypes_cols])} 
            
           FROM `{dataset.project}.{dataset.dataset_id}.{data_type}_variants` as v
    LEFT JOIN (
        SELECT {",".join([f"gt.{f}" for f in genotypes_cols])}, 
               {",".join([f"meta.{f}" for f in meta_cols if f not in genotypes_cols])} 
        FROM `{dataset.project}.{dataset.dataset_id}.{data_type}_genotypes` as gt
        JOIN (
            SELECT {",".join([f for f in meta_cols])} FROM `{dataset.project}.{dataset.dataset_id}.{data_type}_meta`
             ) as meta
        ON gt.s = meta.s
        ) as g
    ON v.idx = g.v
    
    """


def main(args):

    client = bigquery.Client()
    dataset = client.dataset(args.dataset)


    logger.info("Creating all_meta view")
    create_table(client,
                 dataset.table('all_meta'),
                 sql=create_union_query(client,
                                                  [client.get_table(dataset.table('genomes_meta')),
                                                   client.get_table(dataset.table('exomes_meta'))],
                                                  True
                                                  ),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_meta_table_desc()
                 )

    logger.info("Creating exomes view")
    create_table(client,
                 dataset.table('exomes'),
                 sql=get_data_view(client, 'exomes', dataset),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc('exomes')
                 )

    logger.info("Creating genomes view")
    create_table(client,
                 dataset.table('genomes'),
                 sql=get_data_view(client, 'genomes', dataset),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc('genomes')
                 )

    logger.info("Creating all view")
    create_table(client,
                 dataset.table('all'),
                 sql=create_union_query(client, [dataset.table('exomes'), dataset.table('genomes')]),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc()
                 )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--overwrite', help='If set, will overwrite all existing tables.', action='store_true')
    args = parser.parse_args()

    main(args)
import argparse
import sys
from google.cloud import bigquery


def load_parquet_bq(client: bigquery.Client, dataset: str, table: str, parquet_files: str, write_disposition: str):

    dataset_ref = client.dataset(dataset)
    table_ref = dataset_ref.table(table)

    job_config = bigquery.LoadJobConfig()
    job_config.source_format = bigquery.SourceFormat.PARQUET
    job_config.write_disposition = write_disposition
    load_job = client.load_table_from_uri(parquet_files,
                                          table_ref,
                                          job_config=job_config)

    load_job.result()
    print('Successfully loaded {} rows in table {}.{}.'.format(client.get_table(table_ref).num_rows, dataset, table))


def main(args):

    client = bigquery.Client()

    data_types = []
    if args.exomes:
        data_types.append('test_exomes') # TODO: Remove test when ready!
    if args.genomes:
        data_types.append('test_genomes') # TODO: Remove test when ready!

    for data_type in data_types:

        if args.load_metadata:
            load_parquet_bq(client, 'gnomad', f'{data_type}_samples_metadata', args.load_metadata, args.write_disposition)

        if args.load_genotypes:
            load_parquet_bq(client, 'gnomad', f'{data_type}_genotypes', args.load_genotypes, args.load_genotypes)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--load_metadata', help='Loads samples metadata from the specified file(s).')
    parser.add_argument('--load_genotypes', help='Loads genotypes from the specified file(s).')
    parser.add_argument('--create_variants_pairs', help='Creates the variant pairs table from the genotypes table.', action='store_true')
    parser.add_argument('--write_disposition', help='One of WRITE_EMPTY (error if table exists, default), WRITE_APPEND (append to table if exists) or WRITE_TRUNCATE (replace existing table)', default='WRITE_EMPTY')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: Exactly one of --exomes or --genomes must be specified.')


    main(args)

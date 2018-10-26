import argparse
from google.cloud import bigquery


def main(args):

    client = bigquery.Client()

    dataset_ref = client.dataset(args.dataset)
    table_ref = dataset_ref.table(args.table)

    job_config = bigquery.LoadJobConfig()
    job_config.source_format = bigquery.SourceFormat.PARQUET
    job_config.write_disposition = args.write_disposition
    load_job = client.load_table_from_uri(args.parquet_files,
                                          table_ref,
                                          job_config=job_config)

    load_job.result()
    print('Successfully loaded {} rows in table {}.{}.'.format(client.get_table(table_ref).num_rows, args.dataset, args.table))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--parquet_files', help='PARQUET file(s) location on GCS. Supports * patterns.')
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--table', help='Table name.')
    parser.add_argument('--write_disposition', help='One of WRITE_EMPTY (error if table exists, default), WRITE_APPEND (append to table if exists) or WRITE_TRUNCATE (replace existing table)', default='WRITE_EMPTY')
    args = parser.parse_args()

    main(args)

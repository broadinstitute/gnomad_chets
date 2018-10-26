import argparse
from bq_utils import *


def main(args):

    client = bigquery.Client()
    dataset = client.dataset(args.dataset)

    if args.meta_view:
        create_table(client,
                     dataset.table('all_meta'),
                     sql=create_union_query(client,
                                                      [client.get_table(dataset.table('genomes_meta')),
                                                       client.get_table(dataset.table('exomes_meta'))],
                                                      True
                                                      ),
                     view=True,
                     overwrite=args.overwrite
                     )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--meta_view', help='Creates unified exomes/genomes metadata view', action='store_true')
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--overwrite', help='If set, will overwrite all existing tables.', action='store_true')
    args = parser.parse_args()

    main(args)
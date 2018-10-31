import argparse
import bq_import
from bq_table_descriptions import *
from bq_utils import logger
import sys


def main(args):

    def import_data(data_type: str, data: str, description: str):
        parquet_files = f"{args.input_dir}/gnomad_{data_type}_{data}.parquet/*.parquet"
        logger.info(f"Importing {data_type} {data} from {parquet_files}")
        parser = bq_import.get_parser()
        imp_args = parser.parse_args(['--dataset', args.dataset,
                        '--parquet_files', parquet_files,
                        '--table', f'{data_type}_{data}',
                        '--write_disposition', 'WRITE_TRUNCATE' if args.overwrite else 'WRITE_EMPTY',
                        '--description', description])
        bq_import.main(imp_args)

    data_types = (['exomes'] if args.exomes else []) + (['genomes'] if args.genomes else [])
    for data_type in data_types:
        if args.import_meta:
            import_data(data_type, 'meta', get_meta_table_desc(data_type))

        if args.import_variants:
            import_data(data_type, 'variants', get_variants_table_desc(data_type))

        if args.import_genotypes:
            import_data(data_type, 'genotypes', get_genotypes_table_desc(data_type))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Will import exomes data. At least one of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Will import genomes data. At least one of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--import_meta', help='Imports samples metadata.', action='store_true')
    parser.add_argument('--import_variants', help='Imports variants.', action='store_true')
    parser.add_argument('--import_genotypes', help='Imports genotypes.', action='store_true')
    parser.add_argument('--input_dir', help='Input root directory (assumes data was produced with `bq_export.py`). default: gs://gnomad-tmp/bq', default='gs://gnomad-tmp/bq')
    parser.add_argument('--overwrite', help='If set, overwrites existing table.', action='store_true')
    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit("At least one of --exomes or --genomes needs to be specified.")

    main(args)

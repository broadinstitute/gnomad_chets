from gnomad_hail import *
from os.path  import splitext
import hail as hl

def normalize_names(ht: hl.Table) -> hl.Table:
    return ht.rename({x: x.replace(" ", "_").lower() for x in ht.row})


# Imports Sam's LGMD case meta
def import_lgmd_case_meta(input_dir: str) -> hl.Table:
    ht = hl.import_table(f'{input_dir}/LGMD_cases-CMG_LGMD_Cases.tsv', impute=True)
    # Lowercase, remove spaces
    ht = normalize_names(ht)
    ht = ht.drop("")
    ht = ht.key_by(s=ht.individual_id)
    return ht


# Import VCF metrics and project info from Lauren's meta
def import_vcf_meta(input_dir: str) -> hl.Table:
    ht = hl.import_table(f'{input_dir}/LGMD_VCF_Detail_Metrics_V2_01302019.txt')
    ht = normalize_names(ht)

    # Make case/control a bool
    ht = ht.transmute(
        is_case=hl.case()
            .when(ht['case/control_status'] == "Case", True)
            .when(ht['case/control_status'] == "Control", False)
            .or_missing()
    )
    ht = ht.key_by(s=ht.sample_alias)
    return ht


# Import seqr individual metadata
def import_seqr_individual_meta(input_dir: str) -> hl.Table:
    files = hl.utils.hadoop_ls(f'{input_dir}/seqr_exports/')
    def import_file(path):
        ht = hl.import_table(path, quote='"')
        ht = normalize_names(ht)
        return ht.key_by('individual_id')
    hts = [import_file(file['path']) for file in files]
    ht = hts[0].union(*hts[1:])
    return ht.distinct()


def import_myoseq_center_meta(input_dir: str) -> hl.Table:
    ht = hl.import_table(f'{input_dir}/MYOSEQ_Centers.tsv', impute=True)
    ht = normalize_names(ht)
    ht = ht.transmute(is_european=hl.cond(ht['european?'] == "Yes", True, False))
    ht = ht.drop("google_map")
    ht = ht.rename({'': 'latitude', '_1': 'longitude'})
    return ht.key_by('z')


def main(args):

    ht = import_vcf_meta(args.input_dir)
    lgmd_case_meta = import_lgmd_case_meta(args.input_dir)
    seqr_meta = import_seqr_individual_meta(args.input_dir)

    # Give priority to lgmd case meta fields (should be equivalent anyways)
    fields_in_both=set(lgmd_case_meta.row_value).intersection(set(seqr_meta.row_value))
    logger.info("Found the following fields in both seqr and lgmd case meta, using lgmd case meta value: {}".format(fields_in_both))
    seqr_meta = seqr_meta.drop(*fields_in_both)

    ht = ht.annotate(
        **lgmd_case_meta[ht.key],
        **seqr_meta[ht.key]
    )

    centers_meta_ht = import_myoseq_center_meta(args.input_dir)
    ht = ht.key_by(center_id=ht.s[:4])
    ht = ht.annotate(**centers_meta_ht[ht.key])
    ht = ht.key_by('s')

    # Check that affected status matches is_case
    n_diff = ht.aggregate(hl.agg.count_where((ht.affected_status == "Affected") != ht.is_case))
    if n_diff:
        logger.warn(f"Found {n_diff} individuals with different case/control status based on seqr vs project")
    else:
        ht = ht.drop("affected_status")

    ht.write(args.output, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', help='Location of the base location of TSV metadata files to import', default="gs://gnomad/projects/compound_hets/myoseq/original_meta_files")
    parser.add_argument('--output', help='Output HT file location.', default="gs://gnomad/projects/compound_hets/myoseq/original_meta_files/myoseq_meta.ht")
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
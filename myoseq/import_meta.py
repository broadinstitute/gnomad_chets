from gnomad_hail import *
from os.path  import splitext
import hail as hl


def main(args):
    output_file = splitext(args.meta)[0] + ".ht"

    ht = hl.import_table(args.meta, impute=True)

    # Lowercase, remove spaces
    ht = ht.rename({x: x.replace(" ", "_").lower() for x in ht.row})

    # Make case/control a bool
    ht = ht.transmute(
        is_case=hl.case()
            .when(ht['case/control_status'] == "Case", True)
            .when(ht['case/control_status'] == "Control", False)
            .or_missing()
    )

    ht = ht.key_by(s=ht.sample_alias)
    ht.write(output_file, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', help='Location of the TSV metadata file to import', default="gs://gnomad/projects/compound_hets/myoseq/LGMD_VCF_Detail_Metrics_V2_01302019.txt")
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
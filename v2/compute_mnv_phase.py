import hail as hl
from gnomad.utils.liftover import get_liftover_genome
import argparse
from gnomad_chets.compute_phase import compute_phase, flatten_phased_ht, liftover_expr


def main(args):
    # Load MNV data
    mnv_ht = hl.import_table(args.input_csv, impute=True)
    mnv_ht = mnv_ht.key_by(
        locus1=hl.locus(mnv_ht.snv1.split("-")[0], hl.int(mnv_ht.snv1.split("-")[1]), reference_genome='GRCh37'),
        alleles1=[mnv_ht.snv1.split("-")[2], mnv_ht.snv1.split("-")[3]],
        locus2=hl.locus(mnv_ht.snv2.split("-")[0], hl.int(mnv_ht.snv2.split("-")[1]), reference_genome='GRCh37'),
        alleles2=[mnv_ht.snv2.split("-")[2], mnv_ht.snv2.split("-")[3]],
    )

    # Add phase information
    phased_ht = compute_phase(mnv_ht)
    phased_ht = phased_ht.annotate(
        **mnv_ht[phased_ht.key]
    )

    # Flatten and export
    flatten_phased_ht(phased_ht).export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_csv', help='CMG file containing variants to phase.',
                          default='gs://gcp-public-data--gnomad/release/2.1/mnv/gnomad_mnv_coding.tsv')
    parser.add_argument('--out', help='Output TSV file',
                          default='gs://gnomad/projects/compound_hets/gnomad_2.1_mnv_coding_phased.tsv')

    args = parser.parse_args()

    main(args)

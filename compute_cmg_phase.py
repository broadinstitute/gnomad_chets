import hail as hl
from gnomad.utils.liftover import get_liftover_genome
from gnomad.utils.slack import try_slack
import argparse
from compute_phase import compute_phase, flatten_phased_ht, liftover_expr


def load_cmg(cmg_csv: str) -> hl.Table:
    cmg_ht = hl.import_table(cmg_csv, impute=True, delimiter=",", quote='"')

    cmg_ht = cmg_ht.transmute(
        locus1_b38=hl.locus("chr" + hl.str(cmg_ht.chrom_1), cmg_ht.pos_1, reference_genome='GRCh38'),
        alleles1_b38=[cmg_ht.ref_1, cmg_ht.alt_1],
        locus2_b38=hl.locus("chr" + hl.str(cmg_ht.chrom_2), cmg_ht.pos_2, reference_genome='GRCh38'),
        alleles2_b38=[cmg_ht.ref_2, cmg_ht.alt_2]
    )

    liftover_references = get_liftover_genome(cmg_ht.rename({'locus1_b38': 'locus'}))
    lifted_over_variants = hl.sorted(
        hl.array([
            liftover_expr(cmg_ht.locus1_b38, cmg_ht.alleles1_b38, liftover_references[1]),
            liftover_expr(cmg_ht.locus2_b38, cmg_ht.alleles2_b38, liftover_references[1])
        ]),
        lambda x: x.locus
    )

    cmg_ht = cmg_ht.key_by(
        locus1=lifted_over_variants[0].locus,
        alleles1=lifted_over_variants[0].alleles,
        locus2=lifted_over_variants[1].locus,
        alleles2=lifted_over_variants[1].alleles
    )

    return cmg_ht.annotate(
        bad_liftover=(
                hl.is_missing(cmg_ht.locus1) |
                hl.is_missing(cmg_ht.locus2) |
                (cmg_ht.locus1.sequence_context() != cmg_ht.alleles1[0][0]) |
                (cmg_ht.locus2.sequence_context() != cmg_ht.alleles2[0][0])
        )
    )


def main(args):
    # Load CMG data
    cmg_ht = load_cmg(args.cmg)

    # Add phase information
    phased_ht = compute_phase(cmg_ht)

    # Flatten and export
    flatten_phased_ht(phased_ht).export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--cmg', help='CMG file containing variants to phase.',
                          default='gs://gnomad/projects/compound_hets/Feb_2020_CMG_Compound_het_list.csv')
    data_grp.add_argument('--out', help='Output TSV file',
                          default='gs://gnomad/projects/compound_hets/Feb_2020_CMG_Compound_het_list_phased.tsv')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

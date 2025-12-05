import hail as hl
from gnomad.utils.liftover import get_liftover_genome
from phasing import compute_phase, flatten_phased_ht
import argparse
from math import ceil
import logging
from resources import LEAST_CONSEQUENCE, MAX_FREQ

logger = logging.getLogger("compute_phase")
logger.setLevel(logging.INFO)


def get_sorted_variants_expr(
        locus1: hl.expr.LocusExpression,
        alleles1: hl.expr.ArrayExpression,
        locus2: hl.expr.LocusExpression,
        alleles2: hl.expr.ArrayExpression
) -> hl.expr.StructExpression:
    if locus1.dtype.reference_genome.name == 'GRCh37':
        variants = [
            hl.struct(
                locus=locus1,
                alleles=alleles1
            ),
            hl.struct(
                locus=locus2,
                alleles=alleles2
            )
        ]
    else:
        logger.warning("Variants are not on GRCh37; they will be lifted over.")
        _, destination_ref = get_liftover_genome(hl.struct(locus=locus1))
        variants = [
            liftover_expr(
                locus=locus1,
                alleles=alleles1,
                destination_ref=destination_ref
            ),
            liftover_expr(
                locus=locus2,
                alleles=alleles2,
                destination_ref=destination_ref
            )
        ]

    sorted_variants = hl.sorted(variants, key=lambda x: x.locus)
    return hl.struct(
        locus1=sorted_variants[0].locus,
        alleles1=sorted_variants[0].alleles,
        locus2=sorted_variants[1].locus,
        alleles2=sorted_variants[1].alleles
    )


def read_variants_ht(path: str) -> hl.Table:
    variants_ht = hl.read_table(path)

    # Make sure that types match
    assert (
            isinstance(variants_ht.key[0], hl.expr.LocusExpression) &
            (variants_ht.key[1].dtype == hl.tarray(hl.tstr)) &
            isinstance(variants_ht.key[2], hl.expr.LocusExpression) &
            (variants_ht.key[3].dtype == hl.tarray(hl.tstr))
    )

    variants_ht = variants_ht.key_by(
        **get_sorted_variants_expr(
            variants_ht.key[0],
            variants_ht.key[1],
            variants_ht.key[2],
            variants_ht.key[3]
        )
    ).persist()

    return variants_ht


def variants_ht_from_text(variants: str) -> hl.Table:
    v1, v2 = variants.split(",")
    chr1, pos1, ref1, alt1 = v1.strip().split(":")
    chr2, pos2, ref2, alt2 = v2.strip().split(":")

    ref = "GRCh38" if chr1.startswith("chr") else "GRCh37"

    variants_ht = hl.Table.parallelize(
        rows=[
            get_sorted_variants_expr(
                hl.locus(chr1, int(pos1), reference_genome=ref),
                hl.array([ref1, alt1]),
                hl.locus(chr2, int(pos2), reference_genome=ref),
                hl.array([ref2, alt2])
            )
        ],
        key=['locus1', 'alleles1', 'locus2', 'alleles2']
    ).persist()

    logger.info("Found the following variant pair to phase:")
    variants_ht.show()

    return variants_ht


def main(args):
    # Load data
    variants_ht = read_variants_ht(args.ht) if args.ht else variants_ht_from_text(args.variants)

    # Add phase
    phased_ht = compute_phase(variants_ht, args.least_consequence, args.max_freq)

    # Write results
    if args.out.endswith(".ht"):
        phased_ht.write(args.out, overwrite=args.overwrite)
    else:
        phased_ht = flatten_phased_ht(phased_ht)
        phased_ht.export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--ht', help='HT containing variants. Needs to be keyed by locus1, alleles1, locus2, alleles2.')
    data_grp.add_argument('--variants', help='Variants to phase in format chr1:pos1:ref1:alt1,chr2:pos2:ref2:alt2. Note that chromosome needs to start with "chr" for GRCh38 variants')
    parser.add_argument('--least_consequence', help=f'Includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--max_freq', help=f'If specified, maximum global adj AF for genotypes table to emit. (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--out', help="Output file path. Output file format depends on extension (.ht, .tsv or .tsv.gz)")
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()

    main(args)

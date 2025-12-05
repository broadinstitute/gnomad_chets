"""
This module contains functions for computing phase information for variant pairs.

# Phasing models

Reference table for GT counts

| v0/v1 | BB    | Bb    | bb    |
|-------|-------|-------|-------|
| AA    | n00   | n01   | n02   |
| Aa    | n10   | n11   | n12   |
| aa    | n20   | n21   | n22   |

Corresponding indices

| v0/v1 | BB    | Bb    | bb    |
|-------|-------|-------|-------|
| AA    | 0     | 1     | 2     |
| Aa    | 3     | 4     | 5     |
| aa    | 6     | 7     | 8     |
"""

import logging

import hail as hl
from math import ceil

from gnomad.resources.grch37.gnomad import public_release

from chet_utils import vep_genes_expr
from resources import LEAST_CONSEQUENCE, MAX_FREQ, phased_vp_count_ht_path

logger = logging.getLogger("phasing")
logger.setLevel(logging.INFO)


def chet_likelihood_expr(gt_counts, e: float = 1e-6, distance: int = None):
    """
    ### Het model
    | Haplotype | Freq    |
    |-----------|:-------:|
    | aB        | p       |
    | Ab        | q       |
    | ab        | e       |
    | AB        | 1-p-q-e |


    Therefore, we have the following frequencies:

    | v0/v1 |           BB          |           Bb          |       bb      |
    |-------|:---------------------:|:---------------------:|:-------------:|
    | AA    | (1-p-q-e)<sup>2</sup> |     2*(1-p-q-e)*q     | q<sup>2</sup> |
    | Aa    |     2*(1-p-q-e)*p     | 2*(p*q + (1-p-q-e)*e) |     2*q*e     |
    | aa    |      p<sup>2<sup>     |         2*p*e         | e<sup>2</sup> |

    :param gt_counts:
    :param e:
    :param distance:
    :return:
    """
    n = 2 * hl.sum(gt_counts)
    p = (gt_counts[3] + gt_counts[4] + gt_counts[7] + 2 * gt_counts[6]) / n
    q = (gt_counts[1] + gt_counts[4] + gt_counts[5] + 2 * gt_counts[2]) / n
    x = 1 - p - q - e

    # Compute log-likelihoods
    def compute_chet_log_like(n,p,q,x):
        res = (
            hl.cond(
                (p > 0) & (q > 0),
                hl.fold(
                    lambda i, j: i + j[0] * j[1],
                    0,
                    hl.zip(gt_counts,
                           [hl.log10(x) * 2, hl.log10(2 * x * q), hl.log10(q) * 2,
                            hl.log10(2 * x * p), hl.log10(2 * (p * q + x * e)), hl.log10(2 * q * e),
                            hl.log10(p) * 2, hl.log10(2 * p * e), hl.log10(e) * 2]
                           )
                ),
                -1e-31
            )
        )
        # If desired, add distance posterior based on value derived from regression
        if distance is not None:
            res = res + hl.max(-6, hl.log10(0.03 + 0.03 * hl.log(distance - 1)))

        return res

    return hl.bind(compute_chet_log_like, n, p, q, x)


def same_hap_likelihood_expr(gt_counts, e: float = 1e-6, distance: int = None):
    """
    ### Same haplotype model

    | Haplotype | Frequency |
    |-----------|:---------:|
    | aB        |     p     |
    | Ab        |     e     |
    | ab        |     q     |
    | AB        | 1-p-q-e   |

    With: p >= q and p = 0 if in perfect LD.


    Therefore, we have the following frequencies:

    | v0/v1 |           BB          |           Bb          |       bb      |
    |-------|:---------------------:|:---------------------:|:-------------:|
    | AA    | (1-p-q-e)<sup>2</sup> |     2*(1-p-q-e)*e     | e<sup>2</sup> |
    | Ab    |     2*(1-p-q-e)*p     | 2*(p*e + (1-p-q-e)*q) |     2*q*e     |
    | ab    |      p<sup>2<sup>     |         2*p*q         | q<sup>2</sup> |

    :param gt_counts:
    :param e:
    :param distance:
    :return:
    """
    n = 2 * hl.sum(gt_counts)
    f1 = hl.sum(gt_counts[3:6] + 2 * hl.sum(gt_counts[6:])) / n
    f2 = (gt_counts[1] + gt_counts[4] + gt_counts[7] + 2 * (gt_counts[2] + gt_counts[5] + gt_counts[8])) / n
    p = hl.cond(f1 > f2, f1, f2)
    q = (gt_counts[4] + gt_counts[5] + gt_counts[7] + 2 * gt_counts[8]) / n
    x = 1 - p - q - e

    # Compute log-likelihoods
    def compute_same_hap_log_like(n,p,q,x):
        res = (
            hl.cond(
                q > 0,
                hl.fold(
                    lambda i, j: i + j[0] * j[1],
                    0.0,
                    hl.zip(gt_counts,
                           [hl.log10(x) * 2, hl.log10(2 * x * e), hl.log10(e) * 2,
                            hl.log10(2 * x * p), hl.log10(2 * (p * e + x * q)), hl.log10(2 * q * e),
                            hl.log10(p) * 2, hl.log10(2 * p * q), hl.log10(q) * 2]
                           )
                ),
                -1e31  # Very large negative value if no q is present
            )
        )

        # If desired, add distance posterior based on value derived from regression
        if distance is not None:
            res = res + hl.max(-6, hl.log10(0.97 - 0.03 * hl.log(distance + 1)))

        return res

    return hl.bind(compute_same_hap_log_like, n, p, q, x)


def get_em_expr(gt_counts):
    hap_counts = hl.experimental.haplotype_freq_em(gt_counts)
    return hl.bind(
        lambda x: hl.struct(
            hap_counts=x,
            p_chet=(x[1] * x[2]) / (x[0] * x[3] + x[1] * x[2])
        ),
        hap_counts
    )


def get_em_expressions(gt_counts):
    return dict(
        em=hl.struct(
            raw=get_em_expr(gt_counts.raw),
            adj=get_em_expr(gt_counts.adj),
        ),
        em_plus_one=hl.struct(
            raw=get_em_expr(gt_counts.raw + [0, 0, 0, 0, 1, 0, 0, 0, 0]),
            adj=get_em_expr(gt_counts.adj + [0, 0, 0, 0, 1, 0, 0, 0, 0]),
        )
    )


def get_lr_expressions(gt_counts):
    def get_lr_annotation(gt_counts):
        same_hap_likelihood = same_hap_likelihood_expr(gt_counts)
        chet_likelihood = chet_likelihood_expr(gt_counts)
        return hl.bind(
            lambda x, y: hl.struct(
                same_hap_like=x,
                chet_like=y,
                lr_chet=y - x
            ),
            same_hap_likelihood,
            chet_likelihood
        )

    return dict(
        likelihood_model=hl.struct(
            raw=get_lr_annotation(gt_counts.raw),
            adj=get_lr_annotation(gt_counts.adj)
        )
    )


def get_single_het_expressions(gt_counts):
    def get_single_het_expr(gt_counts):
        return hl.bind(
            lambda x:
            hl.cond(x[1] > x[3],
                    (x[1] + x[2]) / hl.sum(hl.range(1, 9).filter(lambda x: x % 3 > 0).map(lambda y: x[y])),
                    (x[3] + x[6]) / hl.sum(x[3:])
                    ),
            gt_counts
        )

    return dict(
        singlet_het_ratio=hl.struct(
            raw=get_single_het_expr(gt_counts.raw),
            adj=get_single_het_expr(gt_counts.adj)
        )
    )


def flatten_gt_counts(gt_counts: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
    """
    Flattens the GT count array into a struct

    :param gt_counts: Array of GT counts
    :return: Struct with GT counts using ref/het/hom items
    """
    # [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb]
    return hl.struct(
        ref_ref=gt_counts[0],
        ref_het=gt_counts[1],
        ref_hom=gt_counts[2],
        het_ref=gt_counts[3],
        het_het=gt_counts[4],
        het_hom=gt_counts[5],
        hom_ref=gt_counts[6],
        hom_het=gt_counts[7],
        hom_hom=gt_counts[8]
    )

def flatten_phased_ht(phased_ht: hl.Table) -> hl.Table:
    """
    Flatten a phased Hail table by restructuring nested phase information into a flat format.
    
    This function takes a Hail table containing phased variant pair information and 
    transforms it into a flattened structure where nested genotype counts and EM 
    (Expectation-Maximization) results are expanded into individual columns. The 
    function also extracts variant coordinates and alleles from the nested locus/allele 
    structures.

    The input ```phased_ht``` should have the following structure:
        - locus1, locus2: Locus fields for the two variants
        - alleles1, alleles2: Array fields containing reference and alternate alleles
        - phase_info: Struct containing nested phase information with:
            - gt_counts: Struct with 'raw' and 'adj' arrays of genotype counts.
            - em: Struct with 'raw' and 'adj' EM results containing p_chet values

    The output ```flattened_ht``` will have the following structure:
        - chrom: Chromosome (from locus1.contig)
        - pos1, pos2: Positions from locus1.position and locus2.position
        - ref1, alt1: Reference and alternate alleles from alleles1
        - ref2, alt2: Reference and alternate alleles from alleles2
        - raw_*: Flattened genotype counts from phase_info.gt_counts.raw
        - adj_*: Flattened genotype counts from phase_info.gt_counts.adj
        - em_p_chet_raw: EM p_chet value from phase_info.em.raw.p_chet
        - em_p_chet_adj: EM p_chet value from phase_info.em.adj.p_chet

        The genotype count fields follow the pattern:
            - ref_ref, ref_het, ref_hom: Reference/heterozygous/homozygous for 
              variant 1, reference for variant 2
            - het_ref, het_het, het_hom: Heterozygous for variant 1, ref/het/hom for 
              variant 2
            - hom_ref, hom_het, hom_hom: Homozygous for variant 1, ref/het/hom for
              variant 2
        
    :param phased_ht: A Hail table containing phased variant pair data.
    :return: The flattened phased HT.
    """
    phased_ht = phased_ht.key_by()

    def flatten_phase_dict(expr: hl.expr.StructExpression) -> hl.expr.StructExpression:
        """
        Flatten the phase information from a dictionary structure into a flat structure.

        :param expr: The phase information struct to flatten.
        :return: The flattened phase information struct.
        """
        return hl.struct(
            raw=flatten_gt_counts(expr.gt_counts.raw),
            adj=flatten_gt_counts(expr.gt_counts.adj),
            em_p_chet_raw=expr.em.raw.p_chet,
            em_p_chet_adj=expr.em.adj.p_chet,
        )

    return phased_ht.transmute(
        chrom=phased_ht.locus1.contig,
        pos1=phased_ht.locus1.position,
        ref1=phased_ht.alleles1[0],
        alt1=phased_ht.alleles1[1],
        pos2=phased_ht.locus2.position,
        ref2=phased_ht.alleles2[0],
        alt2=phased_ht.alleles2[1],
        **{
            k: v for k, v in flatten_phase_dict(phased_ht.phase_info).items()
        }
    ).flatten()


def get_ac_from_gt_counts(gt_counts: hl.expr.ArrayNumericExpression, v1: bool) -> hl.expr.Float32Expression:
    if v1:
        return hl.sum(gt_counts[3:6])+2*hl.sum(gt_counts[6:])
    else:
        return (
                gt_counts[1] + gt_counts[4] + gt_counts[7] +
                2 *(gt_counts[2] + gt_counts[5] + gt_counts[8])
        )


def get_phased_gnomad_ht(
        ht: hl.Table,
        em: bool = True,
        lr: bool = True,
        shr: bool = True
) -> hl.Table:
    expr_fun = []

    if em:
        expr_fun.append(get_em_expressions)

    if lr:
        expr_fun.append(get_lr_expressions)

    if shr:
        expr_fun.append(get_single_het_expressions)

    if not expr_fun:
        raise (Exception("No expressions to annotate"))

    # Support for both exploded or dict versions of gt_counts
    # dict
    if isinstance(ht.gt_counts, hl.expr.DictExpression):
        ht = ht.select(
            phase_info=ht.gt_counts.map_values(
                lambda pop_count: hl.bind(
                    lambda x: hl.struct(
                        gt_counts=x,
                        **{
                            k: v for f in expr_fun for k, v in f(x).items()
                        }
                    ),
                    hl.struct(
                        raw=pop_count.raw.map(lambda y: hl.int32(y)),
                        adj=pop_count.adj.map(lambda z: hl.int32(z))
                    )
                )
            )
        )
    # exploded
    else:
        ht = ht.annotate(
            **{
                k: v for f in expr_fun for k, v in f(ht.gt_counts).items()
            }
        )

    return ht


def explode_phase_info(ht: hl.Table, remove_all_ref: bool = True) -> hl.Table:
    """
    Explode phase information from a dictionary structure into individual rows per 
    genetic ancestry group.
    
    This function transforms a Hail table where phase_info is stored as a dictionary
    (keyed by genetic ancestry group) into a flattened structure where each genetic 
    ancestry group's phase information becomes a separate row. The function also 
    optionally filters out variant pairs that have no alternate allele carriers in 
    any genetic ancestry group.

    The input ```ht``` should have the following structure:
        - phase_info: Dictionary with genetic ancestry group names as keys and phase 
          information as values
        - Each phase_info value contains gt_counts and em fields
    
    The output ```ht``` will have the following structure:
        - Each row represents one genetic ancestry group's phase information for a 
          variant pair
        - pop: Genetic ancestry group identifier extracted from the original dictionary 
          key
        - phase_info: The phase information struct for that genetic ancestry group 
          containing:
            - gt_counts: Struct with 'raw' and 'adj' genotype count arrays
            - em: Struct with 'raw' and 'adj' EM results
        - All other fields from the input table are preserved

    :param ht: The phased HT to explode.
    :param remove_all_ref: If True, filters out variant pairs where all samples across
        all genetic ancestry groups are homozygous reference (i.e., no alternate 
        allele carriers).
        This is determined by checking if the sum of genotype counts for alternate
        alleles (indices 1-8) is greater than 0. Default is True.
    :return: The exploded phased HT.
    """
    ht = ht.transmute(phase_info=hl.array(ht.phase_info))
    ht = ht.explode('phase_info')
    ht = ht.transmute(
        pop=ht.phase_info[0],
        phase_info=ht.phase_info[1]
    )

    if remove_all_ref:
        ht = ht.filter(hl.sum(ht.phase_info.gt_counts.raw[1:]) > 0)

    return ht


# TODO: How to handle one variant absent from gnomAD?
def annotate_unphased_pairs(
    unphased_ht: hl.Table,
    n_variant_pairs: int,
    least_consequence: str = LEAST_CONSEQUENCE,
    max_af: float = MAX_FREQ
):
    """
    Annotate unphased variant pairs with frequency-based phase estimates from gnomAD SNV data.
    
    This function processes variant pairs that were not found in gnomAD phased data and 
    estimates their phase relationships using individual variant frequencies from 
    gnomAD exomes. The function explodes variant pairs into individual variants, 
    annotates them with gnomAD frequency data, and then reconstructs genotype count 
    estimates assuming variants are unlinked (independent).
    
    The function applies several filtering steps:
        1. Filters out variant pairs where at least one variant has AF > max_af
        2. Filters out variant pairs where variants are not in the same gene with at 
           least the specified consequence severity
        3. Estimates genotype counts assuming independence between variants
    

    The output ```unphased_ht``` will have the following structure:
        - locus1, locus2, alleles1, alleles2: Original variant pair identifiers
        - pop: Population identifier (set to 'all')
        - phase_info: Struct containing estimated phase information with:
            - gt_counts: Struct with 'raw' and 'adj' genotype count arrays (9-element 
              arrays following the standard format: [AABB, AABb, AAbb, AaBB, AaBb, 
              Aabb, aaBB, aaBb, aabb])
            - em: Struct with 'raw' and 'adj' EM results from get_em_expr()
        - vep_genes: Array of genes containing the variants (filtered to non-empty 
          entries)

    .. note::

        The genotype count estimation assumes variants are unlinked (independent), 
        which may not be accurate for closely linked variants. The function logs the 
        number of variant pairs excluded due to frequency and VEP filtering criteria.

    :param unphased_ht: A Hail table containing variant pairs that were not found in 
        gnomAD phased data. Expected to have locus1, locus2, alleles1, and alleles2 
        fields.
    :param n_variant_pairs: Total number of variant pairs being processed (used for
        repartitioning optimization).
    :param least_consequence: Minimum VEP consequence severity for variant inclusion.
        Variants must be in the same gene with consequences at least as severe as this
        threshold to be included in the analysis. Default is ```LEAST_CONSEQUENCE```.
    :param max_af: Maximum global adjusted allele frequency threshold. Variant pairs
        where at least one variant has AF > max_af are excluded from the analysis.
        Default is ```MAX_FREQ```.
    :return: The unphased HT with estimated phase information.
    """
    # Explode variant pairs.
    unphased_ht = unphased_ht.annotate(
        las=[
            hl.tuple([unphased_ht.locus1, unphased_ht.alleles1]),
            hl.tuple([unphased_ht.locus2, unphased_ht.alleles2])
        ]
    ).explode('las', name='la')

    unphased_ht = unphased_ht.key_by(
        locus=unphased_ht.la[0],
        alleles=unphased_ht.la[1]
    ).checkpoint(hl.utils.new_temp_file("unphased_ht"))

    # Annotate single variants with gnomAD freq.
    gnomad_ht = public_release('exomes').ht()
    gnomad_ht = gnomad_ht.semi_join(unphased_ht).repartition(
        ceil(n_variant_pairs / 10000),
        shuffle=True
    ).checkpoint(hl.utils.new_temp_file("gnomad_ht"))

    missing_freq = hl.struct(
        AC=0,
        AF=0,
        AN=125748 * 2,  # TODO:set to no missing for now.
        homozygote_count=0
    )

    logger.info(
        f"{gnomad_ht.count()}/{unphased_ht.count()} single variants from the unphased "
        "pairs found in gnomAD."
    )

    gnomad_indexed = gnomad_ht[unphased_ht.key]
    gnomad_freq = gnomad_indexed.freq
    unphased_ht = unphased_ht.annotate(
        adj_freq=hl.or_else(
            gnomad_freq[0],
            missing_freq
        ),
        raw_freq=hl.or_else(
            gnomad_freq[1],
            missing_freq
        ),
        vep_genes=vep_genes_expr(gnomad_indexed.vep, least_consequence),
        max_af_filter=gnomad_indexed.freq[0].AF <= max_af
        # pop_max_freq=hl.or_else(
        #     gnomad_exomes.popmax[0],
        #     missing_freq.annotate(
        #         pop=hl.null(hl.tstr)
        #     )
        # )
    ).checkpoint(hl.utils.new_temp_file("unphased_ht"))

    loci_expr = hl.sorted(
        hl.agg.collect(
            hl.tuple([
                unphased_ht.locus,
                hl.struct(
                    adj_freq=unphased_ht.adj_freq,
                    raw_freq=unphased_ht.raw_freq,
                    # pop_max_freq=unphased_ht.pop_max_freq
                )
            ])
        ),
        lambda x: x[0]  # Sort by locus.
    ).map(lambda x: x[1])   # Get rid of locus.

    vp_freq_expr = hl.struct(
        v1=loci_expr[0],
        v2=loci_expr[1]
    )

    def get_gt_counts(freq: str) -> hl.expr.ArrayExpression:
        """
        Get the genotype counts for a given frequency type under the independence assumption.

        The array is in the following order: [AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB,
        aaBb, aabb]

        This function estimates genotype counts assuming the two variants are 
        independent (unlinked and mutually exclusive). Under this assumption, compound 
        heterozygous (AaBb) and double homozygous alternate (aabb) genotypes are set to 
        zero, as well as mixed heterozygous-homozygous combinations (Aabb, aaBb).

        :param freq: The frequency type to get the genotype counts for ('raw_freq' or 
            'adj_freq').
        :return: The genotype counts array with zeros for AaBb (index 4), Aabb (index 5),
            aaBb (index 7), and aabb (index 8) due to the independence assumption.
        """
        return hl.array([
            hl.min(vp_freq_expr.v1[freq].AN, vp_freq_expr.v2[freq].AN),  # AABB.
            vp_freq_expr.v2[freq].AC - (2 * vp_freq_expr.v2[freq].homozygote_count),  # AABb.
            vp_freq_expr.v2[freq].homozygote_count,  # AAbb.
            vp_freq_expr.v1[freq].AC - (2 * vp_freq_expr.v1[freq].homozygote_count),  # AaBB.
            0,  # AaBb.
            0,  # Aabb.
            vp_freq_expr.v1[freq].homozygote_count,  # aaBB.
            0,  # aaBb.
            0  # aabb.
        ])

    gt_counts_raw_expr = get_gt_counts('raw_freq')
    gt_counts_adj_expr = get_gt_counts('adj_freq')

    # gt_counts_pop_max_expr = get_gt_counts('pop_max_freq')
    unphased_ht = unphased_ht.group_by(
        unphased_ht.locus1,
        unphased_ht.alleles1,
        unphased_ht.locus2,
        unphased_ht.alleles2
    ).aggregate(
        pop='all',  # TODO Add option for multiple pops?
        phase_info=hl.struct(
            gt_counts=hl.struct(
                raw=gt_counts_raw_expr,
                adj=gt_counts_adj_expr
            ),
            em=hl.struct(
                raw=get_em_expr(gt_counts_raw_expr),
                adj=get_em_expr(gt_counts_raw_expr)
            )
        ),
        vep_genes=hl.agg.collect(unphased_ht.vep_genes).filter(lambda x: hl.len(x) > 0),
        max_af_filter=hl.agg.all(unphased_ht.max_af_filter)

        # pop_max_gt_counts_adj=gt_counts_raw_expr,
        # pop_max_em_p_chet_adj=get_em_expr(gt_counts_raw_expr).p_chet,
    ) # .key_by()

    unphased_ht = unphased_ht.transmute(
        vep_filter=(hl.len(unphased_ht.vep_genes) > 1) &
        (hl.len(unphased_ht.vep_genes[0].intersection(unphased_ht.vep_genes[1])) > 0)
    )

    max_af_filtered, vep_filtered = unphased_ht.aggregate([
        hl.agg.count_where(~unphased_ht.max_af_filter),
        hl.agg.count_where(~unphased_ht.vep_filter)
    ])
    if max_af_filtered > 0:
        logger.info(f"{max_af_filtered} variant-pairs excluded because the AF of at least one variant was > {max_af}")
    if vep_filtered > 0:
        logger.info(f"{vep_filtered} variant-pairs excluded because the variants were not found within the same gene with a csq of at least {least_consequence}")

    unphased_ht = unphased_ht.filter(
        unphased_ht.max_af_filter &
        unphased_ht.vep_filter
    )

    return unphased_ht.drop('max_af_filter', 'vep_filter')


def compute_phase(
    variants_ht: hl.Table,
    least_consequence: str = LEAST_CONSEQUENCE,
    max_freq: float = MAX_FREQ
) -> hl.Table:
    """
    Compute phase information for variant pairs by combining gnomAD phased data with single-variant estimates.
    
    This function performs a two-step process to determine phase relationships between 
    variant pairs:
        1. First, it looks up existing phase information from gnomAD phased variant pair 
           data
        2. For variant pairs not found in gnomAD phased data, it computes phase 
           estimates from individual variant frequencies using the 
           ```annotate_unphased_pairs``` function
    
    The function joins the input variant pairs with gnomAD exome phased variant pair 
    counts, explodes the phase information by genetic ancestry group, and combines 
    results from both phased and unphased variant pairs into a unified table.

    The input ```variants_ht``` should have the following structure:
        - locus1, locus2: Locus fields for the two variants
        - alleles1, alleles2: Array fields containing reference and alternate alleles

    The output ```phased_ht``` will have the following structure:
        - All original variant pair fields (locus1, locus2, alleles1, alleles2)
        - pop: Genetic ancestry group identifier (from exploded phase_info)
        - phase_info: Struct containing:
            - gt_counts: Struct with 'raw' and 'adj' genotype count arrays
            - em: Struct with 'raw' and 'adj' EM (Expectation-Maximization) results
        - Additional fields from annotate_unphased_pairs for unphased variant pairs

    :param variants_ht: A Hail table containing variant pairs to analyze.
    :param least_consequence: Minimum VEP consequence severity for variant inclusion.
        Variants with consequences at least as severe as this threshold are included.
        Default is ```LEAST_CONSEQUENCE```.
    :param max_freq: Maximum global adjusted allele frequency threshold for variant 
        inclusion. Variants with frequencies above this threshold are excluded. 
        Default is ```MAX_FREQ```.
    :return: The phased HT.
    """
    n_variant_pairs = variants_ht.count()
    logger.info(f"Looking up phase for {n_variant_pairs} variant pair(s).")

    # Join with gnomad phased variants.
    vp_ht = hl.read_table(phased_vp_count_ht_path('exomes'))
    phased_ht = vp_ht.semi_join(variants_ht)
    n_phased = phased_ht.count()
    
    # Explode phase_info by genetic ancestry group (pop in v2).
    phased_ht = explode_phase_info(phased_ht)
    phased_ht = phased_ht.transmute(
        phase_info=phased_ht.phase_info.select('gt_counts', 'em')
    ).repartition(ceil(n_variant_pairs / 10000), shuffle=True)
    phased_ht = phased_ht.checkpoint(hl.utils.new_temp_file("phased_ht"))

    # If not all pairs had at least one carrier of both, then compute phase estimate 
    # from single variants.
    logger.info(
        f"{n_phased}/{n_variant_pairs} variant pair(s) found with carriers of both in " 
        "gnomAD."
    )

    if n_phased < n_variant_pairs:
        unphased_ht = variants_ht.anti_join(vp_ht)
        unphased_ht = annotate_unphased_pairs(
            unphased_ht,
            n_variant_pairs,
            least_consequence,
            max_freq
        )
        phased_ht = phased_ht.union(
            unphased_ht,
            unify=True
        )

    return phased_ht

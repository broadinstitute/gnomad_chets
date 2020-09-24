import hail as hl

"""
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
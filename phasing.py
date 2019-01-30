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
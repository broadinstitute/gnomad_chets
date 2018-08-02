from gnomad_hail import *
import argparse
import hail as hl


def main(args):
    if args.import_vcf:
        mt = hl.import_vcf(args.import_vcf, call_fields=['GT'], force_bgz=True)

        mt.write("{}.mt".format(args.output))

    if args.aggregate_genotypes:
        mt = hl.read_matrix_table("{}.mt".format(args.output))

        vp_sites_ht = hl.import_table('gs://gnomad/projects/compound_hets/gnomad_exomes_chrom21_variant_pairs.tsv',
                                      types={'chrom1': hl.tstr, 'pos1': hl.tint32, 'pos2': hl.tint32}, delimiter="\s")

        vp_sites_ht = vp_sites_ht.transmute(
            locus1=hl.locus(vp_sites_ht.chrom1, vp_sites_ht.pos1),
            locus2=hl.locus(vp_sites_ht.chrom1, vp_sites_ht.pos2),
            alleles1=[vp_sites_ht.ref1, vp_sites_ht.alt1],
            alleles2=[vp_sites_ht.ref2, vp_sites_ht.alt2]
        )

        # Mark v1 sites to keep them
        # v1_sites_ht = vp_ht.group_by(vp_ht.locus1, vp_ht.lo)
        # v1_sites_mt = hl.MatrixTable.from_rows_table(vp_ht, partition_key='locus1')
        # v1_sites_mt = v1_sites_mt.key_by('locus1','alleles1', partition_key='locus1')
        # mt = mt.annotate_rows(
        #    v1=hl.is_defined(v1_sites_mt[(mt.locus, mt.alleles),:])
        # )

        # Compress genotype data
        n_samples = hl.literal(mt.count_cols())
        mt = mt.select_rows(
            mt.locus,
            mt.alleles,
            missing=hl.agg.collect_as_set(hl.agg.filter(hl.is_missing(mt.GT), mt.s)),
            hets=hl.agg.collect(hl.agg.filter(mt.GT.is_het(), hl.struct(s=mt.s, hap1_alt=mt.GT[0] == 1))),
            hom_alts=hl.agg.collect_as_set(hl.agg.filter(mt.GT.is_hom_var(), mt.s)),
        )
        mt = mt.drop()
        mt = mt.filter_cols(False)

        # Create v1 table
        v1_sites_ht = vp_sites_ht.group_by('locus1', 'alleles1').aggregate(
            v2=hl.agg.collect(hl.struct(locus=vp_sites_ht.locus2, alleles=vp_sites_ht.alleles2)))
        v1_sites_mt = hl.MatrixTable.from_rows_table(v1_sites_ht.key_by('locus1', 'alleles1'), partition_key='locus1')
        mt = mt.annotate_rows(
            **v1_sites_mt[(mt.locus, mt.alleles), :]
        )
        v1_ht = mt.filter_rows(hl.is_defined(mt.v2)).rows()
        v1_ht = v1_ht.explode('v2')
        v1_ht = v1_ht.transmute(
            locus1=v1_ht.locus,
            alleles1=v1_ht.alleles,
            locus2=v1_ht.v2.locus,
            alleles2=v1_ht.v2.alleles,
            missing1=v1_ht.missing,
            hets1=v1_ht.hets,
            hom_alts1=v1_ht.hom_alts,

        )

        # Create v2 table
        v2_sites_mt = hl.MatrixTable.from_rows_table(vp_sites_ht.key_by('locus2', 'alleles2'), partition_key='locus2')
        v2_ht = mt.filter_rows(hl.is_defined(v2_sites_mt[(mt.locus, mt.alleles), :])).rows()
        v2_ht = v2_ht.transmute(
            locus2=v2_ht.locus,
            alleles2=v2_ht.alleles,
            missing2=v2_ht.missing,
            hets2=v2_ht.hets,
            hom_alts2=v2_ht.hom_alts,
        )

        # Join the tables (fingers crossed!!!)
        vp_ht = v1_ht.key_by('locus2', 'alleles2').join(v2_ht.key_by('locus2', 'alleles2'))
        vp_ht = vp_ht.annotate_globals(n_samples=n_samples)
        vp_ht.write("{}.phased.ht".format(args.output))

    if args.compute_results:
        vp_ht = hl.read_table("{}.phased.ht".format(args.output))

        ##TODO: Remove these 3 lines when HT re-created
        mt = hl.read_matrix_table("{}.mt".format(args.output))
        n_samples = mt.count_cols()
        print("Found {} samples in MT.".format(n_samples))
        vp_ht = vp_ht.annotate_globals(n_samples=n_samples)

        # Computes haplotypes
        vp_ht = vp_ht.drop(vp_ht.v2)
        vp_ht = vp_ht.transmute(
            chrom1=vp_ht.locus1.contig,
            pos1=vp_ht.locus1.position,
            chrom2=vp_ht.locus2.contig,
            pos2=vp_ht.locus2.position,
            ref1=vp_ht.alleles1[0],
            alt1=vp_ht.alleles1[1],
            ref2=vp_ht.alleles2[0],
            alt2=vp_ht.alleles2[1],
            # n00=vp_ht.n_samples - hl.len(
            #     hl.set(vp_ht.hets1.map(lambda x: x.s)).union(
            #         hl.set(vp_ht.hets2.map(lambda x: x.s))).union(
            #         vp_ht.missing1.union(
            #             vp_ht.missing2.union(
            #                 vp_ht.hom_alts1.union(
            #                     vp_ht.hom_alts2
            #                 )
            #             )
            #         )
            #     )),
            n10=hl.len(vp_ht.hets1.filter(
                lambda x:
                ~vp_ht.hom_alts2.contains(x.s) &
                ~vp_ht.missing2.contains(x.s) &
                ~vp_ht.hets2.any(lambda y: (x.s == y.s) & (x.hap1_alt == y.hap1_alt))

            )) + hl.sum(
                vp_ht.hom_alts1.map(
                    lambda x:
                    (hl.case()
                     .when(vp_ht.hets2.any(lambda y: x == y.s), 1)
                     .when(vp_ht.hom_alts2.contains(x) |
                           vp_ht.missing2.contains(x), 0)
                     .default(2)
                     )

                )
            ),
            n01=hl.len(vp_ht.hets2.filter(
                lambda x:
                ~vp_ht.hom_alts1.contains(x.s) &
                ~vp_ht.missing1.contains(x.s) &
                ~vp_ht.hets1.any(lambda y: (x.s == y.s) & (x.hap1_alt == y.hap1_alt))
            )) + hl.sum(
                vp_ht.hom_alts2.map(
                    lambda x:
                    (hl.case()
                     .when(vp_ht.hets1.any(lambda y: x == y.s), 1)
                     .when(vp_ht.hom_alts1.contains(x) |
                           vp_ht.missing1.contains(x), 2)
                     .default(2)
                     )

                )
            ),
            n11=(
                    hl.len(vp_ht.hets1.filter(
                        lambda x:
                        vp_ht.hom_alts2.contains(x.s) |
                        vp_ht.hets2.any(lambda y: (x.s == y.s) & (x.hap1_alt == y.hap1_alt)))
                    ) +
                    hl.sum(
                        vp_ht.hom_alts1.map(
                            lambda x:
                            (hl.case()
                             .when(vp_ht.hets2.any(lambda y: y.s == x), 1)
                             .when(vp_ht.hom_alts2.contains(x), 2)
                             .default(0)
                             )
                        )
                    )
            ),
            n_missing=hl.len(vp_ht.missing1.union(vp_ht.missing2))
        )

        vp_ht = vp_ht.annotate(
            n00=(
                    2 * (vp_ht.n_samples - vp_ht.n_missing) -
                    vp_ht.n01 -
                    vp_ht.n10 -
                    vp_ht.n11
            )
        )

        vp_ht = vp_ht.annotate(
            prob_eagle=(vp_ht.n00 * vp_ht.n11) / (vp_ht.n00 * vp_ht.n11 + vp_ht.n10 * vp_ht.n01)
        )

        vp_ht.export("{}.phased.tsv.bgz".format(args.output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--import_vcf', help='import target Eagle VCF and creates phase genotype representation.')
    parser.add_argument('--aggregate_genotypes', help='Creates variant-pair table with aggregated genotypes.',
                        action='store_true')
    parser.add_argument('--compute_results', help='Computes resulting phased haplotype counts.', action='store_true')
    parser.add_argument('--variant_pairs', help="File with variant pairs.")
    parser.add_argument('--output', help='Output prefix.',
                        default='gs://gnomad/projects/compound_hets/eagle/vcf/exomes.chr21.no_singletons.no_fams.split.phased')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

import hail as hl
from resources import LEAST_CONSEQUENCE, MAX_FREQ
from resources import (
    pbt_comparison_full_mt_path,
    pbt_comparison_vp_count_ht_path,
    pbt_comparison_phased_vp_count_ht_path,
    full_mt_path,
    vp_ann_ht_path,
    pbt_phase_count_ht_path
)
from gnomad_qc.v2.resources import get_gnomad_meta, get_gnomad_data
from create_vp_matrix import create_full_vp, create_vp_ann, create_vp_summary
import argparse
import logging
from phasing import get_ac_from_gt_counts, get_phased_gnomad_ht

logger = logging.getLogger("create_gnomad_pbt_comparison_data")

def main(args):

    hl.init(log="/tmp/hail_comp_vp.log")
    data_type = 'exomes' if args.exomes else 'genomes'

    if args.create_full_vp:
        logger.info(f"Generating gnomAD VP MT for PBT VPs, excluding PBT samples.")
        # Load PBT VP MT
        pbt_vp_mt = hl.read_matrix_table(full_mt_path(data_type, True, args.least_consequence, args.max_freq, args.chrom))
        pbt_samples = pbt_vp_mt.cols().key_by('s')

        mt = get_gnomad_data(data_type)
        mt = mt.select_entries(
            GT=hl.or_missing(mt.GT.is_non_ref(), mt.GT),
            PID=mt.PID,
            missing=hl.is_missing(mt.GT),
            adj=mt.adj
        ).select_cols().select_rows()
        meta = get_gnomad_meta('exomes')
        mt = mt.filter_cols(meta[mt.col_key].release & hl.is_missing(pbt_samples[mt.col_key]))
        vp_mt = create_full_vp(
            mt,
            vp_list_ht=pbt_vp_mt.rows(),
            data_type=data_type
        )
        vp_mt.write(
            pbt_comparison_full_mt_path(
                data_type=data_type,
                least_consequence=args.least_consequence,
                max_freq=args.max_freq,
                chrom=args.chrom
            ),
            overwrite=args.overwrite
        )

    if args.create_vp_summary:
        logger.info("Creating VP summary")
        mt = hl.read_matrix_table(
            pbt_comparison_full_mt_path(
                data_type=data_type,
                least_consequence=args.least_consequence,
                max_freq=args.max_freq,
                chrom=args.chrom
            )
        )
        meta = get_gnomad_meta(data_type).select('pop', 'release')
        mt = mt.annotate_cols(**meta[mt.col_key])
        ht = create_vp_summary(mt)
        ht = ht.checkpoint(
            pbt_comparison_vp_count_ht_path(
                data_type=data_type,
                least_consequence=args.least_consequence,
                max_freq=args.max_freq,
                chrom=args.chrom
            ),
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite
        )

        logger.info("Phasing VP summary")
        ht = get_phased_gnomad_ht(ht)
        ht.write(
            pbt_comparison_phased_vp_count_ht_path(
                data_type=data_type,
                least_consequence=args.least_consequence,
                max_freq=args.max_freq,
                chrom=args.chrom
            ),
            overwrite=args.overwrite
        )

    if args.export:
        data_type = 'exomes'
        pbt = hl.read_table(pbt_phase_count_ht_path(data_type, pbt=True))
        pbt_ann = hl.read_table(vp_ann_ht_path(data_type, pbt=True))
        pbt = pbt.annotate(
            **pbt_ann[pbt.locus1, pbt.alleles1, pbt.locus2, pbt.alleles2],
            distance=pbt.locus2.position - pbt.locus1.position
        )

        discordant_within_pop_expr = hl.array(pbt.phase_by_pop).any(
            lambda x: (x[0] != 'all') & (x[1].adj.n_same_hap > 0) & (x[1].adj.n_chet > 0)
        )
        pbt = pbt.annotate(
            phase_by_pop=hl.array(pbt.phase_by_pop),
            discordant_within_pop=discordant_within_pop_expr,
            discordant_between_pops=~discordant_within_pop_expr & (pbt.phase_by_pop['all'].adj.n_same_hap > 0) & (pbt.phase_by_pop['all'].adj.n_chet > 0),
            **pbt_ann[pbt.key],
            distance=pbt.locus2.position - pbt.locus1.position
        )
        pbt = pbt.explode('phase_by_pop')
        pbt = pbt.transmute(
            pop=pbt.phase_by_pop[0],
            **pbt.phase_by_pop[1]
        )
        # Drop RF-filtered sites
        pbt = pbt.filter((hl.len(pbt.filters1) == 0) & (hl.len(pbt.filters2) == 0))
        # Filter sites that are in LCR / Decoy
        pbt = pbt.filter(pbt.lcr1 | pbt.lcr2 | pbt.decoy1 | pbt.decoy2 | pbt.segdup1 | pbt.segdup2, keep=False)

        # Drop sites with inconsistent trio phasing?
        pbt = pbt.filter(~pbt.discordant_within_pop)
        # pbt = pbt.filter(pbt.adj.n_same_hap + pbt.adj.n_chet > 0)
        # Drop sites that are too frequent in a given pop
        pbt = pbt.filter((pbt.freq1[pbt.pop].AF <= 0.05) & (pbt.freq2[pbt.pop].AF <= 0.05))
        # Drop sites that really are het non-ref
        pbt = pbt.filter(pbt.distance > 0)
        # filter to autosomes
        autosomes = hl.parse_locus_interval('1-22')
        pbt = pbt.filter(
            autosomes.contains(pbt.locus1)
        )
        phase_ht = hl.read_table(
            pbt_comparison_phased_vp_count_ht_path(
                data_type=data_type,
                least_consequence=args.least_consequence,
                max_freq=args.max_freq,
                chrom=args.chrom
            )
        )
        pbt = pbt.annotate(
            trio_chet=hl.struct(
                raw=hl.case()
                    .when((pbt.raw.n_same_hap > 0) & (pbt.raw.n_chet == 0), False)
                    .when((pbt.raw.n_same_hap == 0) & (pbt.raw.n_chet > 0), True)
                    .or_missing(),
                adj=hl.case()
                    .when((pbt.adj.n_same_hap > 0) & (pbt.adj.n_chet == 0), False)
                    .when((pbt.adj.n_same_hap == 0) & (pbt.adj.n_chet > 0), True)
                    .or_missing()
            ),
            **phase_ht[pbt.locus1, pbt.alleles1, pbt.locus2, pbt.alleles2].phase_info[pbt.pop]
        )

        pbt = pbt.filter(~hl.is_nan(pbt.gt_counts.adj[0]) & (pbt.pop != 'oth')).key_by()
        rf_features = {
            'snv1': pbt.snv1,
            'snv2': pbt.snv2,
            'cpg1': hl.or_else(pbt.cpg1, False),
            'cpg2': hl.or_else(pbt.cpg2, False),
            'distance': pbt.distance
        }
        rf_features.update({f'n{i}{j}': pbt.gt_counts.adj[(3 * i) + j] for i in [0, 1, 2] for j in [0, 1, 2]})

        ac1 = get_ac_from_gt_counts(pbt.gt_counts.adj, True)
        ac2 = get_ac_from_gt_counts(pbt.gt_counts.adj, False)
        an = 2 * hl.sum(pbt.gt_counts.adj)
        pbt_df = pbt.select(
            locus1=pbt.locus1,
            ref1=pbt.alleles1[0],
            alt1=pbt.alleles1[1],
            locus2=pbt.locus2,
            ref2=pbt.alleles2[0],
            alt2=pbt.alleles2[1],
            pop=pbt.pop,
            trio_chet=pbt.trio_chet.adj,
            em=pbt.em.adj.p_chet,
            singlet_het_ratio=pbt.singlet_het_ratio.adj,
            lr=pbt.likelihood_model.adj,
            AC1=ac1,
            AC2=ac2,
            AF1=ac1 / an,
            AF2=ac2 / an,
            n_var_gnomad=(ac1 > 0) + (ac2 > 0),
            discordant_between_pops=pbt.discordant_between_pops,
            discordant_within_pop=pbt.discordant_within_pop,
            **rf_features
        ).flatten().to_pandas() # NOTE: The serialization to pandas happens because this code comes from a notebook initially

        with hl.utils.hadoop_open('gs://gnomad/projects/compound_hets/pbt_annotated.csv', 'w') as f:
            pbt_df.to_csv(f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    data_grp = parser.add_mutually_exclusive_group(required=True)
    data_grp.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    data_grp.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--create_full_vp', help='Creates the VP MT.', action='store_true')
    parser.add_argument('--create_vp_summary', help='Creates a summarised VP table, with counts in release samples only. If --pbt is specified, then only sites present in PBT samples are used and counts exclude PBT samples.',
                        action='store_true')
    parser.add_argument('--export', help='Export data in CSV format.',
                        action='store_true')
    parser.add_argument('--least_consequence', help=f'Includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: {LEAST_CONSEQUENCE})',
                        default=LEAST_CONSEQUENCE)
    parser.add_argument('--max_freq', help=f'If specified, maximum global adj AF for genotypes table to emit. (default: {MAX_FREQ:.3f})', default=MAX_FREQ, type=float)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--chrom', help='Only run on given chromosome')

    args = parser.parse_args()
    main(args)

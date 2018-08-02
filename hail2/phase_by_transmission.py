from gnomad_hail import *
from gnomad_hail.resources.phasing import *
import hail as hl


# def get_variant_pairs_per_transctipt(mt: hl.MatrixTable, data_type: str, max_freq: float):
#     vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
#     freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
#
#     mt = mt.select_rows(
#             transcript=(
#                 vep_ht[mt.row_key].vep.transcript_consequences
#                 .filter(
#                     lambda tc: (tc.biotype == 'protein_coding') &
#                                (tc.allele_num == mt.a_index))
#                 ).map(
#                     lambda x: x.transcript_id
#                 ),
#             af=hl.float32(hl.find(lambda x: x.meta.get('group') == 'adj', freq[mt.row_key].freq).AF[1])
#         )
#
#     mt = mt.filter_rows(mt.af <= max_freq)
#
#     mt = mt.explode_rows(mt.transcript)
#     mt.describe()
#
#     mt = (
#         mt.group_rows_by(mt.transcript)
#         .aggregate(non_refs=hl.agg.collect(
#             hl.agg.filter(mt.GT.is_non_ref(), hl.tuple(mt.s, mt.locus, mt.alleles)))
#         )
#     )
#
#     mt = mt.annotate_rows(
#         variant_pairs=(
#             mt.non_refs
#             .group_by(lambda x: x[0])
#             .
#         )
#     )

def get_fam_ht(data_type: str):
    fam_ht = hl.import_fam(fam_path(data_type), delimiter='\\t')
    mat_ht = fam_ht.key_by(id=fam_ht.mat_id)
    pat_ht = fam_ht.key_by(id=fam_ht.pat_id)
    fam_ht = fam_ht.select(fam_id=hl.int(fam_ht.fam_id), proband=hl.struct(), is_female=hl.or_missing(fam_ht.is_female, hl.struct()))
    fam_ht = fam_ht.union(mat_ht.select(fam_id=hl.int(mat_ht.fam_id), proband=hl.null(hl.tstruct()), is_female=hl.struct()))
    return fam_ht.union(pat_ht.select(fam_id=hl.int(pat_ht.fam_id), proband=hl.null(hl.tstruct()), is_female=hl.null(hl.tstruct())))


def main(args):
    data_type = 'exomes' if args.exomes else 'genomes'

    if args.phase:
        mt = get_gnomad_data(data_type, split=False)
        ped = hl.Pedigree.read(fam_path(data_type), delimiter='\\t')

        ped_samples = hl.literal(set([s for trio in ped.complete_trios() for s in [trio.s, trio.pat_id, trio.mat_id]]))

        mt = mt.filter_cols(ped_samples.contains(mt.s) & mt.meta.high_quality)
        mt = mt.select_cols().select_rows()
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

        tm = hl.trio_matrix(mt, ped)
        tm = phase_trio_matrix_by_transmission(tm)

        pmt = explode_trio_matrix(tm)
        write_temp_gcs(pmt, pbt_phased_trios_mt_path(data_type, split=False), overwrite=args.overwrite)
        #pmt.write(pbt_phased_trios_mt_path(data_type, split=False), overwrite=args.overwrite)

        pmt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type, split=False))
        pmt = split_multi_dynamic(pmt)
        pmt.write(pbt_phased_trios_mt_path(data_type), overwrite=args.overwrite)

    if args.create_sparse_matrix:
        mt = hl.read_matrix_table(pbt_phased_trios_mt_path(data_type))

        fam_ht = get_fam_ht(data_type)
        mt = mt.annotate_cols(**fam_ht[mt.col_key])
        mt = mt.filter_cols(hl.is_missing(mt.proband))

        meta_ht = get_gnomad_meta(data_type)
        mt = mt.annotate_cols(pop=meta_ht[mt.col_key].pop)

        vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
        freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))

        consequences_to_keep = hl.literal(set(["inframe_insertion", "inframe_deletion", "missense_variant", "synonymous_variant"]))
        missense_consequences = hl.literal(set(["inframe_insertion", "inframe_deletion", "missense_variant"]))

        mt = mt.annotate_globals(
            vep_consequences_to_keep=consequences_to_keep,
            missense_consequences=missense_consequences
        )

        mt = mt.select_rows(
            vep=(
                vep_ht[mt.row_key].vep.transcript_consequences
                    .filter(
                    lambda tc: (tc.biotype == 'protein_coding') &
                               (tc.allele_num == mt.a_index) &
                               ((hl.is_defined(tc.lof) & (tc.lof == "HC")) | tc.consequence_terms.any(lambda c: consequences_to_keep.contains(c)))
                ).map(
                    lambda x: x.select(
                        'transcript_id',
                        lof=hl.or_missing(x.lof == "HC", hl.struct()),
                        missense=hl.or_missing(hl.is_missing(x.lof) & x.consequence_terms.any(lambda c: missense_consequences.contains(c)), hl.struct()),
                        canonical=hl.or_missing(x.canonical == 1, hl.struct())
                    )
                )
            ),
            af=hl.float32(hl.find( lambda x: x.meta.get('group') == 'adj', freq[mt.row_key].freq).AF[1]),
            was_split=hl.or_missing(mt.was_split, hl.struct())
        )

        mt = mt.filter_rows(hl.is_defined(mt.vep) & (hl.len(mt.vep) > 0) & (mt.af > 0) & (mt.af <= 0.01))
        mt = mt.explode_rows(mt.vep)

        mt = mt.transmute_rows(**mt.vep)
        mt = mt.select_cols('pop', 'fam_id', 'is_female')
        mt = mt.select_entries('PBT_GT', 'adj')
        mt.describe()

        mt = create_grouped_sparse_genotypes_matrix_table(
            mt, col_groups=['pop'], row_groups=['transcript_id'],
            call='PBT_GT',
            tmp_dir='gs://gnomad-tmp/pbt_parents_sparse_genotype_matrix'
        )
        write_temp_gcs(mt, pbt_phased_parents_sparse_genotype_matrix_path(data_type), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--phase', help='Phase gnomAD by transmission.',
                        action='store_true')
    parser.add_argument('--create_sparse_matrix', help='Creates per-gene sparse matrix with the phased trios.',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)




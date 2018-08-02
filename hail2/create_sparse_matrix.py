from gnomad_hail import *
from gnomad_hail.resources.phasing import *
import hail as hl
from gnomad_hail.utils.constants import CSQ_CODING_HIGH_IMPACT, CSQ_CODING_MEDIUM_IMPACT, POLYPHEN_PREDICTIONS, SIFT_PREDICTIONS



def main(args):

    data_type = 'exomes' if args.exomes else 'genomes'
    
    mt = get_gnomad_data(data_type, release_samples=True)
    mt = mt.select_cols(pop=mt.meta.pop)
    mt = mt.select_rows('a_index')

    vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep'))
    cadd = hl.read_table(cadd_ht_path)
    freq = hl.read_table(annotations_ht_path(data_type, 'frequencies'))

    consequences_to_keep = hl.literal(set(["inframe_insertion", "inframe_deletion", "missense_variant",]))
    # consequences_index = hl.literal({pred: i for i, pred in enumerate(CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT)})
    # polyphen_pred_index = hl.literal({pred: i for i, pred in enumerate(POLYPHEN_PREDICTIONS)})
    # sift_pred_index = hl.literal({pred: i for i, pred in enumerate(SIFT_PREDICTIONS)})

    mt = mt.annotate_globals(
        vep_consequences_to_keep=consequences_to_keep
        # polyphen_pred=POLYPHEN_PREDICTIONS,
        # sift_pred=SIFT_PREDICTIONS
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
                    #lof=hl.cond(hl.is_defined(x.lof), hl.struct(), hl.null(hl.tstruct())),
                    # polyphen_prediction=polyphen_pred_index[x.polyphen_prediction],
                    # sift_prediction=sift_pred_index[x.sift_prediction],
                    canonical=hl.or_missing(x.canonical == 1, hl.struct())
                )
            )
        ),
        cadd=hl.int(cadd[mt.row_key].phred_score),
        af=hl.float32(hl.find(lambda x: x.meta.get('group') == 'adj', freq[mt.row_key].freq).AF[1])
    )

    mt = mt.filter_rows(hl.is_defined(mt.vep) & (hl.len(mt.vep) > 0) & (mt.af > 0) & (mt.af <= 0.01))
    mt = mt.explode_rows(mt.vep)

    mt = mt.transmute_rows(**mt.vep)
    mt.describe()
    mt = create_grouped_sparse_genotypes_matrix_table(mt, ['transcript_id'], ['pop'])
    write_temp_gcs(grouped_sparse_genotype_matrix_path(data_type), overwrite=args.overwrite)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
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




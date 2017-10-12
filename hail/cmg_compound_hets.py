from utils import *
from resources import *
from hail import *
import argparse
from compound_hets_utils import *

def get_gnomad_phased_kt(kt, data_type, all_variants_vds, va_to_add=None, split_by_pop=False):
    """

    Adds phase, counts and annotations based on gnomAD to a KeyTable of variant pairs.

    :param KeyTable kt: Input KeyTable
    :param str data_type: Which gnomAD data to annotate with. Either "exomes" or "genomes".
    :param VariantDataset all_variants_vds: A VDS containing all variants in the KeyTable along with a `va.gene` annotation to aggregate on
    :param va_to_add: variant annotations to keep from gnomAD release in the form: destination_col: expr
    :param bool split_by_pop: Whether to split EM calculations by population or not
    :return: KeyTable with annotations
    :rtype: KeyTable
    """

    variants_list = list(all_variants_vds.query_variants('variants.collectAsSet()'))
    gnomad = get_gnomad_data(kt.hc, data_type, hardcalls="raw", split=True)
    gnomad = gnomad.filter_variants_list(variants_list)
    gnomad = gnomad.annotate_variants_vds(all_variants_vds, expr="va.gene = vds.gene")

    gnomad = gnomad.annotate_variants_expr('va = select(va, gene)')

    agg_keys = {'v1': 'v1', 'v2': 'v2'}
    chet_annotations = ['genotype_counts', 'haplotype_counts', 'prob_same_haplotype']
    if split_by_pop:
        gnomad = gnomad.annotate_samples_expr('sa = {pop : sa.meta.population}')
        agg_keys['pop'] = '`sa.pop`'
        chet_annotations.append('pop')
    else:
        gnomad = gnomad.annotate_samples_expr('sa = {}')

    gnomad_kt = gnomad.phase_em(['va.gene'], 1,
                                      sa_keys='sa.pop' if split_by_pop else None,
                                      variant_pairs= kt)

    gnomad_kt = aggregate_em_ref(gnomad_kt, agg_keys)

    gnomad_kt = kt.key_by(['v1','v2']).join(gnomad_kt.key_by(['v1','v2']), how="left")
    gnomad_kt = gnomad_kt.rename({
        ann: '{0}_{1}'.format("ge" if data_type == "exomes" else "gg", ann) for ann in chet_annotations
    })

    if va_to_add is not None:
        release = get_gnomad_public_data(kt.hc, data_type, split=True)
        release = release.filter_variants_list(variants_list)
        release_kt = release.variants_table()
        release_kt = release_kt.annotate(['{0} = {1}'.format(ann, expr) for ann,expr in va_to_add.iteritems()])
        release_kt = release_kt.select(['v'] + [ann for ann in va_to_add.keys()])

        for v in ['1', '2']:
            gnomad_kt = gnomad_kt.key_by('v{}'.format(v)).join(
                release_kt.rename({ann: '{0}{1}'.format(ann, v) for ann in va_to_add.keys()})
                , how="left")


    return gnomad_kt


def pre_process_data(kt, variant_expr, gene_expr, fam_expr, other_expr = None):
    """

    Pre-process the data from a KeyTable of family/variants with 2 lines per chet variant to a list of fam/variant pairs

    :param KeyTable kt: KeyTable
    :param str variant_expr: Expression to get a Variant from the KeyTable (e.g. Variant(chrom, pos, ref, alt)
    :param str gene_expr: Expression to get the gene from the KeyTable
    :param str fam_expr: Expression to get the family from the KeyTable
    :param dict of str:str other_expr: Any other per-variant name: expressions to be collected and kept in the resulting KeyTable (e.g. impact: effect => impact1 and impact2)
    :return: KeyTable of variant pairs
    :rtype: KeyTable

    """

    if other_expr is None:
        other_expr = {}

    kt = kt.annotate(['v = {}'.format(variant_expr),
                      'fam = {}'.format(fam_expr),
                      'gene = {}'.format(gene_expr)])
    all_variants_vds = VariantDataset.from_table(kt.select(['v', 'gene']).key_by(['v']))

    kt = kt.aggregate_by_key(['fam = fam', 'gene = gene'],
                             ['variants = v.collect()'] +
                             ['{0}s =  {1}.collect()'.format(ann, old_ann) for ann, old_ann in other_expr.iteritems()])

    kt = kt.filter('variants.length == 2')

    kt = kt.annotate(['v1 = variants[0]',
                      'v2 = variants[1]'] +
                      ['{0}{1} = {0}s[{2}]'.format(ann, num, num-1) for num in [1,2] for ann in other_expr.keys()])
    return kt, all_variants_vds


def main(args):
    hc = HailContext()
    kt = hc.import_table(args.variants_file, impute=True)

    if args.cmg:
        kt, all_variants_vds = pre_process_data(kt,
                                variant_expr='Variant(gDNA.replace("[ >]",":").replace("^chr",""))',
                                gene_expr='Gene',
                                fam_expr='`Family ID`',
                                other_expr={'impact' : '`Variant type`',
                                            'poo': 'POO'})
    else:
        kt, all_variants_vds = pre_process_data(kt,
                                variant_expr='Variant(chrom,pos,ref,alt)',
                                gene_expr='gene',
                                fam_expr='family',
                                other_expr={'impact': 'effect',
                                            'child_gt': 'GT_genotype_0',
                                            'parent_1': 'GT_genotype_0'})

    va_to_add = {'{}_ac_raw': 'va.info.AC_raw',
                 '{}_an_raw': 'va.info.AN_raw',
                 '{}_ac': 'va.info.AC',
                 '{}_an': 'va.info.AN'
                 }

    kt = get_gnomad_phased_kt(kt.key_by(['v1','v2']),
                              "exomes",
                              all_variants_vds,
                              { k.format("ge"):v for k,v in va_to_add.iteritems()},
                              args.split_by_pop)
    kt = get_gnomad_phased_kt(kt.key_by(['v1','v2']),
                              "genomes",
                              all_variants_vds,
                              {k.format("gg"): v for k, v in va_to_add.iteritems()},
                              args.split_by_pop)

    kt.write(args.output + ".kt", overwrite=True)

    kt = hc.read_table(args.output + ".kt")

    kt, count_cols_ge = flatten_counts(kt, gc_ann="ge_genotype_counts", hc_ann="ge_haplotype_counts", out_prefix="ge_")
    kt, count_cols_gg = flatten_counts(kt, gc_ann="gg_genotype_counts", hc_ann="gg_haplotype_counts", out_prefix="gg_")

    kt.export(args.output +  ".kt.tsv.gz")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmg', help='File is in CMG data format.', required=False, action='store_true')
    parser.add_argument('--seqr', help='File is in SEQR data format.', required=False, action='store_true')
    parser.add_argument('--variants_file', help='File containing variants.', required=True)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--split_by_pop', help='Splits data by population when computing EM', required=False, action='store_true')
    parser.add_argument('--output', help='Output prefix', required=True)
    args = parser.parse_args()

    if int(args.cmg) + int(args.seqr) != 1:
        sys.exit("One and only one of --cmg or --seqr should be specified.")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)





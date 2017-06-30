
from utils import *
import argparse


def write_hardcalls(vds, sample_group_filters, output, fam_file = None, overwrite = False, medians = True):
    """

    Writes multi-allelic hardcalls with the following annotations:
     - AC_unrelated
     - variantType
     - nAltAllele
     - Allele-specific statistics and callstats for sample groups (using raw genotypes)
     - Omni, HapMap, 1KG high conf SNVs, Mills

    :param VariantDataSet vds: Full VDS
    :param dict of str:str sample_group_filters: Mapping of sample group name -> expression. E.g. { "release_samples_raw": "sa.meta.keep" }
    :param str output: Output VDS path
    :param str fam_file: Fam pedigree file location
    :param bool overwrite: Whether to overwrite the output if present
    :param bool medians: Whether to compute the allele-specific stats medians
    """

    hc = vds.hc

    if fam_file is not None:
        vds = vds.annotate_samples_table(hail.KeyTable.import_fam(fam_file), root='sa.fam')

    allele_annotations = [
        "va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.nNonRefAlleles).sum()"]
    for group, filter_expr in sample_group_filters.iteritems():
        allele_annotations.extend(
            get_allele_stats_expr("va.stats.%s" % group, medians=medians, samples_filter_expr=filter_expr))

    variant_annotations = [get_variant_type_expr()]
    variant_annotations.append('va.nAltAlleles = v.altAlleles.filter(a => !a.isStar).length')
    for group, filter_expr in sample_group_filters.iteritems():
        if filter_expr:
            filter_expr = '.filter(g => %s)' % filter_expr

        variant_annotations.extend(["va.calldata.%s = gs%s.callStats(g => v)" % (group, filter_expr),
                                    "va.stats.%s.qd = (va.stats.%s.qual / va.stats.%s.nrdp).map(x => [35,x].min)" % (
                                    group, group, group)])

    vds = (
        vds.annotate_variants_vds(hc.read(hapmap_vds_path), expr='va.hapmap = isDefined(vds)')
            .annotate_variants_vds(hc.read(omni_vds_path), expr='va.omni = isDefined(vds)')
            .annotate_variants_vds(hc.read(mills_vds_path), expr='va.mills = isDefined(vds)')
            .annotate_variants_vds(hc.read(kgp_high_conf_snvs_vds_path), 'va.kgp_high_conf = isDefined(vds)')
            .annotate_alleles_expr(allele_annotations)
            .annotate_variants_expr(variant_annotations)
    )
    if args.sites_only:
        vds = vds.drop_samples()
    else:
        vds = vds.hardcalls()

    vds.write(output, overwrite=overwrite)


def write_split_hardcalls(hardcalls_vds, sample_group_filters, output, fam_file = None, overwrite = False, medians = True):
    """
    Takes multi-allelic hardcalls as input and writes the split version, splitting multi-allelic annotations.

    :param VariantDataset hardcalls_vds: Multi-allelic hardcalls VDS
    :param dict of str:str sample_group_filters: Mapping of sample group name -> expression. E.g. { "release_samples_raw": "sa.meta.keep" }
    :param str output: Output VDS path
    :param str fam_file: Fam pedigree file location
    :param bool overwrite: Whether to overwrite the output if present
    :param bool medians: Whether to compute the allele-specific stats medians
    """

    a_ann = ['va.AC_unrelated']
    r_ann = []
    for group in sample_group_filters.keys():
        a_ann.extend([x.format(group) for x in [
            'va.stats.{0}.gq',
            'va.stats.{0}.dp',
            'va.stats.{0}.nrq',
            'va.stats.{0}.ab',
            'va.stats.{0}.best_ab',
            'va.stats.{0}.pab',
            'va.stats.{0}.nrdp',
            'va.stats.{0}.combined_pAB',
            'va.stats.{0}.qual',
            'va.stats.{0}.qd'

        ]])
        r_ann.extend([x.format(group) for x in [
            'va.calldata.{0}.AC',
            'va.calldata.{0}.AF'
        ]])
        if medians:
            a_ann.extend([x.format(group) for x in [
                'va.stats.{0}.gq_median',
                'va.stats.{0}.dp_median',
                'va.stats.{0}.nrq_median',
                'va.stats.{0}.ab_median',
                'va.stats.{0}.pab_median']])

    vds = (
        hardcalls_vds
            .annotate_variants_expr(['va.nonsplit_alleles = v.altAlleles.map(a => a.alt)',
                                     'va.hasStar = v.altAlleles.exists(a => a.isStar)'])
            .split_multi()
            .annotate_variants_expr([
            'va.wasMixed = va.variantType == "mixed"',
            'va.alleleType = if(v.altAllele.isSNP) "snv"'
            '   else if(v.altAllele.isInsertion) "ins"'
            '   else if(v.altAllele.isDeletion) "del"'
            '   else "complex"'])
            .annotate_variants_expr(index_into_arrays(a_based_annotations=a_ann, r_based_annotations=r_ann, drop_ref_ann=True))
    )

    if fam_file is not None:
        vds = vds.tdt(hail.Pedigree.read(fam_file))

    vds.write(output, overwrite=overwrite)


def main(args):

    hc = hail.HailContext(log='/variantqc.log')

    if args.debug:
        logger.setLevel(logging.DEBUG)

    hardcalls_path = args.output + ".raw_hardcalls.vds"
    hardcalls_split_path =  args.output + ".raw_hardcalls.split.vds"

    if args.genomes:
        sample_group_filters = {"all_samples_raw": '',
                                "qc_samples_raw": 'sa.meta.qc_sample || (sa.in_exomes && sa.qc_pass)',
                                "release_samples_raw": 'sa.meta.keep'
                                }
        fam_file = genomes_fam_path
    else:
        sample_group_filters = {"all_samples_raw": '',
                                "qc_samples_raw": 'sa.meta.drop_status == "keep" || '
                                                 '(!isMissing(sa.fam.famID) && !("hard" ~ sa.meta.drop_condense)) || '
                                                 's == "C1975::NA12878" || s == "CHMI_CHMI3_Nex1" || (sa.in_genomes && sa.qc_pass)',
                                "release_samples_raw": 'sa.meta.drop_status == "keep"'
                                }
        fam_file = exomes_fam_path

    # Create hardcalls file with raw annotations
    if args.write_hardcalls:
        vds = add_exomes_sa(hc.read(full_exome_vds_path)) if args.exomes else add_genomes_sa(hc.read(full_genome_vds_path))
        write_hardcalls(vds, sample_group_filters, hardcalls_path, fam_file=fam_file, overwrite=args.overwrite,
                        medians=True)

    if args.write_split_hardcalls:
        hardcalls_vds = hc.read(hardcalls_path)
        write_split_hardcalls(hardcalls_vds, sample_group_filters, hardcalls_split_path, fam_file=fam_file,
                              overwrite=args.overwrite,
                              medians=True)

    if args.slack_channel:
        send_message(args.slack_channel, 'Create hardcalls %s is done processing!' % args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--sites_only', help='Creates a sites vds', action='store_true')
    parser.add_argument('--write_hardcalls', help='Creates a hardcalls vds', action='store_true')
    parser.add_argument('--write_split_hardcalls', help='Creates a split hardcalls vds from the hardcalls vds', action='store_true')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--output', '-o', help='Output prefix', required=True)
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    args = parser.parse_args()

    if args.write_hardcalls:
        if int(args.exomes) + int(args.genomes) != 1:
            sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    main(args)


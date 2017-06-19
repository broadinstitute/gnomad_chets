import hail

hc = hail.HailContext()
#vds = hc.import_vcf("/Users/laurent/tools/hail/src/test/resources/mendel.vcf")
vds = hc.read("/Users/laurent/projects/gnomad/one_chunk/gnomAD.variant_filtered.1415.1k_variants.vds")
print("Read VDS gnomAD.variant_filtered.1415.1k_variants.vds.")
print(vds.variant_schema)

def test_annotation(vds, ann, n = 10):
    vds = vds.annotate_variants_expr([ann])
    print(vds.query_variants('variants.map(v => %s).take(%d)[0]' % (ann.split("=")[0],n))[0])
    return(vds)


def test_qd(vds, out = "/Users/laurent/tmp/qd.out.txt"):
    # exp = ['va.nrdp = gs.filter(g => g.isCalledNonRef).map(g => g.dp).sum()',
    #                                                  # 'va.nrdp2 = gs.filter(g => g.isCalledNonRef || ( g.dp < 15 && g.gq < 6)).map(g => g.dp).sum()',
    #                                                 # 'va.qual2 = -10*gs.map(g => if(g.pl[0] > 3000) 3000 else log10(g.gp[0])).sum()',
    #                                                 # 'va.qual3 = -10*gs.filter(g => !(g.isHomRef && g.gq < 10 && g.dp < 20 && g.ad[1] == 0)).map(g => if(g.pl[0] > 3000) 3000 else log10(g.gp[0])).sum()',
    #                                                  'va.qual4 = -10*gs.filter(g => g.isCalledNonRef).map(g => if(g.pl[0] > 3000) -300 else log10(g.gp[0])).sum()',
    #                                                 'va.combined_pAB = let hetSamples = gs.filter(g => g.isHet).map(g => log(g.pAB())).collect() in orMissing(!hetSamples.isEmpty, -10*log10(pchisqtail(-2*hetSamples.sum(),2*hetSamples.length)))'
    #                                               #   'va.nNonRef = gs.filter(g => g.isCalledNonRef).count()',
    #                                               #   'va.nSamples = gs.filter(g => g.isCalled).count()',
    #                                               #   'va.best_non_ref = gs.filter(g => g.isCalledNonRef).collect().sortBy(g => g.gq)[1:10]',
    #                                               #   'va.first_non_ref = gs.filter(g => g.isCalledNonRef).take(5)',
    #                                               #   'va.gp0_stats = gs.map(g => g.gp[0]).stats()',
    #                                               #   'va.max_pl = gs.map(g => g.pl.max).max()',
    #                                               #   'va.max_ref_pl = gs.filter(g => g.isHomRef).map(g => g.pl.max).max()',
    #                                               #   'va.max_non_ref_pl = gs.filter(g => g.isCalledNonRef).map(g => g.pl.max).max()',
    #                                               #   'va.min_pl = gs.map(g => g.pl.max).min()',
    #                                               #   'va.min_ref_pl = gs.filter(g => g.isHomRef).map(g => g.pl.min).min()',
    #                                               #   'va.min_non_ref_pl = gs.filter(g => g.isCalledNonRef).map(g => g.pl.min).min()'
    #
    #                                                 ]

    def get_new_allele_stats_expr(root="va.stats", samples_filter_expr=''):
        if samples_filter_expr:
            samples_filter_expr = "&& " + samples_filter_expr

        stats = [
            '%s.pab = gs.filter(g => g.isHet %s).map(g => g.pAB()).stats()',
            '%s.nrdp = gs.filter(g => g.isCalledNonRef %s).map(g => g.dp).sum()',
            '%s.qual = -10*gs.filter(g => g.isCalledNonRef %s).map(g => if(g.pl[0] > 3000) -300 else log10(g.gp[0])).sum()',
            '%s.combined_pAB = let hetSamples = gs.filter(g => g.isHet %s).map(g => log(g.pAB())).collect() in orMissing(!hetSamples.isEmpty, -10*log10(pchisqtail(-2*hetSamples.sum(),2*hetSamples.length)))',
            '%s.pab_median = gs.filter(g => g.isHet %s).map(g => g.pAB()).collect().median']

        stats_expr = [x % (root, samples_filter_expr) for x in stats]

        return stats_expr

    vds = vds.filter_multi().annotate_alleles_expr(get_new_allele_stats_expr())
    vds.export_variants(out,
                        'chrom = v.contig,'
                        'pos = v.start,'
                        'ref = v.ref,'
                        'alt = v.alt,'
                        'qual = va.qual,'
                 #       'qual2 = va.qual2,'
                 #       'qual3 = va.qual3,'
                        'qual2 = va.stats.qual[0],'
                        'qd = va.info.QD,'
                        'qd2 = [35,va.stats.qual[0] / va.stats.nrdp[0]].min,'
                        'dp = va.info.DP,'
                        'nrdp = va.stats.nrdp[0],'
                       # 'nrdp2 = va.nrdp2,'
                        'qd_prime = va.qual/va.stats.nrdp[0]'
                        # 'nNonRef = va.nNonRef,'
                        # 'best_non_ref= va.best_non_ref,'
                        # 'gp0_mean = va.gp0_stats.mean,'
                        # 'gp0_stdev = va.gp0_stats.stdev,'
                        # 'gp0_min = va.gp0_stats.min,'
                        # 'gp0_max = va.gp0_stats.max,'
                        # 'max_pl = va.max_pl,'
                        # 'max_ref_pl = va.max_ref_pl,'
                        # 'max_non_ref_pl = va.max_non_ref_pl,'
                        # 'min_pl = va.min_pl,'
                        # 'min_ref_pl = va.min_ref_pl,'
                        # 'min_non_ref_pl = va.min_non_ref_pl,'
                        # 'first_non_ref = va.first_non_ref'
                        )
    return vds

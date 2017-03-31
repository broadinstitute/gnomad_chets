from utils import *


def run_sanity_checks(vds,pops):

    queries = []
    metrics = ['AC','AN','Hom']

    #Filter counts
    queries.append('variants.map(v => va.filters).counter()')

    #Check that raw is always larger than adj
    queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                    '.exists(i => va.info.%s[i] > va.info.%s_raw[i])).count()' %(x,x) for x in metrics])
    queries.append('variants.filter(v => va.info.STAR_AC_raw > va.info.STAR_AC).count()')

    #Check that male + female == total
    queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                    '.exists(i => va.info.%s_Male[i] + va.info.%s_Female[i] != va.info.%s[i])).count()' %(x,x,x) for x in metrics])

    #Check that sum(pops) == total
    for metric in metrics:
        queries.extend(['variants.filter(v => range(v.nAltAlleles)'
                        '.exists(i => %s != va.info.%s[i])).count()' % ( " + ".join(["va.info.%s_%s" %(metric,pop) for pop in pops]), metric)])

    print queries

run_sanity_checks(None,["AFR","oth"])

# test_path = "/Users/laurent/tools/hail/src/test/resources/split_test.vcf"
# tmp = "/Users/laurent/tmp/test2.vds"
#
# hc = hail.HailContext()
#
# vds = hc.import_vcf(test_path)
# vds.write(tmp, overwrite=True)
#
# split = vds.split_multi()
# split = split.annotate_variants_expr("va.split_allele = v.alt")
# split = split.annotate_variants_expr('va.split_allele2 = v.alt + "blah" ')
#
# vds = annotate_non_split_from_split(hc, tmp, split, ["va.split_allele",'va.split_allele2'])
# print(vds.query_variants('variants.map(v => v.nAltAlleles == va.split_allele.length).counter()'))
# print(vds.query_variants('variants.map(v => v.altAlleles.map(a => a.alt)).collect()'))
# print(vds.query_variants('variants.map(v => va.split_allele).collect()'))
# print(vds.query_variants('variants.map(v => v.nAltAlleles == va.split_allele2.length).counter()'))
# print(vds.query_variants('variants.map(v => va.split_allele2).collect()'))







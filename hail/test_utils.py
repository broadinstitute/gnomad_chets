import hail

hc = hail.HailContext()
vds = hc.import_vcf("/Users/laurent/tools/hail/src/test/resources/mendel.vcf")

def test_annotation(vds, ann):
    vds = vds.annotate_variants_expr([ann])
    print(vds.query_variants('variants.map(v => %s).collect()[0]' % ann.split("=")[0])[0])
    return(vds)
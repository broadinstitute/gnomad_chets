from utils import *


def get_kt_sparse_vsm(vds, keys):
    kt = vds.filter_genotypes('!g.isHomRef').genotypes_table()
    kt = kt.annotate('stuff = {v:v, va:va, s:s, sa:sa, gt:g.gt}')  # Dump sa ?
    kt = kt.aggregate_by_key(keys, 'index(stuff.collect(), v)')


def annotate_methylation(vds):
    mkt = vds.hc.read_table(methylation_kt).select(['locus', 'MEAN'])
    vds = vds.annotate_variants_table(mkt, root='va.methylation.value')
    vds = vds.annotate_variants_expr(
        ['va.methylated_cpg = v.altAllele().isTransition() && va.methylation.value >= 0.25'])
    #logger.debug('Number of methylated CpGs: {}'.format(
    #    str(vds.query_variants(['variants.map(v => va.methylated_cpg).counter()'])[0])))
    return vds


def annotate_gene_impact(vds):
    vds = filter_vep_to_canonical_transcripts(vds)
    vds = process_consequences(vds)
    impact = {x: "high" for x in ["transcript_ablation",
                                  "splice_acceptor_variant",
                                  "splice_donor_variant",
                                  "stop_gained",
                                  "frameshift_variant",
                                  "stop_lost"]
              }

    impact.update({x: "medium" for x in ["start_lost",  # new in v81
                                         "initiator_codon_variant",  # deprecated
                                         "transcript_amplification",
                                         "inframe_insertion",
                                         "inframe_deletion",
                                         "missense_variant",
                                         "protein_altering_variant",  # new in v79
                                         "splice_region_variant",
                                         ]})

    impact.update({x: "low" for x in [
        "incomplete_terminal_codon_variant",
        "stop_retained_variant",
        "synonymous_variant"]})

    vds = vds.annotate_global('global.impact', impact, TDict(TString(), TString()))

    vds = vds.annotate_variants_expr([
        'va.impact = if("-LC" ~ va.vep.worst_csq_suffix ) "medium" else if (global.impact.contains(va.vep.worst_csq)) global.impact[va.vep.worst_csq] else NA:String',
        'va.gene = if("-HC" ~ va.vep.worst_csq_suffix) '
        '   va.vep.transcript_consequences.find(x => x.lof == "HC").gene_symbol '
        'else'
        '   va.vep.transcript_consequences.find(x => x.consequence_terms.toSet.contains(va.vep.worst_csq)).gene_symbol',
        'va.alleleType = if(v.altAllele.isSNP) "SNP" else "indel"'
    ])

    #logger.debug("Variant impact counts: {0}, Variants with gene annotations: {1}".format(*vds.query_variants(
    #    ['variants.map(v => va.impact).counter()', 'variants.map(v => isDefined(va.gene)).counter()'])))
    vds = vds.filter_variants_expr('isDefined(va.impact) && isDefined(va.gene)')
    return vds

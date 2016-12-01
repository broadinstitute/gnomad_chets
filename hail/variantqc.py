__author__ = 'konrad'

import random
from utils import *

hc = HailContext()

# Input files
root = 'file:///mnt/lustre/konradk/exac'
meta_path = '%s/super_meta.txt.bgz' % root
vds_path = 'file:///mnt/lustre/lfran/exac2/exacv2.vds'
v1_vds_path = 'file:///mnt/lustre/lfran/exac/exac_all.new.vds'

# QC meta
raw_fam_path = '%s/variantqc/exac2.qctrios.raw.fam' % root
fam_path = '%s/variantqc/exac2.qctrios.fam' % root
qc_sites_path = '%s/variantqc/qc_sites.vds' % root
truth_dir = '%s/variantqc/truth_sets' % root

# All hardcalls
raw_hardcallvds_path = '%s/hardcalls/v2/exacv2.raw.hardcalls.qc.vds' % root
hardcallvds_path = '%s/hardcalls/v2/exacv2.hardcalls.qc.vds' % root
v1_hardcallvds_path = '%s/hardcalls/v1/exacv1.hardcalls.qc.vds' % root

# All splitmulti hardcalls
raw_split_hardcallvds_path = '%s/hardcalls/v2/exacv2.raw.hardcalls.splitmulti.qc.vds' % root
split_hardcallvds_path = '%s/hardcalls/v2/exacv2.hardcalls.splitmulti.qc.vds' % root
v1_split_hardcallvds_path = '%s/hardcalls/v1/exacv1.hardcalls.splitmulti.qc.vds' % root


fam_training_proportion = 0.8


def prepare_files():
    CHROMS = map(str, range(1, 23))
    CHROMS.extend(['X', 'Y'])

    # Prepare some intervals files
    with open('%s/intervals/chrX.txt', 'w') as f:
        f.write('X:1-1000000000')
    with open('%s/intervals/chrY.txt', 'w') as f:
        f.write('Y:1-1000000000')
    with open('%s/intervals/autosomes.txt', 'w') as f:
        f.write('\n'.join(['%s:1-1000000000' % x for x in CHROMS]))

    local_fam_path = fam_path.replace('file://', '')
    with open(local_fam_path) as f:
        g = open(local_fam_path.replace('.fam', '.train.fam'), 'w')
        h = open(local_fam_path.replace('.fam', '.test.fam'), 'w')
        for line in f:
            if random.random() < fam_training_proportion:
                g.write(line)
            else:
                h.write(line)
        g.close()
        h.close()


def write_hardcalls(input_path, output_path, adj=True, partitions=10000):
    """
    hail read -i $v1vds \
    filtergenotypes -c 'g.gq>=20 && g.dp>=10' --keep \
    filtervariants expr -c 'gs.filter(g => g.isCalledNonRef).count() > 0' --keep \
    hardcalls \
    repartition -n 10000 \
    write -o $v1hardcallvds
    """
    if adj:
        return (hc.read(input_path)
                .filter_genotypes('g.gq>=20 && g.dp>=10')
                .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')
                .hardcalls()
                .repartition(partitions)
                .write(output_path)
        )
    else:
        return (hc.read(input_path)
                .hardcalls()
                .repartition(partitions)
                .write(output_path)
        )


def write_split(input_path, output_path):
    return hc.read(input_path)\
        .annotate_variants_expr('va.isMixed = v.nAltAlleles > 1 && v.altAlleles.map(a => if(a.isSNP) 1 else 2).toSet.size > 1')\
        .split_multi().write(output_path)

# Prepare all metadata

# write_hardcalls(vds_path, raw_hardcallvds_path, adj=False)
# write_hardcalls(vds_path, hardcallvds_path)
# write_hardcalls(v1_vds_path, v1_hardcallvds_path)
# write_split(raw_hardcallvds_path, raw_split_hardcallvds_path)
# write_split(hardcallvds_path, split_hardcallvds_path)
# write_split(v1_hardcallvds_path, v1_split_hardcallvds_path)


# Raw QC for determing true positives
hardcallvds = hc.read(raw_split_hardcallvds_path)
output_path = '%s/variantqc/v2_tdt.raw.vds' % root

(hardcallvds.filter_variants_intervals('%s/intervals/autosomes.txt' % root)
 .tdt(fam_path)
 .write(output_path))


hardcallvds = hc.read(split_hardcallvds_path)
rf_vds = hc.read('%s/sites/exac2.sites.RF2.vds' % root)
output_path = '%s/variantqc/v2_variantqc.vds' % root

new_vds = (hardcallvds.filter_variants_intervals('%s/intervals/autosomes.txt' % root)
           .annotate_samples_table(meta_path, 'sample', root='sa.meta', impute=True)
           .annotate_samples_fam(fam_path, root='sa.fam')
           .filter_samples_expr('sa.meta.drop_status == "keep" || !isMissing(sa.fam.famID)')
           .annotate_variants_intervals('%s/intervals/exome_evaluation_regions.v1.intervals' % root, root='va.evaluation_interval')
           .annotate_variants_intervals('%s/intervals/high_coverage.auto.interval_list' % root, root='va.high_coverage_interval')
           .annotate_variants_table('%s/variantqc/v2.lmendel' % root, 'SNP', impute=True, root='va.mendel')
           .annotate_variants_table('%s/variantqc/validatedDN.cut.txt.bgz' % root, 'Variant(CHROM, POSITION.toInt, REF, ALT)', impute=True, code='va.validated_denovo = table.DataSet')
           .annotate_variants_vds(rf_vds, root='va.rf')
           .annotate_variants_expr('va.AC_unrelated = gs.filter(g => g.isCalledNonRef && isMissing(sa.fam.patID)).map(g => g.oneHotAlleles(v)).sum(),'
                                   'va.pass = va.filters.contains("PASS")')
           .tdt(fam_path)
           .filter_samples_expr('sa.meta.drop_status == "keep"')
           .variant_qc()
           .write(output_path)
)

columns = '''chrom = v.contig,
pos = v.start,
ref = v.ref,
alt = v.alt,
evaluation_interval = !isMissing(va.evaluation_interval),
high_coverage_interval = va.high_coverage_interval,
vqslod = va.info.VQSLOD,
type = va.rf.variantType,
rfprob_all = va.rf.RF.probability[2],
pass = va.pass,
is_mixed = va.isMixed,
qd = va.info.QD,
wassplit = va.wasSplit,
mendel_errors = va.mendel.N,
validated_denovo = va.validated_denovo,
ac_unrelated = va.AC_unrelated[1].toInt,
transmitted = va.tdt.nTransmitted,
untransmitted = va.tdt.nUntransmitted,
ac_orig = va.info.AC[va.aIndex - 1],
ab_mean = va.info.ab_stats[va.aIndex - 1].mean,
callrate = va.qc.callRate,
af = va.qc.AF,
ac = va.qc.AC'''

hc.read('%s/variantqc.vds' % root).export_variants('%s/variantqc.txt.bgz' % root, columns)
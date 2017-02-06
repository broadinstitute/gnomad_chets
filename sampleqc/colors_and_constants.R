color_amr = k_amr = '#ED1E24'
color_eur = k_eur = '#6AA5CD'
color_afr = k_afr = '#941494'
color_sas = k_sas = '#FF9912'
color_eas = k_eas = '#108C44'
color_oth = k_oth = '#ABB9B9'
color_mde = k_mde = '#33CC33'
color_asj = k_asj = 'coral'

color_nfe = k_nfe = color_eur
color_fin = k_fin = '#002F6C'

color_syn = k_syn = '#AAAAAA'
color_mis = k_mis = '#FF6103'
color_lof = k_lof = '#9D1309'

color_cpg = '#2E9FFE'
color_ti = '#458B00'
color_tv = '#EA4444'

lof_like = c('frameshift_variant','essential_splice','stop_gained')
mis_like = c('missense_variant','inframe_indel','stop_lost',
             'mature_miRNA_variant','start_lost')
syn_like = c('synonymous_variant','3_prime_UTR_variant','5_prime_UTR_variant',
             'extended_splice','stop_retained_variant','non_coding_transcript_exon_variant',
             'intron_variant','intergenic_variant','regulatory_region_variant')

format_vep_category = function(category_list) {
  return(category_list %>%
           gsub("_"," ", .) %>%
           gsub('stop gained', 'nonsense', .) %>%
           gsub("inframe indel", "in-frame indel", .) %>%
           gsub("initiator codon", "start lost", .) %>%
           gsub(" variant", "", .) %>%
           gsub("transcript exon", "transcript", .) %>%
           gsub(" prime ","'", .) %>%
           gsub("probably damaging", "prob damaging", .) %>%
           gsub("possibly damaging", "poss damaging", .))
}
variant_types = c('transversion', 'non-CpG transition', 'CpG transition')
variant_type_colors = c(color_tv, color_ti, color_cpg)

pops <- c('afr', 'amr', 'eas', 'fin', 'nfe', 'sas')
pop_colors = c('afr' = color_afr,
               'amr' = color_amr,
               'eas' = color_eas,
               'fin' = color_fin,
               'eur' = color_nfe,
               'nfe' = color_nfe,
               'oth' = color_oth,
               'sas' = color_sas,
               'mde' = color_mde,
               'asj' = color_asj,
               'uniform' = 'pink',
               'consanguineous' = 'pink',
               'sas_non_consang' = 'orange',
               'exac' = 'gray',
               'est' = 'black')
pop_names = c('afr' = 'African',
              'amr' = 'Latino',
              'eas' = 'East Asian',
              'fin' = 'Finnish',
              'nfe' = 'European',
              'oth' = 'Other',
              'sas' = 'South Asian',
              'mde' = 'Middle Eastern',
              'asj' = 'Ashkenazi Jewish',
              'uniform' = 'Uniform',
              'sas_non_consang' = 'South Asian (F < 0.05)',
              'consanguineous' = 'South Asian (F > 0.05)',
              'exac' = 'ExAC')

agilent = c("Exome Express", "Standard Exome Sequencing", "Standard Exome Sequencing v2",
            "WHOLE EXOME HYB SEL & SEQ", "HYB SEL & SEQ")
ice = c("Standard Exome Sequencing v3", "Standard Exome Sequencing v4", 
        "Nextera Exome", "G4L WES + Array v1", "Express Human WES (Standard Coverage) v1",
        "Exome Express v2", "Exome Express v3",
        "Express Human WES (Deep Coverage) v1")
ice150 = c("G4L WES + Array v2", "Standard Exome Sequencing v5") # ICE with 150bp reads
wgs = c('PCR Free', 'PCR Plus')

platform_alphas = c('agilent' = 0.8, 'ice' = 0.8, 'ice150' = 0.8, 'nimblegen' = 0.8, 'gnomAD' = 0.8, 'multiple' = 0.5, 'unknown' = 0.5)
platform_colors = c('agilent' = '#F8766D', 'ice' = '#00BA38', 'ice150' = '#619CFF', 'nimblegen' = '#B79F00', 'gnomAD' = '#F564E3', 'multiple' = '#F564E3', 'unknown' = 'gray')
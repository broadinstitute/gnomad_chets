
#### Color constants ####

color_syn = k_syn = '#AAAAAA'
color_mis = k_mis = '#FF6103'
color_lof = k_lof = '#9D1309'

var_type_aliases = c('syn' = 'Synonymous', 'mis' = 'Missense', 'lof' = 'pLoF', 'os' = 'Other splice')
colors = c('synonymous_variant' = color_syn,
           'missense_variant' = color_mis,
           'stop_gained' = color_lof,
           'Synonymous' = color_syn,
           'Missense' = color_mis,
           'synonymous' = color_syn,
           'missense' = color_mis,
           'nonsense' = color_lof,
           'Other splice' = color_os,
           'LoF' = color_lof,
           'pLoF' = color_lof)
mut_cols = c("transversion" = "#EA4444", "non-CpG transition" = "#458B00", "CpG transition" = '#2E9FFE',
             'ACG' = '#422efe', 'CCG' = '#2e73fe', 'GCG' = '#2ec3fe', 'TCG' = '#2efef7')
dataset_colors = c('ExAC' = '#4682B4', 'gnomAD' = '#73AB3D', 
                   'exomes' = '#4682B4', 'genomes' = '#73AB3D')

color_all = 'black'
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

color_cpg = '#2E9FFE'
color_ti = '#458B00'
color_tv = '#EA4444'

lof_like = c('frameshift_variant','essential_splice','stop_gained','splice_donor_variant','splice_acceptor_variant')
mis_like = c('missense_variant','inframe_indel','stop_lost',
             'mature_miRNA_variant','start_lost')
syn_like = c('synonymous_variant','3_prime_UTR_variant','5_prime_UTR_variant', 'splice_region_variant',
             'extended_splice','stop_retained_variant','non_coding_transcript_exon_variant',
             'upstream_gene_variant', 'downstream_gene_variant',
             'intron_variant','intergenic_variant','regulatory_region_variant')

variant_category_colors = c(rep(color_lof, length(lof_like)),
                            rep(color_mis, length(mis_like)),
                            rep(color_syn, length(syn_like)))
names(variant_category_colors) = c(format_vep_category(lof_like),
                                   format_vep_category(mis_like),
                                   format_vep_category(syn_like))

variant_type_colors = c(color_tv, color_ti, color_cpg, color_cpg)
names(variant_type_colors) = c('transversion', 'non-CpG transition', 'CpG transition', 'CpG')

pops <- c('afr', 'amr', 'eas', 'fin', 'nfe', 'sas')
pop_colors = c(
  'all' = color_all,
  'afr' = color_afr,
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
   'est' = 'black',
   'bgr' = '#66C2A5',
   'est' = '#4891D9',
   'nwe' = '#C60C30',
   'seu' = '#3ca021',
   'swe' = 'purple',
   'onf' = color_nfe,
   'kor' = '#4891D9',
   'jpn' = '#BC002D',
   'oea' = color_eas,
   'unk' = color_oth
  )
pop_names = c('oth' = 'Other',
              'all' = 'All',
              'afr' = 'African/African-American',
              'amr' = 'Latino',
              'eas' = 'East Asian',
              'fin' = 'Finnish',
              'eur' = 'European',
              'nfe' = 'European',
              'sas' = 'South Asian',
              'mde' = 'Middle Eastern',
              'asj' = 'Ashkenazi Jewish',
              'uniform' = 'Uniform',
              'sas_non_consang' = 'South Asian (F < 0.05)',
              'consanguineous' = 'South Asian (F > 0.05)',
              'exac' = 'ExAC',
              'bgr' = 'Bulgarian',
              'est' = 'Estonian',
              'nwe' = 'Northern European',
              'seu' = 'Southern European',
              'swe' = 'Swedish',
              'onf' = 'Other Non-Finnish European',
              'kor' = 'Korean',
              'jpn' = 'Japanese',
              'oea' = 'Other East Asian',
              'unk' = 'Unknown')

phase_colors = c(
  'Compound heterozygous' = 'blue',
  'Same haplotype' = 'red',
  'NA' = 'gray'
)

#### Functions ####
get_pr_data <- function(pbt, metric, min_ac = 2, nbins=100){
  
  get_pr_data_per_phase <- function(pbt, metric, nbins, chet){
    
    pbt_ann = pbt %>%
      rename(
        score = !!metric
        ) %>%
      mutate(
        phase_truth = ifelse(rep(chet, dim(pbt)[1]), trio_chet, !trio_chet),
        score = ifelse(rep(chet, dim(pbt)[1]), score, 1-score)
      )
    
    # This total includes variant pairs with 1 or both variants missing from gnomAD
    total_tp = pbt_ann %>%
      filter(!is.na(phase_truth)) %>%
      group_by(pop) %>%
      summarise(
        n = sum(phase_truth)
      ) %$%
      setNames(n, pop)
    
    score_binned_chet = pbt_ann %>%
      filter(!is.na(phase_truth)) %>%
      filter(n_var_gnomad == 2) %>%
      filter(AC1 >= min_ac & AC2 >= min_ac) %>%
      group_by(pop)  %>%
      arrange(desc(score)) %>% 
      mutate(
        cum_tp = cumsum(phase_truth),
        cum_fp = cumsum(!phase_truth),
        bin = ntile(score, nbins)
      ) %>%
      group_by(pop, bin) %>%
      summarise(
        score = max(score),
        cum_tp = last(cum_tp),
        cum_fp = last(cum_fp)
      ) %>%
      mutate(
        precision = cum_tp /(cum_tp + cum_fp),
        recall = cum_tp / total_tp[pop]
      )
  }
  
  return(
    get_pr_data_per_phase(pbt, metric, nbins, T) %>%
    mutate(
      phase="Compound heterozygous"
    ) %>%
    bind_rows(
      get_pr_data_per_phase(pbt, metric, nbins, F) %>%
        mutate(
          phase = "Same haplotype"
        )
    ) %>%
      mutate(metric = metric)
  )
}

get_pr_plot <- function(pbt, metrics, min_ac=2, nbins=100){
  
  print(sprintf("Computing PR data for metric: %s", metrics[1]))
  plot_data = get_pr_data(pbt, metrics[1], min_ac, nbins) %>%
    mutate(metric=metrics[1])
  
  for (metric in metrics[-1]){
    print(sprintf("Computing PR data for metric: %s", metric))
    plot_data %<>%
      bind_rows(
        get_pr_data(pbt, metric, min_ac, nbins) %>%
          mutate(metric=metric)
      )
  }
  
  return(
    plot_data %>%
      ggplot(aes(recall, precision, col=pop)) + 
      scale_color_manual(values=pop_colors) + 
      geom_point() + 
      facet_grid(phase ~ metric)
  )
}

get_pbt_binned_metrics <- function(pbt, metrics, nbins=100){
  
  get_pbt_binned_metric <- function(pbt, metric, nbins){
    return(
      pbt %>%
        rename(metric=!!metric) %>%
        filter(!is.na(metric)) %>%
        mutate(
          bin = ntile(metric, nbins)
        ) %>%
        group_by(bin) %>%
        mutate(
          score=min(metric)
        ) %>%
        summarise(
          score=min(score),
          n=n(),
          n_trio_chet=sum(trio_chet),
          prop_trio_chet=n_trio_chet/n
        ) %>%
        arrange(
          desc(bin)
        ) %>%
        mutate(
          cumul_n_trio_chet=cumsum(n_trio_chet),
          cumul_prop_trio_chet=cumul_n_trio_chet/cumsum(n),
          metric=metric
          )
    )
  }
  
  binned_metrics = get_pbt_binned_metric(pbt, metrics[1], nbins) 
  
  for(metric in metrics[-1]){
    binned_metrics %<>%
      bind_rows(
        get_pbt_binned_metric(pbt, metric, nbins) 
      )
  }
  return(binned_metrics)
}

plot_pbt_binned_metrics <- function(pbt_binned_metrics){
  decile_denum = max(pbt_binned_metrics$bin) / 10
  return(
    pbt_binned_metrics %>%
      mutate(decile = as.factor(ceiling(bin / decile_denum))) %>%
      ggplot(aes(score, prop_trio_chet, col=decile)) + 
      geom_point() + 
      geom_abline(slope=1, intercept = 0) + 
      facet_grid(.~metric) + 
      scale_y_continuous(breaks = seq(0,1,0.1))
  )
}

plot_cumul_pbt_binned_metrics <- function(pbt_binned_metrics){
  decile_denum = max(pbt_binned_metrics$bin) / 10
  return(
    pbt_binned_metrics %>%
      mutate(decile = as.factor(ceiling(bin / decile_denum))) %>%
      ggplot(aes(score, cumul_prop_trio_chet, col=decile)) + geom_point() + facet_grid(.~metric)
  )
}

plot_binned_metrics_with_threshold <- function(pbt_binned_metrics, chet_threshold, same_hap_threshold){
  return(
    pbt_binned_metrics %>%
    left_join(
      chet_threshold %>%
        select(metric, chet_threshold)
    ) %>% 
    left_join(
      same_hap_threshold %>%
        select(metric, same_hap_threshold)
    ) %>%
    mutate(
      phase = factor(
        case_when(
        score < same_hap_threshold ~ 'Same haplotype', 
        score > chet_threshold ~ "Compound heterozygous",
        T ~ 'NA'
      ),
      levels=c("Compound heterozygous", 'Same haplotype', 'NA')
      )
    ) %>% 
    ggplot(aes(score, prop_trio_chet, col=phase)) + 
      geom_point() + 
      scale_color_manual(values=phase_colors) + 
      facet_grid(.~metric) + 
      theme(
        legend.position = "top",
        legend.title = element_blank()
        )
  )
}

add_ac_bin = function(pbt, merge_symetric_bins = TRUE){
  
  return(
    pbt %>%
      mutate(
        .AC1=ifelse(merge_symetric_bins & AC1 > AC2, AC2, AC1),
        .AC2=ifelse(merge_symetric_bins & AC1 > AC2, AC1, AC2)
      ) %>%
      mutate(
        ac1_bin = factor(
          case_when(
            .AC1 == 0 ~ 'AC = 0',
            .AC1 == 1 ~ 'AC = 1',
            .AC1 == 2 ~ 'AC = 2',
            .AC1 < 5 ~ '2 < AC < 5',
            .AC1 < 10 ~ '5 ≤ AC < 10',
            .AC1 < 100 ~ '10 ≤ AC < 100',
            .AC1 < 1000 ~ '100 ≤ AC < 1000',
            .AC1 >= 1000 ~ 'AC ≥ 1000'
          ),
          levels=c(
            'AC = 0',
            'AC = 1',
            'AC = 2',
            '2 < AC < 5',
            '5 ≤ AC < 10',
            '10 ≤ AC < 100',
            '100 ≤ AC < 1000',
            'AC ≥ 1000'
          )
        ),
        ac2_bin = factor(
          case_when(
            .AC2 == 0 ~ 'AC = 0',
            .AC2 == 1 ~ 'AC = 1',
            .AC2 == 2 ~ 'AC = 2',
            .AC2 < 5 ~ '2 < AC < 5',
            .AC2 < 10 ~ '5 ≤ AC < 10',
            .AC2 < 100 ~ '10 ≤ AC < 100',
            .AC2 < 1000 ~ '100 ≤ AC < 1000',
            .AC2 >= 1000 ~ 'AC ≥ 1000'
          ),
          levels=c(
            'AC = 0',
            'AC = 1',
            'AC = 2',
            '2 < AC < 5',
            '5 ≤ AC < 10',
            '10 ≤ AC < 100',
            '100 ≤ AC < 1000',
            'AC ≥ 1000'
          )
        )
      ) %>%
      select(-.AC1, -.AC2)
  )
}

get_pbt_binned_by_ac = function(pbt, merge_symetric_bins = TRUE){
  return(
    pbt %>%
    filter(pop != 'all') %>%
    add_ac_bin(merge_symetric_bins) %>%
    gather(
      key = "metric",
      value = "score",
      em,
      rf_prediction
    ) %>%
    left_join(
      same_hap_threshold %>%
        select(metric, same_hap_threshold)
    ) %>%
    left_join(
      chet_threshold %>%
        select(metric, chet_threshold)
    ) %>%
    group_by(
      metric,
      # pop,
      ac1_bin,
      ac2_bin
    ) %>%
    summarise(
      n = n(),
      n_chet = sum(trio_chet),
      n_same_hap = sum(!trio_chet),
      prop_trio_chet = sum(trio_chet) / n(),
      chet_precision = sum(trio_chet & score > chet_threshold) / sum(score > chet_threshold),
      chet_recall = sum(trio_chet & score > chet_threshold) / sum(trio_chet),
      same_hap_precision = sum(!trio_chet & score < same_hap_threshold) / sum(score < same_hap_threshold),
      same_hap_recall = sum(!trio_chet & score < same_hap_threshold) / sum(!trio_chet)
    )
  )
}


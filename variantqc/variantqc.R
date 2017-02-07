options(stringsAsFactors = F)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pracma)
# source('exac_qc_constants.R')

ggplot2 = function(...) {
  return(ggplot(...) + theme_bw() + 
           theme(panel.border=element_blank(), 
                 legend.key = element_blank(),
                 text = element_text(size=20))
           )
}

read_variant_data = function(fname='', bins=100) {
  if (fname == '') {
    variant_data = read.table(gzfile('data/exome_variantqc.txt.bgz', 'r'), header=T)
  } else {
    variant_data = read.table(fname, header=T)
  }
  
  variant_data$pos_id = paste(variant_data$chrom, formatC(variant_data$pos,width=9,flag='0'), variant_data$ref, variant_data$alt, sep='_')
  variant_data$indel = nchar(variant_data$ref) + nchar(variant_data$alt) != 2
  variant_data$transition = (variant_data$ref == 'A' & variant_data$alt == 'G') | 
    (variant_data$ref == 'G' & variant_data$alt == 'A') | 
    (variant_data$ref == 'C' & variant_data$alt == 'T') | 
    (variant_data$ref == 'T' & variant_data$alt == 'C')
  variant_data$transition[variant_data$indel] = NA
  variant_data$bases_inserted = nchar(variant_data$alt) - nchar(variant_data$ref)
  variant_data$bases_inserted[!variant_data$indel] = NA
  variant_data$inframe = variant_data$bases_inserted %% 3 == 0
  variant_data$insertion = variant_data$bases_inserted > 0
  
  variant_data$evaluation_interval = as.logical(variant_data$evaluation_interval)
  variant_data$high_coverage_interval = as.logical(variant_data$high_coverage_interval)
  variant_data$wassplit = as.logical(variant_data$wassplit)
  variant_data$pass = as.logical(variant_data$pass)
  variant_data$chrpos = paste(variant_data$chrom, variant_data$pos)
  variant_data$model = ifelse(variant_data$indel | variant_data$type == 'mixed', 'indel', 'snp')
  
  variant_data = variant_data %>% filter(alt != '*' & !(chrom %in% c('X', 'Y')) & !is.na(vqslod))
  
  variant_data$vqslod_qd = variant_data$vqslod
  variant_data$vqslod_qd[variant_data$qd < 4] = min(variant_data$vqslod)
  
  if ('rfprob' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(rf_cut = bins - ntile(rfprob, bins))
  }
  if ('rfprob_all' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(rf_all_cut = bins - ntile(rfprob_all, bins))
  }
  if ('rfprob_bal' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(rf_bal_cut = bins - ntile(rfprob_bal, bins))
  }
  if ('ab_mean' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(ab_cut = bins - ntile(ab_mean, bins))
  }
  
  if ('an_qc_raw' %in% colnames(variant_data) & !('callrate' %in% colnames(variant_data))) {
    variant_data$callrate = variant_data$an_qc_raw/max(variant_data$an_qc_raw)
  }
  
  variant_data %>% group_by(model) %>% 
    mutate(vqs_cut = bins - ntile(vqslod, bins),
           vqs_qd_cut = bins - ntile(vqslod_qd, bins)) %>% ungroup
}

get_biallelic_data = function(variant_data) {
  allelic_state = plyr::count(subset(variant_data, select=c(chrpos)))
  biall = subset(allelic_state, freq == 1)
  biallelic_data = subset(variant_data, chrpos %in% biall$chrpos)
}

rebin_data = function(variant_data, rebin) {
  variant_data = variant_data %>% mutate(vqs_cut = rebin - ntile(vqslod, rebin))
  variant_data = variant_data %>% mutate(vqs_qd_cut = rebin - ntile(vqslod_qd, rebin))
  if ('rfprob' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(rf_cut = rebin - ntile(rfprob, rebin))
  }
  if ('rfprob_all' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(rf_all_cut = rebin - ntile(rfprob_all, rebin))
  }
  if ('rfprob_bal' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(rf_bal_cut = rebin - ntile(rfprob_bal, rebin))
  }
  if ('ab_mean' %in% colnames(variant_data)) {
    variant_data = variant_data %>% mutate(ab_cut = rebin - ntile(ab_mean, rebin))
  }
  variant_data
}

plot_singleton_titv = function(variant_data, hc_region=T, rebin=100, cumulative=F, return=F, datasets=c('rf', 'vqs')) {
  process_data = variant_data %>% filter(type == 'snv' & !indel & ac_orig == 1)
  if (hc_region) {
    process_data = process_data %>% filter(callrate > 0.99 & evaluation_interval & high_coverage_interval)
  }

  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(tis = sum(transition), tvs = sum(!transition), titv=tis/tvs,
              vqslod = max(vqslod), rfprob = max(rfprob)) -> bin_data
  
  if (cumulative) {
    bin_data = bin_data %>%
      arrange(cut) %>%
      mutate(tis = cumsum(tis), tvs = cumsum(tvs), titv=tis/tvs)
  }
  
  p = ggplot2(bin_data) + aes(x = cut, y = titv, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
  if (return) { return(bin_data) } else { return(p) }
}

plot_singleton_indel_ratio = function(variant_data, hc_region=T, rebin=100, cumulative=F, return=F, datasets=c('rf', 'vqs')) {
  process_data = variant_data %>% filter(indel & ac_orig == 1)
  if (hc_region) {
    process_data = process_data %>% filter(callrate > 0.99 & evaluation_interval & high_coverage_interval)
  }
  
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(ins = sum(insertion), dels = sum(!insertion), indel_ratio=ins/dels) -> bin_data
  if (cumulative) {
    bin_data = bin_data %>%
      arrange(cut) %>%
      mutate(ins = cumsum(ins), dels = cumsum(dels), indel_ratio=ins/dels)
  }
  
  p = ggplot2(bin_data) + aes(x = cut, y = indel_ratio, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
  if (return) { return(bin_data) } else { return(p) }
}

plot_transmission = function(variant_data, indels=F, rebin=100, return=F, correction=F, datasets=c('rf', 'vqs')) {
  process_data = variant_data %>%
    filter(indel == indels & 
             ac_orig <= 2 & ac_unrelated == 1 &
             transmitted + untransmitted > 0 & 
             callrate > 0.99 &
             !wassplit)
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  pass_cut = get_vqslod_pass_value(process_data)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(tx = sum(transmitted), ntx = sum(untransmitted > 0)) %>%
    arrange(cut) %>%
    mutate(tx = cumsum(tx), ntx = cumsum(ntx),
              transmission_rate = tx/(tx + ntx)) -> bin_data
  
  if (correction) {
    process_data = variant_data %>%
      filter(indel == indels & !wassplit & callrate > 0.99 & ac_orig <= 2)
    process_data %>%
      gather('metric', 'cut', rf_cut, vqs_cut, vqs_qd_cut) %>%
      group_by(metric, cut) %>%
      summarize(singleton = sum(ac_orig == 1), doubleton = sum(ac_orig == 2)) %>%
      arrange(cut) %>%
      mutate(singleton = cumsum(singleton), doubleton = cumsum(doubleton),
                correction_factor = doubleton/singleton) -> correction_data
    bin_data = bin_data %>%
      left_join(select(correction_data, metric, cut, correction_factor)) %>%
      mutate(corrected_transmission_rate = transmission_rate/correction_factor)
  }
  
  new_cut = ifelse(indels, 75, 90)
  bin_data$label = NA
  bin_data$label[bin_data$metric == 'vqslod' & bin_data$cut == pass_cut] = 'PASS'
  if (new) bin_data$label[bin_data$metric == 'vqslod' & bin_data$cut == new_cut] = 'New cutoff'
  p = ggplot2(bin_data) + geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  
  if (correction) {
    p = p + aes(x = cut, y = corrected_transmission_rate, color = metric)
  } else {
    p = p + aes(x = cut, y = transmission_rate, color = metric) +
      geom_hline(yintercept = 0.5)
  }
  print(p)
  if (return) { return(vqs_cumulative) } else { return(p) }
}

plot_singleton_fs_inframe = function(variant_data, rebin=100, return=F, correction=F, datasets=c('rf', 'vqs')) {
  process_data = variant_data %>% 
    filter(indel & ac_orig == 1 &
             callrate > 0.99 & evaluation_interval & high_coverage_interval)
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(fs = sum(!inframe), inframe = sum(inframe), 
              fs_in_ratio=fs/inframe) -> bin_data
  if (correction) {
    process_data %>%
      gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
      group_by(metric, cut) %>%
      summarize(one = sum(abs(bases_inserted) == 1, na.rm=T), 
                two = sum(abs(bases_inserted) == 2, na.rm=T), 
                three = sum(abs(bases_inserted) == 3, na.rm=T),
                correction_factor = one/two) -> correction_data
    # ggplot2(correction_data) + aes(x = vqs_cut, y = one/two) + geom_point()
    bin_data = bin_data %>%
      left_join(select(correction_data, metric, cut, correction_factor)) %>%
      mutate(corrected_fs_ratio = fs_in_ratio/correction_factor)
  }
  
  p = ggplot2(bin_data) + geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  if (correction) {
    p = p + aes(x = cut, y = corrected_fs_ratio, color = metric) 
  } else {
    p = p + aes(x = cut, y = fs_in_ratio, color = metric)
  }
  print(p)
  if (return) { return(bin_data) } else { return(p) }
}

plot_singleton_mendel = function(variant_data, indels=F, rebin=100, return=F, datasets=c('rf', 'vqs')) {
  process_data = filter(variant_data, indel == indels & ac_orig == 1)
  
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(mendel_errors=sum(mendel_errors == 1, na.rm=T), n=n()) -> bin_data
  
  p = ggplot2(bin_data) + aes(x = cut, y = mendel_errors, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
  if (return) { return(bin_data) } else { return(p) }
}

plot_denovo_sensitivity = function(variant_data, indels=F, rebin=100, return=F, datasets=c('rf', 'vqs')) {
  process_data = filter(variant_data, indel == indels & ac_orig == 1)
  
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(denovos = sum(!is.na(validated_denovo)), n=n()) %>%
    arrange(cut) %>%
    mutate(denovos = cumsum(denovos), denovos=denovos/max(denovos),
           n = cumsum(n))-> bin_data
  
  p = ggplot2(bin_data) + aes(x = cut, y = denovos, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
  if (return) { return(bin_data) } else { return(p) }
}

plot_indel_length_distribution = function(variant_data, rebin=100) {
  process_data = filter(variant_data, indel & callrate > 0.99 & ac_orig == 1)
  
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(one = sum(abs(bases_inserted) == 1, na.rm=T), 
              two = sum(abs(bases_inserted) == 2, na.rm=T), 
              n = n()) -> bin_data
  
  p = ggplot2(bin_data) + aes(x = cut, y = one/n, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
}

plot_proportion_singleton = function(variant_data, indels=F, rebin=100) {
  process_data = filter(variant_data, indel == indels)
  
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(singleton = sum(ac_orig == 1, na.rm=T), 
              common = sum(ac_orig >= 100, na.rm=T), 
              n = n()) -> bin_data
  
  p = ggplot2(bin_data) + aes(x = cut, y = singleton/n, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
}

plot_proportion_training = function(variant_data, indels=F, rebin=100, train_type='TP') {
  process_data = filter(variant_data, indel == indels)
  
  if (rebin > 0) process_data = rebin_data(process_data, rebin)
  
  process_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    group_by(metric, cut) %>%
    summarize(training = sum(train == train_type, na.rm=T),
              n = n()) -> bin_data
  
  p = ggplot2(bin_data) + aes(x = cut, y = training/n, color = metric) + 
    geom_point() + guides(size=F) + scale_size_identity() + scale_color_discrete(na.value='black')
  print(p)
}

all_summary_stats = function() {
  variant_data = read_variant_data('exome_variantqc.all.balanced.txt.bgz')
  plot_summary_stats(variant_data, 'all_exome_plots_all_bal.pdf')
  variant_data = read_variant_data('exome_variantqc.hardfilteronly.txt.bgz')
  plot_summary_stats(variant_data, 'all_exome_plots_hf_bal.pdf')
  variant_data = read_variant_data('exome_variantqc.hardfilteronly.unbalanced.txt.bgz')
  plot_summary_stats(variant_data, 'all_exome_plots_hf_unbal.pdf')
  variant_data = read_variant_data('exome_variantqc.unbalanced.txt.bgz')
  plot_summary_stats(variant_data, 'all_exome_plots_unbal.pdf')
  variant_data = read_variant_data('exome_variantqc.hardfilteronly.unbalanced.allsamples.txt.bgz')
  plot_summary_stats(variant_data, 'all_exome_plots_hf_unbal_allsamp.pdf')
  variant_data = read_variant_data('exome_variantqc.hf.rebalanced.txt.bgz')
  plot_summary_stats(variant_data, 'all_exome_plots_hf_rebal.pdf')
}



plot_summary_stats = function(variant_data, out_fname='all_exome_plots.pdf') {
  
  pdf(out_fname, height=4, width=6)
  # Titv/indels
  plot_singleton_titv(variant_data)
  # plot_singleton_titv(variant_data, cumulative=T)
  plot_singleton_indel_ratio(variant_data)
  # plot_singleton_indel_ratio(variant_data, cumulative=T)
  
  # De novos and mendel errors
  plot_denovo_sensitivity(variant_data)
  plot_denovo_sensitivity(variant_data, indels = T)
  plot_singleton_mendel(variant_data)
  plot_singleton_mendel(variant_data, indels=T)
  
  # Fs/in-frame ratio
  plot_singleton_fs_inframe(variant_data)
  plot_singleton_fs_inframe(variant_data, correction = T)
  
  # Misc metrics
  plot_indel_length_distribution(variant_data)
  plot_proportion_singleton(variant_data)
  
  # Allele balance
  # filter(variant_data, ac_orig == 1) %>%
    # ggplot2 + aes(x = ab_median) + geom_histogram() + xlab('% Alt Allele')
  
  # rfprob density
  p = ggplot2(variant_data) + aes(x = rfprob, col = type) + geom_density()
  print(p)
  
  # variant_data %>% group_by(rf_cut) %>% 
    # summarize(prop_mixed = sum(type == 'mixed')/n()) %>%
    # ggplot2 + aes(x = rf_cut, y = prop_mixed) + geom_bar(stat='identity')
  
  dev.off()
}

process_concordance = function(input_data, bins=100) {
  input_data$indel = as.logical(input_data$indel)
  input_data$vqslod_qd = input_data$vqslod
  input_data$vqslod_qd[input_data$qd < 4] = min(input_data$vqslod)
  
  input_data$model = ifelse(input_data$type %in% c('snv', 'multi-snv'), 'snp', 'indel')
  
  input_data$bases_inserted = nchar(input_data$alt) - nchar(input_data$ref)
  
  input_data = input_data %>% mutate(rf_cut = bins - ntile(rfprob, bins))
  input_data %>% group_by(model) %>% 
    mutate(vqs_cut = bins - ntile(vqslod, bins)) %>% ungroup %>% select(-concordance)
}

plot_truth_pr = function(input_data, indels=F, rebin=100, datasets=c('vqs', 'rf'), strict=F) {
  use_data = input_data %>% filter(indel == indels)
  
  total_fns = sum(use_data$called_gt == 'missing' & use_data$truth_gt != 'missing')
  total_fps = sum(use_data$truth_gt == 'missing')
  total_tps = sum(use_data$called_gt != 'missing' & use_data$truth_gt != 'missing')
  # total_tps = sum(use_data$called_gt == use_data$truth_gt)
  # print(paste('Overall precision:', total_tps/(total_tps + total_fps)))
  # print(paste('Overall recall:', total_tps/(total_tps + total_fns)))
  if (indels) {
    
  }
  
  if (rebin > 0) use_data = rebin_data(use_data, rebin)
  
  use_data = use_data %>%
    gather_('metric', 'cut', paste0(datasets, '_cut')) %>%
    mutate(cut = ifelse(!is.na(cut), cut, 101)) %>%
    group_by(metric, cut)
  
  use_data %>%
    summarize(tp = sum(truth_gt != "missing" & called_gt != "missing"), 
           fp = sum(truth_gt == "missing"),
           fn = sum(truth_gt != "missing"),
           n = n(),
           rfprob = max(rfprob)) %>%
    arrange(metric, cut) -> process_data
  
  binned_data = process_data %>%
    mutate(cum_tp = cumsum(tp),
           cum_fp = cumsum(fp),
           cum_n = cumsum(n)) %>%
    arrange(desc(cut)) %>%
    mutate(cum_fn = cumsum(fn),
           precision = cum_tp / (cum_tp + cum_fp),
           recall = cum_tp / (cum_tp + cum_fn))
  
  # TODO: figure out final point recall inflation
  endpoints = data.frame(metric=rep(paste0(datasets, '_cut'), 2),
                         precision=rep(c(0, 1), each=length(datasets)),
                         recall=rep(c(1, 0), each=length(datasets)))
  print(binned_data %>% group_by(metric) %>% 
          select(metric, precision, recall) %>%
          bind_rows(endpoints) %>%
          arrange(recall) %>% 
          summarize(auprc = trapz(recall, precision)))
  
  cut_point = ifelse(indels, 75, 90)
  label_data = subset(binned_data, cut == cut_point)
  title = paste(max(binned_data$cum_n), ifelse(indels, 'Indels', 'SNPs'))
  p = ggplot2(binned_data) + aes(recall, precision, col = metric) + 
    geom_point() + ggtitle(title) + geom_label_repel(aes_(label=cut_point), label_data, min.segment.length=unit(0, 'in'))
  print(p)
  binned_data %>% filter(metric == 'rf_cut') %>%
    ggplot2 + aes(recall, precision, col = rfprob) +   
    geom_point() + ggtitle(title) + scale_color_gradient(low="red") -> p
  print(p)
}

truth_stats = function(type = 'stats') {
  # Synthetic-diploid
  syndip_data = read.table(paste0('data/exome_syndip.', type, '.txt.bgz'), header=T) %>% process_concordance
  # table(select(syndip_data, truth_gt, called_gt))
  plot_truth_pr(syndip_data)
  plot_truth_pr(syndip_data, indels = T)
  
  # 12878
  na12878_data = read.table(paste0('data/exome_na12878.', type, '.txt.bgz'), header=T) %>% process_concordance
  # table(select(na12878_data, truth_gt, called_gt))
  plot_truth_pr(na12878_data)
  plot_truth_pr(na12878_data, indels = T)
}

plot_truth_stats = function(type='stats') {
  pdf(paste0('truth_', type, '.pdf'))
  truth_stats(type)
  dev.off()
}

all_truth_stats = function() {
  plot_truth_stats()
  plot_truth_stats('exome_unbalanced.stats')
  plot_truth_stats('exome_hardfilteronly.unbalanced.stats')
  plot_truth_stats('exome_hardfilteronly.stats')
  plot_truth_stats('exome_hardfilteronly.unbalanced.fsqd.stats')
  plot_truth_stats('exome_hardfilteronly.unbalanced.allsamples.stats')
  plot_truth_stats('exome_hf.rebalanced.stats')
}

  # V1 V2 concordance
  # Samples
#   sample_concordance = read.delim('data/v1_compare_sample_concordance.txt.bgz', header=T)
#   ggplot(sample_concordance) + aes(x = correct/total) + geom_histogram()
#   subset(sample_concordance, correct/total < 0.95)
#   
#   # Variants
#   variant_concordance = read.delim('data/v1_compare_concordance.txt.bgz', header=T)
#   # variant_concordance = read.delim('data/v1_compare_highcov_concordance.txt.bgz', header=T)
#   variant_concordance = process_concordance(variant_concordance)
#   
#   library(jsonlite)
#   library(parallel)
#   # conc_data = mclapply(variant_concordance$concordance, get_stats_l, mc.cores = 4)
#   conc_data = lapply(variant_concordance$concordance, get_stats_l)
#   conc_data = data.frame(t(simplify2array(conc_data)))
#   colnames(conc_data) = c('ref', 'right_call', 'wrong_call', 'missing_left', 'missing_right', 'missing_gt_left', 'missing_gt_right')

get_stats_l = function(x) {
  m = fromJSON(x)
  # 1: No Data (missing variant)
  # 2: No Call (missing genotype call)
  # 3: Hom Ref
  # 4: Heterozygous
  # 5: Hom Var
  ref = m[3,3]
  right_call = m[4,4] + m[5,5]
  wrong_call = m[3,4] + m[4,3] + m[3,5] + m[5,3] + m[4,5] + m[5,4]
  missing_left = m[1,3] + m[1,4] + m[1,5]
  missing_right = m[3,1] + m[4,1] + m[5,1]
  missing_gt_left = m[2,3] + m[2,4] + m[2,5]
  missing_gt_right = m[3,2] + m[4,2] + m[5,2]
  return(c(ref, right_call, wrong_call, missing_left, missing_right, missing_gt_left, missing_gt_right))
}



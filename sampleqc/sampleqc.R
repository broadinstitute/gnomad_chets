options(stringsAsFactors = F)
library(plyr)
library(dplyr)
library(tidyverse)
library(broom)
library(ggrepel)
library(rpart)
library(randomForest)
library(magrittr)
# install.packages(c('plyr', 'dplyr', 'ggplot2', 'magrittr', 'shiny', 'DT', 'tidyverse', 'broom', 'ggrepel', 'randomForest', 'ROCR'))
# devtools::install_github('baptiste/ggflags')

PROBABILITY_CUTOFF = 0.9

source('colors_and_constants.R')

# Platform data
read_missingness_pca_data = function() {
  missingness_pca_data = read.delim('data/missingness_pca_data.txt.gz', header=F)
  colnames(missingness_pca_data) = c('sample_name_in_vcf', paste0('pc', 1:10))
  return(missingness_pca_data)
}
read_platforms_data = function() {
  platforms = read.csv('data/ExAc_Exomes_Since_12Oct2012_June2016.csv', header=T)
  platforms %<>% select(project_or_cohort = Seq.Project, product = Product) %>% distinct
  platform2 = read.csv('data/Products for Konrad.csv', header=T)
  platform2 %<>% unite(product, Product1, Product2, sep='|') %>%
    mutate(product = gsub('\\|$', '', product, perl=T)) %>% 
    select(project_or_cohort = project, product)
  platforms %<>% rbind(platform2)
  platforms$gross_platform = 'unknown'
  platforms$gross_platform[platforms$product %in% agilent] = 'agilent'
  platforms$gross_platform[platforms$product %in% ice] = 'ice'
  platforms$gross_platform[platforms$product %in% ice150] = 'ice150'
  platforms %<>% 
    group_by(project_or_cohort) %>% 
    summarize(gross_platform = paste(unique(gross_platform), collapse='|'), 
              platform = paste(unique(product), collapse='|'))
  platforms$gross_platform[grepl('|', platforms$gross_platform, fixed=T)] = 'multiple'
  platforms %<>% rbind(data.frame(project_or_cohort = 'FINNMETSEQ', gross_platform = 'nimblegen', platform = 'Nimblegen (Finns)'))
  return(platforms)
}
read_samples_platforms = function() {
  samples = read.delim('data/pid_sid.txt', header=T, sep='\t')
  names(samples)[2] = 'sample'
  
  platforms = read_platforms_data()
  samples %<>% left_join(platforms)
  samples$gross_platform[is.na(samples$platform)] = 'unknown'
  return(samples)
}
final_platforms = function() {
  samples = read_samples_platforms()
  missingness_pca = read_missingness_pca_data()
  
  samples %<>% inner_join(missingness_pca)
  
  # Classifying multiple into ice/agilent/ice150/unknown
  library(rpart)
  
  training_data = subset(samples, gross_platform %in% c('ice', 'ice150', 'agilent'))
  training_data$gross_platform %<>% factor
  
  fit = rpart(gross_platform ~ pc1 + pc2 + pc3, training_data)
  
  function() {
    model_data = subset(samples, gross_platform %in% c('ice', 'ice150', 'agilent'))
    model_data$gross_platform %<>% factor
    
    training_data = sample_frac(model_data, 0.9)
    test_data = subset(model_data, !(sample %in% training_data$sample))
    
    fit = rpart(gross_platform ~ pc1 + pc2 + pc3, training_data)
    fit_data = data.frame(predict(fit, test_data), sample_name_in_vcf = test_data$sample_name_in_vcf,
                          actual_platform=test_data$gross_platform)
    fit_data %<>%
      gather(predicted_platform, probability, -sample_name_in_vcf, -actual_platform) %>%
      group_by(sample_name_in_vcf) %>%
      slice(which.max(probability))
    
    sum(fit_data$actual_platform == fit_data$predicted_platform)/nrow(fit_data) # 0.9904
  }
  
  multiple_data = subset(samples, gross_platform == 'multiple')
  
  fit_data = data.frame(predict(fit, multiple_data), sample_name_in_vcf = multiple_data$sample_name_in_vcf)
  fit_data %<>%
    gather(predicted_platform, probability, -sample_name_in_vcf) %>%
    group_by(sample_name_in_vcf) %>%
    slice(which.max(probability))
  
  samples %<>% left_join(fit_data, by='sample_name_in_vcf') %>%
    mutate(gross_platform = ifelse(!is.na(predicted_platform), predicted_platform, gross_platform))
  return(samples)
}

# Step 1: Generate final platforms
function() {
  samples = final_platforms()
  output_file = gzfile('data/final_platforms.tsv.gz', 'w')
  write.table(samples, file=output_file, quote=F, row.names=F, sep='\t')
  close(output_file)
}

# Step 2: assign populations
function() {
  forest_data = final_population()
  table(known = forest_data$known_pop, pred = forest_data$predicted_pop)
  table(forest_data$predicted_pop)
  
  output_file = gzfile('data/final_populations.tsv.gz', 'w')
  write.table(forest_data, file=output_file, quote=F, row.names=F, sep='\t')
  close(output_file)
  
  forest_data %>% filter(source == 'ExAC') %>%
    select(sample, final_pop=predicted_pop) %>% 
    write.table(file='data/pops_exac.txt', quote=F, row.names=F, sep='\t')
  
  forest_data %>% filter(source == 'gnomAD') %>%
    select(sample, final_pop=predicted_pop) %>% 
    write.table(file='data/pops_gnomad.txt', quote=F, row.names=F, sep='\t')
}

# Step 3: final samples
function() {
  final_data = read_metadata(type = 'final')
  output_file = gzfile('data/final_samples.tsv.gz', 'w')
  write.table(final_data$sample, file=output_file, quote=F, row.names=F, col.names=F, sep='\t')
  close(output_file)
}

# Population assignment code
final_population = function() {
  forest_data = get_forest_data()
  
  # Evaluation and plotting
  # evaluate_rf()
  
  forest_data$predicted_pop[forest_data$probability < PROBABILITY_CUTOFF] = 'oth'
  forest_data$predicted_pop[forest_data$predicted_pop == 'eur'] = 'nfe'
  
  forest_data %$% table(source, predicted_pop)
  return(forest_data)
}
evaluate_rf = function(save_plots=T) {
  if (save_plots) pdf('rf_evaluate.pdf', width=6, height=4)
  # Get all data
  data = exac_and_gnomad('all')
  all_known = get_known_samples(data)
  all_known_data = data %>% inner_join(all_known)
  
  # Split into training and testing
  set.seed(42)
  training_subset = all_known_data %>% group_by(known_pop) %>% sample_frac(0.8)
  testing_subset = filter(all_known_data, !(sample %in% training_subset$sample))
  
  # Hold out each projects
  
  # Run model
  test_output = pop_forest(training_subset, testing_subset)
  
  # Evaluate model
  test_output %<>% left_join(testing_subset)
  table(known = test_output$known_pop, pred = test_output$predicted_pop)
  test_results = test_output %$% table(known = known_pop, pred = predicted_pop)
  print(paste0('Error: ', round(100*(1 - sum(diag(test_results))/sum(test_results)), 3), '%'))
  print(test_results)

  test_results = test_output %>% filter(probability > PROBABILITY_CUTOFF) %$% table(known = known_pop, pred = predicted_pop)
  print(paste0('Error: ', round(100*(1 - sum(diag(test_results))/sum(test_results)), 3), '%'))
  print(test_results)
  
  per_pop_sens = table(known = test_output$known_pop, pred = test_output$predicted_pop == test_output$known_pop)
  per_pop_sens = data.frame(prop=per_pop_sens[,2]/(per_pop_sens[,1] + per_pop_sens[,2]))
  print('Sensitivity per pop:')
  print(per_pop_sens)

  per_pop_spec = table(known = test_output$predicted_pop, pred = test_output$predicted_pop == test_output$known_pop)
  per_pop_spec = data.frame(prop=per_pop_spec[,2]/(per_pop_spec[,1] + per_pop_spec[,2]))
  print('Specificity per pop:')
  print(per_pop_spec)
  
  accuracy_by_probability = test_output %>% ungroup %>% group_by(probability) %>%
    summarize(total = n(), right = sum(predicted_pop == known_pop), proportion_right = right/total)
  p1 = ggplot(accuracy_by_probability) + aes(x = probability, y = proportion_right) + geom_point()
  print(p1)
  
  accuracy_by_probability_cumulative = test_output %>% ungroup %>% 
    arrange(desc(probability)) %>% 
    mutate(cumulative_right = cumsum(predicted_pop == known_pop), total = cumsum(!is.na(sample)),
           proportion_right = cumulative_right/total) %>%
    group_by(probability) %>% summarize(proportion_right = last(proportion_right), total=last(total))
  p2 = ggplot(accuracy_by_probability_cumulative) + aes(x = probability, y = proportion_right) + geom_point()
  print(p2)
  
  accuracy_by_probability_cumulative_per_pop_sens = test_output %>% ungroup %>% 
    group_by(known_pop) %>%
    arrange(desc(probability)) %>% 
    mutate(cumulative_right = cumsum(predicted_pop == known_pop), total = cumsum(!is.na(sample)),
           proportion_right = cumulative_right/total) %>%
    group_by(known_pop, probability) %>% summarize(proportion_right = last(proportion_right), total=last(total))
  p3 = ggplot(accuracy_by_probability_cumulative_per_pop_sens) + aes(x = probability, y = proportion_right, col=known_pop) +
    geom_point() + scale_color_manual(values=pop_colors) + facet_grid(known_pop~.)
  print(p3)
  
  accuracy_by_probability_cumulative_per_pop_spec = test_output %>% ungroup %>% 
    group_by(predicted_pop) %>%
    arrange(desc(probability)) %>% 
    mutate(cumulative_right = cumsum(predicted_pop == known_pop), total = cumsum(!is.na(sample)),
           proportion_right = cumulative_right/total) %>%
    group_by(predicted_pop, probability) %>% summarize(proportion_right = last(proportion_right), total=last(total))
  p4 = ggplot(accuracy_by_probability_cumulative_per_pop_spec) + aes(x = probability, y = proportion_right, col=predicted_pop) +
    geom_point() + scale_color_manual(values=pop_colors) + facet_grid(predicted_pop~.)
  print(p4)
  
  library(ROCR)
  
  pred <- prediction(test_output$probability, test_output$predicted_pop == test_output$known_pop)
  perf <- performance(pred, measure = "prec", x.measure = "rec") 
  plot(perf, colorize=T, lwd=3)
  
  table(forest_data$predicted_pop)
  table(known = forest_data$known_pop, pred = forest_data$predicted_pop)
  removed_by_prob = forest_data %>% ungroup %>% arrange(probability) %>%
    mutate(total = cumsum(!is.na(sample))) %>% group_by(probability) %>%
    summarize(total=last(total))
  
  p5 = ggplot(removed_by_prob) + aes(x = probability, y = total) + geom_point() + scale_y_log10() + annotation_logticks()
  print(p5)
  
  if (save_plots) dev.off()
}
get_forest_data = function(separate_estonians=F) {
  data = exac_and_gnomad()
  all_known = get_known_samples(data, separate_estonians=separate_estonians)
  
  # Run random forest
  all_known_data = data %>% inner_join(all_known)
  fit_data = pop_forest(all_known_data, data)
  data %>% left_join(all_known) %>% left_join(fit_data, by='sample')
}
pop_forest = function(training_data, data, ntree=100, seed=42) {
  set.seed(seed)
  forest = randomForest(as.factor(known_pop) ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6,
                        data = training_data,
                        importance = T,
                        ntree = ntree)
  
  fit_data = data.frame(predict(forest, data, type='prob'), sample = data$sample)
  fit_data %>%
    gather(predicted_pop, probability, -sample) %>%
    group_by(sample) %>%
    slice(which.max(probability))
}
read_1kg_pops = function() {
  kg_data = read.delim('pop_assignment/1kg_all_pop.tsv', header=F)
  names(kg_data) = c('sample', 'super_pop', 'pop')
  kg_data %<>% mutate(super_pop=tolower(super_pop), pop=tolower(pop))
  kg_data$super_pop[kg_data$super_pop == 'san'] = 'sas'
  kg_data$super_pop[kg_data$super_pop == 'asn'] = 'eas'
  kg_data$super_pop[kg_data$pop == 'fin'] = 'fin'
  return(kg_data)
}
get_known_samples = function(data, separate_estonians=F, europe_only=F) {
  # Europeans
  german = filter(data, project_or_cohort == 'C1708') %>% select(sample) %>% mutate(known_pop='de')
  icr = filter(data, project_or_cohort %in% c('ICR1000', 'ICR142')) %>% select(sample) %>% mutate(known_pop='gb')
  atvb = filter(data, project_or_cohort == 'C1017') %>% select(sample) %>% mutate(known_pop='it')
  regicor = filter(data, project_or_cohort == 'C1568') %>% select(sample) %>% mutate(known_pop='es')
  bulgarian_trios = filter(data, project_or_cohort %in% c('Bulgarian_Trios', 'C533', 'C821', 'C952')) %>% select(sample) %>% mutate(known_pop='bg')
  eur = rbind(german, icr, atvb, regicor, bulgarian_trios)
  
  est = filter(data, project_or_cohort %in% c('G89634', 'G94980')) %>% select(sample) %>% mutate(known_pop='ee')
  
  # Finns from Mitja
  finn_metadata = read.delim('pop_assignment/99percent_finns_plus_AD_IBD_NFID.tsv.gz', header=T)
  colnames(finn_metadata) = tolower(colnames(finn_metadata))
  finn_samples = subset(finn_metadata, percent_finnish > 0.99) %>% select(sample_name_in_vcf)
  finns = subset(data, sample_name_in_vcf %in% finn_samples$sample_name_in_vcf) %>% select(sample) %>% mutate(known_pop='fi')
  
  if (europe_only) {
    # leicester = filter(data, project_or_cohort == 'C1830') %>% select(sample) %>% mutate(known_pop='gb')
    return(distinct(rbind(eur, est, finns)))
  }
  # If not Europe, combine them all
  eur$known_pop = 'eur'
  est$known_pop = if (separate_estonians) 'est' else 'eur'
  finns$known_pop = 'fin'
  
  # 1kg - 2187
  kg_pops = read_1kg_pops()
  kg_pops %<>% subset(sample %in% data$sample) %>% select(sample, known_pop=super_pop)
  
  # MDE = Walsh, ciliopathies c('C871', 'C1441')
  # MDE - 963
  # mdes = subset(data, project_or_cohort %in% c('C871', 'C1441')) %>% select(sample)
  # mdes$known_pop = 'mde'
  
  # AJs - 2726 samples
  ajs = read.table('pop_assignment/aj.9_ids')
  colnames(ajs) = 'sample'
  ajs %<>% subset(sample %in% data$sample)
  ajs$known_pop = 'asj'
  
  # PROMIS
  sas = subset(data, grepl('PROMIS', description)) %>% select(sample)
  sas$known_pop = 'sas'
  
  # T2D SIGMA
  sigma = filter(data, grepl('SIGMA', description) & exac_version == 'ExACv1') %>% select(sample)
  sigma$known_pop = 'amr'
  
  # African-Americans
  aa_t2d = c('C773', 'C1002') # T2D-GENES AA cohorts
  aa_jhs = 'C1567' 
  aa_biome = 'C1956'
  afr_cohorts = c(aa_t2d, aa_jhs, aa_biome)
  afr = filter(data, project_or_cohort %in% afr_cohorts) %>% select(sample)
  afr$known_pop = 'afr'
  
  # EAS
  taiwanese_trios = c('C1397', 'C1443', 'C1506', 'C1867', 'C978')
  eas_t2d = 'C774'
  singapore = 'C1940'
  hkg = 'C1980'
  korean = 'C1982'
  eas_cohorts = c(taiwanese_trios, eas_t2d, singapore, hkg, korean)
  eas = filter(data, project_or_cohort %in% eas_cohorts) %>% select(sample)
  eas$known_pop = 'eas'
  
  all_known = rbind(kg_pops, ajs, sas, finns, sigma, eur, afr, eas, est)
  return(distinct(all_known))
}

# Read final assignments
read_platform_data = function() {
  return(read.delim('data/final_platforms.tsv.gz', header=T))
}
read_population_data = function() {
  return(read.delim('data/final_populations.tsv.gz', header=T))
}

# Final metadata functions
read_metadata = function(type='final') {
  metadata = read.delim('data/super_meta.txt.gz', header=T, sep='\t')
  # metadata$sample %<>% gsub(' ', '_', .)
  
  metadata$project_or_cohort = metadata$pid
  metadata$project = metadata$project_or_cohort
  metadata$project[grepl('ESP_', metadata$project_or_cohort) || grepl('ESP', metadata$description)] = 'ESP'
  metadata$project[grepl('1kg_', metadata$project_or_cohort) || grepl('1KG', metadata$description)] = '1kg'
  metadata$project[grepl('TCGA', metadata$description)] = 'TCGA'
  metadata$project[grepl('GTEx', metadata$description)] = 'GTEx'
  
  library(assertthat)
  assert_that(199558 == nrow(metadata)) # cp1
  assert_that(139082 == sum(!grepl('hard', metadata$drop_condense))) # cp2
  assert_that(136204 == sum(!grepl('hard', metadata$drop_condense) & !grepl('duplicate', metadata$drop_condense) )) # cp3
  assert_that(127970 == sum(metadata$drop_status == 'keep' | metadata$drop_condense == 'pop' | metadata$drop_condense == 'related_gno' | metadata$drop_condense == 'pop,related_gno')) # cp4
  assert_that(123136 == sum(metadata$drop_status == 'keep')) # final
  
  # QC samples
  # assert_that(127896 == sum((metadata$sample %in% read.table('exac2.qctrios.fam')$V2 & metadata$description %in% c('Bulgarian_Trios', 'Bulgarian_Trios_External', 'TaiTrios_Nexome_G4L', 'Taiwanese_Trios')) | metadata$drop_status == 'keep'))
  # subset(metadata, (sample %in% read.table('exac2.qctrios.fam')$V2 & description %in% c('Bulgarian_Trios', 'Bulgarian_Trios_External', 'TaiTrios_Nexome_G4L', 'Taiwanese_Trios')) | drop_status == 'keep') 
  # %>% select(sample) 
  # %>% write.table(quote=F, row.names=F, file='release_plus_trios.txt')
  
  if (type == 'final') {
    metadata %<>% subset(drop_status == 'keep')
  } else if (type == 'permission') {
    metadata %<>% subset(permission == 'YES')
  } else if (type == 'cp5') { # 131250
    metadata %<>% filter(drop_status == 'keep' | drop_condense == 'pop')
  } else if (type == 'cp5_mde') { # 132158
    metadata %<>% filter((drop_status == 'keep' | drop_condense == 'pop') | 
                           (project_or_cohort %in% c('C871', 'C1441') & drop_reasons == 'NR'))
  } else if (type == 'cp4') { # 132248
    metadata %<>% filter(drop_status == 'keep' | drop_condense == 'pop' | drop_reasons == 'related_gno')
  } else if (type == 'cp4_mde') { # 133161
    metadata %<>% filter((drop_status == 'keep' | drop_condense == 'pop' | drop_reasons == 'related_gno') |
                             (project_or_cohort %in% c('C871', 'C1441') & (drop_reasons == 'NR' | drop_reasons == 'NR,related_gno')))
  }
  metadata$gross_platform[is.na(metadata$gross_platform)] = 'unknown'
  
  # Fun fact: we have 5 different sample identifiers
  # sample_names = select(metadata, underscored_sample, sample, sample_name, sample_name_in_vcf, old_sample_name, old_sample_name_in_vcf)
  # sapply(colnames(sample_names), function(x) {
  #   sapply(colnames(sample_names), function(y) {
  #     all(sample_names[,x] == sample_names[,y])
  #   })
  # })
  
  return(metadata)
}
read_exac_gnomad_pca_data = function() {
  data = read.delim('data/gnomad.pca.txt.gz', header=T, sep='\t')
  colnames(data) = tolower(colnames(data))
  return(data)
}
read_gnomad_metadata = function() {
  data = read.delim('data/gnomad.final.all_meta.txt.gz', header=T)
  colnames(data) = tolower(colnames(data))
  data$project = data$project_or_cohort
  data$combined_sample = paste0('genome_', gsub(' ', '_', data$sample))
  data %>% 
    mutate(permission = ifelse(releasable, 'YES', 'NO'), ex_in = 'internal') %>%
    select(sample, combined_sample, sex, callrate:rinsertiondeletion, project_or_cohort, project, pi, platform = product_simplified, 
           permission, description = research_project, keep, population = final_pop)
}

exac_and_gnomad = function(type='all') {
  data = read_metadata(type)
  data$combined_sample = paste0('exome_', gsub(' ', '_', data$sample))
  data$keep = data$drop_status == 'keep'
  gnomad_meta = read_gnomad_metadata()
  pca_data = read_exac_gnomad_pca_data()
  
  return(data %>% bind_rows(gnomad_meta) %>%
           inner_join(pca_data, by = c('combined_sample' = 'sample')) %>%
           mutate(source = ifelse(grepl('genome_', combined_sample), 'gnomAD', 'ExAC')))
}

final_gnomad_meta = function(write=F) {
  shiny_data = final_population()
  shiny_data$gross_platform[shiny_data$source == 'gnomAD'] = 'gnomAD'
  
  platform_data = read_platform_data()
  colnames(platform_data)[7:16] %<>% paste0('missingness_', .)
  platform_data$combined_sample = paste0('exome_', gsub(' ', '_', platform_data$sample))
  
  shiny_data %<>% left_join(select(platform_data, combined_sample, missingness_pc1:missingness_pc10))
  
  shiny_data %<>%
    mutate(overall_platform = factor(gross_platform, levels = c('unknown', 'multiple', 'nimblegen', 'ice150', 'ice', 'agilent', 'gnomAD'))) %>%
    mutate(predicted_pop = factor(predicted_pop, levels = c('oth', 'nfe', 'amr', 'sas', 'fin', 'eas', 'afr', 'asj'))) %>%
    arrange(overall_platform)
  if (write) write.table(shiny_data, 'gnomAD_super_super_meta.txt', quote=F, row.names=F, sep='\t')
  shiny_data
}



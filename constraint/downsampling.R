methylation_bins = 5
color_cpg = '#2E9FFE'
color_ti = '#458B00'
color_tv = '#EA4444'
variant_types = c('transversion', 'non-CpG transition', 'CpG transition')
variant_type_colors = c(color_tv, color_ti, color_cpg)
variant_type_fill = scale_fill_manual(name='Variant', values=c("CpG transition" = color_cpg, "non-CpG transition" = color_ti, 'transversion' = color_tv))
variant_type_legend = scale_color_manual(name='Variant', values=c("CpG transition" = color_cpg, "non-CpG transition" = color_ti, 'transversion' = color_tv))

prop_observed = tbl_df(read.table('data/proportion_observed.txt.bgz', header=T))
prop_observed_long = prop_observed %>%
  dplyr::select(-sum_coverage, -proportion_observed) %>%
  gather(property, num, -context, -ref, -alt, -methylation_level, -annotation, -possible_variants) %>%
  separate(property, c('property', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 123135, as.numeric(n)),
         num = ifelse(is.na(num), 0, num),
         transition = (ref == 'A' & alt == 'G') | (ref == 'G' & alt == 'A') | (ref == 'C' & alt == 'T') | (ref == 'T' & alt == 'C'),
         cpg = (ref == 'C' & alt == 'T' & substr(context, 2, 3) == 'CG') | (ref == 'G' & alt == 'A' & substr(context, 1, 2) == 'CG'),
         methylation = round(methylation_level/methylation_bins),
         variant_type = case_when(
           !transition ~ 'transversion',
           cpg ~ 'CpG transition',
           !cpg & transition ~ 'non-CpG transition',
           TRUE ~ 'wut'
         ),
         variant_type_methylation = ifelse(variant_type == 'CpG transition', paste(variant_type, methylation), variant_type))
  

prop_observed_long %>% 
  filter(annotation == 'synonymous_variant' & n == 123136) %>%
  group_by(variant_type, methylation) %>%
  summarize(observed = sum(num)) -> total_observed
  
prop_observed_long %>%
  filter(annotation == 'synonymous_variant' & n != 123135) %>%
  group_by(variant_type, methylation, variant_type_methylation, n) %>%
  left_join(total_observed) %>%
  summarize(possible_variants = sum(possible_variants),
            num = sum(num),
            proportion_observed = num/possible_variants,
            proportion_observed_thus_far = num/mean(observed)) %>%
  ungroup %>%
  mutate(variant_type_methylation = factor(variant_type_methylation, levels = names(variant_type_values))) %>%
  arrange(variant_type_methylation) %>%
  group_by(variant_type_methylation) %>%
  mutate(
    added_variants = ifelse(!is.na(lag(num)), num - lag(num), num),
    s_n = sqrt(n),
    added_s_n = ifelse(!is.na(lag(s_n)), s_n - lag(s_n), s_n),
    added_n = ifelse(!is.na(lag(n)), n - lag(n), n),
    added_variants_per_sqrt_n = added_variants/added_s_n,
    added_variants_per_n = added_variants/added_n
         ) -> collapsed_prop_observed

ggplot(collapsed_prop_observed) + 
  aes(x = n, y = proportion_observed) +
  scale_x_sqrt(limits=c(1000, 120000)) +
  geom_line(aes(color = variant_type_methylation)) + theme_classic()

plot_prop_observed_methylation_points = function() {
  ggplot(collapsed_prop_observed) + 
    # aes(x = n, y = proportion_observed_thus_far, col=paste(variant_type, methylation)) +
    aes(x = n, y = proportion_observed) +
    # geom_line() + 
    scale_x_sqrt(limits=c(1000, 120000)) +
    geom_line(aes(color = variant_type)) +
    scale_color_manual(name='Variant', values=c("CpG transition" = '#FFFFFF', "non-CpG transition" = color_ti, 'transversion' = color_tv)) +
    geom_point(data=subset(collapsed_prop_observed, variant_type == 'CpG transition'), size=2, shape=21, stroke=0, aes(fill = methylation)) +
    scale_fill_gradientn(name='Mean Methylation\nof CpG transition', #trans = "reverse", 
                         # colors=c("#9FFFFF", '#619CFF', "#3A5E9A")) +
                         colors=c("#3A5E9A", '#619CFF', "#9FFFFF")) + #132B43->56B1F7, 3A5E9A->9FFFFF
    guides(
      color = guide_legend(order=1),
      fill = guide_colorbar(order=2)
    ) + theme_classic()
}

plot_prop_observed_methylation_lines = function() {
  library(RColorBrewer)
  variant_type_values = c()
  x = max(collapsed_prop_observed$methylation, na.rm=T)
  for (i in seq(x, 0, -1)) {
    variant_type_values[paste('CpG transition', i)] = colorRampPalette(brewer.pal(9, 'Blues'))(x+2)[i+2]
  }
  variant_type_values["CpG transition"] = color_cpg
  variant_type_values["non-CpG transition"] = color_ti
  variant_type_values['transversion'] = color_tv
  
  collapsed_prop_observed %>%
    ggplot() + 
    aes(x = n, y = proportion_observed) +
    scale_x_sqrt(limits=c(1000, 120000)) +
    geom_line(aes(color = variant_type_methylation)) +
    scale_color_manual(name='Variant (methylation level)', values=variant_type_values) + theme_classic()
}

derivative_plot = function() {
  ggplot(collapsed_prop_observed) + 
    aes(x = s_n, y = proportion_observed) +
    geom_line(aes(color = variant_type_methylation)) + theme_classic()
  
  ggplot(collapsed_prop_observed) + 
    aes(x = n, y = num) +
    geom_line(aes(color = variant_type_methylation)) + theme_classic()
  collapsed_prop_observed %>%
    ggplot + theme_classic() + 
    # aes(x = s_n, y = added_variants_per_sqrt_n) +
    aes(x = n, y = added_variants_per_n) +
    # scale_x_sqrt(limits=c(1000, 120000)) + 
    # scale_x_log10() + scale_y_log10() +
    # geom_line(aes(color = variant_type_methylation, fill=variant_type_methylation)) 
    geom_smooth(aes(color = variant_type_methylation, fill=variant_type_methylation))
  
  collapsed_prop_observed %>%
    ggplot + theme_classic() + 
    aes(x = s_n, y = added_variants_per_sqrt_n/possible_variants) +
    # scale_x_sqrt(limits=c(1000, 120000)) + 
    # scale_x_log10() + scale_y_log10() +
    # geom_line(aes(color = variant_type_methylation, fill=variant_type_methylation)) 
    geom_smooth(aes(color = variant_type_methylation, fill=variant_type_methylation))
}

collapsed_prop_observed %>%
  filter(proportion_observed > 0.23) %>%
  summarize(at_20=min(n)) -> got_to_20

prop_observed_long %>%
  filter(annotation == 'synonymous_variant') %>%
  left_join(got_to_20) %>%
  mutate(road_to_20 = n/at_20) %>%
  filter(road_to_20 < 1) %>%
  group_by(variant_type_methylation, road_to_20) %>%
  summarize(possible_variants = sum(possible_variants),
            num = sum(num),
            proportion_observed = num/possible_variants) %>%
  ggplot + aes(x = road_to_20, y = proportion_observed, col=variant_type_methylation) +
  geom_line() + theme_classic()
  # scale_color_manual(values = variant_type_colors) + 
  scale_color_manual(name='Variant (methylation level)', values=variant_type_values) + 


syn_exons = tbl_df(read.table('data/syn.txt.bgz', header=T)) %>% left_join(genes, by=c('transcript' = 'ensembl_transcript_id'))
syn_exons_long = syn_exons %>%
  dplyr::select(transcript, exon, variant_count:median_coverage, external_gene_name) %>%
  gather(property, num, -transcript, -exon, -external_gene_name, -aggregate_mutation_rate) %>%
  separate(property, c('property', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 123135, as.numeric(n)))
syn_exons_long %>% 
  group_by(property, n) %>% 
  summarise(num=sum(num)) -> all_exons
all_exons %>%
  filter(property == 'variant_count') %>%
  ggplot + aes(x = n, y = num, col=property) + geom_line()

syn_exons_long %>%
  filter(property == 'variant_count') %>%
  group_by(n) %>%
  do(tidy(lm(num ~ aggregate_mutation_rate, data=.))) -> exon_lms
ggplot(exon_lms) + aes(x = n, y = estimate, col = term) + geom_line()

syn_depth = tbl_df(read.table('data/syn_depth_explore.txt.bgz', header=T))
syn_depth_long = syn_depth %>%
  gather(property, num, -median_coverage) %>%
  separate(property, c('property', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 124000, as.numeric(n)))
syn_depth_long %>%
  filter(property == 'expected' & median_coverage > 1 & median_coverage < 50) %>%
  inner_join(syn_depth_long %>% filter(property == 'observed'), by=c('median_coverage', 'n')) %>%
  rename(expected = num.x, observed = num.y) %>%
  mutate(obs_to_exp = observed/expected) %>%
  group_by(n) %>%
  do(tidy(lm(obs_to_exp ~ log(median_coverage), data=.))) -> depth_lms
ggplot(depth_lms) + aes(x = n, y = estimate) + geom_line() + facet_grid(term ~ ., scales = 'free')

constraint_data = tbl_df(read.table('data/constraint.txt.bgz', header=T))
table(constraint_data$variant_count == constraint_data$variant_count_n123136)

constraint_long = constraint_data %>%
  gather(property, num, -annotation, -transcript) %>%
  separate(property, c('property', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 123135, as.numeric(n)))

constraint_long %>% 
  group_by(annotation, property, n) %>% 
  summarise(num=sum(num)) -> all_transcripts

subset(all_transcripts, annotation == 'synonymous_variant' & (property %in% c('variant_count', 'expected_variant_count_adj'))) %>% 
  ggplot + aes(x=n, y=num, col=property) + geom_line(size=1) + xlab('Number of individuals') + ylab('Variant Count')
ggslackr(channels='#constraint', file='observed and expected by sample size ')

all_transcripts %>%
  ungroup %>%
  filter(annotation %in% c('synonymous_variant', 'missense_variant', 'stop_gained') & property %in% c('variant_count', 'expected_variant_count_adj')) %>% 
  mutate(property = factor(property, levels=c('variant_count', 'expected_variant_count_adj'))) %>%
  ggplot + aes(x=n, y=num, col=annotation, linetype=property) + geom_line(size=1) + xlab('Number of individuals') + ylab('Variant Count')

library(ggjoy)
obs_exp = constraint_long %>% filter(property == 'variant_count') %>%
  inner_join(constraint_long %>% filter(property == 'expected_variant_count_adj'), by=c('transcript', 'annotation', 'n'))

obs_exp %<>%
  rename(observed = num.x, expected = num.y) %>%
  dplyr::select(-property.x, -property.y)

obs_exp %>%
  mutate(obs_exp = observed/expected) %>%
  filter(n %% 10000 == 0) %>%
  filter(annotation %in% c('synonymous_variant', 'missense_variant', 'stop_gained')) -> temp

temp %>%
  ggplot + aes(x = obs_exp, y = as.factor(n), fill=annotation) +
  scale_fill_manual(values=colors) + geom_joy(alpha=0.5) + xlim(0, 2) +
  xlab('Observed/Expected') + ylab('Number of individuals')

temp %>%
  group_by(n, annotation) %>%
  filter(annotation %in% c('synonymous_variant', 'missense_variant')) %>%
  summarize(mean_obs_exp = mean(obs_exp), sd = 1.96*sd(obs_exp, na.rm=T)/sqrt(n())) %>%
  # summarize(obs_exp = sum(observed)/sum(expected)) %>%
  # ggplot + aes(x = n, y = mean_obs_exp, ymin = mean_obs_exp - sd, ymax = mean_obs_exp + sd, fill=annotation) + 
  # scale_fill_manual(values=colors) + geom_ribbon() + geom_line(size=1) +
  ggplot + aes(x = n, y = sd, color=annotation) + 
  scale_color_manual(values=colors) + geom_line(size=1) +
  ylab('Observed/Expected') + xlab('Number of individuals')

filter(temp, annotation %in% c('synonymous_variant')) %>%
  ggplot + aes(x = obs_exp, y = as.factor(n), fill=annotation) + xlim(c(0, 2)) + 
  geom_joy(alpha=0.5) + xlab('Observed/Expected') + ylab('Number of individuals') +
  scale_fill_manual(values=colors)
  


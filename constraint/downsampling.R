syn_exons = tbl_df(read.table('data/syn.txt.bgz', header=T)) %>% left_join(genes, by=c('transcript' = 'ensembl_transcript_id'))
syn_exons_long = syn_exons %>%
  dplyr::select(transcript, exon, variant_count:median_coverage, external_gene_name) %>%
  gather(type, num, -transcript, -exon, -external_gene_name, -aggregate_mutation_rate) %>%
  separate(type, c('type', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 123135, as.numeric(n)))
syn_exons_long %>% 
  group_by(type, n) %>% 
  summarise(num=sum(num)) -> all_exons
all_exons %>%
  filter(type == 'variant_count') %>%
  ggplot + aes(x = n, y = num, col=type) + geom_line()

syn_exons_long %>%
  filter(type == 'variant_count') %>%
  group_by(n) %>%
  do(tidy(lm(num ~ aggregate_mutation_rate, data=.))) -> exon_lms
ggplot(exon_lms) + aes(x = n, y = estimate, col = term) + geom_line()

syn_depth = tbl_df(read.table('data/syn_depth_explore.txt.bgz', header=T))
syn_depth_long = syn_depth %>%
  gather(type, num, -median_coverage) %>%
  separate(type, c('type', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 124000, as.numeric(n)))
syn_depth_long %>%
  filter(type == 'expected' & median_coverage > 1 & median_coverage < 50) %>%
  inner_join(syn_depth_long %>% filter(type == 'observed'), by=c('median_coverage', 'n')) %>%
  rename(expected = num.x, observed = num.y) %>%
  mutate(obs_to_exp = observed/expected) %>%
  group_by(n) %>%
  do(tidy(lm(obs_to_exp ~ log(median_coverage), data=.))) -> depth_lms
ggplot(depth_lms) + aes(x = n, y = estimate) + geom_line() + facet_grid(term ~ ., scales = 'free')

constraint_data = tbl_df(read.table('data/constraint.txt.bgz', header=T))
table(constraint_data$variant_count == constraint_data$variant_count_n123136)

constraint_long = constraint_data %>%
  gather(type, num, -annotation, -transcript) %>%
  separate(type, c('type', 'n'), '_n') %>%
  mutate(n = ifelse(is.na(n), 123135, as.numeric(n)))

constraint_long %>% 
  group_by(annotation, type, n) %>% 
  summarise(num=sum(num)) -> all_transcripts

subset(all_transcripts, annotation == 'synonymous_variant' & (type %in% c('variant_count', 'expected_variant_count_adj'))) %>% 
  ggplot + aes(x=n, y=num, col=type) + geom_line(size=1) + xlab('Number of individuals') + ylab('Variant Count')
ggslackr(channels='#constraint', file='observed and expected by sample size ')

all_transcripts %>%
  ungroup %>%
  filter(annotation %in% c('synonymous_variant', 'missense_variant', 'stop_gained') & type %in% c('variant_count', 'expected_variant_count_adj')) %>% 
  mutate(type = factor(type, levels=c('variant_count', 'expected_variant_count_adj'))) %>%
  ggplot + aes(x=n, y=num, col=annotation, linetype=type) + geom_line(size=1) + xlab('Number of individuals') + ylab('Variant Count')

library(ggjoy)
obs_exp = constraint_long %>% filter(type == 'variant_count') %>%
  inner_join(constraint_long %>% filter(type == 'expected_variant_count_adj'), by=c('transcript', 'annotation', 'n'))

obs_exp %<>%
  rename(observed = num.x, expected = num.y) %>%
  dplyr::select(-type.x, -type.y)

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
  


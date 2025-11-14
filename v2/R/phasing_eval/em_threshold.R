
# Some stats
pbt %>%
  group_by(
    has_em=!is.na(em),
    AC1_gt1=AC1>0,
    AC2_gt1=AC2>0
  ) %>%
  summarise(
    n=n()
  )


#### Evaluation of intervals ####

pbt_by_interval = pbt %>%
  filter(!discordant_between_pops) %>%
  group_by(
    pop,
    locus1.contig
  ) %>%
  filter(AC1>5 & AC2>5) %>%
  mutate(
    bin=ntile(locus1.position, 1000)
  ) %>%
  group_by(
    pop,
    locus1.contig,
    bin
  ) %>%
  summarise(
    n=n(),
    pos=min(locus1.position),
    prop_chet=sum(trio_chet)/n(),
    n_bad_em=sum(em > 0.95 & !trio_chet)
  )

pbt_by_interval %>%
  ggplot(aes( prop_chet)) + geom_histogram()+ facet_grid(pop ~ ., scales="free")

pbt_by_interval2 = pbt %>%
  filter(!discordant_between_pops) %>%
  group_by(
    pop,
    locus1.contig
  ) %>%
  # filter(AC1>5 & AC2>5) %>%
  mutate(
    bin=ntile(locus1.position, 1000)
  ) %>%
  group_by(
    pop,
    locus1.contig,
    bin
  ) %>%
  summarise(
    n=n(),
    pos=min(locus1.position),
    prop_chet=sum(trio_chet)/n(),
    n_bad_em=sum(em > 0.95 & !trio_chet)
  )

#### PR plots ####

# pr_data = get_pr_data(pbt, 'em')
get_pr_plot(
  pbt %>%
    filter(locus1.contig > 19),
  c("em", "rf_prediction"),
  min_ac=2
)

get_pr_plot(
  pbt,
  c("em", "rf_prediction"),
  min_ac=2
)

#### Picking cutoff ####

pbt_binned_metrics =get_pbt_binned_metrics(
  pbt %>%
    filter(pop != 'all') %>%
    # filter(AC1>1 & AC2>1) %>%
    filter(!discordant_between_pops),
  c("em", "rf_prediction"),
  nbins=100
)
plot_pbt_binned_metrics(pbt_binned_metrics)
plot_cumul_pbt_binned_metrics(pbt_binned_metrics)

pbt_binned_metrics =get_pbt_binned_metrics(
  pbt %>%
    filter(pop != 'all') %>%
    filter(AC1>5 & AC2>5) %>%
    filter(!discordant_between_pops),
  c("em", "rf_prediction")
)
plot_pbt_binned_metrics(pbt_binned_metrics)

pbt_binned_metrics =get_pbt_binned_metrics(
  pbt %>%
    filter(pop != 'all') %>%
    filter(AC1>1 & AC2>1) %>%
    filter((AC1 < 3 & AC2 < 10) | (AC1 < 10 & AC2 < 3)) %>%
    filter(!discordant_between_pops),
  c("em", "rf_prediction")
)
plot_pbt_binned_metrics(pbt_binned_metrics)

pbt_binned_metrics =get_pbt_binned_metrics(
  pbt %>%
    filter(pop != 'all') %>%
    filter(AC1>1 & AC2>1) %>%
    filter((AC1 >= 3 | AC2 >= 10) & (AC1 >= 10 | AC2 >= 3)) %>%
    filter(!discordant_between_pops),
  c("em", "rf_prediction")
)
plot_pbt_binned_metrics(pbt_binned_metrics)

# Same hap threshold
same_hap_threshold = pbt_binned_metrics %>% 
  group_by(metric) %>%
  filter(prop_trio_chet < 0.1) %>%
  summarise(
    same_hap_threshold = max(score),
    overall_prop_trio_chet = mean(prop_trio_chet),
    nbins = n()
  )
same_hap_threshold

# Compound het threshold
chet_threshold = pbt_binned_metrics %>% 
  group_by(metric) %>%
  filter(prop_trio_chet > 0.9) %>%
  summarise(
    chet_threshold = min(score),
    overall_prop_trio_chet = mean(prop_trio_chet),
    nbins = n()
  )
chet_threshold

plot_binned_metrics_with_threshold(
  pbt_binned_metrics,
  chet_threshold,
  same_hap_threshold
)
same_hap_threshold
#### Metrics w.r.t. AC ####

# Max 

# Create the bins
pbt_by_ac = get_pbt_binned_by_ac(pbt)

  # mutate( # Allows for better scales in plots
  #   chet_recall = ifelse(chet_recall == 0, NA_real_, chet_recall),
  #   same_hap_recall = ifelse(same_hap_recall == 0, NA_real_, same_hap_recall)
  # )

# Prop trio chet
pbt_by_ac %>%
  # filter(pop == "all") %>%
  filter(n > 29) %>%
  filter(metric == "em") %>% # No need to double-count
  ggplot(aes(ac1_bin, ac2_bin, fill=log10(n))) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle=90))

# Prop trio chet
pbt_by_ac %>%
  # filter(pop == "all") %>%
  filter(n > 29) %>%
  filter(metric == "em") %>% # No need to double-count
  ggplot(aes(ac1_bin, ac2_bin, fill=prop_trio_chet)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle=90))


# precision and recall
pbt_by_ac %>%
  filter(as.numeric(ac1_bin) > 1 & as.numeric(ac2_bin) > 1) %>%
  filter(n > 29) %>%
  ggplot(aes(ac1_bin, ac2_bin, fill=chet_precision)) + 
  geom_tile() +
  facet_grid(.~metric) + 
  theme(axis.text.x = element_text(angle=90))


pbt_by_ac %>%
  filter(as.numeric(ac1_bin) > 1 & as.numeric(ac2_bin) > 1) %>%
  # filter(ac1_bin != "AC = 0" & ac2_bin != "AC = 0") %>% 
  filter(n > 29) %>%
  ggplot(aes(ac1_bin, ac2_bin, fill=chet_recall)) + 
  geom_tile() +
  facet_grid(.~metric) + 
  theme(axis.text.x = element_text(angle=90))

pbt_by_ac %>%
  filter(as.numeric(ac1_bin) > 1 & as.numeric(ac2_bin) > 1) %>%
  filter(n > 29) %>%
  ggplot(aes(ac1_bin, ac2_bin, fill=same_hap_precision)) + 
  geom_tile() +
  facet_grid(.~metric) + 
  theme(axis.text.x = element_text(angle=90))

pbt_by_ac %>%
  filter(as.numeric(ac1_bin) > 1 & as.numeric(ac2_bin) > 1) %>%
  filter(n > 29) %>%
  ggplot(aes(ac1_bin, ac2_bin, fill=same_hap_recall)) + 
  geom_tile() +
  facet_grid(.~metric) + 
  theme(axis.text.x = element_text(angle=90))

#### Very low ACs ####

get_pr_plot(
  pbt %>%
    filter(AC1 == 1 & AC2 == 1),
  c("em", "rf_prediction")
)

pbt %>%
  filter(AC1 == 1 & AC2 == 1) %>%
  group_by(
    pop,
    trio_chet
  ) %>%
  summarise(
    n=n()
  )




pbt_binned %>% View()
pbt_binned %>% 
  filter(
    score > 0.5 & prop_trio_chet < 0.75
  )

x = pbt %>%
  filter(pop == 'all') %>%
  filter(AC1>1 & AC2>1) %>%
  filter(!discordant_between_pops) %>%
  mutate(
    bin = ntile(em, 100)
  ) %>%
  filter(bin==97)

x %<>%
  mutate(
    pos_bin = ntile(locus1.position, 50)
  ) %>%
  group_by(pos_bin) %>%
  summarise(
    n=n(),
    contig = min(locus1.contig),
    pos = min(locus1.position),
    prop_chet=sum(trio_chet)/n()   
  ) 

x %>%
  ggplot(aes(pos, prop_chet)) + geom_point()


x %>% View()

table(x$trio_chet)



em_binned_pbt = pbt %>%
  filter(pop != 'all') %>%
  # filter(AC1>1 & AC2>1) %>%
  filter(!discordant_between_pops) %>%
  group_by(pop) %>%
  mutate(
    bin = ntile(em, 100)
  ) %>%
  group_by(pop, bin) %>%
  summarise(
    n=n(),
    score=min(em),
    prop_trio_chet=sum(trio_chet)/n()
  ) 
em_binned_pbt %>%
ggplot(aes(score, prop_trio_chet, color=pop)) + geom_point()

pbt %>%
  filter(pop != 'all') %>%
  filter(AC1>1 & AC2>1) %>%
  filter(!discordant_between_pops) %>%
  group_by( pop) %>%
  mutate(
    bin = ntile(rf_prediction, 100)
  ) %>%
  group_by(pop, bin) %>%
  summarise(
    n=n(),
    score=min(rf_prediction),
    prop_trio_chet=sum(trio_chet)/n()
  ) %>%
  ggplot(aes(score, prop_trio_chet, color=pop)) + geom_point()

pbt %>%
  ggplot(aes(em)) + geom_density()

pbt %>%
  ggplot(aes(rf_prediction)) + geom_density()

# Plots above suggest thresholds:
same_hap_threshold = 0.0294 # Correspond to start of bins with >90% correct same hap phasing
chet_threshold = 0.7 # Corresponds to start of bins with >90% correct chet phasing



rf_pbt = pbt %>%
  filter(pop != 'all') %>%
  filter(AC1>0 & AC2>0) %>%
  mutate(
    bin = ntile(rf_prediction, 100)
  ) %>%
  group_by(bin) %>%
  summarise(
    n=n(),
    score=min(rf_prediction),
    prop_trio_chet=sum(trio_chet)/n()
  ) 

rf_pbt %>%
  ggplot(aes(score, prop_trio_chet)) + geom_point()

pbt %>%
  ggplot(aes(rf_prediction)) + geom_density()

#### TEST ####
bad_em = pbt %>%
  filter(em > 0.95) %>%
  filter(!trio_chet)

bad_em_by_ac = get_pbt_binned_by_ac(bad_em)
bad_em_by_ac %>%
  ggplot(aes(ac1_bin, ac2_bin, fill=n)) + geom_tile()


best_em_by_ac = pbt %>%
  filter(em > 0.95) %>%
  get_pbt_binned_by_ac()
  
best_em_by_ac %>%
  # filter(pop == "all") %>%
  filter(n > 29) %>%
  ggplot(aes(ac1_bin, ac2_bin, fill=chet_precision)) + 
  geom_tile() +
  facet_grid(.~metric) + 
  theme(axis.text.x = element_text(angle=90))

# Look at RF score distribution for chets with AC < 3
pbt %>%
  filter((AC1 < 3 & AC2 < 5) | (AC2 < 3 & AC1 < 5)) %>%
  filter(trio_chet) %>%
  ggplot(aes(rf_prediction)) + geom_histogram()

pbt %>%
  filter((AC1 < 3 & AC2 < 5) | (AC2 < 3 & AC1 < 5)) %>%
  filter(trio_chet) %>%
  filter(em > 0.95) %>%
  ggplot(aes(rf_prediction)) + geom_histogram()

test = glm(
  trio_chet ~ distance + cpg1 + cpg2 + discordant_between_pops + snv1 + snv2,
  pbt %>%
    filter(em > 0.95),
  family = binomial
)

pbt %>% 
  filter(em == 1.0 & !trio_chet & AC1 > 100 & AC2 > 100) %>% 
  View()
  


  
  group_by(ac1_bin, ac2_bin) %>% 
  summarise(n=n()) %>% View()


  
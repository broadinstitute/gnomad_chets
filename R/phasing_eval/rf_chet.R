library(randomForest)

# Load PBT data

pbt = read_csv("pbt_annotated.csv") 
pbt %<>% select(-X1)
pbt %<>% filter(!is.na(trio_chet))

table(pbt$locus1.contig)

# Train RF model holding off chrom 20-22
pbt %<>% mutate(
  trio_chet_f = as.factor(trio_chet),
  log_distance = log(distance)
  )
rf = randomForest(
  trio_chet ~ n00 + n01 + n02 + n10 + n11 + n12 + n20 + n21 + n22 + log_distance + cpg1 + cpg2 + snv1 + snv2,
  data=pbt %>%
    filter(locus1.contig < 20),
  sampsize=750,
  strata=trio_chet_f,
  importance=T,
  ntree=200
)

# Evaluate on chrom 20-22
pbt_eval = pbt %>%
  filter(locus1.contig >= 20) 
pbt_eval$prediction = predict(object=rf, pbt_eval)
pbt_eval %>%
  group_by(
    trio_chet_f
  ) %>%
  summarise(
    false=sum(prediction<0.5),
    true=sum(prediction>0.5),
    class.error=sum((prediction > 0.5) != trio_chet_f) /n()
  )

# rf2 = randomForest(
#   trio_chet_f ~ n00 + n01 + n02 + n10 + n11 + n12 + n20 + n21 + n22 + distance + cpg1 + cpg2,
#   data=pbt %>%
#     filter(locus1.contig < 20),
#   sampsize=750,
#   strata=trio_chet_f,
#   importance=T,
#   ntree=100
# )

pbt$rf_prediction = predict(rf, pbt)

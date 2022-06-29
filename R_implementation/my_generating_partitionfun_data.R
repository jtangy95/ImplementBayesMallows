library(dplyr)
load("/Users/changtaeyeong///Desktop/BayesMallowsRankModel/BayesMallows/data-raw/footrule_cardinalities.Rdata", verbose = T)
# exact sequence for cardinalities available given footrule distance
partition_function_data <- tibble(
  n_items = seq_along(seq2),
  metric = rep("footrule", length(seq2)),
  values = seq2
)
rm(seq2)

load("~/Desktop/BayesMallowsRankModel/BayesMallows/data-raw/spearman_cardinalities.Rdata", verbose = T)
# exact sequence for cardinalities available given spearman distance
partition_function_data <- tibble(
  n_items = seq_along(seq2),
  metric = rep("spearman", length(seq2)),
  values = seq2
) %>%
  bind_rows(partition_function_data)
rm(seq2)

## set working directory to use relative path
setwd("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/R_implementation")

save(partition_function_data, file = "./my_partition_function_data.RData")



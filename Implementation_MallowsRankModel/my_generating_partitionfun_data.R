library(dplyr)
load("~/Desktop/BayesMallowsRankModel/BayesMallows/data-raw/footrule_cardinalities.Rdata", verbose = T)
partition_function_data <- tibble(
  n_items = seq_along(seq2),
  metric = rep("footrule", length(seq2)),
  values = seq2
)
rm(seq2)
load("~/Desktop/BayesMallowsRankModel/BayesMallows/data-raw/spearman_cardinalities.Rdata", verbose = T)
partition_function_data <- tibble(
  n_items = seq_along(seq2),
  metric = rep("spearman", length(seq2)),
  values = seq2
) %>%
  bind_rows(partition_function_data)
rm(seq2)


save(partition_function_data, file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_partition_function_data.RData")


library(BayesMallows)
data(package = "BayesMallows")
data("potato_true_ranking")
data("potato_visual")
data("potato_weighing")

potatoes_rank = list("visual_assess" = potato_visual, "weighing_assess" = potato_weighing, "true_rank" = "potato_true_ranking")

save(potatoes_rank, file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_potatoes.RData")

data("beach_preferences")

save(beach_preferences, file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_beach.RData")

data("sushi_rankings")

save(sushi_rankings, file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_sushi.RData")

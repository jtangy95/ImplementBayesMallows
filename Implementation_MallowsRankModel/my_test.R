source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_partition_function.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_tidy_mcmc.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_mallows.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_plot.BayesMallows.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_consensus.R")

load(file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_potatoes.RData", verbose = T)


fit = compute_mallows(potatoes_rank$visual_assess, verbose = F, nmc = 100000)
fit_is = compute_mallows(potatoes_rank$visual_assess, verbose = F, nmc = 100000, importance_sampling = T)

fit$rho
fit$alpha
fit$rho_acceptance
fit$alpha_acceptance

fit_is$rho
fit_is$alpha
fit_is$rho_acceptance
fit_is$alpha_acceptance

fit$burnin = 1000
plot(fit, parameter = "alpha")
plot(fit, parameter = "rho", items = c(1,3,5,7,9,11))
plot(fit, parameter = "rho", items = c("P1", "P10"))
compute_consensus(fit, parameter = "rho" , type = "CP")
compute_consensus(fit, parameter = "rho" , type = "MAP")

load(file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_sushi.RData", verbose = T)

fit = compute_mallows(sushi_rankings, nmc = 100000, verbose = F)
fit_is = compute_mallows(sushi_rankings, nmc = 100000, verbose = F, importance_sampling = T)

fit$rho
fit$alpha
fit$rho_acceptance
fit$alpha_acceptance

fit_is$rho
fit_is$alpha
fit_is$rho_acceptance
fit_is$alpha_acceptance

fit_is$burnin = 1000
plot(fit_is, parameter = "alpha")
plot(fit_is, parameter = "rho", items = c(1,3,5,7))
plot(fit_is, parameter = "rho", items = c("shrimp", "egg", "squid"))
compute_consensus(fit_is, parameter = "rho" , type = "CP")
compute_consensus(fit_is, parameter = "rho" , type = "MAP")


source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_mallows_partial.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_tidy_mcmc_partial.R")
library(dplyr)

create_miss = function(rankings, fix = NULL, lambda = NULL){
  N = nrow(rankings) 
  n = ncol(rankings)
  if(!is.null(fix)){
    for(j in 1:N){
      rankings[j, which(rankings[j,] %in% seq(from = n, to = n - (fix - 1), by = -1) )] = NA
    }
  } else if(!is.null(lambda)){
    miss = extraDistr::rtpois(N, lambda, a = 0, b= n)
    for(j in 1:N){
      rankings[j, which(rankings[j,] %in% seq(from = n, to = n - (miss[j] - 1), by = -1) )] = NA
    }
  } 
  return(rankings)
}

potato = potatoes_rank$visual_assess
potato_miss1 = create_miss(potatoes_rank$visual_assess, fix = 5)
potato_miss2 = create_miss(potatoes_rank$visual_assess, lambda = 10)

fit = compute_mallows_partial(potato_miss1, verbose =F, nmc = 100000, save_aug = T)
fit$burnin = 1000
plot(fit, parameter = "rho", items = 1:10)
compute_consensus(fit, parameter = "rho", type = "CP")
fit$augmented_data %>% filter(assessor == 3, iteration == 100)
fit$aug_acceptance
compute_consensus(fit, parameter = "Rtilde", type = "CP", assessors = 2)
compute_consensus(fit, parameter = "Rtilde", type = "MAP", assessors = 2)
compute_consensus(fit, parameter = "Rtilde", type = "CP", assessors = c(2,3))
compute_consensus(fit, parameter = "Rtilde", type = "MAP", assessors = c(2,3))



source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_topk.R")
fit = compute_mallows_partial(potato_miss2, verbose =F, nmc = 100000, save_aug = T)
fit$burnin = 1000
predict_top_k(fit, k = 10) %>% filter(assessor ==1)
predict_top_k(fit, k = 10) %>% filter(assessor ==10)
compute_consensus(fit, parameter = "rho", type = "CP")
plot_top_k(fit, k=10)

sushi_miss = create_miss(sushi_rankings, fix = 5)
fit = compute_mallows_partial(sushi_miss, verbose = T, nmc = 10000, save_aug = T, aug_thinning = 10)
fit$burnin = 1000
predict_top_k(fit, k = 7) %>% filter(assessor ==1)
predict_top_k(fit, k = 7) %>% filter(assessor ==10)
compute_consensus(fit, parameter = "rho", type = "CP")

source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_mallows_cluster.R")

fit = compute_mallows_cluster(potato, verbose = T, nmc = 1000, n_clusters = 2)
fit$rho[,,1]
fit$rho[,,100]
fit$rho_acceptance
fit$alpha_acceptance
fit$cluster_probs[, 900]
fit$cluster_assignment[, 10:30]
fit$cluster_assignment[, 80:100]

fit = compute_mallows_cluster(sushi_rankings, nmc = 10000, verbose = T, n_clusters = 3)
fit$rho[,,1]
fit$rho[,,100]
fit$rho_acceptance
fit$alpha_acceptance
fit$cluster_probs[, 900]
fit$cluster_assignment[1:40, 10:30]
fit$cluster_assignment[1:40, 80:100]
fit$cluster_assignment[1:40, 1000:1020]
fit$cluster_assignment[1:40, 3000:3020]
fit$cluster_assignment[1:40, 9000:9020]

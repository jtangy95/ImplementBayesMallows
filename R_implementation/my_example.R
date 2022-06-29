## set working directory to use relative path
setwd("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/R_implementation")

source("./my_compute_mallows_combined.R")
source("./my_plot.BayesMallows_cluster.R")
source("./my_compute_consensus.R")
source("./my_topk.R")
load(file = "./dataset_potatoes.RData", verbose = T)
load(file = "./dataset_sushi.RData", verbose = T)
load(file = "./dataset_beach.RData", verbose = T)
library(dplyr)

fit = compute_mallows(potatoes_rank$visual_assess, verbose = T, nmc = 10000)
fit_is = compute_mallows(potatoes_rank$visual_assess, verbose = T, nmc = 10000, importance_sampling = T)

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


fit = compute_mallows(sushi_rankings, nmc = 10000, verbose = T)
fit_is = compute_mallows(sushi_rankings, nmc = 10000, verbose = T, importance_sampling = T)

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


fit = compute_mallows(sushi_rankings, nmc = 10000, verbose = T, n_clusters = 5, clus_thin = 10, rho_thinning = 10)

fit$rho %>% filter(iteration == 1, cluster == "Cluster 1")
fit$rho %>% filter(iteration == 1, cluster == "Cluster 2")
fit$rho %>% filter(iteration == 1, cluster == "Cluster 3")

fit$rho %>% filter(iteration == 501, cluster == "Cluster 1")
fit$rho %>% filter(iteration == 501, cluster == "Cluster 2")
fit$rho %>% filter(iteration == 501, cluster == "Cluster 3")

fit$rho_acceptance
fit$alpha_acceptance

fit$cluster_probs
fit$cluster_probs %>% filter(iteration > 500 & iteration < 530)
fit$cluster_assignment 
fit$cluster_assignment %>% filter(iteration > 10 & iteration < 100) %>% filter(assessor == 1)
fit$cluster_assignment %>% filter(iteration > 10 & iteration < 100) %>% filter(assessor == 1000)
fit$cluster_assignment %>% filter(iteration > 10 & iteration < 100) %>% filter(assessor == 3000)
fit$cluster_assignment %>% filter(iteration > 1000 & iteration < 1100) %>% filter(assessor == 1)
fit$cluster_assignment %>% filter(iteration > 1000 & iteration < 1100) %>% filter(assessor == 1000)

fit$within_cluster_distance
fit$within_cluster_distance %>% filter(iteration == 10000)

fit$burnin = 1000
plot(fit, parameter = "alpha")
plot(fit, parameter = "rho", items = c(1,3,5))
plot(fit, parameter = "rho", items = c("shrimp", "egg", "squid"))
plot(fit, parameter = "cluster_probs")
plot(fit, parameter = "cluster_assignment")

compute_consensus(fit, type = "CP", parameter = "rho") %>% filter(cluster == "Cluster 3")
compute_consensus(fit, type = "MAP", parameter = "rho") %>% filter(cluster == "Cluster 3")


fit = compute_mallows(preferences = beach_preferences, nmc = 1e4, verbose = T,
                      leap_size = 2, alpha_prop_sd = 0.1, lambda = 0.1, alpha_jump =  100, save_aug = T, aug_thinning = 100)
fit$rho
fit$alpha
fit$rho_acceptance
fit$alpha_acceptance
fit$save_aug
fit$burnin = 1e3
compute_consensus(fit, parameter = "rho", type = "CP")
fit$augmented_data %>% filter(assessor == 3, iteration == 1001)
fit$augmented_data %>% filter(assessor == 3, iteration == 2001)
fit$aug_acceptance
compute_consensus(fit, parameter = "Rtilde", type = "CP", assessors = 2)
compute_consensus(fit, parameter = "Rtilde", type = "MAP", assessors = 5)
compute_consensus(fit, parameter = "Rtilde", type = "CP", assessors = c(2,3))
compute_consensus(fit, parameter = "Rtilde", type = "MAP", assessors = c(2,3))
predict_top_k(fit, k = 3) %>% filter(assessor ==1)
predict_top_k(fit, k = 3) %>% filter(assessor ==10)
plot_top_k(fit, k=3, assessors = 1:20)
plot_top_k(fit, k=3, assessors = c(1,3,5,7,9))

source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/BayesMallows/R/compute_mallows_mixtures.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/BayesMallows/R/plot_elbow.R")
library(parallel)
(nCores=detectCores()-1)
cl<-makeCluster(nCores)
mixtures = compute_mallows_mixtures(n_clusters=3:5, rankings=sushi_rankings, nmc=10000, 
                                    rho_thinning = 10, clus_thin = 10, include_wcd=T)
plot_elbow(mixtures, burnin=1000)



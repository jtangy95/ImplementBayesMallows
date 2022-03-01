library(Rcpp)
library(RcppArmadillo)
sourceCpp('/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_sample_mallows.cpp')

validate_permutation <- function(vec){
  return(all(sort(vec) == seq_along(vec)))
}

sample_mallows <- function(rho0, alpha0, n_samples,
                           leap_size = max(1L, floor(length(rho0)/5)),
                           metric = "footrule",
                           burnin = 1000,
                           thinning = 1000){
  if(!(validate_permutation(rho0) && !any(is.na(rho0)) )){
    stop("rho0 must be a proper ranking with no missing values.")
  }
  
  if(n_samples <= 0){
    stop("n_samples must be positive.")
  }
  
  n_items <- length(rho0)
  
  samples <- t(rmallows(
    rho0 = rho0,
    alpha0 = alpha0,
    n_samples = n_samples,
    burnin = burnin,
    thinning = thinning,
    leap_size = leap_size,
    metric = metric
  ))
   
  colnames(samples) <- paste("Item", 1:n_items)
  
  return(samples)
}

sample_mallows_mixed <- function(rho0, alpha0, n_samples, n_clusters = 1L, 
                           leap_size = max(1L, floor(nrow(rho0)/5)),
                           psi = n_samples / n_clusters,
                           metric = "footrule",
                           burnin = 1000,
                           thinning = 1000){
  if(!(all(apply(rho0, 2, validate_permutation)))){
    stop("each column vector of rho0 must be a proper ranking")
  }
  if(any(is.na(rho0))){
    stop("rho0 must have no missing values.")
  }
  
  if(n_samples <= 0){
    stop("n_samples must be positive.")
  }
  if(ncol(rho0)!= n_clusters){
    stop("number of columns of rho0 should be each to n_clusters")
  }
  if(length(alpha0)!= n_clusters){
    stop("number of elements of alpha0 should be each to n_clusters")
  }
  
  n_items <- nrow(rho0)
  
  tau = rgamma(n_clusters, psi)
  tau = tau / sum(tau)
  z = sample(1:n_clusters, n_samples, replace = T, prob = tau)
  Nc = as.numeric(table(z))
  cumNc = c(0, cumsum(Nc))
  samples <- matrix(0, nrow = n_samples, ncol = n_items)
  for(i in 1:n_clusters){
    cluster_sample <- t(rmallows(
      rho0 = rho0[,i],
      alpha0 = alpha0[i],
      n_samples = Nc[i],
      burnin = burnin,
      thinning = thinning,
      leap_size = leap_size,
      metric = metric
    ))
    samples[(cumNc[i]+1):cumNc[i+1],] = cluster_sample
  }
  
  colnames(samples) <- paste("Item", 1:n_items)
  
  return(samples)
}

# n_item = 10
# rho0 = sample(1:n_item, n_item)
# n_samples = 40
# alpha0 = 1
# samples = sample_mallows(rho0, alpha0,n_samples)
# samples
# 
# n_clusters = 5
# 
# rho0 = matrix(0, ncol= n_clusters, nrow = n_item)
# alpha0 = rep(0, n_clusters)
# for(i in 1:n_clusters){
#   rho0[, i] = sample(1:n_item, n_item)
#   alpha0[i] = rexp(n=1 ,rate = 3)
# }
# samples_cluster = sample_mallows_mixed(rho0, alpha0, n_samples, n_clusters)
# samples_cluster

create_topk = function(rankings, k){
  N = nrow(rankings) 
  n = ncol(rankings)

  for(j in 1:N){
    rankings[j, which(rankings[j,] %in% (k+1):n )] = NA
  }

  return(rankings)
}

# samples
# top5samples = create_topk(samples, 5)
# top5samples


library(Rcpp)
library(RcppArmadillo)

# validate_permutation prepared also for incomplete ranking
validate_permutation <- function(vec){
  if(!any(is.na(vec))){
    # complete ranking
    return(all(sort(vec) == seq_along(vec)))
  } else if(all(is.na(vec))){
    # totally missing rank
    return(TRUE)
  } else {
    # partial ranking
    return(all(vec[!is.na(vec)] <= length(vec)) &&
             all(vec[!is.na(vec)] >= 1) && !any(duplicated(vec[!is.na(vec)])))
  }
}

## set working directory to use relative path
setwd("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/R_implementation")

sourceCpp('./my_run_mcmc_partial.cpp')

compute_mallows_partial <- function(rankings ,
                            metric = "footrule",
                            nmc = 10000L,
                            leap_size = max(1L, floor(ncol(rankings) / 5)),
                            rho_init = NULL,
                            rho_thinning = 1L,
                            alpha_prop_sd = 0.1,
                            alpha_init = 1,
                            alpha_jump = 1L,
                            lambda = 0.001,
                            alpha_max = 1e6,
                            logz_estimate = NULL,
                            importance_sampling = FALSE,
                            is_nmc = 10000L,
                            aug_thinning = 1L,
                            save_aug = FALSE,
                            verbose = FALSE,
                            validate_rankings = TRUE,
                            seed = NULL
){
  
  if(!is.null(seed)) set.seed(seed)
  
  # Check that at most one of rankings and preferences is set
  if(is.null(rankings)){
    stop("Rankings must be provided.")
  }
  
  if(nmc <= 0) stop("nmc must be strictly positive")
  
  # Check that we do not jump over all alphas
  if(alpha_jump >= nmc) stop("alpha_jump must be strictly smaller than nmc")
  
  # Check that we do not jump over all rhos or augmented ranks
  if(rho_thinning >= nmc) stop("rho_thinning must be strictly smaller than nmc")
  if(aug_thinning >= nmc) stop("aug_thinning must be strictly smaller than nmc")
  
  if(lambda <= 0) stop("exponential rate parameter lambda must be strictly positive")
  
  # Check that all rows of rankings are proper permutations
  if(!is.null(rankings) && validate_rankings && !all(apply(rankings, 1, validate_permutation))){
    stop("invalid permutations provided in rankings matrix")
  }
  
  # Find the number of items
  n_items <- ncol(rankings)
  
  
  if(!is.null(rho_init)) {
    if(!validate_permutation(rho_init)) stop("rho_init must be a proper permutation")
    if(any(is.na(rho_init))) stop("rho_init cannot have missing values")
    if(length(rho_init) != n_items) stop("rho_init must have the same number of items as implied by rankings or preferences")
    rho_init <- matrix(rho_init, ncol = 1)
  }
  
  exact_available = ifelse(  (metric == "footrule" && n_items > 50)||(metric == "spearman" && n_items > 14) , FALSE, TRUE )
  if(verbose){
    print("Before running Metropolis-Hastings algorithm , yielding partition function values...")
  }
  
  logz_list <- prepare_partition_function(importance_sampling, exact_available, logz_estimate, metric, n_items, is_nmc)
  
  
  # Transpose rankings to get samples along columns, since we typically want
  # to extract one sample at a time. armadillo is column major, just like rankings
  fit <- run_mcmc_partial(rankings = t(rankings),
                  nmc = nmc,
                  cardinalities = logz_list$cardinalities,
                  logz_estimate = logz_list$logz_estimate,
                  rho_init = rho_init,
                  metric = metric,
                  lambda = lambda,
                  alpha_max = alpha_max,
                  leap_size = leap_size,
                  alpha_prop_sd = alpha_prop_sd,
                  alpha_init = alpha_init,
                  alpha_jump = alpha_jump,
                  rho_thinning = rho_thinning,
                  aug_thinning = aug_thinning,
                  save_aug = save_aug,
                  verbose = verbose
  )
  
  if(verbose){
    print("Metropolis-Hastings algorithm completed. Post-processing data.")
  }
  
  # Change rho_acceptance and alpha_acceptance from frequency to relative frequency
  fit$rho_acceptance <- fit$rho_acceptance/nmc
  fit$alpha_acceptance <- fit$alpha_acceptance/nmc
  fit$aug_acceptance <- fit$aug_acceptance/nmc
  
  # Add some arguments
  fit$metric <- metric
  fit$lambda <- lambda
  fit$nmc <- nmc
  fit$alpha_jump <- alpha_jump
  fit$rho_thinning <- rho_thinning
  fit$leap_size <- leap_size
  fit$alpha_prop_sd <- alpha_prop_sd
  fit$save_aug <- save_aug
  fit$aug_thinning <- aug_thinning
  
  # Add names of item
  if(!is.null(colnames(rankings))) {
    fit$items <- colnames(rankings)
  } else {
    fit$items <- paste("Item", seq(from = 1, to = nrow(fit$rho), by = 1))
  }
  
  fit <- tidy_mcmc_partial(fit)
  
  # Add class attribute
  class(fit) <- "myBayesMallows"
  
  return(fit)
  
}



library(Rcpp)
library(RcppArmadillo)

## set working directory to use relative path
setwd("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/R_implementation")

sourceCpp('./my_importance_sampling.cpp')
load("./my_partition_function_data.RData", verbose = T)

# library(parallel)
# cl <- makeCluster(detectCores() - 1)
# stopCluster(cl)

estimate_partition_function <- function(alpha_vector, n_items, metric,
                                        nmc=10000L, degree = 10,  
                                        seed = NULL){
  
  stopifnot(degree < length(alpha_vector))
  
  # if(!is.null(cl)){
  #   # Split nmc into each cluster
  #   nmc_vec <- rep(floor(nmc / length(cl)), length(cl))
  #   i <- 1
  #   while(sum(nmc_vec) != nmc){
  #     nmc_vec[i] <- nmc_vec[i] + 1
  #     i = i+1
  #     if(i == length(cl)) break()
  #   }
  #   parallel::clusterExport(cl, c("alpha_vector", "n_items", "metric", "seed"),
  #                           envir = environment())
  #   estimates <- parallel::parLapply(cl, nmc_vec, function(x){
  #     if(!is.null(seed)) set.seed(seed)
  #     compute_importance_sampling_estimate(alpha_vector = alpha_vector, n_items = n_items,
  #                                          metric = metric, nmc = x)
  #   })
  #   log_z <- purrr::map2_dbl(estimates[[1]], estimates[[2]], mean)
  # } else {
  
  if(!is.null(seed)) set.seed(seed)
  log_z <- as.vector(
    compute_importance_sampling_estimate(
      alpha_vector = alpha_vector, n_items = n_items,
      metric = metric, nmc = nmc))
  
  # Compute the estimate at each discrete alpha value
  estimate <- dplyr::tibble(
    alpha = alpha_vector,
    log_z = log_z
  )
  
  # Fit a regression model
  form <- stats::as.formula(paste("log_z ~ ",
                                  paste("I( alpha^", seq(from = 1, to = degree, by = 1), ")",
                                        collapse = "+")))
  stats::lm(form, data = estimate)$coefficients
  
}


# library(dplyr)

prepare_partition_function <- function(importance_sampling, exact_available, logz_estimate, metric, n_items, is_nmc){
  # First, has the user supplied an estimate?
  if(!is.null(logz_estimate)){
    return(list(cardinalities = NULL, logz_estimate = logz_estimate))
  }
  
  # Second, does user wants importance_sampling estimates ? OR is an exact seq. for cardinalities unavailable?
  if(importance_sampling || !exact_available){
    alpha_vector <- seq(from = 0.01, to = 10, length = 100)
    logz_estimate = estimate_partition_function(alpha_vector, n_items, metric, nmc = is_nmc)
    return(list(cardinalities = NULL, logz_estimate = logz_estimate))
  }
  
  # Third, is an exact sequence for cardinalities available ?
  if(exact_available){
    met = metric ; n = n_items 
    relevant_params <- unlist(dplyr::filter(partition_function_data, metric == met, n_items == n)$values)
    # relevant_params <- unlist(partition_function_data %>%
    #                             dplyr::filter(metric == met , n_items == n) %>% .$values)
    return(list(cardinalities = relevant_params, logz_estimate = NULL))
  }
  
}


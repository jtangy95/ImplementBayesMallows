tidy_mcmc_partial <- function(fit){
  
  fit <- tidy_rho(fit)
  fit <- tidy_alpha(fit)
  # fit <- tidy_cluster_assignment(fit)
  # fit <- tidy_cluster_probabilities(fit)
  # fit <- tidy_wcd(fit)
  fit <- tidy_augmented_data(fit)
  fit <- tidy_augmentation_acceptance(fit)
  # fit <- tidy_error_probability(fit)
  
  return(fit)
}



tidy_rho <- function(fit){
  # Tidy rho
  n_items = nrow(fit$rho)
  n_stored = ncol(fit$rho)
  # Item1, Item2, Item3, ...., Item1, Item2, Item3
  # Iteration1, Iteration1, ..., Iteration1, Iteration1, Iteration1, Iteration2
  
  value <- as.vector(fit$rho)
  
  item <- rep(fit$items, times = n_stored)
  item <- factor(item, levels = fit$items)
  
  iteration <- rep(seq(from = 1, to = n_stored * fit$rho_thinning, by = fit$rho_thinning),
                   each = n_items)
  
  # Store the final rho as a tibble
  fit$rho <- dplyr::tibble(
    item = item,
    iteration = iteration,
    value = value
  )
  
  return(fit)
}

tidy_alpha <- function(fit){
  # Tidy alpha
  n_stored <- length(fit$alpha)
  
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  value <- as.vector(fit$alpha)
  
  iteration <- seq(from = 1, to = n_stored * fit$alpha_jump, by = fit$alpha_jump)
  
  fit$alpha <- dplyr::tibble(
    iteration = iteration,
    value = value
  )
  
  return(fit)
}

tidy_augmented_data <- function(fit){
  # Tidy augmented data, or delete
  if(fit$save_aug){
    
    # (n_items, n_assessors, n_stored_iterations)
    augdata_dims <- dim(fit$augmented_data)
    
    # Item1, Item2, ..., Item1, Item2, ..., Item1, Item2, ..., Item1, Item2
    # Assessor1, Assessor1, ..., Assessor2, Assessor2, ... Assessor1, Assessor1, ..., Assessor2, Assessor2
    # Iteration1, Iteration1, ..., Iteration1, Iteration1, ..., Iteration2, Iteration2, ... Iteration2, Iteration2
    value <- c(fit$augmented_data)
    
    item <- rep(fit$items, times = augdata_dims[2] * augdata_dims[3])
    item <- factor(item, levels = fit$items)
    
    assessor <- rep(seq(from = 1, to = augdata_dims[2], by = 1), each = augdata_dims[1],
                    times = augdata_dims[3])
    
    iteration <- rep(seq(from = 1, to = augdata_dims[3] * fit$aug_thinning, by = fit$aug_thinning),
                     each = augdata_dims[1] * augdata_dims[2])
    
    fit$augmented_data <- dplyr::tibble(
      iteration = iteration,
      assessor = assessor,
      item = item,
      value = value
    )
  } else {
    fit$augmented_data <- NULL
  }
  
  
  
  return(fit)
}

tidy_augmentation_acceptance <- function(fit){
  # Augmentation acceptance
  
  if(fit$any_missing){
    fit$aug_acceptance <- dplyr::tibble(acceptance_rate = c(fit$aug_acceptance))
    fit$aug_acceptance <- dplyr::mutate(fit$aug_acceptance,
                                        assessor = dplyr::row_number())
    fit$aug_acceptance <- dplyr::select(fit$aug_acceptance,
                                        assessor, acceptance_rate)
    
  } else {
    fit$aug_acceptance <- NULL
  }
  return(fit)
}

tidy_mcmc <- function(fit){
  
  fit <- tidy_rho(fit)
  fit <- tidy_alpha(fit)
  # fit <- tidy_cluster_assignment(fit)
  # fit <- tidy_cluster_probabilities(fit)
  # fit <- tidy_wcd(fit)
  # fit <- tidy_augmented_data(fit)
  # fit <- tidy_augmentation_acceptance(fit)
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


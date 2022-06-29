tidy_mcmc <- function(fit){
  
  fit <- tidy_rho(fit)
  fit <- tidy_alpha(fit)
  fit <- tidy_cluster_assignment(fit)
  fit <- tidy_cluster_probabilities(fit)
  fit <- tidy_wcd(fit)
  fit <- tidy_augmented_data(fit)
  fit <- tidy_augmentation_acceptance(fit)
  
  return(fit)
}



tidy_rho <- function(fit){
  # Tidy rho
  n_items = nrow(fit$rho)
  n_clusters = ncol(fit$rho)
  n_stored = dim(fit$rho)[3]
  # Item1, Item2, Item3, ...., Item1, Item2, Item3
  # Cluster1, Cluster1, Cluster1, ..., Cluster2, Cluster2, Cluster2
  # Iteration1, Iteration1, ..., Iteration1, Iteration1, Iteration1, Iteration2
  
  value <- as.vector(fit$rho)
  
  item <- rep(fit$items, times = n_clusters * n_stored)
  item <- factor(item, levels = fit$items)
  
  cluster <- rep(1:n_clusters, each = n_items, times = n_stored)
  cluster <- cluster <- factor(paste("Cluster", cluster),
                               levels = paste("Cluster", sort(unique(cluster))))
  
  iteration <- rep(seq(from = 1, to = n_stored * fit$rho_thinning, by = fit$rho_thinning),
                   each = n_items * n_clusters)
  
  # Store the final rho as a tibble
  fit$rho <- dplyr::tibble(
    item = item,
    cluster = cluster,
    iteration = iteration,
    value = value
  )
  
  return(fit)
}

tidy_alpha <- function(fit){
  # Tidy alpha
  n_clusters <- nrow(fit$alpha)
  n_stored <- ncol(fit$alpha)
  
  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  value <- as.vector(fit$alpha)
  
  cluster <- rep(1:n_clusters, times = n_stored)
  cluster <- cluster <- factor(paste("Cluster", cluster),
                               levels = paste("Cluster", sort(unique(cluster))))
  
  iteration <- rep(seq(from = 1, to = n_stored * fit$alpha_jump, by = fit$alpha_jump), 
                   each = n_clusters)
  
  fit$alpha <- dplyr::tibble(
    cluster = cluster,
    iteration = iteration,
    value = value
  )
  
  return(fit)
}

tidy_cluster_probabilities <- function(fit){
  
  # Tidy cluster probabilities
  if(fit$n_clusters > 1){
    n_clusters <- nrow(fit$cluster_probs)
    n_stored <- ncol(fit$cluster_probs)
    value <- as.vector(fit$cluster_probs)
    
    # Cluster1, Cluster2, ..., Cluster1, Cluster2
    # Iteration1, Iteration1, ..., Iteration2, Iteration2
    cluster <- rep(1:n_clusters, times = n_stored)
    cluster <- cluster <- factor(paste("Cluster", cluster),
                                 levels = paste("Cluster", sort(unique(cluster))))
    
    iteration <- rep(seq(from = 1, to = n_stored * fit$clus_thin, by = fit$clus_thin),
                     each = n_clusters)
    
    fit$cluster_probs <- dplyr::tibble(
      cluster = cluster,
      iteration = iteration,
      value = value
    )
  } else {
    fit$cluster_probs <- NULL
  }
  
  return(fit)
}

tidy_cluster_assignment <- function(fit){
  
  # Tidy cluster assignment
  
  if(fit$n_clusters > 1){
    n_assessors <- nrow(fit$cluster_assignment)
    n_stored <- ncol(fit$cluster_assignment)
    value <- paste("Cluster", as.vector(fit$cluster_assignment))
    
    # Assessor1, Assessor2, ..., Assessor1, Assessor2
    # Iteration1, Iteration1, ..., Iteration2, Iteration2
    
    assessor <- rep(1:n_assessors, times = n_stored)
    iteration <- rep(seq(from = 1, to = n_stored * fit$clus_thin, by = fit$clus_thin),
                     each = n_assessors)
    
    fit$cluster_assignment <- dplyr::tibble(
      assessor = assessor,
      iteration = iteration,
      value = value
    )
  } else {
    fit$cluster_assignment <- NULL
  }
  
  return(fit)
}

tidy_wcd <- function(fit){
  # Tidy the within-cluster distances, or delete the empty matrix
  
  
  if(fit$include_wcd){
    n_clusters <- nrow(fit$within_cluster_distance)
    n_stored <- ncol(fit$within_cluster_distance)
    value <- as.vector(fit$within_cluster_distance)
    
    # Cluster1, Cluster2, ..., Cluster1, Cluster2
    # Iteration1, Iteration1, ..., Iteration2, Iteration2
    cluster <- rep(1:n_clusters, times = n_stored)
    cluster <- factor(paste("Cluster", cluster),
                      levels = paste("Cluster", sort(unique(cluster))))
    
    iteration <- rep(1:n_stored, each = n_clusters)
    
    fit$within_cluster_distance <- dplyr::tibble(
      cluster = cluster,
      iteration = iteration,
      value = value
    )
    
  } else {
    fit$within_cluster_distance <- NULL
  }
  
  return(fit)
}

tidy_augmented_data <- function(fit){
  # Tidy augmented data, or delete
  if(fit$save_aug){
    
    # (n_items, n_assessors, n_stored_iterations)
    n_items = nrow(fit$augmented_data)
    n_assessors = ncol(fit$augmented_data)
    n_stored = dim(fit$augmented_data)[3]
    
    # Item1, Item2, ..., Item1, Item2, ..., Item1, Item2, ..., Item1, Item2
    # Assessor1, Assessor1, ..., Assessor2, Assessor2, ... Assessor1, Assessor1, ..., Assessor2, Assessor2
    # Iteration1, Iteration1, ..., Iteration1, Iteration1, ..., Iteration2, Iteration2, ... Iteration2, Iteration2
    value <- c(fit$augmented_data)
    
    item <- rep(fit$items, times = n_assessors * n_stored)
    item <- factor(item, levels = fit$items)
    
    assessor <- rep(1:n_assessors, each = n_items, times = n_stored)
    
    iteration <- rep(seq(from = 1, to = n_stored * fit$aug_thinning, by = fit$aug_thinning),
                     each = n_items * n_assessors)
    
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
  
  if(fit$partial_rank|| fit$preference_compare){
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


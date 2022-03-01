compute_consensus <- function(model_fit, type = "CP", burnin = model_fit$burnin,
                              parameter = "rho", assessors = 1L){
  
  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)
  
  stopifnot(class(model_fit) == "myBayesMallows")
  
  if(parameter == "Rtilde" && !inherits(model_fit$augmented_data, "data.frame")){
    stop(paste("For augmented ranks, please refit",
               "model with option 'save_aug = TRUE'."))
  }
  
  if(parameter == "rho"){
    # Filter out the pre-burnin iterations
    
    df = dplyr::filter(model_fit$rho, iteration > burnin)
    
    # Find the problem dimensions
    n_rows <- nrow(dplyr::distinct(df, item))
    
    # Check that there are rows.
    stopifnot(n_rows > 0)
    
    # Check that the number of rows are consistent with the information in
    # the model object
    stopifnot(model_fit$n_items == n_rows)
    
    df <- if(type == "CP"){
      .compute_cp_consensus_cluster(df)
    } else if(type == "MAP"){
      .compute_map_consensus_cluster(df)
    }
    
    
  } else if(parameter == "Rtilde"){
    # Filter out the pre-burnin iterations and get the right assessors
    df <- dplyr::filter(model_fit$augmented_data, iteration > burnin, assessor %in% assessors)
    
    # Find the problem dimensions
    n_rows <- nrow(dplyr::distinct(df, assessor, item))
    
    # Check that there are rows.
    stopifnot(n_rows > 0)
    
    # Check that the number of rows are consistent with the information in
    # the model object
    stopifnot(length(assessors) * model_fit$n_items == n_rows)
    
    # Treat assessors as clusters
    df <- dplyr::rename(df, cluster = "assessor")
    
    df <- if(type == "CP"){
      .compute_cp_consensus_cluster(df)
    } else if(type == "MAP"){
      .compute_map_consensus_cluster(df)
    }
    
    if("cluster" %in% names(df)){
      df <- dplyr::rename(df, assessor = "cluster")
    }
    
  }
  
  return(df)
  
}
.compute_cp_consensus_cluster <- function(df){
  
  # Convert items and cluster to character, since factor levels are not needed in this case
  df <- dplyr::mutate(df, across(c(item, cluster) , as.character))
  
  # Group by item, cluster, and value
  df <- dplyr::group_by(df, item, cluster, value)
  
  # Find the count of each unique combination (value, item, cluster)
  df <- dplyr::count(df)
  
  # Arrange according to value, per item and cluster
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, item, cluster)
  df <- dplyr::arrange(df, value, .by_group = TRUE)
  
  # Find the cumulative probability, by dividing by the total
  # count in (item, cluster) and the summing cumulatively
  df <- dplyr::mutate(df, cumprob = cumsum(n/sum(n)))
  
  # Find the CP consensus per cluster, using the find_cpc_cluster function
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, cluster)
  df <- find_cpc_cluster(df)
  
  # If there is only one cluster, we drop the cluster column
  if(length(unique(df$cluster)) == 1){
    df <- dplyr::select(df, -cluster)
  }
  
  return(df)
}

.compute_cp_consensus <- function(df){
  
  # Convert items to character, since factor levels are not needed in this case

  df <- dplyr::mutate(df, across(c(item) , as.character))
  
  # Group by item and value
  df <- dplyr::group_by(df, item, value)
  
  # Find the count of each unique combination (value, item)
  df <- dplyr::count(df)
  
  # Arrange according to value, per item
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, item)
  df <- dplyr::arrange(df, value, .by_group = TRUE)
  
  # Find the cumulative probability, by dividing by the total
  # count in (item, cluster) and the summing cumulatively
  df <- dplyr::mutate(df, cumprob = cumsum(n/sum(n)))
  
  # Find the CP consensus per cluster, using the find_cpc function
  df <- dplyr::ungroup(df)
  df <- find_cpc(df)

  
  return(df)
}

find_cpc_cluster <- function(group_df){
  # Declare the result dataframe before adding rows to it
  result <- dplyr::tibble(
    cluster = character(),
    ranking = numeric(),
    item = character(),
    cumprob = numeric()
  )
  n_items <- max(group_df$value)
  
  for(i in seq(from = 1, to = n_items, by = 1)){
    # Filter out the relevant rows
    tmp_df <- dplyr::filter(group_df, value == i)
    
    # Remove items in result
    tmp_df <- dplyr::anti_join(tmp_df, result, by = c("cluster", "item"))
    tmp_df <- tidyr::replace_na(tmp_df, list(cumprob =0))
    
    # Keep the max only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    tmp_df <- dplyr::filter(tmp_df, cumprob == max(cumprob) )
    
    # Add the ranking
    tmp_df <- dplyr::mutate(tmp_df, ranking = i)
    
    # Select the columns we want to keep, and put them in result
    result <- dplyr::bind_rows(result,
                               dplyr::select(tmp_df, cluster, ranking, item, cumprob))
    
    if(nrow(result) == n_items){
      break
    }
    
  }
  result <- dplyr::group_by(result, cluster)
  result <- dplyr::arrange(result, ranking, .by_group = TRUE)
  result <- dplyr::ungroup(result)
  
  return(result)
}


# Internal function for finding CP consensus.
find_cpc <- function(group_df){
  # Declare the result dataframe before adding rows to it
  result <- dplyr::tibble(
    ranking = numeric(),
    item = character(),
    cumprob = numeric()
  )
  n_items <- max(group_df$value)
  
  for(i in seq(from = 1, to = n_items, by = 1)){
    # Filter out the relevant rows
    tmp_df <- dplyr::filter(group_df, value == i)
    
    # Remove items already in result
    tmp_df <- dplyr::anti_join(tmp_df, result, by = "item")
    tmp_df <- tidyr::replace_na(tmp_df, list(cumprob =0))
    
    # Keep the max only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    tmp_df <- dplyr::filter(tmp_df, cumprob == max(cumprob) )
    
    # Add the ranking
    tmp_df <- dplyr::mutate(tmp_df, ranking = i)
    
    # Select the columns we want to keep, and put them in result
    result <- dplyr::bind_rows(result,
                               dplyr::select(tmp_df, ranking, item, cumprob))
    
    if(nrow(result) == n_items){
      break
    }
    
  }
  return(result)
}


.compute_map_consensus_cluster <- function(df){
  
  # Store the total number of iterations after burnin
  n_samples <- length(unique(df$iteration))
  
  # Spread to get items along columns
  df <- tidyr::pivot_wider(df, names_from = "item", values_from = "value")
  
  # Group by everything except iteration, and count the unique combinations
  df <- dplyr::group_by(df, across(!iteration))
  df <- dplyr::count(df)
  df <- dplyr::ungroup(df)
  # Keep only the maximum per cluster
  df <- dplyr::group_by(df, cluster)
  df <- dplyr::mutate(df, n_max = max(n))
  df <- dplyr::filter(df, n == n_max)
  df <- dplyr::ungroup(df)
  
  # Compute the probability
  df <- dplyr::mutate(df, prob = n / n_samples)
  df <- dplyr::select(df, -c(n_max, n))
  
  # Now collect one set of ranks per cluster
  df <- tidyr::pivot_longer(df, cols = -c("cluster", "prob"), names_to = "item", values_to = "ranking")
  
  # Sort according to cluster and ranking
  df <- dplyr::arrange(df, cluster, ranking)
  df <- dplyr::select(df, cluster, ranking, item, prob)
  
  if(length(unique(df$cluster)) == 1){
    df <- dplyr::select(df, -cluster)
  }
  
  
  
  return(df)
  
}


.compute_map_consensus <- function(df){
  
  # Store the total number of iterations after burnin
  n_samples <- length(unique(df$iteration))
  
  # Spread to get items along columns
  df <- tidyr::pivot_wider(df, names_from = "item", values_from = "value")
  
  # Group by everything except iteration, and count the unique combinations
  df <- dplyr::group_by(df, across(!iteration))
  df <- dplyr::count(df)
  df <- dplyr::ungroup(df)
  # Keep only the maximum per cluster
  df <- dplyr::mutate(df, n_max = max(n))
  df <- dplyr::filter(df, n == n_max)
  
  # Compute the probability
  df <- dplyr::mutate(df, prob = n / n_samples)
  df <- dplyr::select(df, -c(n_max, n))
  
  # Now collect one set of ranks per cluster
  df <- tidyr::pivot_longer(df, cols = -c("prob"), names_to = "item", values_to = "ranking")
  
  # Sort according to cluster and ranking
  df <- dplyr::arrange(df, ranking)
  df <- dplyr::select(df, ranking, item, prob)
  # df <- dplyr::relocate(df, ranking)
  # df <- dplyr::relocate(df, prob, .after = last_col())
  
  return(df)
  
}

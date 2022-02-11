assign_cluster <- function(model_fit, burnin = model_fit$burnin, soft = TRUE, expand = FALSE){
  
  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)
  
  df <- dplyr::filter(model_fit$cluster_assignment, iteration > burnin)
  
  # Compute the probability of each iteration
  df <- dplyr::group_by(df, assessor)
  df <- dplyr::mutate(df, probability = 1 / dplyr::n())
  df <- dplyr::ungroup(df)
  
  # Compute the probability of each cluster, per assessor
  df <- dplyr::group_by(df, assessor, value)
  df <- dplyr::summarise(df, probability = sum(probability), .groups = "drop")
  df <- dplyr::rename(df, cluster = value)
  
  if(expand){
    df <- tidyr::complete(
      dplyr::group_by(df, assessor),
      cluster = unique(df$cluster),
      fill = list(probability = 0)
    )
    df <- dplyr::ungroup(df)
  }
  
  # Compute the MAP estimate per assessor
  map <- dplyr::group_by(df, assessor)
  map <- dplyr::mutate(map, max_prob = max(probability))
  map <- dplyr::filter(map, probability == max_prob)
  
  # Deal with the possible case of ties
  map <- dplyr::filter(map, dplyr::row_number() == 1)
  map <- dplyr::ungroup(map)
  map <- dplyr::select(map, -c(probability, max_prob))
  map <- dplyr::rename(map, MAP_cluster = cluster)
  
  # Join map back onto df
  df <- dplyr::inner_join(df, map, by = "assessor")
  
  if(!soft){
    df <- dplyr::filter(df, cluster == MAP_cluster)
    df <- dplyr::select(df, -cluster)
  }
  
  return(df)
}
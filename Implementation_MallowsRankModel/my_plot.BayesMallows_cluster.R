source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_assign_cluster.R")

plot.myBayesMallows <- function(x, burnin = x$burnin, parameter = "alpha", items = NULL, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.
  
  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  if(x$nmc <= burnin) stop("nmc must be <= burnin")
  
  stopifnot(parameter %in% c("alpha", "rho", "cluster_probs", "cluster_assignment"))
  
  if(parameter == "alpha") {
    
    df = dplyr::filter(x$alpha, iteration > burnin)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ylab("Posterior density")
    
    if(x$n_clusters > 1){
      p <- p + ggplot2::facet_wrap(~ cluster, scales = "free_x")
    }
    
    return(p)
    
  } else if(parameter == "rho") {
    
    if(is.null(items) && x$n_items > 5){
      message("Items not provided by user. Picking 5 at random.")
      items <- sample.int(x$n_items, 5)
    } else if (is.null(items) && x$n_items > 0) {
      items <- 1:x$n_items
    }
    
    if(!is.character(items)){
      items <- x$items[items]
    } else{
      if(!all(items %in% x$items)){
        stop("Irrelevant items are given as input for plotting rho")
      }
    }
    
    df <- dplyr::filter(x$rho, iteration > burnin, item %in% items)
    
    # Compute the density, rather than the count, since the latter
    # depends on the number of Monte Carlo samples
    df <- dplyr::group_by(df, cluster,item, value)
    df <- dplyr::summarize(df, freq = dplyr::n(), .groups = "drop_last")
    df <- dplyr::mutate(df, prop = freq / sum(freq))
    
    # Finally create the plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = prop)) +
      ggplot2::geom_col() +
      ggplot2::scale_x_continuous(breaks = 1:x$n_items) +
      ggplot2::xlab("rank") +
      ggplot2::ylab("Posterior probability") 
    if(x$n_clusters ==1){
      p <- p + ggplot2::facet_wrap(~ item)
    } else{
      p <- p + ggplot2::facet_wrap(~ cluster + item)
    }
      
    return(p)
  } else if(parameter == "cluster_probs"){
    
    if(x$n_clusters == 1){
      stop("If n_clusters is equal to 1 then plotting cluster_probs is meaningless")
    }
    
    df <- dplyr::filter(x$cluster_probs, iteration > burnin)
    
    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(tau[c])) +
      ggplot2::ylab("Posterior density") +
      ggplot2::facet_wrap(~ cluster)
    
  } else if(parameter == "cluster_assignment"){
    
    if(x$n_clusters == 1){
      stop("If n_clusters is equal to 1 then plotting cluster_assignment is meaningless")
    }
    
    # First get one cluster per assessor, and sort these
    df <- assign_cluster(x, burnin = burnin, soft = FALSE, expand = FALSE)
    df <- dplyr::arrange(df, MAP_cluster)
    assessor_order <- dplyr::pull(df, assessor)
    
    # Next, rerun with 'soft = TRUE' option to get the probability of all clusters
    df <- assign_cluster(x, burnin = burnin, soft = TRUE, expand = TRUE)
    # Then order the assessors according to assessor_order
    df <- dplyr::mutate(df, assessor = factor(assessor, levels = assessor_order))
    
    # Now make a plot
    ggplot2::ggplot(df, ggplot2::aes(assessor, cluster)) +
      ggplot2::geom_tile(ggplot2::aes(fill = probability)) +
      ggplot2::scale_fill_continuous(low = "yellow", high = "red") +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      ) +
      ggplot2::xlab(paste0("Assessors (", min(assessor_order), " - ", max(assessor_order), ")"))
    
  }
}

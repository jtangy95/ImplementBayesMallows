predict_top_k <- function(model_fit, burnin = model_fit$burnin,
                          k = 3){
  
  validate_top_k(model_fit, burnin)
  .predict_top_k(model_fit, burnin, k)
}



.predict_top_k <- function(model_fit, burnin, k){
  
  rankings <- dplyr::filter(model_fit$augmented_data, iteration > burnin, value <= k)
  n_samples <- length(unique(rankings$iteration))
  rankings <- dplyr::mutate(rankings, across(c(item), as.character))
  rankings <- dplyr::group_by(rankings, assessor, item)
  rankings <- dplyr::summarize(rankings, prob = dplyr::n()/n_samples, .groups = "drop")
  
  # this find all missing possible combinations of assessors and items and fill it with zero
  rankings <- tidyr::complete(
    dplyr::group_by(rankings, assessor),
    item = model_fit$items,
    fill = list(prob = 0)
  )
  rankings <- dplyr::arrange(rankings, desc(prob), .by_group = TRUE)
  
  return(rankings)
}


validate_top_k <- function(model_fit, burnin){
  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)
  
  if(!exists("augmented_data", model_fit)){
    stop("model_fit must have element augmented_data. Please set save_aug = TRUE
         in compute_mallows in order to create a top-k plot.")
  }
}

# 
# plot_top_k <- function(model_fit, burnin = model_fit$burnin,
#                        k = 3,
#                        rel_widths = c(1, 10)){
#   
#   validate_top_k(model_fit, burnin)
#   
#   # Extract post burn-in rows with value <= k
#   rho <- dplyr::filter(model_fit$rho, iteration > burnin, value <= k)
#   n_samples <- length(unique(rho$iteration))
#   # Factors are not needed in this case
#   rho <- dplyr::mutate(rho, across(c(item), as.character))
#   rho <- dplyr::group_by(rho, item)
#   rho <- dplyr::summarize(rho, prob = dplyr::n()/n_samples, .groups = "drop")
#   
#   # Find the complete set of items
#   rho <- tidyr::complete(
#     dplyr::group_by(rho),
#     item = model_fit$items,
#     fill = list(prob = 0))
#   rho <- dplyr::ungroup(rho)
#   
#   # Sort the items according to probability
#   # This factor level works for plotting the heat-map below
#   item_ordering <- rev(compute_consensus(model_fit, type = "CP", burnin = burnin)$item)
#   rho <- dplyr::mutate(rho, item = factor(.data$item, levels = item_ordering))
#   
#   # # Trick to make the plot look nicer
#   # if(model_fit$n_clusters == 1){
#   #   rho <- dplyr::mutate(rho, cluster = "")
#   # }
#   
#   rankings <- .predict_top_k(model_fit, burnin = burnin, k = k)
#   
#   # Sorting the items according to their probability in rho
#   # This factor level works for plotting the heatmap below
#   rankings <- dplyr::mutate(rankings, item = factor(item, levels = item_ordering))
#   
#   assessor_plot <- ggplot2::ggplot(rankings, ggplot2::aes(assessor,item)) +
#     ggplot2::geom_tile(ggplot2::aes(fill = prob), colour = "skyblue") +
#     ggplot2::scale_fill_continuous(low = "yellow", high = "red") +
#     ggplot2::xlab("Assessor") +
#     ggplot2::scale_x_continuous(breaks = 1:model_fit$n_assessors) +
#     ggplot2::theme(
#       legend.title = ggplot2::element_blank(),
#       axis.title.y = ggplot2::element_blank(),
#       axis.text.y = ggplot2::element_blank(),
#       axis.ticks.y = ggplot2::element_blank()
#     )
#   
#   rho_plot <- ggplot2::ggplot(rho, ggplot2::aes("", item)) +
#     ggplot2::geom_tile(ggplot2::aes(fill = prob), colour = "skyblue") +
#     ggplot2::scale_fill_continuous(low = "yellow", high = "red") +
#     ggplot2::ylab("Item") +
#     ggplot2::xlab(expression(rho)) +
#     ggplot2::theme(legend.position = "none")
#   
#   # if(model_fit$n_clusters > 1){
#   #   rho_plot <- rho_plot + ggplot2::facet_wrap(~ .data$cluster)
#   # }
#   
#   # rel_widths stands for relative widths of plots
#   cowplot::plot_grid(rho_plot, assessor_plot, rel_widths = rel_widths)
# }



plot_top_k <- function(model_fit, burnin = model_fit$burnin,
                       k = 3, assessors = NULL,
                       rel_widths = c(1, 9)){
  
  validate_top_k(model_fit, burnin)
  
  # Extract post burn-in rows with value <= k
  rho <- dplyr::filter(model_fit$rho, iteration > burnin, value <= k)
  n_samples <- length(unique(rho$iteration))
  # Factors are not needed in this case
  rho <- dplyr::mutate(rho, across(c(item), as.character))
  rho <- dplyr::group_by(rho, item)
  rho <- dplyr::summarize(rho, prob = dplyr::n()/n_samples, .groups = "drop")
  
  # Find the complete set of items
  rho <- tidyr::complete(
    dplyr::group_by(rho),
    item = model_fit$items,
    fill = list(prob = 0))
  rho <- dplyr::ungroup(rho)
  
  # Sort the items according to probability
  # This factor level works for plotting the heat-map below
  item_ordering <- rev(compute_consensus(model_fit, type = "CP", burnin = burnin)$item)
  rho <- dplyr::mutate(rho, item = factor(.data$item, levels = item_ordering))
  
  # # Trick to make the plot look nicer
  # if(model_fit$n_clusters == 1){
  #   rho <- dplyr::mutate(rho, cluster = "")
  # }
  
  rankings <- .predict_top_k(model_fit, burnin = burnin, k = k)
  
  # Sorting the items according to their probability in rho
  # This factor level works for plotting the heatmap below
  rankings <- dplyr::mutate(rankings, item = factor(item, levels = item_ordering))
  if(!is.null(assessors)){
    rankings <- dplyr::filter(rankings, assessor %in% assessors)
  }
  
  assessor_plot <- ggplot2::ggplot(rankings, ggplot2::aes(assessor,item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = prob), colour = "skyblue") +
    ggplot2::scale_fill_continuous(low = "yellow", high = "red") +
    ggplot2::xlab("Assessor") +
    ggplot2::scale_x_continuous(breaks = 1:model_fit$n_assessors) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  
  rho_plot <- ggplot2::ggplot(rho, ggplot2::aes("", item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = prob), colour = "skyblue") +
    ggplot2::scale_fill_continuous(low = "yellow", high = "red") +
    ggplot2::ylab("Item") +
    ggplot2::xlab(expression(rho)) +
    ggplot2::theme(legend.position = "none")
  
  if(model_fit$n_clusters > 1){
    rho_plot <- rho_plot + ggplot2::facet_wrap(~ cluster)
  }
  
  # rel_widths stands for relative widths of plots
  cowplot::plot_grid(rho_plot, assessor_plot, rel_widths = rel_widths)
}

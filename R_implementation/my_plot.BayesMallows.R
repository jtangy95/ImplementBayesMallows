plot.myBayesMallows <- function(x, burnin = x$burnin, parameter = "alpha", items = NULL, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.
  
  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  if(x$nmc <= burnin) stop("nmc must be <= burnin")
  
  stopifnot(parameter %in% c("alpha", "rho"))
  
  if(parameter == "alpha") {
    
    df = dplyr::filter(x$alpha, iteration > burnin)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ylab("Posterior density")
    
    
    return(p)
    
  } else if(parameter == "rho") {
    
    if(is.null(items) && x$n_items > 5){
      message("Items not provided by user. Picking 5 at random.")
      items <- sample.int(x$n_items, 5)
    } else if (is.null(items) && x$n_items > 0) {
      items <- seq.int(from = 1, to = x$n_items)
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
    df <- dplyr::group_by(df, item, value)
    suppressMessages(df <- dplyr::summarize(df, freq = dplyr::n()))
    df <- dplyr::mutate(df, prop = freq / sum(freq))
    
    # Finally create the plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = prop)) +
      ggplot2::geom_col() +
      ggplot2::scale_x_continuous(breaks = 1:x$n_items) +
      ggplot2::xlab("rank") +
      ggplot2::ylab("Posterior probability") +
      ggplot2::facet_wrap(~ item)
    
    return(p)
  }  
}

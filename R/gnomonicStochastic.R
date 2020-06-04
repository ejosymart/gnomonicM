#' Gnomonic stochastic
#'
#' Estimate natural mortality based on gnomonic interval approach with different distribution in fecundity.
#'
#' @param nInterval a numeric value that represents the number of gnomonic intervals.
#' @param eggDuration a numeric value with the eggstage duration- first gnomonic interval (days).
#' @param addInfo a numeric vector with additional information related to the observed duration of the others gnomonic intervals. Write \code{addInfo = NULL} if you do not provide additional information.
#' @param longevity a numeric value indicating the lifespan of the species (days).
#' @param fecundity a numeric value indicating the mean or the mode of the fecundity as the number of eggs produced for a female if a normal or triangular distribution is assumed, respectively.
#' @param sd_fecundity a numeric value indicating the standard deviation of fecundity if a normal distribution is assumed.
#' @param min_fecundity a numeric value indicating the minimum range of fecundity if a uniform or triangle distribution is assumed.
#' @param max_fecundity a numeric value indicating the maximum range of fecundity if a uniform or triangle distribution is assumed.
#' @param distr a character string indicating the distribution to be applied: \code{"uniform"}, \code{"triangle"} or \code{"normal"}
#' @param a_init a numeric value indicating the initial parameter related to the proportionality constant optimized by iterative solution via univariate (1-dim.) minimization.
#' @param niter a single numeric value representing the number of iterations.
#' @param seed a single value, interpreted as an integer.
#' @return A list of class 'gnomosBoot'.
#'
#' \code{a} the proportionality constant.
#'
#' \code{G} the 'n' iter values of constant proportion of the overall natural death rate.
#'
#' \code{mean_G} the mean of constant proportion of the overall natural death rate,
#'
#' \code{M} a dataframe with the M values for each gnomonic intervals for each 'n' iteration.
#'
#' \code{fecundity} the 'n' iter values of fecundity based on the distribution assumed.
#'
#' \code{results} a dataframe with the duration ("interval_duration_day"), mean, confidence interval and standard deviation of natural mortality ("M_lower", "M", "M_upper", "M_sd") for each gnomonic interval.
#' @details Estimate natural mortality (M) based on gnomonic interval approach .
#'
#' The argument \code{nInterval} is NULL by default. If you have -at least- one observed value for the duration of the other gnomonic intervals
#' you should provide this as a vector which length must be nInterval - 1, for example \code{addInfo = c(3, NA, NA, NA, NA, NA)}) for a \code{nInterval = 7}.
#'
#' The argument \code{fecundity} requires a character string indicating the name of the distribution of fecundity values to be used in the
#' analysis (i.e. \code{fecundity = "uniform"}).
#'
#' The argument \code{niter} requires a number which is related with the number of observations. If length(n) > 1, the length is taken to be the number required.
#' can be calculated from each bootstrap sample (median and confidence intervals).
#' @examples
#' #The values are based on Caddy (1996).
#' modelBoot <- gnomonicStochastic(nInterval = 7, eggDuration = 2, addInfo = NULL, longevity = 365,
#' distr = "uniform", min_fecundity = 100000, max_fecundity = 300000, niter = 999, a_init = 2)
#'
#' # 'niter' parameters:
#' modelBoot$a
#' modelBoot$G
#' modelBoot$mean_G
#' modelBoot$M
#' modelBoot$fecundity
#' modelBoot$results
#' @export
gnomonicStochastic <- function(nInterval, eggDuration, addInfo = NULL, longevity,
                               fecundity = NULL, sd_fecundity = NULL, min_fecundity = NULL, max_fecundity = NULL,
                               distr = "uniform", a_init,  niter = 999, seed = 7388){

  set.seed(seed)

  if(any(distr == "uniform") & any(is.null(min_fecundity),
                                   is.null(max_fecundity)))
    stop("HEY! 'uniform' distribution requires 'fecundity', 'min_fecundity' and 'max_fecundity' values.")

  if(any(distr == "uniform") & any(min_fecundity >= max_fecundity))
    stop("HEY! 'max_fecundity' must be greater than 'min_fecundity' value. Review those inputs!")


  if(any(distr == "triangle") & any(is.null(min_fecundity),
                                    is.null(max_fecundity),
                                    is.null(fecundity)))
    stop("HEY! 'triangle' distribution requires 'fecundity', 'min_fecundity' and 'max_fecundity' values.")

  if(any(distr == "triangle") & any(min_fecundity >= max_fecundity))
    stop("HEY! 'max_fecundity' must be greater than 'min_fecundity' value. Review those inputs!")

  if(any(distr == "triangle") & any(min_fecundity >= fecundity))
    stop("HEY! 'fecundity' must be greater than 'min_fecundity' value. Review those inputs!")

  if(any(distr == "triangle") & any(fecundity >= max_fecundity))
    stop("HEY! 'max_fecundity' must be greater than 'fecundity' value. Review those inputs!")


  if(any(distr == "normal") & any(is.null(sd_fecundity),
                                  is.null(fecundity)))
    stop("HEY! 'normal' distribution requires 'fecundity', 'sd_fecundity' value. Review this value!")


  if(is.null(addInfo)){

    cat("--------------------------------------------------------", "\n\n")
    cat('You are only considering the egg stage duration =', eggDuration, "\n\n")
    cat("--------------------------------------------------------", "\n\n")

    output <- .noAddInfo(nInterval = nInterval,
                         eggDuration = eggDuration,
                         longevity = longevity,
                         a_init = a_init)

  }else{

    if(length(addInfo) != nInterval-1) stop('The length of addInfo vector must be equal to nInterval-1')

    output <- .AddInfo(nInterval = nInterval,
                       eggDuration = eggDuration,
                       longevity = longevity,
                       a_init = a_init,
                       addInfo = addInfo)

  }


  if(any(distr == "uniform")){
    print("You are using a 'uniform distribution' for fecundity")
    fec <- runif(n = niter, min = min_fecundity, max = max_fecundity)
  }

  if(any(distr == "triangle")){
    print("You are using a 'triangular distribution' for fecundity")
    fec <- rtriangle(n = niter, a = min_fecundity, b = max_fecundity, c = fecundity)
  }

  if(any(distr == "normal")){
    print("You are using a 'normal distribution' for fecundity")
    fec <- rnorm(n = niter, mean = fecundity, sd = sd_fecundity)
  }

  G   <- -log((2/fec)^(1/nInterval))

  M   <- vector()
  for(j in seq_len(niter)){
    m <- G[j]/output$delta
    M <- rbind(M, m)
    M <- data.frame(M)
  }
  colnames(M) <- paste0("Gnomonic_", seq_len(nInterval))
  rownames(M) <- paste0("Iteration_", seq_len(niter))

  M_mean <- apply(M, 2, mean)
  M_IC   <- apply(M, 2, quantile, probs = c(0.025, 0.975))
  M_sd   <- apply(M, 2, sd)


  tab <- data.frame(Gnonomic_interval = seq_len(nInterval),
                    duration          = round(output$delta, 3),
                    total_duration    = round(cumsum(output$delta), 0),
                    M_lower           = as.numeric(M_IC[1,]),
                    M                 = as.numeric(M_mean),
                    M_upper           = as.numeric(M_IC[2,]),
                    M_sd              = as.numeric(M_sd))

  data <- list(a = output$a, G = G, mean_G = mean(G), M = M, fecundity = fec, results = tab)

  class(data) <- c("gnomosBoot", class(data))

  return(data)
}



#' Print method for gnomosBoot class
#'
#' @param x an object class 'gnomosBoot'.
#' @param \dots Additional arguments to the print method.
#' @return The values of the proportionality constant (a), constant proportion of the overall natural death rate (G) and
#' a data.frame with the duration and natural mortality for each gnomonic interval.
#' @examples
#' #The values are based on Caddy (1996).
#' modelBoot <- gnomonicStochastic(nInterval = 7, eggDuration = 2, addInfo = NULL, longevity = 365,
#' distr = "uniform", min_fecundity = 100000, max_fecundity = 300000, niter = 50, a_init = 2)
#'
#' print(modelBoot)
#' @export
#' @method print gnomosBoot
print.gnomosBoot <- function(x, ...){

  if(!inherits(x, "gnomosBoot"))
    stop("Use only with 'gnomosBoot' objects.")

  data <- x
  cat('Proportionality constant (alpha) =', data$a, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('Mean value of constant proportion of the overall natural death rate (G) =', data$mean_G, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('M values for each gnomonic interval for each iteration:', "\n\n")
  print(data$M)
  cat("--------------------------------------------------------", "\n\n")
  cat('Main results of gnomonic method:', "\n")
  print(data$results)
  return(invisible())
}



#' Plot method for gnomosBoot class
#'
#' @param x an object class 'gnomosBoot'.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col color for the boxplot of M value for each gnomonic intervals.
#' @param boxwex a scale factor to be applied to all boxes in order to make the boxes narrower or wider.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' modelBoot <- gnomonicStochastic(nInterval = 7, eggDuration = 2, addInfo = NULL, longevity = 365,
#' distr = "uniform", min_fecundity = 100000, max_fecundity = 300000, niter = 1000, a_init = 2)
#'
#' plot(modelBoot)
#' @export
#' @method plot gnomosBoot
plot.gnomosBoot <- function(x, xlab = "Gnomonic intervals", ylab = expression(paste(bar(M), " (day"^"-1",")")),
                            col = "lightgrey", boxwex = 0.25, ...){

  if(!inherits(x, "gnomosBoot"))
    stop("Use only with 'gnomosBoot' objects.")

  data <- x

  par(mar = c(4, 6, 4, 1))
  boxplot(data$M, boxwex = boxwex, xlab = xlab,
          ylab = ylab, col = col, axes = FALSE, ...)
  axis(1, seq(from = 1, to = nrow(data$results), by = 1))
  axis(2, las = 2)
  box()

  return(invisible(NULL))
}

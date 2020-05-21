# gnomonicM package: Estimate Natural Mortality for different life --------

#' @importFrom minqa newuoa
#' @importFrom triangle rtriangle
#' @importFrom grDevices colors
#' @importFrom graphics axis box boxplot plot points hist par abline
#' @importFrom stats quantile rnorm runif sd
#' @importFrom utils installed.packages
#' @title Estimate Natural Mortality for different life stages.
#'
#' @description Estimate natural mortality (M) throughout the life history for organisms, mainly fish and invertebrates, based on gnomonic interval approach. It includes estimation of duration of each gnomonic interval (life stage) and the constant probability of death.
#' @name gnomonicM-package
#' @aliases gnomonicM-package gnomonicM
#' @docType package
#' @author Josymar Torrejon-Magallanes <ejosymart@@gmail.com>
#' @details Package: gnomonicM
#' @details Type: Package
#' @details The natural mortality (M) estimation throughout different life stages is based on the gnomonic approach (Caddy, 1991, 1996), but we include new features.
#'
#' In the gnomonic model, the estimation of \eqn{M_{i}} for each gnomonic interval \eqn{\Delta_{i}} requires -at least- information about:
#' (i) number of development stages throughout the life cycle i \eqn{in} \eqn{1, …n}.
#' (ii) the duration of the first life stage corresponding to first gnomonic interval (\eqn{\Delta_{1}}, egg stage),
#' (iii) the mean lifetime fecundity \eqn{MLF}, and
#' (iv) the lifespan of the species. As additional information, the duration of the other developments stages or gnomonic intervals (larval, juvenile, adults) could be provided.
#'
#' According to Caddy (1996) and Martinez-Aguilar (2005), whatever the function of \eqn{M} with age, knowing the number of individuals (\eqn{N})
#' at the beginning of the year and divide a year into specific number \eqn{i = 1, 2, 3, ..., n} of smaller intervals:
#'
#' \deqn{N_{i} = MLF  for i = 1}
#' \deqn{N_{i} = N_{i-1} e^-(M_{i} * \Delta_{i}) for i > 1}
#'
#' where:
#'
#' \eqn{M_{i}} is the average value for natural mortality rate, that integrates the declining death rate through the short time interval duration \eqn{\Delta_{i}}. The \eqn{N_{i}}
#' is the survivors from previous interval, only for the first interval (\eqn{\Delta_{1}}) is assumed that the numbers of haching eggs (initial population) is equivalent to
#' the mean lifetime fecundity (\eqn{MLF}).
#'
#' The duration of first gnomonic interval \eqn{\Delta_{1}} is equal to the time elapsed after the moment of hatching \eqn{t_{1}}.
#' The duration of the subsequent gnomonic intervals (\eqn{i > 1}) are estimated following:
#'
#' \deqn{\Delta_{i} = \Delta_{1}  \alpha  (\alpha + 1)^{i-2}}
#'
#' where,
#' \eqn{\Delta_{i}}: Duration of the gnomonic interval when \eqn{i > 1}.
#' \eqn{\Delta_{1}}: Duration of the first gnomonic interval \eqn{t_{1}}.
#' \eqn{\alpha}: Proportionality constant.
#' \eqn{i}: \eqn{i^{th}} gnomonic interval.
#'
#' The \eqn{M_{i}} is estimated as follows:
#'
#' \deqn{M_{i} = \frac{G}{\Delta_{i, i-1}}
#'
#' where \eqn{G} is the constant proportion of the overall natural death rate. The \eqn{G} value is calculated so that the number  of individuals surviving to the last
#' gnomonic time-interval is \eqn{N_{n} = 2} following the assumptiion of stable population replacement (Caddy, 1996; Martinez-Aguilar, 2005).
#' The new equation for \eqn{G} is expressed:
#'
#' \deqn{G} = -ln((\frac{2}{MLF})^{\frac{1}{n}})
#'
#' The final solution is to estimate the proportionality constant (\alpha) parameter by iterative solution via univariate (1-dim.) minimization.
#'
#' @references Caddy JF (1991). Death rates and time intervals: is there an alternative to the constant natural mortality axiom? Rev Fish Biol Fish 1:109–138. doi: 10.1007/BF00157581.
#' @references Caddy JF (1996). Modelling natural mortality with age in short-lived invertebrate populations: definition of a strategy of gnomonic time division. Aquat Living Resour 9:197–207. doi: 10.1051/alr:1996023.
#' @references Martínez-Aguilar S, Arreguín-Sánchez F, Morales-Bojórquez E (2005). Natural mortality and life history stage duration of Pacific sardine (Sardinops caeruleus) based on gnomonic time divisions. Fish Res 71:103–114. doi: 10.1016/j.fishres.2004.04.008.
#' @references Torrejon-Magallanes J, Arreguín-Sánchez F, Morales-Bojórquez E (in prep). Natural mortality estimation through the life history of Pacific chub mackerel (Scomber japonicus): a gnomonic approach
#' @concept gnomonic
#' @concept natural mortality
#' @concept fecundity
#' @concept lifespan
#' @examples
#' #See examples for functions gnomonic() and gnomonicStochastic().

NULL
#' Gnomonic deterministic
#'
#' Estimate natural mortality based on gnomonic interval approach.
#' @param nInterval a numeric value that represents the number of gnomonic intervals.
#' @param eggDuration a numeric value with the eggstage duration- first gnomonic interval (days).
#' @param addInfo a numeric vector with additional information related to the observed duration of the others gnomonic intervals. Write \code{addInfo = NULL} if you do not provide additional information.
#' @param longevity a numeric value indicating the lifespan of the species (days).
#' @param fecundity a numeric value indicating the mean lifetime fecundity as the number of eggs produced for a female.
#' @param a_init a numeric value indicating the initial parameter related to the proportionality constant which will be optimized.
#' @return A list of class 'gnomos'.
#'
#' \code{a} the proportionality constant.
#'
#' \code{G} the constant proportion of the overall natural death rate.
#'
#' \code{results} a dataframe with the duration ("interval_duration_day") and natural mortality ("M_day" and "M_year") for each gnomonic interval.
#' @details Estimate natural mortality (M) based on gnomonic interval approach.
#'
#' The argument \code{nInterval} is NULL by default. If you have -at least- one observed value for the duration of the other gnomonic intervals
#' you should provide this as a vector which length must be nInterval - 1, for example \code{addInfo = c(3, NA, NA, NA, NA, NA)}) for a \code{nInterval = 7}.
#'
#' @references Caddy JF (1996). Modelling natural mortality with age in short-lived invertebrate populations: definition of a strategy of gnomonic time division. Aquat Living Resour 9:197–207. doi: 10.1051/alr:1996023.
#' @exportClass gnomos
#' @examples
#' #The values are based on Caddy (1996).
#' model <- gnomonic(nInterval = 7, eggDuration = 2, addInfo = NULL,
#' longevity = 365, fecundity = 200000, a_init = 2)
#'
#' model
#' model$a
#' model$G
#' model$results
#'
#' #Additional information for the duration of the second gnomonic intervals.
#' model <- gnomonic(nInterval = 7, eggDuration = 2, addInfo = c(3, NA, NA, NA, NA, NA),
#' longevity = 365, fecundity = 200000, a_init = 2)
#'
#' model
#' model$a
#' model$G
#' model$results
#' @export
gnomonic <- function(nInterval, eggDuration, addInfo = NULL,
                     longevity, fecundity, a_init) {

  if(!is.null(addInfo)){

    if(length(addInfo) != nInterval-1) stop('The length of addInfo vector must be equal to nInterval-1')

    minimize <- function(param, ...){
      d    <- c(eggDuration, addInfo)
      for(i in 2:nInterval){
        if(!is.na(d[i])) next
        d[i] <- d[1]*param*(param+1)^(seq_len(nInterval)[i]-2)
      }
      min <- abs(longevity - sum(d))
      return(min)
    }

    a   <- newuoa(par = a_init, fn = minimize)$par

    delta    <- c(eggDuration, addInfo)
    for(i in 2:nInterval){
      if(!is.na(delta[i])) next
      delta[i] <- delta[1]*a*(a+1)^(seq_len(nInterval)[i]-2)
    }

  }


  if(is.null(addInfo)){

    minimize <- function(param, ...){
      d    <- numeric(nInterval)
      d[1] <- eggDuration
      for(i in 2:nInterval){
        d[i] <- d[1]*param*(param+1)^(seq_len(nInterval)[i]-2)
      }
      min <- abs(longevity - sum(d))
      return(min)
    }

    a   <- newuoa(par = a_init, fn = minimize)$par

    delta    <- numeric(nInterval)
    delta[1] <- eggDuration
    for(i in 2:nInterval){
      delta[i] <- delta[1]*a*(a+1)^(seq_len(nInterval)[i]-2)
    }

  }


  G <- -log((2/fecundity)^(1/nInterval))
  M <- G/delta

  abundance    <- numeric(nInterval)
  for(i in seq_len(nInterval)){
    abundance[i]  <- fecundity*(exp(-G))^(i)
  }

  tab <- data.frame(Gnonomic_interval     = seq_len(nInterval),
                    interval_duration_day = round(delta, 3),
                    total_duration        = round(cumsum(delta), 0),
                    M_day                 = round(M, 3),
                    M_year                = round(M*365, 3),
                    No_Surv               = round(abundance, 0))

  data <- list(a = a, G = G, results = tab)

  class(data) <- c("gnomos", class(data))

  return(data)
}


#' Print method for gnomos class
#'
#' @param x an object class 'gnomos'.
#' @param \dots Additional arguments to the print method.
#' @return The values of the proportionality constant (a), constant proportion of the overall natural death rate (G) and
#' a data.frame with the duration and natural mortality for each gnomonic interval.
#'#' @examples
#' model <- gnomonic(nInterval = 7, eggDuration = 2, addInfo = NULL,
#' longevity = 365, fecundity = 200000, a_init = 2)
#'
#' print(model)
#' @export
#' @method print gnomos
print.gnomos <- function(x, ...){
  if (!inherits(x, "gnomos"))
    stop("Use only with 'gnomos' objects")

  data <- x
  cat('The value of proportionality constant (alpha) =', data$a, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('The value of constant proportion of the overall natural death rate (G) =', data$G, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('Main results of gnomonic methods:', "\n\n")
  print(data$results)
  return(invisible())
}


#' Plot method for gnomos class
#'
#' @param x an object class 'gnomos'.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param bg a background color for the points.
#' @param pch the character indicating the type of plotting.
#' @param cex character expansion in the regression.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' model <- gnomonic(nInterval = 7, eggDuration = 2, addInfo = NULL,
#' longevity = 365, fecundity = 200000, a_init = 2)
#'
#' plot(model)
#' @export
#' @method plot gnomos
plot.gnomos <- function(x, xlab = "Gnomonic intervals", ylab = expression(paste("M (day"^"-1", ")")),
                        bg = "lightgrey", cex = 1.75, pch = 21, ...){

  if (!inherits(x, "gnomos"))
    stop("Use only with 'gnomos' objects")

  data <- x

  par(mar = c(4, 6, 4, 1))
  plot(data$results$M_day, type = "b", cex = cex, pch = pch, lwd = 3, ylim = c(0, 1.1*max(data$results$M_day)),
       xlab = xlab , ylab = ylab, bg = bg , axes = FALSE, ...)
  axis(1, seq(from = 1, to = nrow(data$results), by = 1))
  axis(2, las = 2)
  box()

  return(invisible(NULL))
}





#' Gnomonic stochastic
#'
#' Estimate natural mortality based on gnomonic interval approach with different distribution in fecundity.
#'
#' @param nInterval a numeric value that represents the number of gnomonic intervals.
#' @param eggDuration a numeric value with the eggstage duration- first gnomonic interval (days).
#' @param addInfo a numeric vector with additional information related to the observed duration of the others gnomonic intervals. Write \code{addInfo = NULL} if you do not provide additional information.
#' @param longevity a numeric value indicating the lifespan of the species (days).
#' @param fecundity a numeric value indicating the mean lifetime fecundity as the number of eggs produced for a female if a normal or triangle distribution is assumed.
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
#' @exportClass gnomosBoot
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


  if(!is.null(addInfo)){

    if(length(addInfo) != nInterval-1) stop('The length of addInfo vector must be equal to nInterval-1')

    minimize <- function(param, ...){
      d    <- c(eggDuration, addInfo)
      for(i in 2:nInterval){
        if(!is.na(d[i])) next
        d[i] <- d[1]*param*(param+1)^(seq_len(nInterval)[i]-2)
      }
      min <- abs(longevity - sum(d))
      return(min)
    }

    a   <- newuoa(par = a_init, fn = minimize)$par

    delta    <- c(eggDuration, addInfo)
    for(i in 2:nInterval){
      if(!is.na(delta[i])) next
      delta[i] <- delta[1]*a*(a+1)^(seq_len(nInterval)[i]-2)
    }

  }


  if(is.null(addInfo)){

    minimize <- function(param, ...){
      d    <- numeric(nInterval)
      d[1] <- eggDuration
      for(i in 2:nInterval){
        d[i] <- d[1]*param*(param+1)^(seq_len(nInterval)[i]-2)
      }
      min <- abs(longevity - sum(d))
      return(min)
    }

    a   <- newuoa(par = a_init, fn = minimize)$par

    delta    <- numeric(nInterval)
    delta[1] <- eggDuration
    for(i in 2:nInterval){
      delta[i] <- delta[1]*a*(a+1)^(seq_len(nInterval)[i]-2)
    }
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
    m <- G[j]/delta
    M <- rbind(M, m)
    M <- data.frame(M)
  }
  colnames(M) <- paste0("Gnomonic_", seq_len(nInterval))
  rownames(M) <- paste0("Iteration_", seq_len(niter))

  M_mean <- apply(M, 2, mean)
  M_IC   <- apply(M, 2, quantile, probs = c(0.025, 0.975))
  M_sd   <- apply(M, 2, sd)


  tab <- data.frame(Gnonomic_interval = seq_len(nInterval),
                    duration          = round(delta, 3),
                    total_duration    = round(cumsum(delta), 0),
                    M_lower           = as.numeric(M_IC[1,]),
                    M                 = as.numeric(M_mean),
                    M_upper           = as.numeric(M_IC[2,]),
                    M_sd              = as.numeric(M_sd))

  data <- list(a = a, G = G, mean_G = mean(G), M = M, fecundity = fec, results = tab)

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
  if (!inherits(x, "gnomosBoot"))
    stop("Use only with 'gnomosBoot' objects")

  data <- x
  cat('The value of proportionality constant (alpha) =', data$a, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('The mean value of constant proportion of the overall natural death rate (G) =', data$mean_G, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('The value of M for each gnomonic interval for each iteration:', "\n\n")
  print(data$M)
  cat("--------------------------------------------------------", "\n\n")
  cat('Main results of gnomonic methods:', "\n")
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

  if (!inherits(x, "gnomosBoot"))
    stop("Use only with 'gnomosBoot' objects")

  data <- x

  par(mar = c(4, 6, 4, 1))
  boxplot(data$M, boxwex = boxwex, xlab = xlab,
          ylab = ylab, col = col, axes = FALSE, ...)
  axis(1, seq(from = 1, to = nrow(data$results), by = 1))
  axis(2, las = 2)
  box()

  return(invisible(NULL))
}

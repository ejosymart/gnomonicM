# gnomonicM package: Estimate Natural Mortality for Different Life Stages --------

#' @importFrom minqa newuoa
#' @importFrom triangle rtriangle
#' @importFrom grDevices colors
#' @importFrom kableExtra kbl kable_styling
#' @importFrom graphics axis box boxplot points hist par abline
#' @importFrom stats quantile rnorm runif sd
#' @importFrom utils installed.packages
#' @title Estimate Natural Mortality for Different Life Stages.
#'
#' @description Estimate natural mortality (M) throughout the life history for organisms, mainly fish and invertebrates, based on gnomonic interval approach. It includes estimation of duration of each gnomonic interval (life stage) and the constant probability of death (G).
#' @name gnomonicM-package
#' @aliases gnomonicM-package gnomonicM
#' @docType package
#' @author Josymar Torrejon-Magallanes <ejosymart@@gmail.com>
#' @details Package: gnomonicM
#' @details Type: Package
#' @details The natural mortality (M) estimation throughout different life stages is based on the gnomonic approach (Caddy, 1991, 1996), including new features in this package-version.
#'
#' In the gnomonic model, the estimation of \eqn{M_{i}} for each gnomonic interval \eqn{\Delta_{i}} requires -at least- information about:
#' (i) the number of development stages throughout the life cycle \eqn{i} \eqn{in} \eqn{1,2,3, …n}.
#' (ii) the duration of the first life stage corresponding to first gnomonic interval (\eqn{\Delta_{1}}, egg stage),
#' (iii) the mean lifetime fecundity \eqn{MLF}, and
#' (iv) the longevity of the species. As additional information, the duration of the other developments stages or gnomonic intervals (larval, juvenile, adults) could be provided.
#'
#' According to Caddy (1996) and Martinez-Aguilar (2005), the gnomonic method is supported by a negative exponential function, where the independent variable is
#' \eqn{\Delta_{i}} representing the number of gnomonic intervals from \eqn{i} \eqn{in} \eqn{1,2,3, …n}, the equation is expressed as follows:
#'
#' \deqn{N_{i} = MLF \cdot e^{-(M_{i} \cdot \Delta_{i})};      for i = 1}
#'
#' \deqn{N_{i} = N_{i-1} \cdot e^{-(M_{i} \cdot \Delta_{i})};       for i > 1}
#'
#' where:
#'
#' \eqn{M_{i}} is the average value for natural mortality rate, that integrates the declining death rate through the short time interval duration \eqn{\Delta_{i}}. The \eqn{N_{i}}
#' is the survivors from previous interval, only for the first interval (\eqn{\Delta_{1}}) is assumed that the numbers of hatching eggs (initial population) is equivalent to
#' the mean lifetime fecundity (\eqn{MLF}).
#'
#' The duration of first gnomonic interval \eqn{\Delta_{1}} is equal to the time elapsed after the moment of hatching \eqn{t_{1}}.
#' The duration of the subsequent gnomonic intervals (\eqn{i > 1}) are estimated following:
#'
#' \deqn{\Delta_{i} = \Delta_{1} \cdot \alpha  (\alpha + 1)^{i-2}}
#'
#' where,
#'
#' \eqn{\Delta_{i}}: Duration of the gnomonic interval when \eqn{i > 1}.
#'
#' \eqn{\Delta_{1}}: Duration of the first gnomonic interval \eqn{t_{1}}.
#'
#' \eqn{\alpha}: Proportionality constant.
#'
#' \eqn{i}: \eqn{i^{th}} gnomonic interval.
#'
#' The \eqn{M_{i}} is estimated as follows:
#'
#' \deqn{M_{i} = \frac{G}{\Delta_{i, i-1}}}
#'
#' where \eqn{G} is the constant proportion of the overall natural death rate. The \eqn{G} value is calculated so that the number  of individuals surviving to the last
#' gnomonic time-interval is \eqn{N_{n} = 2} following the assumption of stable population replacement (Caddy, 1996; Martinez-Aguilar, 2005).
#' The new equation for \eqn{G} is expressed:
#'
#' \deqn{G = -ln((\frac{2}{MLF})^{\frac{1}{n}})}
#'
#' The final solution is to estimate the proportionality constant (\eqn{\alpha}) parameter by iterative solution via univariate (1-dim.) minimization.
#'
#' @references Caddy JF (1991). Death rates and time intervals: is there an alternative to the constant natural mortality axiom? Rev Fish Biol Fish 1:109–138. doi:10.1007/BF00157581.
#' @references Caddy JF (1996). Modelling natural mortality with age in short-lived invertebrate populations: definition of a strategy of gnomonic time division. Aquat Living Resour 9:197–207. doi:10.1051/alr:1996023.
#' @references Martínez-Aguilar S, Arreguín-Sánchez F, Morales-Bojórquez E (2005). Natural mortality and life history stage duration of Pacific sardine (Sardinops caeruleus) based on gnomonic time divisions. Fish Res 71:103–114. doi:10.1016/j.fishres.2004.04.008.
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
#' @param eggDuration a numeric value with the egg stage (first gnomonic interval) duration in days.
#' @param addInfo a numeric vector with additional information related to the observed duration of the other gnomonic intervals different than the first interval (egg stage duration). Write \code{addInfo = NULL} if you do not provide additional information.
#' @param longevity a numeric value indicating the lifespan of the species in days.
#' @param fecundity a numeric value indicating the mean lifetime fecundity (MLF) as the number of eggs produced for a female.
#' @param a_init a numeric value indicating the initial parameter related to the proportionality  optimized by iterative solution via univariate (1-dim.) minimization. \code{a_init = 2} as default value.
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
                     longevity, fecundity, a_init = 2){


  if(is.null(addInfo)){

    cat("--------------------------------------------------------", "\n\n")
    cat('No additional information. You are only considering the egg stage duration =', eggDuration, "\n\n")
    cat("--------------------------------------------------------", "\n\n")

    output <- .noAddInfo(nInterval   = nInterval,
                         eggDuration = eggDuration,
                         longevity   = longevity,
                         a_init      = a_init)

  }else{

    if(length(addInfo) != nInterval-1) stop('The length of addInfo vector must be equal to nInterval-1')

    output <- .AddInfo(nInterval   = nInterval,
                       eggDuration = eggDuration,
                       longevity   = longevity,
                       a_init      = a_init,
                       addInfo     = addInfo)

  }

  G <- -log((2/fecundity)^(1/nInterval))
  M <- G/output$delta

  abundance    <- numeric(nInterval)
  for(i in seq_len(nInterval)){
    abundance[i]  <- fecundity*(exp(-G))^(i)
  }

  tab <- data.frame(Gnomonic_interval     = seq_len(nInterval),
                    interval_duration_day = round(output$delta, 3),
                    total_duration        = round(cumsum(output$delta), 0),
                    M_day                 = round(M, 3),
                    M_year                = round(M*365, 3),
                    No_Surv               = round(abundance, 0))

  data <- list(a = output$a, G = G, results = tab)

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

  if(!inherits(x, "gnomos"))
    stop("Use only with 'gnomos' objects.")

  data <- x
  cat('Proportionality constant (alpha) =', data$a, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('Constant proportion of the overall natural death rate (G) =', data$G, "\n\n")
  cat("--------------------------------------------------------", "\n\n")
  cat('Main results of gnomonic method:', "\n\n")
  print(data$results)
  return(invisible())
}


#' Plot method for gnomos class
#'
#' @param x an object class 'gnomos'.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param bg a background color for the points.
#' @param cex character expansion in the regression.
#' @param pch the character indicating the type of plotting.
#' @param dayUnits TRUE by default, to show the M values in 1/day unit. FALSE to show the M values in 1/year units.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' model <- gnomonic(nInterval = 7, eggDuration = 2, addInfo = NULL,
#' longevity = 365, fecundity = 200000, a_init = 2)
#'
#' plot(model)
#' @export
#' @method plot gnomos
plot.gnomos <- function(x, xlab = "Gnomonic intervals", ylab = NULL,
                        bg = "lightgrey", cex = 1.75, pch = 21, dayUnits = TRUE, ...){

  if(!inherits(x, "gnomos"))
    stop("Use only with 'gnomos' objects.")

  data <- x

  if(isTRUE(dayUnits)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar = c(4, 6, 1, 1))
    plot(data$results$M_day, type = "b", cex = cex, pch = pch, lwd = 3, ylim = c(0, 1.1*max(data$results$M_day)),
         xlab = xlab , ylab = expression(paste("M (day"^"-1", ")")), bg = bg , axes = FALSE, ...)
    axis(1, seq(from = 1, to = nrow(data$results), by = 1))
    axis(2, las = 2)
    box()
  }else{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar = c(4, 6, 1, 1))
    plot(data$results$M_year, type = "b", cex = cex, pch = pch, lwd = 3, ylim = c(0, 1.1*max(data$results$M_year)),
         xlab = xlab , ylab = expression(paste("M (year"^"-1", ")")), bg = bg , axes = FALSE, ...)
    axis(1, seq(from = 1, to = nrow(data$results), by = 1))
    axis(2, las = 2)
    box()
  }

  return(invisible(NULL))
}

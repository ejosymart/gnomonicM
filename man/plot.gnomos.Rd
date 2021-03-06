% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnomonic-main.R
\name{plot.gnomos}
\alias{plot.gnomos}
\title{Plot method for gnomos class}
\usage{
\method{plot}{gnomos}(
  x,
  xlab = "Gnomonic intervals",
  ylab = NULL,
  bg = "lightgrey",
  cex = 1.75,
  pch = 21,
  dayUnits = TRUE,
  ...
)
}
\arguments{
\item{x}{an object class 'gnomos'.}

\item{xlab}{a title for the x axis.}

\item{ylab}{a title for the y axis.}

\item{bg}{a background color for the points.}

\item{cex}{character expansion in the regression.}

\item{pch}{the character indicating the type of plotting.}

\item{dayUnits}{TRUE by default, to show the M values in 1/day unit. FALSE to show the M values in 1/year units.}

\item{\dots}{Additional arguments to the plot method.}
}
\description{
Plot method for gnomos class
}
\examples{
model <- gnomonic(nInterval = 7, eggDuration = 2, addInfo = NULL,
longevity = 365, fecundity = 200000, a_init = 2)

plot(model)
}

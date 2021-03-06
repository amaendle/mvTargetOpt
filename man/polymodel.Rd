% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvTargetOpt.R
\name{polymodel}
\alias{polymodel}
\title{polymodel}
\usage{
polymodel(dat, degree = 1, weights = NULL)
}
\arguments{
\item{dat}{list containing the data for predictors \code{dat$x} and descriptors \code{dat$y}}

\item{degree}{integer, degree of polynomial regression}

\item{weights}{vector of weights for the data provided in \code{dat}}
}
\value{
linear model
}
\description{
internal function, returns polynomial regression model of degree \code{degree} for the predictor and descriptor data given in \code{dat}. If \code{weights} are provided, the weighted regression model is returned.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvTargetOpt.R
\name{pred.ptch}
\alias{pred.ptch}
\title{pred.ptch}
\usage{
pred.ptch(
  shift,
  pcnr,
  pls1,
  mknormweights,
  dat1d.PC12,
  mindeg,
  tgmean_norm,
  gr2retlimit,
  wfun,
  sto
)
}
\arguments{
\item{shift}{numeric vector}

\item{pcnr}{integer vector}

\item{pls1}{result of PLS1-function}

\item{mknormweights}{boolean}

\item{dat1d.PC12}{matrix}

\item{mindeg}{integer}

\item{tgmean_norm}{numeric vector}

\item{gr2retlimit}{boolean}

\item{wfun}{function}

\item{sto}{boolean}
}
\value{

}
\description{
Internal function. Extends the function \code{opt.onestep}.
}

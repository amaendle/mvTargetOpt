% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvTargetOpt.R
\name{autosolve}
\alias{autosolve}
\title{autosolve}
\usage{
autosolve(
  startx,
  tgmean,
  tgerr,
  reps = 25,
  maxit = 10,
  reality = foo,
  xeps = 0.01,
  pplot = FALSE,
  pcnr = c(1, 2),
  maxarea = NULL,
  useweights = TRUE,
  mknormweights = F,
  gr2retlimit = T,
  mindeg = 0,
  sequential = F,
  tgpcawg = 1,
  yweights = F,
  datlim = NULL,
  knearest = NULL,
  tgdim = 1,
  ylast = NULL,
  sto = T,
  mod.sd = NULL,
  ...
)
}
\arguments{
\item{startx}{numeric matrix, start values}

\item{tgmean}{numeric v4ector, target value vector}

\item{tgerr}{numeric vector, defines acceptable error range}

\item{reps}{integer, number of repeated measurements}

\item{maxit}{integer, maximum number of iterations}

\item{reality}{function, real model}

\item{xeps}{nunmeric, smallest (reasonably) distinguishable epsilon}

\item{pplot}{boolean, diagnostic plots}

\item{pcnr}{integer vector, defines which principal directions will be considered}

\item{maxarea}{numeric matrix, area range, which will be explored}

\item{useweights}{boolean}

\item{mknormweights}{boolean}

\item{gr2retlimit}{boolean}

\item{mindeg}{integer, minimal degree of order of polynomial model}

\item{sequential}{boolean}

\item{tgpcawg}{numeric}

\item{yweights}{boolean}

\item{datlim}{NULL or integer}

\item{knearest}{integer}

\item{tgdim}{integer}

\item{ylast}{integer}

\item{sto}{boolean}

\item{mod.sd}{numeric}

\item{...}{}
}
\value{
data.frame
}
\description{
Applies the iterative optimization algorithm suggested in this package.
To perform the simulations automatically, the true model has to be specified.
When in practice the true model is unknown, use \code{opt.onestep} function instead to get a candidate for the unknown optimum in a single iteration.
}
\examples{
# example 1: 2x2 MModell
tfoo <- function(x) {
  x1 <- x[,1]
  x2 <- x[,2]
  return( data.frame( y1=0.8*x1 - 1.2*x2,
                      y2=0.8*x1 - 1.2*abs(x2)^0.25
 ) ) }
dstgmean <- tfoo(cbind(1.35,1.4))
tgerr <- c(0.025,0.025)
xeps <- 0.01
startx <- expand.grid(x1=c(-1,3),x2=c(-1,3))
set.seed(123)
autosolve(startx,dstgmean,tgerr,reps=7,maxit=10,tfoo, xeps, F, pcnr=c(1,2), mod.sd=0.2)

# example 2, 3x3 model
set.seed(123)
startx <- data.frame(x1=runif(4,-4,4),x2=runif(4,-4,4),x3=runif(4,-4,4))
tfoo <- function(x) {
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  return( data.frame( y1=0.8*x1 - 1.2*x2,
                      y2=0.8*x1 - 1.2*abs(x2)^0.25,
                      y3=0.8*x1 - 0.6*abs(x2)^0.25 + x3
 ) ) }
tgmean <- tfoo(cbind(1.35,1.4,1.5))#tfoo(cbind(0.35,0.4,0.5))
tgerr <- c(0.5,0.5,0.5)
tmp<-autosolve(startx,tgmean,tgerr*0.125,reps=4,maxit=6,tfoo, xeps=0.01, F, pcnr=c(1,2), mod.sd=0.2)

}

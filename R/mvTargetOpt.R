#' mknorm
#'
#' internal function, standardizes a vector \code{x} by vector of means \code{means} and vector of standard deviations \code{sds}.
#'
#' @param x input vector
#' @param means vector of means
#' @param sds vector of standard deviations
#'
#' @return standardized vector
#' @export=FALSE
mknorm <- function(x, means=NULL, sds=NULL) {
  myx <- as.matrix(x)
  if (is.null(means)) {
    xmeans <- colMeans(myx)
  } else {
    xmeans <- means }
  if (is.null(sds)) {
    xsds <- apply(myx,2,sd)
  } else {
    xsds <- sds }
  xsds[is.na(xsds)] <- 0
  if (sum(xsds==0)>0) {
    warning("zero variances in mknorm")
    xsds[xsds==0] <- 1
  }
 # str(xmeans);print(xmeans);
  myx <- sweep(myx,2,as.numeric(xmeans),"-")
  myx <- sweep(myx,2,as.numeric(xsds),"/")
  return(myx)
}

#' mkreg
#'
#' internal function, inverse of \code{mknorm}. Back-transformation of a standardized vector.
#'
#' @param x input vector (standardized)
#' @param means vector of means
#' @param sds vector of standard deviations
#'
#' @return back-transformed vector
#' @export=FALSE
mkreg <- function(x, means, sds) {
  myx <- as.matrix(x)
  myx <- sweep(myx,2,sds,"*")
  myx <- sweep(myx,2,means,"+")
  return(myx)
}

#' mvdistance
#'
#' internal function, computes the Euclidean (\code{euclid=T}) or Manhattan distances (\code{euclid=F}) between vector \code{y} and the rows of matrix \code{xs}.
#'
#' @param xs matrix, each row representing a vector
#' @param y vector
#' @param euclid \code{TRUE} for Euclidean distance, otherwise MAnhattan distance
#'
#' @return vector of the distances
#' @export=FALSE
mvdistance <- function(xs, y, euclid=F) { # summennorm, manhattan-metrik, x--werte als zeilen in xs
  xs<-as.matrix(xs)
  rs <- sweep(xs,2,y,"-")
  rs <- abs(rs)
  if (euclid==T) rs<- rs^2
  rs<- rowSums(rs)
  if (euclid==T) rs <- sqrt(rs)
  return(rs)
}

#' degByBIC
#'
#' internal function, determines the polynomial model order \code{<=maxorder} for the data \code{dat} which minimizes the BIC information criterion
#'
#' @param dat list containing the data for predictors \code{dat$x} and descriptors \code{dat$y}
#' @param maxorder integer, maximal order of polynomial model
#' @param weights vector of weights for the data in \code{dat}
#' @param mindeg integer,  minimal order of polynomial model
#'
#' @return integer, recommended degree for polynomial model
#' @export=FALSE
degByBIC <- function(dat, maxorder=5, weights=NULL, mindeg=0) {
  currBIC <- numeric()
  if (min(maxorder,length(unique(dat$x))-1)<=mindeg) stop(paste("horror in degbybic:",min(maxorder,length(unique(dat$x))-1)))
 # print(paste("mindeg",mindeg))
#  print(paste("maxarchdeg",min(maxorder,length(unique(dat$x))-1)))
 # print(paste("weights",weights))
  for (i in mindeg:min(maxorder,length(unique(dat$x))-1)) {
    pymo <- NULL
    tryCatch(pymo <- polymodel(dat,i, weights=weights), error=function(e) { warning("error in degByBIC::polymodel") })
  #  print(paste("dat",dat))
  #  print(paste("i",i))
  #  print(paste("weights",weights))
  #  print(paste("pymo",pymo))
    if (is.null(pymo)) break
    currBIC <- c(currBIC, BIC(pymo) ) #currBIC <- c(currBIC, BIC(polymodel(dat,i, weights=weights)) )
  #  print(paste("currBIC",currBIC))
    if (i>mindeg) if (currBIC[length(currBIC)] > currBIC[length(currBIC)-1]) break
  }
#  print(paste("currBIC",currBIC))
  if (is.null(currBIC)) stop("fatal error in degbybic - currBIC is null")
  if (which.min(currBIC)-1+mindeg<mindeg) stop("mindeg: this should better not happen")
  return(which.min(currBIC)-1+mindeg)  #return(nnet::which.is.max(-currBIC)-1)
}

#' polymodel
#'
#' internal function, returns polynomial regression model of degree \code{degree} for the predictor and descriptor data given in \code{dat}. If \code{weights} are provided, the weighted regression model is returned.
#'
#' @param dat list containing the data for predictors \code{dat$x} and descriptors \code{dat$y}
#' @param degree integer, degree of polynomial regression
#' @param weights vector of weights for the data provided in \code{dat}
#'
#' @return linear model
#' @export=FALSE
polymodel <- function(dat, degree=1, weights=NULL) {
  if (sum(is.infinite(weights))>0) {
    stop("polymodel: infinite weights")
  }
  if (degree<1) {
    warning("polymodel < 1")
  }
  if (degree<0) {
    stop("polymodel < 0")
  }
  if (degree>length(unique(dat$x))-1) {
    stop("polymodel: degree > unique")
  }
  if (degree==0) {
    return( lm(y ~ 1, data=dat, weights=weights) )
  } else { #if (length(unique(dat$x))==0) {print(dat) ; stop("polymodel: 0 verschiedene Punkte?!?")}
            #print(dat*weights) ; print(weights)
    return( lm(y ~ poly(x, degree, raw=FALSE), data=dat, weights=weights) )
  }
}

#' polymodel
#'
#' internal function, returns polynomial regression model of degree \code{degree} for the predictor and descriptor data given in \code{dat}. If \code{weights} are provided, the weighted regression model is returned.
#'
#' like \code{polymodel}, but here \code{degree} has no default
#'
#' @param dat list containing the data for predictors \code{dat$x} and descriptors \code{dat$y}
#' @param degree integer, degree of polynomial regression
#' @param weights vector of weights for the data provided in \code{dat}
#'
#' @return linear model
#' @export=FALSE
polymodel2 <- function(dat, degree, weights=NULL) {
  if (sum(is.infinite(weights))>0) { ###
    stop("polymodel2: infinite weights") ###
  } ###
  if (degree<1) {
    warning("polymodel2 < 1")
  }
  if (degree<0) {
    stop("polymodel2 < 0")
  }
  if (degree>length(unique(dat$x))-1) {
    stop("polymodel2: degree > unique")
  }
  if (degree==0) {
    return( lm(y ~ 1, data=dat, weights=weights) )
  } else {
    return( lm(y ~ poly(x, degree, raw=TRUE), data=dat, weights=weights) )
  }
}

getroots2 <- function(dat, degree, target, limit=40, weights=NULL, retlimit=FALSE) {
  if (!is.null(weights)) if (sum(is.na(weights))>0) stop("NA values in weights for getroots2()")
  #linear model          # wie plotlm, nur anders
  fm <- polymodel2(dat, degree, weights=weights)

  #find roots
  coeff <- coefficients(fm)
  if (0<sum(is.na(coeff))) {
    coeff[is.na(coeff)] <- 0 ## NA durch 0 ersetzen, bloedes error handling, vorübergehende notloesung
    warning("getroots: seriously, NA was chosen as coefficient in polymodel2 before. replaced by 0")
  }

  coeff[1]<-coeff[1]-target

  # special case degree 0
  if (degree==0) {
    if (identical(coeff[1],0)) {
      warning("getroots2: all points are cutpoints, degree 0")
    } else {
      warning("getroots2: äöü no cutpoints, degree 0")
    }
    #return(data.frame(x = numeric())) # leeres ergebnis zurück
    return(as.numeric(dat$x[which.min(abs(dat$y-target))]) )#return(data.frame(x = as.numeric(dat$x[which.min(abs(dat$y-target))]) ))
  } else {

    roots <- polyroot( coeff )
    roots <- Re(roots)[abs(Im(roots)) < 1]#1e-1]#1e-6]   ############ACHTUNG, nur experimentel auskommentiert
    if (retlimit==TRUE) {
      if (sum(roots < -abs(limit))+sum(roots > abs(limit))>0) warning(paste("getroots2: roots truncated times",sum(roots < -abs(limit))+sum(roots > abs(limit))))
      roots[roots < -abs(limit)] <- -abs(limit)
      roots[roots > abs(limit)] <- abs(limit)
    } else {
    roots <- roots[abs(roots)<=limit] }
    #    if (dim(roots[abs(roots)<=limit])[1]>0)
    # if (dim(roots)[1]<1)
    #  roots<-sign(roots)*limit*(0.5+0.5*runif(1))
    rootdta <- data.frame(x = roots)
    # rootdta$y <- predict(fm, newdata=rootdta)
    return(rootdta$x)
  }
}


#' PIhcheck
#'
#' internal function, determinas a prediction interval for given linear model \code{model} with confidence level \code{alpha}. currently unused.
#'
#' @param model linear model
#' @param alph numeric, confidence level
#'
#' @return vector, upper and lower prediction intervall
#' @export=FALSE
PIhcheck <- function(model, alph=0.05) {
  if (sum(summary(model)$residuals^2)==0) {
    warning("perfect fit in PIhcheck")
    return(c(0,0))
  } else {
    return(c(2*qnorm(1-alph/2)*sqrt(sum(summary(model)$residuals^2)/qchisq(1-alph/2,df=model$df)), 2*qnorm(1-alph/2)*sqrt(sum(summary(model)$residuals^2)/qchisq(alph/2,df=model$df)) ) )
  }
}

#' sumtheothers
#'
#' internal helper function which is used by \code{euclw} .
#'
#' @param x vector
#'
#' @return vector
#' @export =FALSE
sumtheothers <- function(x) {
  n<-length(x)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(x[-i])
  }
  return(result)
}

#' sumthelowers
#'
#' internal helper function which is used by \code{euclw} .
#'
#' @param x numeric vector
#'
#' @return numeric vector
#' @export=FALSE
sumthelowers <- function(x) {
  n<-length(x)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(x[0:(i-1)])
  }
  return(result)
}

#' euclw
#'
#' Takes a matrix/data.frame \code{x} with each coloumn having the scores corresponding to  a principal component.
#' Computes weights as needed in 1d regression models,
#' i.e. returns for each of the 1d-projections on the PC the Euclidean distances from the n-dimensional point in the n-space.
#'
#' @param x matrix, data.frame
#' @param normalize boolean, \code{TRUE} for normalized weights.
#' @param sto boolean, if \code{TRUE} ignore higher principal component coordinates.
#'
#' @return matrix
#' @export=FALSE
euclw <- function(x, normalize=T, sto=T) {
  if(normalize==2) { ##
    x <- apply(x,2,mknorm) ##
  } ##
  if (sto==T) {
    rs <- t(sqrt(apply(x^2,1,sumtheothers)))
  } else {
    rs <- t(sqrt(apply(x^2,1,sumthelowers)))
  }
  dimnames(rs)[[2]] <- dimnames(x)[[2]]
  if (anyNA(rs)) warning("euclw: NAs in rs")
  if(normalize==T) {
    nrs <- apply(rs,2,mknorm)
    if (anyNA(nrs)) warning("euclw: NAs in nrs")
    nrs
    } else rs
}
#' #' euclw.old
#' #'
#' #' Obsolete, has problems wih precision due to \code{rowsums - x^2}. New version: \code{euclw}
#' #'
#' #' Takes a matrix/data.frame \code{x} with each coloumn having the scores corresponding to  a principal component.
#' #' Computes weights as needed in 1d regression models,
#' #' i.e. returns for each of the 1d-projections on the PC the Euclidean distances from the n-dimensional point in the n-space.
#' #'
#' #' @param x matrix, data.frame
#' #' @param normalize boolean, \code{TRUE} for normalized weights.
#' #' @param sto boolean, if \code{TRUE} ignore higher principal component coordinates.
#' #'
#' #' @return matrix
#' #' @export=FALSE
#' euclw.old <- function(x, normalize=T) {
#'
#'   t2 <- rowSums( x^2 )
#'   t1 <- (x^2)*(-1)
#'   rs <- sqrt(sweep(t1,1,t2,"+"))
#'   if(normalize==T) apply(rs,2,mknorm) else rs
#' }

#' tgpca
#'
#' Function which performs a principal component analysis (PCA) on the descriptor variable data (in the target space) given by \code{dat},
#' In order to choose a certain direction through the target point for the projections,
#' \code{wg} has to be set to 1 -- then the target point is chosen as center for the PCA.
#' If \code{wg} lies between 0 and 1, pseudo observations at the target point are created such that a ratio of \code{wg}
#' of the observations are pseudo observations.
#' Then \code{prcomp} is applied  to the standardized data and pseudo data.
#'
#' @param dat matrix, data.frame
#' @param tgmean numeric vector, optional
#' @param tgerr numeric vector, optional
#' @param wg numeric,  weight for the target value. If wg equals 1 or 2 then the pca is performed with the target value as center
#' @param wfun function, weight function
#' @param mknormweights unused
#' @param yweights  boolean, use weights?
#' @param ylast integer or NULL, if integer, ignore observations older than the last \code{ylast} evaluation points.
#'
#' @return returns the results of the pca and some extra stuff
#' @export=FALSE
#'
#' @examples tgpca(matrix(rnorm(20),ncol=2))
tgpca <- function(dat, tgmean=NULL, tgerr=NULL, wg=1, wfun, mknormweights, yweights=F, ylast=NULL) {    #wg=0.901
  if (dim(dat)[2]!=length(tgmean)) stop(paste0("Dimension error in tgpca: ",dim(dat)[2], " and ",length(tgmean)))
  if (!is.null(tgerr)) {
    if (dim(dat)[2]!=length(tgerr)) stop(paste0("tgerr-Dimension error in tgpca: ",dim(dat)[2], " and ",length(tgerr)))
  }

  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }


  dat.orig<-dat
  nrw <- nrow(dat)
  if (!is.null(ylast)) dat <- dat[max(1,(nrw-ylast+1)):nrw,, drop=FALSE]
  # add ratio wg of pseudo observations
  n1 <- dim(dat)[1]
  n2 <- round(wg*n1/(1-wg))
  if(is.infinite(n2)|(n2<0)) {
    n2<-0
    #wg<-1
  }

  if (yweights==T) { if (is.null(tgmean)) stop("tgmean is NULL in tgpca: if (yweights==T)")
    #yw <- ftweights(  wfun(euclw(sweep(allobs[1:n1,],2,as.numeric(tgmean)), mknormweights))  )
    yw <- ftweights( wfun(sqrt(rowSums(sweep(dat,2,as.numeric(tgmean))^2))) )
    dat <- sqrt(yw)*dat
  }
  if (!is.null(tgmean) & !anyNA(tgmean)  & wg > 0 & wg<1) { # & wg != 0
    pseudobs <- rep.row(as.numeric(tgmean),n2) #apply(t(as.numeric(tgmean)),2,rep,n2)
    if (!is.null(colnames(dat))) colnames(pseudobs)<-colnames(dat)
    allobs <- rbind(pseudobs, dat)
  } else allobs <- dat

  rs<-NULL

  if (wg==1) {
    aobs.mean <- colMeans(allobs)
    aobs.sd <- apply(allobs,2,sd)
    if (dim(allobs)[1]==1) {
      aobs.sd[] <- 0
    }
    if (is.null(tgmean) | anyNA(tgmean)) stop("tgmean set incorrectly in tgpca")
  #  print(paste("aobsmean is:", aobs.mean))
  #  print(paste("tgmean is:", tgmean))
    if (sum(aobs.sd==0)>0) {
      warning("zero variances in tgpca/if(wg==1)")
      aobs.sd[aobs.sd==0] <- 1
    }
    rs$pca <- prcomp(mknorm(allobs), center=(as.numeric(tgmean)-aobs.mean)/aobs.sd , scale=F) # wenn tgmean nicht null
  } else if (wg==2) { #actually this should give the same result as above
    if (is.null(tgmean) | anyNA(tgmean)) stop("tgmean nicht korrekt gegeben in tgpca")
  #  print("allobs"); print(allobs); print("normiert:"); print(mknorm(allobs, means=tgmean))
    rs$pca <- prcomp(mknorm(allobs, means=tgmean), center=F , scale=F)
  } else {
    rs$pca <- prcomp(mknorm(allobs))
  }
  # remove pseudo observations from x
  if (wg<1&0<wg) rs$pca$x <- rs$pca$x[-(1:n2),]
  if(!is.null(tgmean)) {
    if (wg>=1) { # hint: for wg>=1 it should hold pcatg==0
      rs$pcatg <- mknorm(t(as.numeric(tgmean)), tgmean, apply(allobs,2,sd)) %*% rs$pca$rotation
    } else {
      rs$pcatg <- mknorm(t(as.numeric(tgmean)), colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
    }
  } else rs$pcatg <- rep(NA,dim(dat)[2])
  if(!is.null(tgerr)&!anyNA(tgerr)) {
    # print(paste("mult,wg=",wg))
    # print(mknorm(t(as.numeric(tgerr)), 0*colMeans(allobs), apply(allobs,2,sd)))
    # print(rs$pca$rotation)
    rs$pcatgerr <- mknorm(t(as.numeric(tgerr)), 0*colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
  } else {
    namedNA <-t(rep(NA,dim(dat)[2]))
    dimnames(namedNA)[[2]] <- paste0("PC",1:dim(dat)[2])
    rs$pcatgerr <- namedNA #rep(NA,dim(dat)[2]) # because of: myypca$pcatgerr[, "PC1"]
  }
  # add passive obs again
  if (!is.null(ylast)) {
    if (wg>=1) {
      rs$pca$x <- mknorm(dat.orig, tgmean, apply(allobs,2,sd)) %*% rs$pca$rotation
    } else {
      rs$pca$x <- mknorm(dat.orig, colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
    }
  }
  rs$allobsmean <- colMeans(allobs)  # could save it above to make the code less redundant
  rs$allobssd <- apply(allobs,2,sd)   # "
  return(rs)
}

# #' tgpcabeforeweight1
# #'
# #' unused, old version of tgpca before the implementation of special case weight=1
# #'
# #' @param dat
# #' @param tgmean
# #' @param tgerr
# #' @param wg
# #'
# #' @return
# #' @export
# #'
# #' @examples
# tgpcabeforeweight1 <- function(dat, tgmean=NULL, tgerr=NULL, wg=0.901) {   #unelegant (2x): t(as.numeric(tgmean))
#   if (dim(dat)[2]!=length(tgmean)) stop(paste0("Dimension error in tgpca: ",dim(dat)[2], " and ",length(tgmean)))
#   rep.row<-function(x,n){
#     matrix(rep(x,each=n),nrow=n)
#   }
#   # Anteil wg an Pseudobeobachtungen hinzufügen
#   n2 <- round(wg*dim(dat)[1]/(1-wg))
#   if (!is.null(tgmean) & !anyNA(tgmean)  & wg != 0) {
#     pseudobs <- rep.row(as.numeric(tgmean),n2) #apply(t(as.numeric(tgmean)),2,rep,n2)
#     if (!is.null(colnames(dat))) colnames(pseudobs)<-colnames(dat)
#     allobs <- rbind(pseudobs, dat)
#   } else allobs <- dat
#
#   rs<-NULL
#   rs$pca <- prcomp(mknorm(allobs))
#   # pseudos aus x entfernen
#   rs$pca$x <- rs$pca$x[-(1:n2),]
#   #rs$pcatg <- ifelse(is.null(tgmean),NA,mknorm(t(as.numeric(tgmean)), colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation)
#   if(!is.null(tgmean)) {
#     rs$pcatg<-mknorm(t(as.numeric(tgmean)), colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
#   } else rs$pcatg <- rep(NA,dim(dat)[2])
#   #rs$pcatgerr <- ifelse(is.null(tgerr)|anyNA(tgerr),rep(NA,dim(dat)[2]),mknorm(t(as.numeric(tgerr)), 0*colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation)
#   if(!is.null(tgerr)&!anyNA(tgerr)) {
#     rs$pcatgerr <- mknorm(t(as.numeric(tgerr)), 0*colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
#   } else {
#     namedNA <-t(rep(NA,dim(dat)[2]))
#     dimnames(namedNA)[[2]] <- paste0("PC",1:dim(dat)[2])
#     rs$pcatgerr <- namedNA #rep(NA,dim(dat)[2]) # wegen: myypca$pcatgerr[, "PC1"]
#   }
#   rs$allobsmean <- colMeans(allobs)
#   rs$allobssd <- apply(allobs,2,sd)
#   return(rs)
# }

#' mimascores
#'
#' Internal function
#'
#' @param maxarea matrix
#' @param mwgs numeric
#'
#' @return matrix
#' @export=FALSE
mimascores <- function(maxarea, mwgs) {
  if (sum(apply(maxarea,1,diff)<0)>0) stop("mimascores: maximal parameter values must be larger than minimal parameter values")
  mi <- maxarea[,1,drop=F]
  ma <- maxarea[,2,drop=F]

  mimat <- as.vector(mi)*mwgs
  mamat <- as.vector(ma)*mwgs

  mins <- mimat
  mins[mimat>mamat] <- mamat[mimat>mamat]

  maxs <- mamat
  maxs[mimat>mamat] <- mimat[mimat>mamat]

  mmin <- apply(mins,2,max)
  mmax <- apply(maxs,2,min)

  rs <- t(rbind(mmin,mmax))
  return(rs)
}
#' mimascores2
#'
#' Internal function
#'
#' @param maxarea matrix
#' @param mwgs numeric
#'
#' @return matrix
#' @export=FALSE
mimascores2 <- function(maxarea, mwgs) {
  if (sum(apply(maxarea,1,diff)<0)>0) stop("mimascores2: maximal parameter values must be larger than minimal parameter values")
  mi <- maxarea[,1,drop=F]
  ma <- maxarea[,2,drop=F]
  iwgs <- solve(mwgs)

  mimat <- sweep(iwgs,2,as.vector(mi),function(x,a) a/x )
  mamat <- sweep(iwgs,2,as.vector(ma),function(x,a) a/x )

  mins <- mimat
  mins[mimat>mamat] <- mamat[mimat>mamat]

  maxs <- mamat
  maxs[mimat>mamat] <- mimat[mimat>mamat]

  mmin <- apply(mins,1,max)
  mmax <- apply(maxs,1,min)

  rs <- t(rbind(mmin,mmax))
  return(rs)
}
#' ptchoose
#'
#' Internal function. Continues search in unexplored space, when no solution can be found.
#'
#' @param ptchoice numeric, 1 to continue search outside the explored parameters
#' @param dat1d data.frame, predictor data must be in \code{dat1d$x}
#' @param tgmean_norm numeric vector
#' @param maxarea NULL or matrix
#' @param xmeans numeric vector
#' @param xsds numeric vector
#' @param pls1 result of the PLS1 function
#' @param jjj integer, counts how often this function has been called
#'
#' @return matrix
#' @export=FALSE
ptchoose <- function(ptchoice=ptchoice, dat1d, tgmean_norm, maxarea, xmeans, xsds, pls1, jjj) {# jjj only needed for mimascores in pt3
  #        # choose new coordinate in 1d space
  #        cat("1d points avail.:",paste(" ",sort(c(0, unique(dat1d[[1]])))))
  #        cat("1d points avail.:",paste(" ",diff(sort(c(0, unique(dat1d[[1]]))))))
  #        #all new possible measure points:
  #        #cat("meas at:",paste(" ", head( sort(c(0, unique(dat1d[[1]]))) ,-1) + diff(sort(c(0, unique(dat1d[[1]]))))*0.5))
  #        cat("meas at:",paste(" ", head( sort(c(0, unique(dat1d[[1]]))) ,-1) + diff(sort(c(0, unique(dat1d[[1]]))))*0.5))
  #        cat("ext points",
  #            min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2,
  #            " and ",
  #            max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2)
  if (ptchoice==1) { # continue search outside, not so useful, when limits for the parameters have been set
    #message(paste0("Schritt ","i",": modell without roots in any direction: PC ",jjj,". Continue new measurement outrside the already explored interval."))
    if (min(abs(dat1d$y[dat1d$x==min(dat1d$x)]-tgmean_norm)) < min(abs(dat1d$y[dat1d$x==max(dat1d$x)]-tgmean_norm))) {
      newx1d <- min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2
    } else {
      newx1d <- max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    }
  } else if (ptchoice==2) { # 0 better outside of sort?
    newx1d <- c((head( sort(c(0, unique(dat1d[[1]]))) ,-1) + diff(sort(c(0, unique(dat1d[[1]]))))*0.5)[which.max(diff(sort(c(0, unique(dat1d[[1]])))))],
                min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2,
                max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    )
  } else if (ptchoice==20) { # as above, without 0,  i.e. continue search inside (1 point) und am Randand on the outside (2points)
    newx1d <- c((head( sort(c(unique(dat1d[[1]]))) ,-1) + diff(sort(c(unique(dat1d[[1]]))))*0.5)[which.max(diff(sort(c(unique(dat1d[[1]])))))],
                min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2,
                max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    )
  } else if (ptchoice==200) { # aas above, but search inside (1 Punkt) and outside (1Punkt)
    if (min(abs(dat1d$y[dat1d$x==min(dat1d$x)]-tgmean_norm)) < min(abs(dat1d$y[dat1d$x==max(dat1d$x)]-tgmean_norm))) {
      tmp <- min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2
    } else {
      tmp <- max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    }
    newx1d <- c((head( sort(c(unique(dat1d[[1]]))) ,-1) + diff(sort(c(unique(dat1d[[1]]))))*0.5)[which.max(diff(sort(c(unique(dat1d[[1]])))))],
                tmp)
  } else if (ptchoice==3) {# choose only 1 new point, maximum distance of coordinate
    if (!is.null(maxarea)) { # if maxarea is set...
      xcols <- grep("x.", names(xmeans), value = TRUE)
      mima <- mimascores2(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols])), pls1$mod.wgs)[jjj,]
      upts <- unique(dat1d[[1]])  # [[1]]?
      upts <- upts[which(upts>mima[1])] # in the interval?
      upts <- upts[which(upts<mima[2])]
      upts <- sort(c(upts, mima)) # mima-limits, sorted
      measureatpts <- head(upts ,-1) + diff(upts)*0.5
      # if limit is given, choose as distance 2*(smallestobservedxvalue-minpossiblexval) or
      #                                           2*(biggestobservedxvalue-maxpossiblexval)
      if (length(upts)<2) stop("impossible: length(upts)<2 for ptchoice3")
      if (length(upts)<3) {
        warning("unhandled: length(upts)<3 for ptchoice3")
        choosept <- 1 # mal probieren, müsste passen
      } else choosept <- which.max(diff(upts)* c(2,rep(1,length(upts)-3),2)  )
      newx1d <- measureatpts[choosept]

      if (anyNA(newx1d)){ # in case of NAs -> stop
        xcols <- grep("x.", names(xmeans), value = TRUE)
        ttt<-unique(dat1d[[1]])
        ttt <- ttt[which(ttt>mima[1])]
        ttt <- ttt[which(ttt<mima[2])]
        plot(pls1$x.scores%*%solve(pls1$mod.wgs),xlim=c(-9,9),ylim=c(-9,9))
        points(pls1$x.scores[,1,drop=F]%*%solve(pls1$mod.wgs)[1,,drop=F],col="green")
        points(matrix(mima)%*%solve(pls1$mod.wgs)[1,,drop=F],col="red")
        plot(mkreg(pls1$x.scores%*%solve(pls1$mod.wgs),xmeans[xcols],xsds[xcols]),xlim=c(-30,30),ylim=c(-30,30))
        points(mkreg(pls1$x.scores[,1,drop=F]%*%solve(pls1$mod.wgs)[1,,drop=F],xmeans[xcols],xsds[xcols]),col="green")
        points(mkreg(matrix(mima)%*%solve(pls1$mod.wgs)[1,,drop=F],xmeans[xcols],xsds[xcols]),col="red")
        rm(xcols)
        stop("newx1d hat NA in prediction")
      }
      warning(paste("ptchoose: chosenpt in dir",jjj)); #print(newx1d)
    } else { # no limits? Then continue as in ptchoice=1
      if (min(abs(dat1d$y[dat1d$x==min(dat1d$x)]-tgmean_norm)) < min(abs(dat1d$y[dat1d$x==max(dat1d$x)]-tgmean_norm))) {
        newx1d <- min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2
      } else {
        newx1d <- max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
      }
    }
  } # end: if (ptchoice==3)
  return(newx1d)
}
#' ftweights
#'
#' Internal function (for \code{pred.patch}, \code{ored.solution}). Replaces infinite numbers (weights) by large finite values. (clipping)
#'
#' @param w1 numeric vector
#'
#' @return numeric vector
#' @export=FALSE
ftweights <- function(w1) {
  if (0<length(w1)) {
    if (0<sum(is.infinite(w1))) warning("ftweights: infinite weights were clipped")
    w1[is.infinite(w1)]<-1e+306    #1.797693e+308
    w1[w1>1e+306 ]<-1e+306
    return(w1)
  } else w1
}

#' uniqP
#'
#' Internal function.
#' Returns a subset of the  unique points of \code{newx} whgich are not
#' (up to an epsilon \coed{xeps}) identical to the points in \code{datx}.
#'
#' @param newx matrix
#' @param xeps numeric
#' @param datx matrix
#'
#' @return matrix
#' @export=FALSE
uniqP <- function(newx, xeps, datx) {
  if ((!is.null(newx))&(!is.null(unlist(newx)))) {
    dstncs <- matrix(nrow = dim(newx)[1],ncol=dim(datx)[1])
    if (dim(newx)[1]>0)   #nur, wenn neue Punkte vorhanden:
      for (cnt in 1:dim(newx)[1]) {
        dstncs[cnt,] <- mvdistance(datx, newx[cnt,])
      }
    newx <- newx[ suppressWarnings(apply(dstncs,1,min))>xeps,,drop=F]
  }
  return(newx)
}
#' repna
#'
#' Internal function. Replaces NAs in a vector \code{x} by 0.
#'
#' @param x numeric vector
#'
#' @return numeric vector
#' @export=FALSE
repna <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}
#' nclose2mean
#'
#' Internal function. Returns a subset of size \code{n} of the points in \code{datx} which are closest to \code{center}.
#'
#' @param datx matrix
#' @param center numeric vector
#' @param n integer
#'
#' @return matrix
#' @export=FALSE
nclose2mean <- function(datx, center, n=1) {
  # if (!dim(datx)[2]=dim(center)[2]) stop("wrong dimensions in nclose2mean")
   nwas <- nrow(datx)
  if (!is.null(datx)) {
    dstncs <- numeric()
    if (dim(datx)[1]>0)   #íf there are new points:
      for (cnt in 1:dim(datx)[1]) {
        dstncs[cnt] <- mvdistance(datx[cnt,,drop=F],center)
      }
    ind <- sort.int(dstncs,index.return = T)$ix
    if (0<length(ind)) {
      ind <- ind[1:min(n,length(ind))]
      datx <- datx[ind,,drop=F]
    }
  }
  nis <- nrow(datx)
 # message(paste("nwas",nwas,"nis",nis))
  return(datx)
}

#' pred.ptch
#'
#' Internal function. Extends the function \code{pred.solution}.
#'
#' @param shift numeric vector
#' @param pcnr integer vector
#' @param pls1 result of PLS1-function
#' @param mknormweights boolean
#' @param dat1d.PC12 matrix
#' @param mindeg integer
#' @param tgmean_norm numeric vector
#' @param gr2retlimit boolean
#' @param wfun function
#' @param sto boolean
#'
#' @return
#' @export=FALSE
pred.ptch <- function(shift, pcnr, pls1, mknormweights, dat1d.PC12, mindeg, tgmean_norm, gr2retlimit, wfun, sto) { #pcnr, princ comps that have to be dealt with
  if (!is.numeric(shift)) stop("pred.ptch: shift must be numeric")
  xdim <- length(shift)

  pjjj <- setdiff(1:length(shift),pcnr)

  newxs <- list()
  newxs[[1]] <- matrix(shift,nrow=1,ncol=xdim)

  for (jjj in pcnr) {
    lind <- length(newxs)
    newxs[[lind+1]] <-  matrix(NA,nrow=0,ncol=xdim)

    for (z in 1:nrow(newxs[[lind]])) {
      w <-ftweights( wfun(euclw(sweep(pls1$x.scores, 2, unlist(repna(newxs[[lind]][z,]))), mknormweights, sto)))
      degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=w[,jjj], mindeg=mindeg)
      newroots <- getroots2(dat1d.PC12[[jjj]],
                            degjjj,
                            tgmean_norm,
                            limit=1*1.6*diff(range(dat1d.PC12[[jjj]]$x)), #!!!!!!!!
                            weights=w[,jjj], retlimit = gr2retlimit)
      if (length(newroots)>0) {
        updxs <- repna(newxs[[lind]][z,])
        updxs <- matrix(updxs, ncol=xdim, nrow=length(newroots), byrow=T)
        updxs[,jjj] <- newroots
        newxs[[lind+1]] <- rbind(as.matrix(newxs[[lind+1]]),
                                 as.matrix(updxs) )
      } else {
        updxs <- newxs[[lind]][z,]
        #warning("neue nullkoordinaten mit NA eingefügt")
        newxs[[lind+1]] <- rbind(as.matrix(newxs[[lind+1]]),
                                 t(updxs) ) # oder na statt 0?
      }
      names(newxs) <- 1:xdim
    }
  }
}

#' pred.solution
#'
#' Performs one iteration of the implemented optimization procedure.
#' Suggest points for future measurements in order to find a solution that reaches the desired target value.
#' One or more suggested points are returned.
#'
#'
#' @param dat data.frame, data set
#' @param tgmean numerig vector, target mean
#' @param tgerr NULL or numeric vector, maximally accepted deviation from target value
#' @param xeps numeric, smallest delta
#' @param pcnr integer vector,  numbers the principal components that shall be considered
#' @param maxarea matrix or NULL, maximal area that can/should be explored
#' @param ptchoice integer
#' @param useweights boolean
#' @param mknormweights boolean
#' @param allpts boolean
#' @param gr2retlimit boolean
#' @param bpcenter boolean
#' @param mindeg integer, minimal degree for polynomial mdel
#' @param wfun function, weight function
#' @param sequential boolean
#' @param ptchng boolean
#' @param nptc integer
#' @param tgpcawg numeric
#' @param betterweights boolean
#' @param yweights boolean
#' @param datlim integer or NULL
#' @param knearest integer, if specified only the \code{knearest} nearest observations to the target value are considered
#' @param tgdim integer
#' @param ylast integer, if a positive integer is defined, observations from the last \code{ylast} iterations are used only
#' @param sto boolean
#' @param ...
#'
#' @return matrix with recommended points (process parameters for future measurements) in each line
#' @export TRUE
#' @examples
#' library(mvTargetOpt)
#' # Let there be the following true model tfoo:
#' tfoo <- function(x) {
#'   x1 <- x[,1]
#'   x2 <- x[,2]
#'   return( data.frame(y1=0.8*x1 - 1.2*x2,
#'                      y2=0.8*x1 - 1.2*abs(x2)^0.25) )
#' }
#' # assume we have measurements dat taken at startx:
#' startx <- expand.grid(x1=c(0,6,7,8,9,10,11,12),x2=c(-1,3,4,5,6,7))
#' dat <- cbind(startx,tfoo(startx))
#' # assume we want to find process parameters close to tgmean
#' tgmean <- tfoo(cbind(0.35,0.4))
#' # make a guess for a solution based on the specified parameters:
#' pred.solution(dat, tgmean=tgmean, pcnr=1:2)
#' # the same, but adding the error interval, and enable variance control
#' pred.solution(dat, tgmean=tgmean,tgerr=c(0.2,0.2), pcnr=1:2)
pred.solution <- function(dat,tgmean,tgerr=NULL,xeps=0.001,pcnr, maxarea=NULL, ptchoice=1, useweights=TRUE, mknormweights=F, allpts=F, gr2retlimit=TRUE,bpcenter=F,mindeg=0,wfun=function(x) {(1/x)^2}, sequential=F, ptchng=F, nptc=0,tgpcawg=1,betterweights=F,yweights=F,datlim=NULL,knearest=NULL,tgdim=1,ylast=NULL,sto=T,...) {
  debug<-F
  # nptc: how often has ptchoice been performed
  # check data structures for tgmean, tgerr, maxarea
  tgmean<-matrix(tgmean,nrow=1)
  if (!is.null(tgerr)) {
    tgerr <- matrix(tgerr,nrow=1)
  } else tgerr <- matrix(NA,nrow=1,ncol=length(tgmean)) #tgerr<- ifelse(!is.null(tgerr), matrix(tgerr,nrow=1), matrix(NA,nrow=1,ncol=length(tgmean)))

  stepi <- 1 + max(dat$nri)
  if (!is.null(datlim)) {
    dat <- dat[dat$nri >= (max(dat$nri)-datlim) ,] #tail(dat, datlim)
  } else if (!is.null(knearest)&bpcenter!=T) { # k nearest to last observation # !=T, i.e. will be executed also for bpcenter=2
    # alternative: add limit at bpcenter  wfun(pmax(euclw(sweep(pls1$x.scores,2,minpoint), mknormweights),0.0000001))
    ndat <- mknorm(dat$x)
    lobs <- ndat[nrow(ndat),] # ftweights(wfun(sqrt(rowSums(myypca$pca$x^2))))
    idcs <- sort(mvdistance(ndat, lobs, TRUE), index.return=T)$ix[1:min(knearest,nrow(ndat))]  # distance to last observation#sort and select
    dat <- dat[idcs ,]
  }

  if (is.matrix(maxarea)) {
    xdim <- dim(dat$x)[2]#dim(dat)[2]-dim(tgmean)[2]
    if (dim(maxarea)[1] != xdim |  dim(maxarea)[2] != 2)
      stop("wrong dimensions of maxarea in pred.solution")
    rm(xdim)
  } else if (!is.null(maxarea)) stop("maxarea must be a matrix or NULL in pred.solution")

  # init PI; Ydim, Xdim; xmeans, xsds
  PI <- data.frame(pi.l=numeric(), pi.r=numeric(), pc=integer(), nr=integer())
  Ydim <- length(tgmean)
  if (is.matrix(dat)||is.null(dat$y)||is.null(dat$x)) dat <- data.frame(x=I(as.matrix(dat[,1:(dim(dat)[2]-Ydim)])), y=I(as.matrix(dat[,-(1:(dim(dat)[2]-Ydim))])))
  Xdim <- dim(dat$x)[2]
  modeg <- NULL # # Container for used degree for polynomial model

  #WARNING if  DIMENSION and tgmean don't match
  #dimension of x and y
  if (is.vector(dat$y)) {Ydimcheck <- 1} else {Ydimcheck <- dim(dat$y)[2]}
  if (Ydimcheck!=Ydim) stop("pred.solution: dimension of dat$y does not match dimension of tgmean")
  rm(Ydimcheck)

  # DIMENSION REDUCTION X
  xmeans <- colMeans(dat)
  xsds <- apply(dat,2,sd)
  # DIMENSION REDUCTION Y (tgpca)
  if (Ydim>1) {
    myypca <- tgpca(dat$y,tgmean,tgerr, wg=tgpcawg, wfun=wfun, mknormweights=mknormweights, yweights=yweights,ylast=ylast)
    y1d <- myypca$pca$x[,tgdim,drop=F]
    tgmean_norm <- myypca$pcatg[,tgdim]
    tgerr_norm <- myypca$pcatgerr[,tgdim]
    if (betterweights==T) {
      ywghts <- ftweights(wfun(sqrt(rowSums(myypca$pca$x[,-tgdim,drop=F]^2))))
    } else if (betterweights==2) {
      ywghts <- ftweights(wfun(sqrt(rowSums(myypca$pca$x^2))))
    } else if (betterweights==3) {
      ywghtsvar <- rowSums(myypca$pca$x[,-tgdim,drop=F]^2)
    } else ywghts<-1
  } else {  betterweights=F
            ywghts<-1
            ywghtsvar <- 1
    y1d <- dat$y
    tgmean_norm <- tgmean[[1]] ##[[]]
    if (!anyNA(tgerr)) tgerr_norm <- tgerr[[1]] else tgerr_norm <- NA
  }
  # PLS
  doCrosVal <- !(nrow(dat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL - contact Gaston Sanchez for bug report?
  tryCatch( { pls1 = plsdepot::plsreg1(dat$x, y1d, comps = Xdim, crosval=doCrosVal) }
            , error= function(e) { print("plsreg1 throws error:"); print(e); stop("plsreg1 failed in pred.solution")})
  #pls1 = plsdepot::plsreg1(dat$x, y1d, comps = Xdim)
  plsXdim <- dim(pls1$x.scores)[2]
  # treat single dimensions for each principal componant in the list separately
  dat1d.PC12 <- list()
  for (jjj in 1:plsXdim) {
    dat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,paste0("t",jjj)], y=as.numeric(y1d))
  }
  # list shifted 1d point plus shift (for additive==T with weights)
  newx1d.PC12.ext <- list()
  # for each principal component optimize one-dimensional model
  newx1d.PC12 <- list()

  #determine weights ####  bpcenter=T
  if (useweights==TRUE) {
    if (bpcenter==F ) {
      weightsall <- wfun(euclw(pls1$x.scores, mknormweights, sto))
      weightsvar <- (euclw(pls1$x.scores, mknormweights, sto))^2
    } else {
      if (is.vector(dat$y)) {tmp1 <- abs(dat$y-tgmean[[1]]); stop("dat$y is vector in pred.solution")} else {
        tmp1 <- sweep(dat$y,2,as.numeric(tgmean)) # or dat$y?
        tmp1 <- sqrt(rowSums( tmp1^2 ))
      }
      #find minimal distance
      minpoint <- pls1$x.scores[which.min(tmp1),]
      keepidcs <- rep(1,nrow(pls1$x.scores))
      if (!is.null(knearest)&bpcenter==T) {
        keepidcs <- rep(0,nrow(pls1$x.scores))
        idcs <- sort(mvdistance(pls1$x.scores, minpoint, TRUE), index.return=T)$ix[1:min(knearest,nrow(pls1$x.scores))]  # distances to last observation#sort and select
        keepidcs[idcs] <- 1
      }
      weightsall <- keepidcs*wfun(pmax(euclw(sweep(pls1$x.scores,2,minpoint), mknormweights,sto),0.0000001))
      weightsvar <- (euclw(sweep(pls1$x.scores,2,minpoint), mknormweights,sto))^2
    }
  } else {
    weightsall=NULL
    weightsvar <- NULL
  }
  weightsall <- ftweights(weightsall)

  pjjj <- NULL # remember past jjj's
  newx1d.seq <- list()
  newx1d.ext <- matrix(NA,nrow=0,ncol=plsXdim) #? or Xdim
  if (length(pcnr)>plsXdim) {
    pcnr <- pcnr[1:plsXdim]
    warning(paste("pcnr replaced by",pcnr))
  }
  for (iseq in 1:length(pcnr)) newx1d.seq[[iseq]] <- matrix(NA,nrow=0,ncol=iseq)
  for (jjj in pcnr) {                          # don't stop() here!
    if (sum(!is.na(weightsall[,jjj]))==0) .debug<<-T
    if (sum(!is.na(weightsall[,jjj]))==0) warning("no non-NAs in pred.solution for weightsall[,jjj] (with mknormweights=T?) - degByBIC will probably fail")

    if (sequential==F) {
      if (betterweights==3) {
        degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=ftweights(1/(ywghtsvar+weightsvar[,jjj])), mindeg=mindeg)
      } else {
        degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=weightsall[,jjj], mindeg=mindeg)
      }
    modeg <- c(modeg, degjjj)
    if (betterweights==F) {
    newx1d.PC12[[jjj]] <- getroots2(dat1d.PC12[[jjj]],
                                    degjjj,
                                    tgmean_norm,
                                    limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)), weights=weightsall[,jjj],retlimit = gr2retlimit)#weights)
    } else {
      if (betterweights==3) {
        newx1d.PC12[[jjj]] <- getroots2(dat1d.PC12[[jjj]],
                                        degjjj,
                                        tgmean_norm,
                                        limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)), weights=ftweights(1/(ywghtsvar+weightsvar[,jjj])),retlimit = gr2retlimit)#weights)
      } else {
      newx1d.PC12[[jjj]] <- getroots2(dat1d.PC12[[jjj]],
                                      degjjj,
                                      tgmean_norm,
                                      limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)), weights=(ywghts+weightsall[,jjj])/2,retlimit = gr2retlimit)#weights)
      }
    }
    } else  newx1d.PC12[[jjj]] <- NULL
    if (sequential==T) { # above with assumption bpcenter=F, useweights=T?
      # remember past jjj's:
      pjjj <- c(pjjj,jjj)

      if (length(pjjj)==1) {
        w <- ftweights(wfun(euclw(pls1$x.scores, mknormweights,sto)))
        wvar <- euclw(pls1$x.scores, mknormweights,sto)^2
        if (length(dim(w))<2) {
          print(w);
          print(wfun(euclw(pls1$x.scores, mknormweights,sto)));
          print(euclw(pls1$x.scores, mknormweights,sto));
          print(pls1$x.scores);
          stop("dim von w < 2") }
        if (betterweights==3) {
          degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=ftweights(1/(ywghtsvar+wvar[,jjj])), mindeg=mindeg)
        } else {
          degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=w[,jjj], mindeg=mindeg)
        }
        modeg <- c(modeg, degjjj)
      #  if (degjjj==0) stop("degjjj=0 - f2")
        if (betterweights==F) {
        newx1d.seq[[1]] <- data.frame(getroots2(dat1d.PC12[[jjj]],
                                        degjjj,
                                        tgmean_norm,
                                        limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)),
                                        weights=w[,jjj], retlimit = gr2retlimit)    )
        } else {
          if (betterweights==3) {
            newx1d.seq[[1]] <- data.frame(getroots2(dat1d.PC12[[jjj]],
                                                    degjjj,
                                                    tgmean_norm,
                                                    limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)),
                                                    weights=ftweights(1/(ywghtsvar+wvar[,jjj])), retlimit = gr2retlimit)    )
          } else {
          newx1d.seq[[1]] <- data.frame(getroots2(dat1d.PC12[[jjj]],
                                                  degjjj,
                                                  tgmean_norm,
                                                  limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)),
                                                  weights=(ywghts+w[,jjj])/2, retlimit = gr2retlimit)    )
          }
        }#weights=sqrt(ywghts)*sqrt(w[,jjj])
        names(newx1d.seq[[1]]) <- pjjj
      } else { #message(paste("Xdim",Xdim,"lpcnr",pcnr))
        if (nrow(newx1d.seq[[length(pjjj)-1]])>0) { message("i:|")
        for (z in 1:nrow(newx1d.seq[[length(pjjj)-1]])) { # improve each past point with additional principal component
          #weight for data points minus shift
          shift <- list(); for (ishift in 1:Xdim) shift[[ishift]] <- 0 # init shift = 0
          for (ishiftagain in head(pjjj,-1)) {
             #coloumn names are used -- have to be updated/added here
            colnames(newx1d.seq[[length(pjjj)-1]]) <- head(pjjj,-1)
            shift[[ishiftagain]] <- repna(newx1d.seq[[length(pjjj)-1]][z,paste(ishiftagain)]) }
                                                                                                  # check code above...
                                                                    # attention: in the following unlist und x.scores
          w <-ftweights( wfun(euclw(sweep(pls1$x.scores, 2, unlist(shift)), mknormweights,sto)))
          wvar <- euclw(sweep(pls1$x.scores, 2, unlist(shift)), mknormweights,sto)^2
          #if (!exists("w")) stop("w missing in pred")
          if (betterweights==3) {
            degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=ftweights(1/(ywghtsvar+wvar[,jjj])), mindeg=mindeg)
          } else {
            degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=w[,jjj], mindeg=mindeg)
          }
          #modeg <- c(modeg, degjjj)
         # if (degjjj==0) stop("degjjj=0 - f3")
          if (betterweights==F) {
          newroots <- getroots2(dat1d.PC12[[jjj]],
                                             degjjj,
                                             tgmean_norm,
                                             limit=1*1.6*diff(range(dat1d.PC12[[jjj]]$x)), #!!!!!!!!
                                             weights=w[,jjj], retlimit = gr2retlimit)
          } else {
            if (betterweights==3) {
              newroots <- getroots2(dat1d.PC12[[jjj]],
                                    degjjj,
                                    tgmean_norm,
                                    limit=1*1.6*diff(range(dat1d.PC12[[jjj]]$x)),
                                    weights=ftweights(1/(ywghtsvar+wvar[,jjj])), retlimit = gr2retlimit)
            } else {
            newroots <- getroots2(dat1d.PC12[[jjj]],
                                  degjjj,
                                  tgmean_norm,
                                  limit=1*1.6*diff(range(dat1d.PC12[[jjj]]$x)),
                                  weights=(ywghts+w[,jjj])/2, retlimit = gr2retlimit)
            }
          }

          if (length(newroots)>0) {
            #print("new points")
            # print(newx1d.seq[[length(pjjj)-1]])
            # print(newroots)
            # print(w[,jjj]/sum(w[,jjj]))
            newx1d.seq[[length(pjjj)]] <- rbind(as.matrix(newx1d.seq[[length(pjjj)]]),
                                                as.matrix(cbind(newx1d.seq[[length(pjjj)-1]][z,,drop=F],data.frame(  newroots  ) )) )
          } else {
             warning("added new zero-coordinates")
            # print(cbind(newx1d.seq[[length(pjjj)-1]][z,],0 ))
            newx1d.seq[[length(pjjj)]] <- rbind(as.matrix(newx1d.seq[[length(pjjj)]]),
                                                as.matrix(cbind(newx1d.seq[[length(pjjj)-1]][z,,drop=F],NA )) ) # or na instead of 0?
          }
          names(newx1d.seq) <- pjjj
        }
        } else { message("i:nopts")
          newx1d.seq[[length(pjjj)]] <- matrix(NA,nrow=0,ncol=length(pjjj))
        }
      }
    }

    #Varianzkriterium? tgerr ode tgerr_norm
    if (betterweights==3) {
      PI.mima <- t(PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]],weights=ftweights(1/(ywghtsvar+weightsvar[,jjj])), mindeg=mindeg),
                                      weights=weightsall[,jjj])))
    } else {
    PI.mima <- t(PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]],weights=weightsall[,jjj], mindeg=mindeg),weights=weightsall[,jjj])))
    }
 #   message("Prädiktions-Intervallbreite für PC",jjj," in [", PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)))[1],",",PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)))[2],"]")
    if(!anyNA(tgerr_norm)) {   # here no betterweights3 in following line...
      if (PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]],weights=weightsall[,jjj], mindeg=mindeg),weights=weightsall[,jjj]))[1]>2*tgerr_norm & !anyNA(tgerr)) {
        is2high <- 1
       # #if (vartoobig==-1) {vartoobig<-i}
       # warning("Varianz zu groß. Prädiktions-Intervallbreite von [", PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg))),"]>",2*tgerr_norm)
       # #     warning("Varianz zu groß. Prädiktions-Intervallbreite: ", 2*qnorm(1-0.05/2)*sqrt(ssw(dat)),">",2*tgerr)
      } else is2high <- 0
    } else {warning("pred.sol: there is NA in tgerr_norm in prediction function"); is2high<-NA}
    PI <- rbind(PI, cbind(PI.mima,jjj,is2high))
  }
  isprediction <- list()
  # choose points...:
  if (sequential==T) {
    tt <- na.omit(newx1d.seq[[length(pjjj)]])
    newx1d.seq[[length(pjjj)]] <- tt
  }
  for (jjj in pcnr) {
    # ((length(newx1d.PC12[[jjj]])==0 &sequential==F)|(sequential==T)) {   # no solution, choose point
                                                # &jjj==pcnr[1] # error?: & length(newx1d.seq[[length(pjjj)]])==0))
      if (sequential==T) {if (!nrow(newx1d.seq[[length(pjjj)]])==0) {isprediction[[jjj]] <- TRUE;# message(paste("!:",newx1d.seq[[length(pjjj)]]));message(paste("ohoh:",str(newx1d.seq[[length(pjjj)]])));
      break
                          } else message(paste("info:",length(newx1d.seq[[length(pjjj)]]), "; " , nrow(newx1d.seq[[length(pjjj)]])))
        }
      if (sequential==F) {
        if (!length(newx1d.PC12[[jjj]])==0) {isprediction[[jjj]] <- TRUE; next}
      }
    message(paste0("step ",stepi," pcnr", jjj,": no roots for chosen model. measure new point in unexplored space."))

       if (ptchng==T) {
         ## only in first iteration
         if (jjj==pcnr[1]) {
           jind  <- (nptc %% plsXdim)+1
           print(paste("ptchoose",jjj,"pc:", jind, "nptc:", nptc ))
           #how to consider pcnr?
           newx1d.PC12[[jind]] <- ptchoose(ptchoice=ptchoice, dat1d=dat1d.PC12[[jind]],
                                          tgmean_norm=tgmean_norm, maxarea=maxarea,
                                          xmeans=xmeans, xsds=xsds, pls1=pls1, jjj=jind)
           isprediction[[jind]] <- FALSE
           #nptc <- nptc + 1;
           rm(jind)
         }
       } else {
         isprediction[[jjj]] <- FALSE
         newx1d.PC12[[jjj]] <- ptchoose(ptchoice=ptchoice, dat1d=dat1d.PC12[[jjj]],
                                        tgmean_norm=tgmean_norm, maxarea=maxarea,
                                        xmeans=xmeans, xsds=xsds, pls1=pls1, jjj=jjj) # choose in both/all directions
       }
    ourshift <- rep(NA, plsXdim)
    if( length(newx1d.PC12[[jjj]])>1) stop("unhandled situation; pred.sol")
    ourshift[jjj] <- newx1d.PC12[[jjj]]
    newx1d.ext <- rbind(newx1d.ext,
                        pred.ptch(shift=ourshift, pcnr=pcnr, pls1=pls1, mknormweights=mknormweights, dat1d.PC12=dat1d.PC12, mindeg=mindeg, tgmean_norm=tgmean_norm, gr2retlimit=gr2retlimit, wfun=wfun, sto=sto) )
  }

  improvePT <- T #!!
  debugtt=F
  if (improvePT==T) { if(sequential==T) {
   # print(length(pjjj)) ; print(pjjj) ; print(newx1d.seq)
    if (nrow(newx1d.seq[[length(pjjj)]])==0) {
      warning("improvett")
      newx1d.seq[[length(pjjj)]] <- rbind(newx1d.seq[[length(pjjj)]],na.omit(newx1d.ext))
      if (nrow(newx1d.seq[[length(pjjj)]])==0) message("immernoch0")
    }
  }}

  # project 1D to xdim-D: newx1d -> newx
  xcols <- grep("x.", names(xmeans), value = TRUE)
  newx.PC12 <- list()
  #if (additive==TRUE) {
    # have to consider shifts below in back transformation
 # } else { # end if additive
  #   #plsXdim?
    tpcnr <- 1:plsXdim # pcnr # Xdim
    if (sequential==T) {
      if (dim(pls1$mod.wgs)[1]==dim(pls1$mod.wgs)[2]) {
        newx.seq <- mkreg(as.matrix(newx1d.seq[[length(pjjj)]])%*%solve(pls1$mod.wgs)[pcnr,],xmeans[xcols],xsds[xcols])
      } else newx.seq <- mkreg(as.matrix(newx1d.seq[[length(pjjj)]])%*%MASS::ginv(pls1$mod.wgs)[pcnr,],xmeans[xcols],xsds[xcols])
    }
    if (length(pjjj)>0) message(paste("nNewxsq1d",nrow(newx1d.seq[[length(pjjj)]]))) # message only if sequential or length of pjjj >0

    for (jjj in tpcnr) {
      if (length(newx1d.PC12)>=jjj) ##
      {
        if (dim(pls1$mod.wgs)[1]==dim(pls1$mod.wgs)[2]) {
        newx.PC12[[jjj]] <- mkreg(matrix(newx1d.PC12[[jjj]])%*%solve(pls1$mod.wgs)[jjj,],xmeans[xcols],xsds[xcols])
        } else newx.PC12[[jjj]] <- mkreg(matrix(newx1d.PC12[[jjj]])%*%MASS::ginv(pls1$mod.wgs)[jjj,],xmeans[xcols],xsds[xcols])
      } else newx.PC12[[jjj]] <- list(NULL) ##
    }
#  }# end if (additive=TRUE)
  ###if new points are not up to an epsilon identical: uniqP
  for (jjj in tpcnr) {
    if (is.null(unlist(newx.PC12[[jjj]]))) { ## for the rest of the loop ignore tpcnr == jjj
      tpcnr <- setdiff(tpcnr,jjj)
    }
    newx.PC12[[jjj]] <- uniqP(newx.PC12[[jjj]], xeps, dat$x)
  }
  if (sequential==T) newx.seq <- uniqP(newx.seq, xeps, dat$x)

  newx<-NULL
  is.prediction <- numeric(0)
  #message(paste("DATMEAN",colMeans(dat$x)))

  for (jjj in tpcnr) {
    newx <- rbind(newx,newx.PC12[[jjj]] )
    is.prediction <- c(is.prediction, rep(isprediction[[jjj]],nrow(newx.PC12[[jjj]])))
  }
  if (sequential==T) {
    newx.seq <- nclose2mean(newx.seq,colMeans(dat$x))
    newx <- rbind(newx,newx.seq )
    is.prediction <- c(is.prediction, rep(5,nrow(newx.seq)))
  }
  #if (debugtt==T) browser() #
  #if (debugtt==T) is.prediction="debug"
  return(list(newx=newx, dat1d.PC12=dat1d.PC12, newx1d.PC12=newx1d.PC12,  is.prediction = is.prediction, isprediction=isprediction, PI=PI, modeg=modeg))#cbind(PI, is.prediction, unlist(isprediction))))
}

#' plotexp
#'
#' Function which plots...
#'
#' @param dat data.frame, with dat$x matrix of proces parameter data
#' @param target numeric vector
#' @param tgerr numeric vector
#' @param prediction numeric
#' @param model function
#' @param limit0 boolean
#' @param reality function
#'
#' @return plot
#' @export=TRUE
plotexp <- function(dat, target, tgerr, prediction=NULL, model=NULL, limit0=FALSE, reality=NULL) {
  fm <- model

  #cf. code block further below
  preds<-data.frame(x=NULL,y=NULL) #for limit0...
  if (!is.null(prediction)) {
    # roots
    rootdta <- data.frame(x = prediction)
    rootdta$y <- predict(fm, newdata=rootdta)
    preds<-rootdta
  }

  if (!is.null(fm)) {
    #plot regression curve
    #prd <- data.frame(x= seq(1,7,length.out=250))   #data.frame(x = dat$xr)    # xr<-seq(1,7,length.out=n*rpn)
    prd <- data.frame(x= seq(min(0, dat$x,preds$x)*1.4,max(dat$x,preds$x)*1.4,length.out=250))
    #err <- predict(fm, newdata = prd, se.fit = TRUE)
    tmp <- predict.lm(fm, newdata = prd, interval="prediction")
    prd$lci <- tmp[,"lwr"]
    prd$fit <- tmp[,"fit"]
    prd$uci <- tmp[,"upr"]
    rm(tmp)
  }

  myplot <- ggplot() +
    theme_bw()

  #plot reality
  if (!is.null(reality)) {
    myplot <- myplot +
      geom_line(aes(x=x,y=y),data=data.frame(x=seq(1,7,length.out=250), y=reality(seq(1,7,length.out=250))), colour = "grey", size=1.1)
  }

  if (!is.null(fm)) {
    myplot <- myplot +
      geom_ribbon(data=prd, aes(ymin = lci, ymax = uci, x=x), fill="yellow", alpha=0.5, stat = "identity")
  }

  myplot <- myplot +
    geom_point(data = dat, aes(x = x, y = y), color="black", shape=3) +
    labs(x = "Process parameter (p)", y="Descriptor (d)") +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=target-tgerr, ymax=target+tgerr, alpha=0.5, fill="green") +
    geom_hline(yintercept = target)

  if (!is.null(fm)) {
    myplot <- myplot +
      geom_line(aes(x=x,y=fit) ,data=prd)
  }

  preds<-data.frame(x=NULL,y=NULL) #for limit0...
  if (!is.null(prediction)) {
    rootdta <- data.frame(x = prediction)
    rootdta$y <- predict(fm, newdata=rootdta)
    myplot <- myplot +
      geom_point(data = rootdta, aes(x = x, y = y), color="red")
    preds<-rootdta
  }

  if (limit0) {
    myplot <- myplot + coord_cartesian(ylim = c(min(0, dat$y,target)*1.05, max(dat$y,target)*1.05), xlim=c(min(dat$x,preds$x)*1.05, max(dat$x,preds$x)*1.05))
  }
  return(myplot)
}

#' newexp
#'
#' Creates new samples from model function. Used for simulations in the function autosolve as a helper function.
#'
#' @param n integer, number of repeated measurements
#' @param xpos numeric matrix,  coordinates at which the measurement takes place
#' @param foo function, model for the real relationship
#' @param sd , standard deviation
#'
#' @return a matrix containing the samples
#' @export=FALSE
#'
#' @examples
#' #not to be used
newexp <- function(n=25, xpos, foo, sd=0.001) {
  # n Messungen an dieser Stelle durchführen
  #xn<-rep(xpos,n)
  xn <- apply(xpos, 2, rep, times=n)
  if (is.vector(xn)) xn <- t(xn)
  #yn <- foo(xn) + rnorm(n*length(xpos),0,1)
  yn <- foo(xn)
  yn <- yn + rnorm(prod(dim(yn)),0,sd) #########!!!!!!!!!!!!!!!!! so wenig varianz
  return(data.frame(x=I(as.matrix(xn)),y=I(as.matrix(yn))))
}

#' newexp2
#'
#' Creates new samples from model function. Used for simulations in the function autosolve as a helper function.
#'
#' @param n integer, number of repeated measurements
#' @param xpos numeric matrix, coordinates at which the measurement takes place
#' @param foo function, model for the real relationship
#' @param sd numeric, standard deviation
#'
#' @return a matrix containing the samples
#' @export=FALSE
#'
#' @examples
#' #not to be used
newexp2 <- function (n = 25, xpos, foo, sd = 0.001)
{
  xn <- apply(xpos, 2, rep, times = n)
  if (is.vector(xn))
    xn <- t(xn)
  yn <- foo(xn)
  yn <- yn + rnorm(prod(dim(yn)), 0, sd)
  #rs <- list(x = as.matrix(xn), y = as.matrix(yn))
  #class(rs) <- "data.frame"
  rsdim <- dim(as.matrix(xn))[1]
  rs <- data.frame(x = numeric(rsdim), y=numeric(rsdim))
  rs$x <- as.matrix(xn)
  rownames(rs$x)<-rep(NA,dim(rs$x)[1])
  rs$y <- as.matrix(yn)
  return(rs)
}
#' autosolve
#'
#' Applies the iterative optimization algorithm suggested in this package.
#' To perform the simulations automatically, the true model has to be specified.
#' When in practice the true model is unknown, use \code{pred.solution} function instead to get a candidate for the unknown optimum in a single iteration.
#'
#' @param startx numeric matrix, start values
#' @param tgmean numeric v4ector, target value vector
#' @param tgerr numeric vector, defines acceptable error range
#' @param reps integer, number of repeated measurements
#' @param maxit integer, maximum number of iterations
#' @param reality function, real model
#' @param xeps nunmeric, smallest (reasonably) distinguishable epsilon
#' @param pplot boolean, do plots?, pplot=TRUE is deprecated
#' @param pcnr integer vector, defines which principal directions will be considered
#' @param maxarea numeric matrix, area range, which will be explored
#' @param useweights boolean
#' @param mknormweights boolean
#' @param gr2retlimit boolean
#' @param mindeg integer, minimal degree of order of polynomial model
#' @param sequential boolean
#' @param tgpcawg numeric
#' @param yweights boolean
#' @param datlim NULL or integer
#' @param knearest integer
#' @param tgdim integer
#' @param ylast integer
#' @param sto boolean
#' @param mod.sd numeric
#' @param ...
#'
#' @return data.frame
#' @export=TRUE
#'
#' @examples
#' # Erzeuge 3x3 Modell für Simulation
#' startx <- mydata
#' tfoo <- function(x) {
#'   x1 <- x[,1]
#'   x2 <- x[,2]
#'   # x3 <- x[,3]
#'   return( data.frame( y1=0.8*x1 - 1.2*x2,
#'                       y2=0.8*x1 - 1.2*abs(x2)^0.25#,
#'                      #  y3=0.8*x1 - 0.6*abs(x2)^0.25 + x3
#'  ) ) }
#' tgmean <- tfoo(cbind(1.35,1.4))#,1.5))#tfoo(cbind(0.35,0.4,0.5))
#' tgerr <- c(0.2,0.2)#,0.5)
#' #reps <- 5
#' #maxit <-10
#' #reality <- tfoo
#' #xeps <- 0.01
#' #pplot <- TRUE
#' tmp<-autosolve(startx,tgmean,tgerr*0.125,reps=7,maxit=10,tfoo, xeps=0.01, F, pcnr=c(1,2))
#' #nwexample
#' dstgmean <- tfoo(cbind(1.35,1.4))#,1.5))#tfoo(cbind(0.35,0.4,0.5))
#' tgerr <- c(0.2,0.2)#,0.5)
#' reps <- 5
#' maxit <-10
#' reality <- tfoo
#' xeps <- 0.01
#' pplot <- TRUE
#' startx <- expand.grid(x1=c(-1,3),x2=c(-1,3))
#' set.seed(123)
#' autosolve(startx,tgmean,tgerr*0.125,reps=7,maxit=10,tfoo, xeps, F, pcnr=c(1,2))
autosolve <- function(startx, tgmean, tgerr, reps=25, maxit=10, reality=foo, xeps=0.01, pplot=TRUE, pcnr=c(1,2),maxarea=NULL, useweights=TRUE, mknormweights=F,gr2retlimit=T, mindeg=0,sequential=F,tgpcawg=1,yweights=F,datlim=NULL, knearest=NULL,tgdim=1,ylast=NULL,sto=T,mod.sd=NULL,...) {
  tgmean<-matrix(tgmean,nrow=1)
  tgerr<- matrix(tgerr,nrow=1)
  # dimensions
  Xdim <- dim(startx)[2]
  Ydim <- dim(tgmean)[2]
  dat <- data.frame(x=numeric(),y=numeric(),nri=numeric())
  newx <- startx
  vartoobig <- -1;  solfound <- -1;  nomore <- -1
  modeg <- matrix(NA, ncol=2, nrow=0)
  for (i in 1:maxit){
    # if i>0 then change startx: getroots(polymodel(data3,degByBIC(data3, mindeg=mindeg)),tgmean)
    if (i==1) {
      PI <- data.frame(pi.l=numeric(), pi.r=numeric(), pc=integer(), nr=integer(), nri=integer())
      ispred <- F
      is.pred <- NULL
    }
    if (i>1) {
      flag <- TRUE
      tryCatch( { rs <- pred.solution(dat=dat,tgmean=tgmean,tgerr=tgerr,xeps=xeps,pcnr=pcnr,maxarea=maxarea, useweights=useweights, mknormweights=mknormweights,gr2retlimit=gr2retlimit,mindeg=mindeg,sequential=sequential,nptc=sum(!ispred),tgpcawg=tgpcawg,yweights=yweights,datlim=datlim,knearest=knearest,tgdim=1+((i-1-1) %% tgdim), ylast=ylast,sto=sto,...) }
                , error= function(e) { print("break out of for - pred.solution throws error:"); print(e); flag<<-FALSE}) #,additive=additive
      if (!flag) break
      #rs <- pred.solution(dat=dat,tgmean=tgmean,tgerr=tgerr,xeps=xeps,pcnr=pcnr,additive=additive,maxarea=maxarea, useweights=useweights, mknormweights=mknormweights,gr2retlimit=gr2retlimit,mindeg=mindeg,sequential=sequential,nptc=sum(!ispred),tgpcawg=tgpcawg,yweights=yweights,datlim=datlim,...)
      newx <- rs$newx
      # check for NAs
      if (anyNA(newx)){
        print(tail(dat))
        print(newx)
        print(tgmean)
        print(tgerr)
        print(rs)
        stop("NAs in datn")
      }
      dat1d.PC12 <- rs$dat1d.PC12
      newx1d.PC12 <- rs$newx1d.PC12
      rs$PI$nri<-i
      PI <- rbind(PI, rs$PI)
      ispred <- c(ispred,unlist(rs$isprediction))
      is.pred <- c(is.pred,rs$is.prediction)
      modeg <- rbind(modeg, cbind(rs$modeg, i))
      rm(rs)
    }

    if (length(newx)==0) {
      if (nomore==-1) {nomore<-i}
      message(paste0("step ",i,": no new points"))
      break
    }

    datn <- newexp2(n=reps,xpos=newx, foo=reality, sd=mod.sd)
    datn$nri <- i

    # check for NAs
    if (anyNA(datn$y[,"y1",drop=F])){
      print(datn)
      stop("datn has NAs")
    }
    rownames(datn$y) <- rep(NA,dim(datn$y)[1])
    newdat <- rbind(dat,datn) ## !dimnames of datn$y
    if (pplot==TRUE) {
      if (i>1) {
        #where appear the new points in the old coordinate system?
        xcols <- grep("x.", names(xmeans), value = TRUE)
        tmpscores <- mknorm(newdat$x,xmeans[xcols],xsds[xcols])%*%pls1$mod.wgs
        ##### y?, need y with old pca coordinates##############################
        tmpscoresy <- mknorm(newdat$y,myypca$allobsmean,myypca$allobssd) %*% myypca$pca$rotation
        #1d-y?
        newdat1d.PC12 <- list()
        for (jjj in pcnr) {
          newdat1d.PC12[[jjj]] <- data.frame(x=tmpscores[,jjj], y=tmpscoresy[,1]) #newdat1d <-   data.frame(x=tmpscores[,1], y=newdat$y[,"y1"])
          print(plotexp(dat=newdat1d.PC12[[jjj]], target=myypca$pcatg[,"PC1"], myypca$pcatgerr[,"PC1"], prediction=NULL, model=polymodel(dat1d.PC12[[jjj]], degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)), limit0=TRUE) )
        }
      }
      # DIMENSION REDUCTION (for the following plot)
      xmeans <- colMeans(newdat)
      xsds <- apply(newdat,2,sd)
      ### DIMENSION REDUCTION Y; pca for explained y
      if (Ydim>=2) {
        ypca <- tgpca(newdat$y,tgmean,tgerr,wg=tgpcawg)  #besser t(as.numeric(tgmean)), simplify code, cf below

        doCrosVal <- !(nrow(newdat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL
        tryCatch( { pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2], crosval=doCrosVal) }
                  , error= function(e) { print("plsreg1 throws error:"); print(e); stop("plsreg1 failed in line 1572")})
        #pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2]) #oder comps=2 genug hier?
      } else {
        if (anyNA(newdat$y[,"y1",drop=F])){
          print(newdat)
          stop("NAs in newdat$y[,'y1',drop=F]")
        }
        doCrosVal <- !(nrow(newdat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL
        tryCatch( { pls1 = plsdepot::plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2], crosval=doCrosVal) }
                  , error= function(e) { print("plsreg1 throws error:"); print(e); stop("plsreg1 failed in line 1582")})
        #pls1 = plsdepot::plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2])
      }
      newdat1d.PC12<-list()
      for (jjj in pcnr) {
        if (Ydim>=2) {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=ypca$pca$x[,"PC1"])
        } else {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=newdat$y[,"y1"])
        }
      }
      if (Ydim>=2) {
        tgmean_norm <- ypca$pcatg[,"PC1"]
        tgerr_norm <- ypca$pcatgerr[,"PC1"]
      } else {
        tgmean_norm <- mknorm(tgmean, xmeans[c("y")], xsds[c("y")]) # check
        tgerr_norm <-  mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y")]) # check
      }

      #plot: on origginal scale
      # plot(mkreg(pls1$x.scores[,1,drop=F]%*%solve(pls1$mod.wgs)[1,],xmeans[c("x.x1","x.x2")],xsds[c("x.x1","x.x2")]))
      # points(dat$x,col="red")
      # points(mkreg(cbind(pls1$x.scores[,1,drop=F],.21)%*%solve(pls1$mod.wgs)[1:2,],xmeans[c("x.x1","x.x2")],xsds[c("x.x1","x.x2")]), col="green")
      # points(mkreg(cbind(pls1$x.scores[,1,drop=F],-.21)%*%solve(pls1$mod.wgs)[1:2,],xmeans[c("x.x1","x.x2")],xsds[c("x.x1","x.x2")]), col="green")
      for (jjj in pcnr) {
        print(plotexp(dat=newdat1d.PC12[[jjj]],
                      target=tgmean_norm,#t(as.numeric(tgmean)),
                      tgerr_norm, #tgerr,
                      prediction=getroots2(newdat1d.PC12[[jjj]],
                                           degByBIC(newdat1d.PC12[[jjj]], mindeg=mindeg),
                                           tgmean_norm, retlimit = gr2retlimit),# t(as.numeric(tgmean))),
                      model=polymodel(newdat1d.PC12[[jjj]], degByBIC(newdat1d.PC12[[jjj]], mindeg=mindeg)), limit0=T ) )
      }
      #plot using plotly
      # plot_ly(x = newdat$x[,"x1"], y = newdat$x[,"x2"], z = newdat$y[,"y2"], marker = list(color = newdat$y[,"y1"], colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
      #     add_markers() %>%
      #     layout(scene = list(xaxis = list(title = 'Prozessparameter 1'),
      #                        yaxis = list(title = 'Prozessparameter 2'),
      #                        zaxis = list(title = 'Deskriptor 1')),
      #            annotations = list(
      #              x = 1.13,
      #              y = 1.05,
      #              text = 'Wert von Deskriptor bzw. Abstand zu Zielwert',
      #              xref = 'paper',
      #              yref = 'paper',
      #              showarrow = FALSE
      #  ))
      #rgl::plot3d(cbind(newdat$x, newdat$y[,"y1"]), col="red", size=3)
      #rgl::plot3d(cbind(newdat$x, newdat$y[,"y2"]), col="red", size=3)
    }

    if (i>1) {
      #success?
      # DIMENSION REDUCTION
      xmeans <- colMeans(newdat)
      xsds <- apply(newdat,2,sd)
      if (Ydim>=2) {
        ypca <- tgpca(newdat$y,tgmean,tgerr,wg=tgpcawg)
        doCrosVal <- !(nrow(dat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL - contact Gaston Sanchez for Bug report?
        pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2], crosval=doCrosVal) #oder comps=2? genug? #pls1 = plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2])
        #pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2])
      } else {
        if (anyNA(newdat$y[,"y1",drop=F])){
          print(tail(newdat,200))
          print(xmeans)
          print(xsds)
          stop("NAs in anyNA(newdat$y[,'y1',drop=F]")
        }

        doCrosVal <- !(nrow(newdat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL
        tryCatch( { pls1 = plsdepot::plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2], crosval=doCrosVal) }
                  , error= function(e) { print("plsreg1 throws error:"); print(e); stop("plsreg1 failed near line 1546")})
        #pls1 = plsdepot::plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2])
      }
            #stopping critereon
      newdat1d.PC12 <- list()
      for (jjj in pcnr) {
        if (Ydim>=2) {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=ypca$pca$x[,"PC1"]) #data.frame(x=pls1$x.scores[,"t1"], y=newdat$y[,"y1"])
        } else {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=newdat$y[,"y1"]) #data.frame(x=pls1$x.scores[,"t1"], y=newdat$y[,"y1"])
        }
      }
      if (Ydim>=2) {
        tgmean_norm <- ypca$pcatg # mknorm(tgmean, xmeans[c("y")], xsds[c("y")])#mknorm(tgmean, xmeans[c("y.y1","y.y2")], xsds[c("y.y1","y.y2")])
        tgerr_norm <- ypca$pcatgerr # mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y")])#mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y.y1","y.y2")])
      } else {
        tgmean_norm <- mknorm(tgmean, xmeans[c("y")], xsds[c("y")]) # check
        tgerr_norm <-  mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y")]) # check
      }

      for (jjj in pcnr) {
        if (is.null(newx1d.PC12)) break # ok???
        if (jjj>length(newx1d.PC12)) break # ok???
        if (is.null(newx1d.PC12[[jjj]])) next # ok???
        if (degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)==degByBIC(newdat1d.PC12[[jjj]], mindeg=mindeg)) {
          checkpred <- as.data.frame(predict.lm(polymodel(newdat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)), newdata = data.frame(x=newx1d.PC12[[jjj]]), interval="prediction") ) #,row.names=1:(dim(dat)[1]+dim(datn)[1])
          if (any(checkpred$lwr[checkpred$upr<tgmean_norm[1]+tgerr_norm[1]] > tgmean_norm[1]-tgerr_norm[1])) {
            if (solfound==-1) {solfound<-i}
            message(paste0("Solution according to model for PC",jjj))

            dat <- newdat # or return newdat below
            break
          }
          rm(checkpred)
        }
      }
    }
    dat <- newdat
  }
  erfolgwk0 <- suppressWarnings( max(pnorm(tgmean[[1]]+tgerr[[1]], mean=reality(newx)[,1], sd=1) - pnorm(tgmean[[1]]-tgerr[[1]], mean=reality(newx)[,1], sd=1)) )
  return(list(data=dat, data1d=newdat1d.PC12, #roots=roots,
              solfound=solfound, nomorepoints=nomore, variancetoobig=vartoobig, messreihen=max(dat$nri), messungen=dim(dat)[1], erfolgwk=erfolgwk0, PI=PI,  is.pred = is.pred, ispred=ispred, modeg=modeg#, pls1st=pls1st
  ))
}

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
  if (sum(xsds==0)>0) {
    warning("zero variances in mknorm")
    xsds[xsds==0] <- 1
  }
 # str(xmeans);print(xmeans);
  myx <- sweep(myx,2,as.numeric(xmeans),"-")
  myx <- sweep(myx,2,as.numeric(xsds),"/")
  return(myx)
}

mkreg <- function(x, means, sds) {
  myx <- as.matrix(x)
  myx <- sweep(myx,2,sds,"*")
  myx <- sweep(myx,2,means,"+")
  return(myx)
}

mvdistance <- function(xs, y, euclid=F) { # summennorm, manhattan-metrik, x--werte als zeilen in xs
  xs<-as.matrix(xs)
  rs <- sweep(xs,2,y,"-")
  rs <- abs(rs)
  if (euclid==T) rs<- rs^2
  rs<- rowSums(rs)
  if (euclid==T) rs <- sqrt(rs)
  return(rs)
}

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
  if (is.null(currBIC)) stop("fatal error in degbybic - currbic is null")
  if (which.min(currBIC)-1+mindeg<mindeg) stop("mindeg: the horror")
  return(which.min(currBIC)-1+mindeg)  #return(nnet::which.is.max(-currBIC)-1)#?????????????????????
}
# degByBIC <- function(dat, maxorder=5, weights=NULL) {
#   currBIC <- numeric()
#   for (i in 1:min(maxorder,length(unique(dat$x))-1)) {
#     currBIC <- c(currBIC, BIC(polymodel(dat,i, weights=weights)) )
#     if (i>1) if (currBIC[length(currBIC)] > currBIC[length(currBIC)-1]) break
#   }
#   return(nnet::which.is.max(-currBIC))
# }

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
# polymodel <- function(dat, degree=1, weights=NULL) {
#   if (degree<1) {
#     warning("polymodel < 1")
#   }
#   if (degree>length(unique(dat$x))-1) {
#     stop("polymodel: degree > unique")
#   }
#   return( lm(y ~ poly(x, degree, raw=FALSE), data=dat, weights=weights) )
# }

polymodel2 <- function(dat, degree, weights=NULL) {
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
# polymodel2 <- function(dat, degree, weights=NULL) {
#   if (degree<1) {
#     warning("polymodel2 < 1")
#   }
#   if (degree>length(unique(dat$x))-1) {
#     stop("polymodel2: degree > unique")
#   }
#   return( lm(y ~ poly(x, degree, raw=TRUE), data=dat, weights=weights) )
# }

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
# getroots2 <- function(dat, degree, target, limit=40, weights=NULL, retlimit=FALSE) {
#   if (!is.null(weights)) if (sum(is.na(weights))>0) stop("NA values in weights for getroots2()")
#   #getroots <- function(dat,target) {
#   # wie plotlm, nur anders
#   #linear model
#   fm <- polymodel2(dat, degree, weights=weights)
#
#   #find roots
#   coeff <- coefficients(fm)
#   coeff[1]<-coeff[1]-target
#   roots <- polyroot( coeff )
#   roots <- Re(roots)[abs(Im(roots)) < 1]#1e-1]#1e-6]   ############ACHTUNG, nur experimentel auskommentiert
#   if (retlimit==TRUE) {
#     roots[roots < -abs(limit)] <- -abs(limit)
#     roots[roots > abs(limit)] <- abs(limit)
#   } else {
#     roots <- roots[abs(roots)<=limit] }
#   #    if (dim(roots[abs(roots)<=limit])[1]>0)
#   # if (dim(roots)[1]<1)
#   #  roots<-sign(roots)*limit*(0.5+0.5*runif(1))
#   rootdta <- data.frame(x = roots)
#   # rootdta$y <- predict(fm, newdata=rootdta)
#   return(rootdta$x)
# }

PIhcheck <- function(model, alph=0.05) {
  if (sum(summary(model)$residuals^2)==0) {
    warning("perfect fit in PIhcheck")
    return(c(0,0))
  } else {
    return(c(2*qnorm(1-alph/2)*sqrt(sum(summary(model)$residuals^2)/qchisq(1-alph/2,df=model$df)), 2*qnorm(1-alph/2)*sqrt(sum(summary(model)$residuals^2)/qchisq(alph/2,df=model$df)) ) )
  }
}

# helper function to write a more exact euclw function
sumtheothers <- function(x) {
  n<-length(x)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(x[-i])
  }
  return(result)
}
sumthelowers <- function(x) {
  n<-length(x)
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(x[0:(i-1)])
  }
  return(result)
}
euclw <- function(x, normalize=T, sto=T) {
  # takes a matrix/data.frame with each coloumn having the scores corresp to  a principal component
  # computes weights as needed in 1d regression models:
  # i.e. computes for each of the 1d-projections on the PC the eucl distances from the n-dimensional point in the n-space
 # t2 <- rowSums( x^2 )
 # t1 <- (x^2)*(-1)
 # rs <- sqrt(sweep(t1,1,t2,"+"))
  if(normalize==2) {
    x <- apply(x,2,mknorm)
  }
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
euclw.old <- function(x, normalize=T) {
  # takes a matrix/data.frame with each coloumn having the scores corresp to  a principal component
  # computes weights as needed in 1d regression models:
  # i.e. computes for each of the 1d-projections on the PC the eucl distances from the n-dimensional point in the n-space
  #has problems wih precision due to rowsums - x^2
  t2 <- rowSums( x^2 )
  t1 <- (x^2)*(-1)
  rs <- sqrt(sweep(t1,1,t2,"+"))
  if(normalize==T) apply(rs,2,mknorm) else rs
}

#' PCA
#'
#' used for dimension reduction of the target space
#' in order to choose a certain direction through the target point for the projections wg has to be set to 1, then the target point is chosen as center for the pca.
#' If wg between 0 and 1 then pseudo observations at the target point are created such that a ratio of wg of the observations are pseudo observations.
#' uses prcomp on standardized data and pseudo data.
#'
#' @param dat data set
#' @param tgmean optional
#' @param tgerr optional
#' @param wg weight of the target value. If wg equals 1 or 2 then the pca is performed with the target value as center
#'
#' @return returns the results of the pca and some extra stuff
#' @export=FALSE
#'
#' @examples tgpca(matrix(rnorm(20),ncol=2))
tgpca <- function(dat, tgmean=NULL, tgerr=NULL, wg=0.901, wfun, mknormweights, yweights=F, ylast=NULL) {   #unelegant (2x): t(as.numeric(tgmean))
  if (dim(dat)[2]!=length(tgmean)) stop(paste0("Dimension error in tgpca: ",dim(dat)[2], " and ",length(tgmean)))
  if (!is.null(tgerr)) {
    if (dim(dat)[2]!=length(tgerr)) stop(paste0("tgerr-Dimension error in tgpca: ",dim(dat)[2], " and ",length(tgmean)))
  }

  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }


  dat.orig<-dat
  nrw <- nrow(dat)
  if (!is.null(ylast)) dat <- dat[max(1,(nrw-ylast+1)):nrw,]

  # Anteil wg an Pseudobeobachtungen hinzufügen
  n1 <- dim(dat)[1]
  n2 <- round(wg*n1/(1-wg))
  if(is.infinite(n2)|(n2<0)) {
    n2<-0
    #wg<-1
  }

  #yweights = T
  if (yweights==T) { if (is.null(tgmean)) stop("tgmean null in tgpca yw")
    #yw <- ftweights(  wfun(euclw(sweep(allobs[1:n1,],2,as.numeric(tgmean)), mknormweights))  )
    yw <- ftweights( wfun(sqrt(rowSums(sweep(dat,2,as.numeric(tgmean))^2))) )
    # print(paste("sqrtyw",sqrt(yw)))
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
    if (is.null(tgmean) | anyNA(tgmean)) stop("tgmean nicht korrekt gegeben in tgpca")
  #  print(paste("aobsmean is:", aobs.mean))
  #  print(paste("tgmean is:", tgmean))
   # str(tgmean); str(aobs.mean)
    #print("allobs"); print(allobs)
    if (sum(aobs.sd==0)>0) {
      warning("zero variances in tgpca/if(wg==1)")
      aobs.sd[aobs.sd==0] <- 1
    }
    rs$pca <- prcomp(mknorm(allobs), center=(as.numeric(tgmean)-aobs.mean)/aobs.sd , scale=F) # wenn tgmean nicht null
  } else if (wg==2) { #variante von oben, sollte vom ergebnis gleich sein
    if (is.null(tgmean) | anyNA(tgmean)) stop("tgmean nicht korrekt gegeben in tgpca")
  #  print("allobs"); print(allobs); print("normiert:"); print(mknorm(allobs, means=tgmean))
  #  print("tgmean"); print(tgmean); str(tgmean)
    rs$pca <- prcomp(mknorm(allobs, means=tgmean), center=F , scale=F)
  } else {
    rs$pca <- prcomp(mknorm(allobs))
  }
  # pseudos aus x entfernen
  if (wg<1&0<wg) rs$pca$x <- rs$pca$x[-(1:n2),]
  if(!is.null(tgmean)) {
    if (wg>=1) { # hinweis: bei wg>=1 sollte pcatg==0 gelten
      rs$pcatg <- mknorm(t(as.numeric(tgmean)), tgmean, apply(allobs,2,sd)) %*% rs$pca$rotation
    # } else if (wg==2) {
    #   rs$pcatg <- mknorm(t(as.numeric(tgmean)), tgmean, apply(allobs,2,sd)) %*% rs$pca$rotation
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
    rs$pcatgerr <- namedNA #rep(NA,dim(dat)[2]) # wegen: myypca$pcatgerr[, "PC1"]
  }
  # die passiven beobachtungen wieder hinzufügen
  if (!is.null(ylast)) {
    if (wg>=1) {
      rs$pca$x <- mknorm(dat.orig, tgmean, apply(allobs,2,sd)) %*% rs$pca$rotation
    } else {
      rs$pca$x <- mknorm(dat.orig, colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
    }
  }
  rs$allobsmean <- colMeans(allobs)  # kann man auch oben speichern und redundanz vermeiden
  rs$allobssd <- apply(allobs,2,sd)   # ebenso
  return(rs)
}
tgpcabeforeweight1 <- function(dat, tgmean=NULL, tgerr=NULL, wg=0.901) {   #unelegant (2x): t(as.numeric(tgmean))
  if (dim(dat)[2]!=length(tgmean)) stop(paste0("Dimension error in tgpca: ",dim(dat)[2], " and ",length(tgmean)))
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  # Anteil wg an Pseudobeobachtungen hinzufügen
  n2 <- round(wg*dim(dat)[1]/(1-wg))
  if (!is.null(tgmean) & !anyNA(tgmean)  & wg != 0) {
    pseudobs <- rep.row(as.numeric(tgmean),n2) #apply(t(as.numeric(tgmean)),2,rep,n2)
    if (!is.null(colnames(dat))) colnames(pseudobs)<-colnames(dat)
    allobs <- rbind(pseudobs, dat)
  } else allobs <- dat

  rs<-NULL
  rs$pca <- prcomp(mknorm(allobs))
  # pseudos aus x entfernen
  rs$pca$x <- rs$pca$x[-(1:n2),]
  #rs$pcatg <- ifelse(is.null(tgmean),NA,mknorm(t(as.numeric(tgmean)), colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation)
  if(!is.null(tgmean)) {
    rs$pcatg<-mknorm(t(as.numeric(tgmean)), colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
  } else rs$pcatg <- rep(NA,dim(dat)[2])
  #rs$pcatgerr <- ifelse(is.null(tgerr)|anyNA(tgerr),rep(NA,dim(dat)[2]),mknorm(t(as.numeric(tgerr)), 0*colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation)
  if(!is.null(tgerr)&!anyNA(tgerr)) {
    rs$pcatgerr <- mknorm(t(as.numeric(tgerr)), 0*colMeans(allobs), apply(allobs,2,sd)) %*% rs$pca$rotation
  } else {
    namedNA <-t(rep(NA,dim(dat)[2]))
    dimnames(namedNA)[[2]] <- paste0("PC",1:dim(dat)[2])
    rs$pcatgerr <- namedNA #rep(NA,dim(dat)[2]) # wegen: myypca$pcatgerr[, "PC1"]
  }
  rs$allobsmean <- colMeans(allobs)
  rs$allobssd <- apply(allobs,2,sd)
  return(rs)
}

# getminmaxscores <- function(maxarea, mwgs) {
#   vals <- list()
#   for (i in 1:dim(mwgs)[1]) { # for every parameter p[i]
#     vals[[i]] <- t(maxarea)[,i,drop=F]%*%mwgs[i,,drop=F]
#   }
#   mins <- list()
#   for (i in 1:dim(mwgs)[1]) { # for every parameter p[i]
#     mins[[i]] <- apply(vals[[i]],2,min)
#   }
#   mins <- matrix(unlist(mins),byrow=T, ncol=dim(mwgs)[1])
#   mmin <- apply(mins,2,max)
#   maxs <- list()
#   for (i in 1:dim(mwgs)[1]) { # for every parameter p[i]
#     maxs[[i]] <- apply(vals[[i]],2,max)
#   }
#   maxs <- matrix(unlist(maxs),byrow=T, ncol=dim(mwgs)[1])
#   mmax <- apply(maxs,2,min)
#   rs <- t(rbind(mmin,mmax))
#   dimnames(rs)[1] <- dimnames(mwgs)[2]
#   return(rs)
# }
mimascores <- function(maxarea, mwgs) {
  if (sum(apply(maxarea,1,diff)<0)>0) stop("mimascores: maximle Parameterwerte müssen größer als minimale Parameterwerte sein")
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
mimascores2 <- function(maxarea, mwgs) {
  if (sum(apply(maxarea,1,diff)<0)>0) stop("mimascores2: maximle Parameterwerte müssen größer als minimale Parameterwerte sein")
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
ptchoose <- function(ptchoice=ptchoice, dat1d, tgmean_norm, maxarea, xmeans, xsds, pls1, jjj) {# brauche jjj nur für mimascores in pt3
  # dat1d=dat1d.PC12[[jjj]] # return newx1d.PC12[[jjj]]=newx1d
  #        # hier neuen punkt in richtung wählen - bräuchten alle projektionen auf hauptkomponente um größte lücke zu finden
  #        cat("1d points avail.:",paste(" ",sort(c(0, unique(dat1d[[1]])))))
  #        cat("1d points avail.:",paste(" ",diff(sort(c(0, unique(dat1d[[1]]))))))
  #        #alle neuen möglichen messpunkte:
  #        #cat("meas at:",paste(" ", head( sort(c(0, unique(dat1d[[1]]))) ,-1) + diff(sort(c(0, unique(dat1d[[1]]))))*0.5))
  #        cat("meas at:",paste(" ", head( sort(c(0, unique(dat1d[[1]]))) ,-1) + diff(sort(c(0, unique(dat1d[[1]]))))*0.5))
  #        cat("ext points",
  #            min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2,
  #            " and ",
  #            max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2)
  if (ptchoice==1) { # punkt aussen suchen, macht nur begrenzt sinn, wenn Grenzen für Parameter existieren
    #message(paste0("Schritt ","i",": gewähltes Modell findet keine Nullstellen in Richtung PC ",jjj,". Messe an neuem Punkt außerhalb des bisherig betrachteten Bereichs."))
    if (min(abs(dat1d$y[dat1d$x==min(dat1d$x)]-tgmean_norm)) < min(abs(dat1d$y[dat1d$x==max(dat1d$x)]-tgmean_norm))) {
      newx1d <- min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2
    } else {
      newx1d <- max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    }
  } else if (ptchoice==2) { # die 0 macht innerhalb von sort wenig sinn, oder?
    newx1d <- c((head( sort(c(0, unique(dat1d[[1]]))) ,-1) + diff(sort(c(0, unique(dat1d[[1]]))))*0.5)[which.max(diff(sort(c(0, unique(dat1d[[1]])))))],
                min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2,
                max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    )
  } else if (ptchoice==20) { # wie oben, aber ohne 0,  d.h. suche innen (1 Punkt) und am Rand (2Punkte)
    newx1d <- c((head( sort(c(unique(dat1d[[1]]))) ,-1) + diff(sort(c(unique(dat1d[[1]]))))*0.5)[which.max(diff(sort(c(unique(dat1d[[1]])))))],
                min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2,
                max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    )
  } else if (ptchoice==200) { # aber innen (1 Punkt) und am Rand (1Punkt)
    if (min(abs(dat1d$y[dat1d$x==min(dat1d$x)]-tgmean_norm)) < min(abs(dat1d$y[dat1d$x==max(dat1d$x)]-tgmean_norm))) {
      tmp <- min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2
    } else {
      tmp <- max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
    }
    newx1d <- c((head( sort(c(unique(dat1d[[1]]))) ,-1) + diff(sort(c(unique(dat1d[[1]]))))*0.5)[which.max(diff(sort(c(unique(dat1d[[1]])))))],
                tmp)
  } else if (ptchoice==3) {# nur 1 punkt für max abstand
    if (!is.null(maxarea)) { # wenn maxarea gegeben...
      xcols <- grep("x.", names(xmeans), value = TRUE)
      mima <- mimascores2(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols])), pls1$mod.wgs)[jjj,]
      upts <- unique(dat1d[[1]])  # wofür ist die [[1]]?
      upts <- upts[which(upts>mima[1])] # im intervall?
      upts <- upts[which(upts<mima[2])]
      upts <- sort(c(upts, mima)) # mit mima-Grenzen und sortiert
      measureatpts <- head(upts ,-1) + diff(upts)*0.5
      # wenn rand gegeben, dann wähle als abstand  2*(kleinstervorkommenderxwert-minpossiblexval) bzw
      #                                           2*(groestervorkommenderxwert-maxpossiblexval)
      if (length(upts)<2) stop("impossible: length(upts)<2 für ptchoice3")
      if (length(upts)<3) {
        warning("unhandled: length(upts)<3 für ptchoice3")
        choosept <- 1 # mal probieren, müsste passen
      } else choosept <- which.max(diff(upts)* c(2,rep(1,length(upts)-3),2)  )
      newx1d <- measureatpts[choosept]

      if (anyNA(newx1d)){ # nur im Fehlerfall...
        xcols <- grep("x.", names(xmeans), value = TRUE)
        print(maxarea)
        #print(xcols)
        #print(xmeans[xcols])
        #print(xsds[xcols])
        print(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))
        print(mima)
        print(newx1d)
        ttt<-unique(dat1d[[1]])
        # print(ttt)
        ttt <- ttt[which(ttt>mima[1])] # im intervall?
        #  print(ttt)
        ttt <- ttt[which(ttt<mima[2])]
        print(ttt)
        print(upts)
        # print(measureatpts)
        print(choosept)
        #print(plot(pls1$x.scores%*%solve(pls1$mod.wgs)))
        plot(pls1$x.scores%*%solve(pls1$mod.wgs),xlim=c(-9,9),ylim=c(-9,9))
        points(pls1$x.scores[,1,drop=F]%*%solve(pls1$mod.wgs)[1,,drop=F],col="green")
        points(matrix(mima)%*%solve(pls1$mod.wgs)[1,,drop=F],col="red")
        plot(mkreg(pls1$x.scores%*%solve(pls1$mod.wgs),xmeans[xcols],xsds[xcols]),xlim=c(-30,30),ylim=c(-30,30))
        points(mkreg(pls1$x.scores[,1,drop=F]%*%solve(pls1$mod.wgs)[1,,drop=F],xmeans[xcols],xsds[xcols]),col="green")
        points(mkreg(matrix(mima)%*%solve(pls1$mod.wgs)[1,,drop=F],xmeans[xcols],xsds[xcols]),col="red")
        rm(xcols)
        stop("newx1d hat NA in prediction")
      }
      warning(paste("ptchoose: chosenpt in dir",jjj)); print(newx1d) #; message("mimawas:"); print(mima)
    #   message("xmeans:") ; print(xmeans[xcols])
    #   message("xsds:") ; print(xsds[xcols])
    #   message("maxarea norm:") ; print(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols])))
    # print( c(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))[1,1],0)  )
    # print(  pls1$mod.wgs  )
    # print(  t(c(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))[1,1],0))%*% pls1$mod.wgs    )
    # print(  pls1$mod.wgs%*%matrix(c(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))[1,1],0))     )
    # print(  pls1$mod.wgs%*%matrix(c(t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))[1,2],0))     ) ### wgs*pt oder pt*wgs? überprüfen!
    # print(  pls1$mod.wgs%*%matrix(c(0,t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))[2,1]))     )
    # print(  pls1$mod.wgs%*%matrix(c(0,t(mknorm(t(maxarea),xmeans[xcols],xsds[xcols]))[2,2]))     )
    } else { # wenn keine grenzen dann wie ptchoice1
      if (min(abs(dat1d$y[dat1d$x==min(dat1d$x)]-tgmean_norm)) < min(abs(dat1d$y[dat1d$x==max(dat1d$x)]-tgmean_norm))) {
        newx1d <- min(dat1d$x) - (max(dat1d$x) - min(dat1d$x))/2
      } else {
        newx1d <- max(dat1d$x) + (max(dat1d$x) - min(dat1d$x))/2
      }
    }
  } # end: if (ptchoice==3)
  return(newx1d)
}
ftweights <- function(w1) {
  if (0<length(w1)) {
    if (0<sum(is.infinite(w1))) warning("ftweights: infinite weights were clipped")
    w1[is.infinite(w1)]<-1e+306    #1.797693e+308
    w1[w1>1e+306 ]<-1e+306
    return(w1)
  } else w1
}
uniqP <- function(newx, xeps, datx) {
  # neue Punkte um mehr als xeps verschieden?

  # Obere Variante war:
  #   dstncs.PC12 <- list()
  #   for (jjj in tpcnr) {   ### weshalb unten unlist?? hä??
  #     if (is.null(unlist(newx.PC12[[jjj]]))) { ## für den rest der schleife tpcnr == jjj ignorieren
  #       tpcnr <- setdiff(tpcnr,jjj)
  #     } else {     ###nur bestimmen, wenn neue Punkte um mehr als xeps verschieden
  #       dstncs.PC12[[jjj]] <- matrix(nrow = dim(newx.PC12[[jjj]])[1],ncol=dim(dat$x)[1])
  #       if (dim(newx.PC12[[jjj]])[1]>0)   #nur, wenn neu Punkte in PCjjj:
  #         for (cnt in 1:dim(newx.PC12[[jjj]])[1]) {
  #           dstncs.PC12[[jjj]][cnt,] <- mvdistance(dat$x, newx.PC12[[jjj]][cnt,])
  #         }
  #     }
  #   }
  #   for (jjj in tpcnr) {                                         ## warum min? weil abstand zum nächstgelegenen punkt
  #     newx.PC12[[jjj]] <- newx.PC12[[jjj]][ suppressWarnings(apply(dstncs.PC12[[jjj]],1,min))>xeps ,,drop=F]
  #   }

  # Untere Variante war:
  # if (!is.null(newx.seq)) {
  #   #nur bestimmen, wenn neue Punkte um mehr als xeps verschieden
  #   dstncs.seq <- matrix(nrow = dim(newx.seq)[1],ncol=dim(dat$x)[1])
  #   if (dim(newx.seq)[1]>0)   #nur, wenn neu Punkte:
  #     for (cnt in 1:dim(newx.seq)[1]) {
  #       dstncs.seq[cnt,] <- mvdistance(dat$x, newx.seq[cnt,])
  #     }
  #   newx.seq <- newx.seq[ suppressWarnings(apply(dstncs.seq,1,min))>xeps,,drop=F]
  # }
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
repna <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}
nclose2mean <- function(datx, center, n=1) {
  # if (!dim(datx)[2]=dim(center)[2]) stop("wrong dimensions in nclose2mean")
   nwas <- nrow(datx)
  if (!is.null(datx)) {
    dstncs <- numeric()
    if (dim(datx)[1]>0)   #nur, wenn neue Punkte vorhanden:
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

pred.ptch <- function(shift, pcnr, pls1, mknormweights, dat1d.PC12, mindeg, tgmean_norm, gr2retlimit, wfun, sto) { #pcnr, die noch zu bearbeiten sind
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

#' PLSolve.1step prediction / suggestion for new measurements
#'
#' try to solve equation sytem iteratively using PLS regression.
#' PLSolve.1step returns points for new measurement/evaluation of the model equation
#' PLSolve.auto iterative algorithm which tries to autoatically solve the equation. For this it is necessary to specify the true model
#' Suggest points for future measurements in order to find a solution that reaches the desired target value iteratively
#' Eventuell weiter aufspilitten: Varianzabbruch und Delta-Kriterium
#'
#' @param dat data set
#' @param tgmean target mean
#' @param tgerr maximal target error - should be optional
#' @param xeps smallest delta
#' @param pplot do you want to plot?
#' @param pcnr principal components considered
#' @param additive dont use this
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
pred.solution <- function(dat,tgmean,tgerr=NULL,xeps=0.001,pcnr,additive=FALSE, maxarea=NULL, ptchoice=1, useweights=TRUE, mknormweights=F, allpts=F, gr2retlimit=TRUE,bpcenter=F,mindeg=0,wfun=function(x) {(1/x)^2}, sequential=F, ptchng=F, nptc=0,tgpcawg=1,betterweights=F,yweights=F,datlim=NULL,knearest=NULL,tgdim=1,ylast=NULL,sto=T,...) {
  # nptc: how often has ptchoice been performed
  # Datenstrukturen für tgmean, tgerr, maxarea handlen
  tgmean<-matrix(tgmean,nrow=1)
  if (!is.null(tgerr)) {
    tgerr <- matrix(tgerr,nrow=1)
  } else tgerr <- matrix(NA,nrow=1,ncol=length(tgmean)) #tgerr<- ifelse(!is.null(tgerr), matrix(tgerr,nrow=1), matrix(NA,nrow=1,ncol=length(tgmean)))

  stepi <- 1 + max(dat$nri)
 # print(tail(dat, 5))
  #if (!is.null(datlim)) dat <- tail(dat, datlim)  # dat[1:datlim,]
  if (!is.null(datlim)) {
    dat <- dat[dat$nri >= (max(dat$nri)-datlim) ,] #tail(dat, datlim)
  } else if (!is.null(knearest)&bpcenter!=T) { # k nächse zu letzter beobachtung # !=T, d.h. wird ausgeführt für bpcenter=2
    # alternative bei bpcenter ein limit einfügen viell besser  wfun(pmax(euclw(sweep(pls1$x.scores,2,minpoint), mknormweights),0.0000001))
    # letzte beobachtung row.sum  ftweights(wfun(sqrt(rowSums(myypca$pca$x^2))))
    ndat <- mknorm(dat$x)
    lobs <- ndat[nrow(ndat),] # ftweights(wfun(sqrt(rowSums(myypca$pca$x^2))))
    idcs <- sort(mvdistance(ndat, lobs, TRUE), index.return=T)$ix[1:min(knearest,nrow(ndat))]  # abstände zu letzter beobachtung#sortieren und auswählen
    #wfun(pmax(euclw(sweep(ndat,2,lobs), mknormweights),0.0000001))
    dat <- dat[idcs ,]
    #print("NotNow")
  }

  if (is.matrix(maxarea)) {
    #print(dat)
    xdim <- dim(dat$x)[2]#dim(dat)[2]-dim(tgmean)[2]
    if (dim(maxarea)[1] != xdim |  dim(maxarea)[2] != 2)
      stop("wrong dimensions of maxarea in pred.solution")
    rm(xdim)
  } else if (!is.null(maxarea)) stop("maxarea must be matrix or NULL in pred.solution")

  # initialisiere PI; Ydim, Xdim; xmeans, xsds
  PI <- data.frame(pi.l=numeric(), pi.r=numeric(), pc=integer(), nr=integer())
  Ydim <- length(tgmean)
  if (is.matrix(dat)||is.null(dat$y)||is.null(dat$x)) dat <- data.frame(x=I(as.matrix(dat[,1:(dim(dat)[2]-Ydim)])), y=I(as.matrix(dat[,-(1:(dim(dat)[2]-Ydim))])))
  Xdim <- dim(dat$x)[2]   #  if (is.vector(dat$x)) {Xdim <- 1} else {Xdim <- dim(dat$x)[2]}
  modeg <- NULL # # Container für benutzte Modelle/Polynomgrade

  #WARNUNG wenn DIMENSION nich zu tgmean passt
  #dimension von x und y
  if (is.vector(dat$y)) {Ydimcheck <- 1} else {Ydimcheck <- dim(dat$y)[2]}
  if (Ydimcheck!=Ydim) stop("pred.solution: dimension of dat$y does not fit dimension of tgmean")
  rm(Ydimcheck)

  # DIMENSIONSREDUKTION X
  xmeans <- colMeans(dat)
  xsds <- apply(dat,2,sd)
  # DIMENSIONSREDUKTION Y # tg pca
  if (Ydim>1) {
    myypca <- tgpca(dat$y,tgmean,tgerr, wg=tgpcawg, wfun=wfun, mknormweights=mknormweights, yweights=yweights,ylast=ylast)
    y1d <- myypca$pca$x[,tgdim,drop=F] #y1d <- myypca$pca$x[,"PC1",drop=F] # tgdim
    tgmean_norm <- myypca$pcatg[,tgdim] #tgmean_norm <- myypca$pcatg[,"PC1"]
    tgerr_norm <- myypca$pcatgerr[,tgdim] #tgerr_norm <- myypca$pcatgerr[,"PC1"]
    #betterweights=T
    if (betterweights==T) {
      ywghts <- ftweights(wfun(sqrt(rowSums(myypca$pca$x[,-tgdim,drop=F]^2)))) #ywghts <- ftweights(wfun(sqrt(rowSums(myypca$pca$x[,-1,drop=F]^2))))
    } else if (betterweights==2) {
      ywghts <- ftweights(wfun(sqrt(rowSums(myypca$pca$x^2))))
    } else if (betterweights==3) {
      ywghtsvar <- rowSums(myypca$pca$x[,-tgdim,drop=F]^2)
    } else ywghts<-1
  } else {  betterweights=F
            ywghts<-1
            ywghtsvar <- 1
    y1d <- dat$y
    tgmean_norm <- tgmean[[1]] #???? warum liste
    if (!anyNA(tgerr)) tgerr_norm <- tgerr[[1]] else tgerr_norm <- NA #### ACHTUNG geändert: tgerr[[1]]+3 else
  }
  # PLS
  doCrosVal <- !(nrow(dat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL - contact Gaston Sanchez for Bug report?
  tryCatch( { pls1 = plsdepot::plsreg1(dat$x, y1d, comps = Xdim, crosval=doCrosVal) }
            , error= function(e) { print("plsreg1 throws error:"); print(e); stop("plsreg1 failed in pred.solution")})
  #pls1 = plsdepot::plsreg1(dat$x, y1d, comps = Xdim)
  plsXdim <- dim(pls1$x.scores)[2]
  # Daten eindimensional für jede HK im Liste betrachten
  dat1d.PC12 <- list()
  for (jjj in 1:plsXdim) {
    dat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,paste0("t",jjj)], y=as.numeric(y1d))
  }
  # extended version - nur shift-liste # 1d punkt plus shift korrektur (für additive mit weights)
  newx1d.PC12.ext <- list()
  # für jede HK eindimensionale Vorhersage treffen
  newx1d.PC12 <- list()

  #Gewichte bestimmen #### bei bpcenter=T - warum zentrieren mit bestem punkt? Käse?
  if (useweights==TRUE) {
    if (bpcenter==F ) {
      weightsall <- wfun(euclw(pls1$x.scores, mknormweights, sto)) #weightsall <- (1/euclw(pls1$x.scores, mknormweights))^2 # ^1 oder ^2?
      weightsvar <- (euclw(pls1$x.scores, mknormweights, sto))^2
    } else { #sweep(matrix(0:2, nrow = 3, ncol=2), 2, c(0,10))
      #stop("folgende Zeile korrigieren!")
      if (is.vector(dat$y)) {tmp1 <- abs(dat$y-tgmean[[1]]); stop("dat$y is vector in pred.solution")} else {
        tmp1 <- sweep(dat$y,2,as.numeric(tgmean)) # ACHTUNG!!! soll das nicht dat$y sein??????
        tmp1 <- sqrt(rowSums( tmp1^2 ))
      }
      #finde min abstand
      minpoint <- pls1$x.scores[which.min(tmp1),]
      keepidcs <- rep(1,nrow(pls1$x.scores))
      if (!is.null(knearest)&bpcenter==T) { #&bpcenter==T
        keepidcs <- rep(0,nrow(pls1$x.scores))
        idcs <- sort(mvdistance(pls1$x.scores, minpoint, TRUE), index.return=T)$ix[1:min(knearest,nrow(pls1$x.scores))]  # abstände zu letzter beobachtung#sortieren und auswählen
        keepidcs[idcs] <- 1
       # print(keepidcs)
      }
      weightsall <- keepidcs*wfun(pmax(euclw(sweep(pls1$x.scores,2,minpoint), mknormweights,sto),0.0000001))  # weightsall <- (1/pmax(euclw(sweep(pls1$x.scores,2,minpoint), mknormweights),0.001))^2
      weightsvar <- (euclw(sweep(pls1$x.scores,2,minpoint), mknormweights,sto))^2
    }
  } else {
    weightsall=NULL
    weightsvar <- NULL
  }
  weightsall <- ftweights(weightsall)

  pjjj <- NULL # remember past jjj's
  newx1d.seq <- list()
  newx1d.ext <- matrix(NA,nrow=0,ncol=plsXdim) #?? oder Xdim
  #newx1d.seq[[max(pcnr, length(pcnr))]] <- list(NULL)
  if (length(pcnr)>plsXdim) {
    pcnr <- pcnr[1:plsXdim]
    warning(paste("pcnr replaced by",pcnr))
  }
  for (iseq in 1:length(pcnr)) newx1d.seq[[iseq]] <- matrix(NA,nrow=0,ncol=iseq)
  for (jjj in pcnr) {                          # kein stop ! sonst irgendwas passiert!
   # browser()
    if (sum(!is.na(weightsall[,jjj]))==0) warning("no non-NAs in pred.solution for weightsall[,jjj] (with mknormweights=T?) - degByBIC will probably fail in next line")

    if (sequential==F) {
      if (betterweights==3) {
        degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=ftweights(1/(ywghtsvar+weightsvar[,jjj])), mindeg=mindeg)
      } else {
        degjjj <- degByBIC(dat1d.PC12[[jjj]], weights=weightsall[,jjj], mindeg=mindeg)
      }
    modeg <- c(modeg, degjjj)
    #if (degjjj==0) stop("degjjj=0 - f1")
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
    }#sqrt(ywghts)*sqrt(weightsall[,jjj])
    } else  newx1d.PC12[[jjj]] <- NULL # macht das was kaputt??
    if (sequential==T) { # wie wäre obiges dann? assume bpcenter=F, useweights=T
      # remember past jjj's:
      pjjj <- c(pjjj,jjj)

      if (length(pjjj)==1) { # (is.null(pjjj))
        w <- ftweights(wfun(euclw(pls1$x.scores, mknormweights,sto)))
        wvar <- euclw(pls1$x.scores, mknormweights,sto)^2
        if (length(dim(w))<2) { print(w); print(wfun(euclw(pls1$x.scores, mknormweights,sto))); print(euclw(pls1$x.scores, mknormweights,sto)); print(pls1$x.scores); stop("dim von w < 2") }
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
        for (z in 1:nrow(newx1d.seq[[length(pjjj)-1]])) { # jeder bisherige punkt soll verbessert werden
          #gewicht für datenpunkte minus shift
          shift <- list(); for (ishift in 1:Xdim) shift[[ishift]] <- 0 # shift initialisieren mit 0
          for (ishiftagain in head(pjjj,-1)) {
            #browser() #spaltennamen werden benutzt und müssen hier spätestens aktualisiert/ergänzt werden
            colnames(newx1d.seq[[length(pjjj)-1]]) <- head(pjjj,-1)
            shift[[ishiftagain]] <- repna(newx1d.seq[[length(pjjj)-1]][z,paste(ishiftagain)]) }
                                                                                                  # achtung oben, hk, nicht allgemeingültig so...
                                                                    # achtung, unten genaugenommen nicht ok: unlist und x.scores
         # print(shift);print(unlist(shift))
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
          }#weights=sqrt(ywghts)*sqrt(w[,jjj])

          if (length(newroots)>0) {
            #print("neue punkte")
            # print(newx1d.seq[[length(pjjj)-1]])
            # print(newroots)
            # #print(w)
            # print(w[,jjj]/sum(w[,jjj]))
           # browser()
            newx1d.seq[[length(pjjj)]] <- rbind(as.matrix(newx1d.seq[[length(pjjj)]]),
                                                as.matrix(cbind(newx1d.seq[[length(pjjj)-1]][z,,drop=F],data.frame(  newroots  ) )) )
            #print("rbind over")
          } else {
             warning("neue nullkoordinaten mit NA eingefügt")
            # print(cbind(newx1d.seq[[length(pjjj)-1]][z,],0 ))
        #    browser()
            newx1d.seq[[length(pjjj)]] <- rbind(as.matrix(newx1d.seq[[length(pjjj)]]),
                                                as.matrix(cbind(newx1d.seq[[length(pjjj)-1]][z,,drop=F],NA )) ) # oder na statt 0?
          }
          names(newx1d.seq) <- pjjj
        }
        } else { message("i:nopts")
          newx1d.seq[[length(pjjj)]] <- matrix(NA,nrow=0,ncol=length(pjjj)) #newx1d.seq[[length(pjjj)-1]] # korrigiert das den F?
        }
      }
    }
    #print(warnings())

    if (additive==TRUE) {
       #
    }#end if additive
    #Varianzkriterium? tgerr ode tgerr_norm
    if (betterweights==3) {
    #  browser()
      PI.mima <- t(PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]],weights=ftweights(1/(ywghtsvar+weightsvar[,jjj])), mindeg=mindeg),
                                      weights=weightsall[,jjj])))
    } else {
    #  browser()
    PI.mima <- t(PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]],weights=weightsall[,jjj], mindeg=mindeg),weights=weightsall[,jjj])))
    }
 #   message("Prädiktions-Intervallbreite für PC",jjj," in [", PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)))[1],",",PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)))[2],"]")
    if(!anyNA(tgerr_norm)) {   # here no betterweights3 in following line...
   #   browser()
      if (PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]],weights=weightsall[,jjj], mindeg=mindeg),weights=weightsall[,jjj]))[1]>2*tgerr_norm & !anyNA(tgerr)) {
        is2high <- 1
       # #if (vartoobig==-1) {vartoobig<-i}
       # warning("Varianz zu groß. Prädiktions-Intervallbreite von [", PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg))),"]>",2*tgerr_norm)
       # #     warning("Varianz zu groß. Prädiktions-Intervallbreite: ", 2*qnorm(1-0.05/2)*sqrt(ssw(dat)),">",2*tgerr)
      } else is2high <- 0
    } else {warning("pred.sol: there is NA in tgerr_norm in prediction function"); is2high<-NA}
 #   browser()
    PI <- rbind(PI, cbind(PI.mima,jjj,is2high))
  }
  isprediction <- list()
  # Punkte nach "System" wählen:
  if (sequential==T) {
    tt <- na.omit(newx1d.seq[[length(pjjj)]])  # alternativ repna...?
    #if (nrow(tt)>0) newx1d.seq[[length(pjjj)]] <- tt # lieber nicht omit? nur wenn trotzdem noch punkte da sind

    newx1d.seq[[length(pjjj)]] <- tt # oder doch auf jeden fall?
  }
  for (jjj in pcnr) {
    # ((length(newx1d.PC12[[jjj]])==0 &sequential==F)|(sequential==T)) {   # keine Vorhersage, Punkt "systematisch" wählen
                                                # &jjj==pcnr[1] # fehler?: & length(newx1d.seq[[length(pjjj)]])==0))
      if (sequential==T) {if (!nrow(newx1d.seq[[length(pjjj)]])==0) {isprediction[[jjj]] <- TRUE;# message(paste("ohoh:",newx1d.seq[[length(pjjj)]]));message(paste("ohoh:",str(newx1d.seq[[length(pjjj)]])));
      break
                          } else message(paste("gut:",length(newx1d.seq[[length(pjjj)]]), "und" , nrow(newx1d.seq[[length(pjjj)]])))
        }
      if (sequential==F) {
        if (!length(newx1d.PC12[[jjj]])==0) {isprediction[[jjj]] <- TRUE; next}
      }
    message(paste0("Schritt ",stepi," pcnr", jjj,": gewähltes Modell findet keine Nullstellen in keiner Richtung. Messe an neuem Punkt außerhalb des bisherig betrachteten Bereichs."))

     #  print(paste("ptchoose",jjj,"pc", pcnr[(nptc %% Xdim)+1] )) # Xdim statt length(pcnr)
      # print(str(dat1d.PC12))
      # print(dat1d.PC12)
       if (ptchng==T) {
         ## nur wenn erster durchgang
         if (jjj==pcnr[1]) {
           #jind <- pcnr[(nptc %% Xdim)+1] # Xdim statt length(pcnr)
         #  browser() #plsXdim?
           jind  <- (nptc %% plsXdim)+1
           print(paste("ptchoose",jjj,"pc:", jind, "nptc:", nptc ))
           #wie soll pcnr brücksichtigt werden?
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
                                        xmeans=xmeans, xsds=xsds, pls1=pls1, jjj=jjj) # wählt ggf in beiden/allen richtungen
       }
 #   browser() # plsXdim?
    ourshift <- rep(NA, plsXdim)
    if( length(newx1d.PC12[[jjj]])>1) stop("nicht gehandlet pred.sol")
    ourshift[jjj] <- newx1d.PC12[[jjj]]
    newx1d.ext <- rbind(newx1d.ext,
                        pred.ptch(shift=ourshift, pcnr=pcnr, pls1=pls1, mknormweights=mknormweights, dat1d.PC12=dat1d.PC12, mindeg=mindeg, tgmean_norm=tgmean_norm, gr2retlimit=gr2retlimit, wfun=wfun, sto=sto) )
  }

  improvePT <- T #!!!!!!!!!!!!!!!##############
  if (improvePT==T) { if(sequential==T) {
    #print("vorhert")
   # print(length(pjjj)) ; print(pjjj) ; print(newx1d.seq)
    if (nrow(newx1d.seq[[length(pjjj)]])==0) { warning("improvett"); print("improvett")
      newx1d.seq[[length(pjjj)]] <- rbind(newx1d.seq[[length(pjjj)]],na.omit(newx1d.ext))
      # tt <- rbind(tt,newx1d.ext)
      # tt <- na.omit(tt)
      if (nrow(newx1d.seq[[length(pjjj)]])==0) message("immernoch0")
    }
   # print("NAC")
  }}

  # Vom 1D in den xdim-D: newx1d -> newx
  xcols <- grep("x.", names(xmeans), value = TRUE)
  newx.PC12 <- list()
  if (additive==TRUE) {
    # wir müssen hier die neuen shifts berücksichtige unten beim rücktransformieren
  } else { # end if additive
  #  browser() #plsXdim?
    tpcnr <- 1:plsXdim # pcnr # Xdim
    if (sequential==T) {
      if (dim(pls1$mod.wgs)[1]==dim(pls1$mod.wgs)[2]) {
        newx.seq <- mkreg(as.matrix(newx1d.seq[[length(pjjj)]])%*%solve(pls1$mod.wgs)[pcnr,],xmeans[xcols],xsds[xcols])
      } else newx.seq <- mkreg(as.matrix(newx1d.seq[[length(pjjj)]])%*%MASS::ginv(pls1$mod.wgs)[pcnr,],xmeans[xcols],xsds[xcols])
    }
   # print("vorher")
    if (length(pjjj)>0) message(paste("nNewxsq1d",nrow(newx1d.seq[[length(pjjj)]]))) # ausgabe nur wenn sequential bzw pjjj-länge >0
   # print("nachher")

  #  browser()
    for (jjj in tpcnr) {#warning("nr3")
      if (length(newx1d.PC12)>=jjj) # passiert (nicht) wenn sequential
      {
        if (dim(pls1$mod.wgs)[1]==dim(pls1$mod.wgs)[2]) {
        newx.PC12[[jjj]] <- mkreg(matrix(newx1d.PC12[[jjj]])%*%solve(pls1$mod.wgs)[jjj,],xmeans[xcols],xsds[xcols])
        } else newx.PC12[[jjj]] <- mkreg(matrix(newx1d.PC12[[jjj]])%*%MASS::ginv(pls1$mod.wgs)[jjj,],xmeans[xcols],xsds[xcols])
      } else newx.PC12[[jjj]] <- list(NULL) #macht das was kaputt?
    }
  }
  ###nur bestimmen, wenn neue Punkte um mehr als xeps verschieden: uniqP
  for (jjj in tpcnr) {   ### weshalb unten unlist?? hä??
    if (is.null(unlist(newx.PC12[[jjj]]))) { ## für den rest der schleife tpcnr == jjj ignorieren
      tpcnr <- setdiff(tpcnr,jjj)
    }
    newx.PC12[[jjj]] <- uniqP(newx.PC12[[jjj]], xeps, dat$x)
  }
  if (sequential==T) newx.seq <- uniqP(newx.seq, xeps, dat$x)

  newx<-NULL
  is.prediction <- numeric(0)
  message(paste("DATMEAN",colMeans(dat$x)))

  for (jjj in tpcnr) {
    newx <- rbind(newx,newx.PC12[[jjj]] )
    is.prediction <- c(is.prediction, rep(isprediction[[jjj]],nrow(newx.PC12[[jjj]])))
  }
  if (sequential==T) {
    newx.seq <- nclose2mean(newx.seq,colMeans(dat$x))
    newx <- rbind(newx,newx.seq )
    is.prediction <- c(is.prediction, rep(5,nrow(newx.seq)))
    #message(paste("seq"))
    #print(newx.seq)
  }
  return(list(newx=newx, dat1d.PC12=dat1d.PC12, newx1d.PC12=newx1d.PC12,  is.prediction = is.prediction, isprediction=isprediction, PI=PI, modeg=modeg))#cbind(PI, is.prediction, unlist(isprediction))))
}

plotexp <- function(dat, target, tgerr, prediction=NULL, model=NULL, limit0=FALSE, reality=NULL) {
  fm <- model

  # todo: plotbereich von 1 bis 7 automatisch bestimmen lassen?

  #block doppelt, ier ohne plot, unten mit plot
  preds<-data.frame(x=NULL,y=NULL) #für limit0...
  if (!is.null(prediction)) {
    # roots
    rootdta <- data.frame(x = prediction)
    rootdta$y <- predict(fm, newdata=rootdta)
    #   myplot <- myplot +
    #    geom_point(data = rootdta, aes(x = x, y = y), color="red")
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

  #realität plotten?
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
    labs(x = "Prozessparameter (p)", y="Deskriptor (d)") +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=target-tgerr, ymax=target+tgerr, alpha=0.5, fill="green") +
    geom_hline(yintercept = target)

  if (!is.null(fm)) {
    myplot <- myplot +
      geom_line(aes(x=x,y=fit) ,data=prd)
  }

  preds<-data.frame(x=NULL,y=NULL) #für limit0...
  if (!is.null(prediction)) {
    # roots
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

pred.solution.old <- function(dat,tgmean,tgerr=NULL,xeps=0.001,pcnr,additive=FALSE, maxarea=NULL, ptchoice=1, useweights=TRUE) {
  tgmean<-matrix(tgmean,nrow=1)
  if (!is.null(tgerr)) {
    tgerr <- matrix(tgerr,nrow=1)
  } else tgerr <- matrix(NA,nrow=1,ncol=length(tgmean))
  #tgerr<- ifelse(!is.null(tgerr), matrix(tgerr,nrow=1), matrix(NA,nrow=1,ncol=length(tgmean)))
  #WARNUNG wenn DIMENSION nich tzu tgmean passt??
  #dimension von x und y
#  if (is.vector(dat$y)) {Ydim <- 1} else {Ydim <- dim(dat$y)[2]}
#  if (is.vector(dat$x)) {Xdim <- 1} else {Xdim <- dim(dat$x)[2]}

  if (is.matrix(maxarea)) {
    if (dim(maxarea)[1] != dim(tgmean)[2] |  dim(maxarea)[2] != 2)
      stop("wrong dimensions of maxarea")
  }

  Ydim <- length(tgmean)
  if (is.matrix(dat)||is.null(dat$y)||is.null(dat$x)) dat <- data.frame(x=I(as.matrix(dat[,1:(dim(dat)[2]-Ydim)])), y=I(as.matrix(dat[,-(1:(dim(dat)[2]-Ydim))])))
  Xdim <- dim(dat$x)[2]
  # DIMENSIONSREDUKTION X
  xmeans <- colMeans(dat)
  xsds <- apply(dat,2,sd)
  # DIMENSIONSREDUKTION Y # tg pca
  if (Ydim>1) {
    myypca <- tgpca(dat$y,tgmean,tgerr)
    y1d <- myypca$pca$x[,"PC1",drop=F]
    tgmean_norm <- myypca$pcatg[,"PC1"]
    tgerr_norm <- myypca$pcatgerr[,"PC1"]
  } else {
    y1d <- dat$y
    tgmean_norm <- tgmean[[1]]
    #tgerr_norm <- ifelse(!is.na(tgerr), tgerr[[1]], NA)
    if (!anyNA(tgerr)) tgerr_norm <- tgerr[[1]]+3 else tgerr_norm <- NA
  }
  # PLS
  pls1 = plsdepot::plsreg1(dat$x, y1d, comps = Xdim)
  # Daten eindimensional für jede HK betrachten
  dat1d.PC12 <- list()
  for (jjj in pcnr) {
    #browser()
    dat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,paste0("t",jjj)], y=as.numeric(y1d))
  }
  # für jede HK eindimensionale Vorhersage treffen
  newx1d.PC12 <- list()
  for (jjj in pcnr) {
    if (useweights==TRUE) {
      # get weights for regression - might work more easy with list operation
      weights <- rep(0,nrow(dat1d.PC12[[jjj]]))
      for (kkk in setdiff(pcnr,jjj)) {
        weights <- weights+(dat1d.PC12[[kkk]]$x)^2
      }
      weights <- sqrt(weights) # eucl distance
    } else {
      weights<-NULL }
    if (!is.null(weights)) weights <- (1/mknorm(weights))^2
    newx1d.PC12[[jjj]] <- getroots2(dat1d.PC12[[jjj]],
                                    degByBIC(dat1d.PC12[[jjj]]),
                                    tgmean_norm,
                                    limit=1.6*diff(range(dat1d.PC12[[jjj]]$x)), weights=weights)
    if (additive==TRUE) {
      for (ipnts in 1:nrow(newx1d.PC12[[jjj]])) {

      }
      # dann punkt verbessern in richtung dimension lll
      # hier fehlt schleife. für jeden punkt
      ipcnr <- setdiff(pcnr,jjj)
      shift <- list(); for (ishift in pcnr) shift[[ishift]] <- 0
      shift[[jjj]] <- wohin#############
      for (lll in setdiff(pcnr,jjj)) {
        #iterativ
        ipcnr <- setdiff(ipcnr,lll)
        if (useweights==TRUE) { #über abstände
          weights <- rep(0,nrow(dat1d.PC12[[lll]]))
          for (kkk in setdiff(pcnr,lll)) {
            weights <- weights+(dat1d.PC12[[kkk]]$x-shift[kkk])^2
          }
          # weights distances from dim lll at predicted point
          #addinf for doím lll
        } else {
          weights<-NULL }


      }
    }
    #Varianzkriterium? tgerr ode tgerr_norm
    message("Prädiktions-Intervallbreite für PC",jjj," in [", PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]])))[1],",",PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]])))[2],"]")
    if(!anyNA(tgerr_norm)) {
      if (PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]])))[1]>2*tgerr_norm & !is.na(tgerr)) {
        ############QQQQQQQQQQ!!!!!!!!!!!!!!!!!: the condition has length > 1 and only the first element will be used
        #if (vartoobig==-1) {vartoobig<-i}
        warning("Varianz zu groß. Prädiktions-Intervallbreite von [", PIhcheck(polymodel(dat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]]))),"]>",2*tgerr_norm)
        #     warning("Varianz zu groß. Prädiktions-Intervallbreite: ", 2*qnorm(1-0.05/2)*sqrt(ssw(dat)),">",2*tgerr)
      }
    }
  }
  isprediction <- list()
  # Alle Richtungen erzwingen? Ja!
  if (FALSE) { # dont force solutions in all PC directions
    # if (length(unlist(newx1d.PC12))==0 ) {
    #   message(paste0("Schritt ",i,": gewähltes Modell findet keine Nullstellen in keiner Richtung. Messe an neuem Punkt außerhalb des bisherig betrachteten Bereichs."))
    #   for (jjj in pcnr) {
    #     if (min(abs(dat1d.PC12[[jjj]]$y[dat1d.PC12[[jjj]]$x==min(dat1d.PC12[[jjj]]$x)]-tgmean_norm)) < min(abs(dat1d.PC12[[jjj]]$y[dat1d.PC12[[jjj]]$x==max(dat1d.PC12[[jjj]]$x)]-tgmean_norm))) {
    #       newx1d.PC12[[jjj]] <- min(dat1d.PC12[[jjj]]$x) - (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2
    #     }
    #     else {
    #       newx1d.PC12[[jjj]] <- max(dat1d.PC12[[jjj]]$x) + (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2
    #     }
    #   }
    # }
  } else {
    for (jjj in pcnr) {
      if (length(newx1d.PC12[[jjj]])==0 ) {
        isprediction[[jjj]] <- FALSE
        # hier neuen punkt in richtung wählen - bräuchten alle projektionen auf hauptkomponente um größte lücke zu finden
        cat("1d points avail.:",paste(" ",sort(c(0, unique(dat1d.PC12[[jjj]][[1]])))))
        cat("1d points avail.:",paste(" ",diff(sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))))))
        #alle neuen möglichen messpunkte:
        #cat("meas at:",paste(" ", head( sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))) ,-1) + diff(sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))))*0.5))
        cat("meas at:",paste(" ", head( sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))) ,-1) + diff(sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))))*0.5))
        cat("ext points",
            min(dat1d.PC12[[jjj]]$x) - (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2,
            " and ",
            max(dat1d.PC12[[jjj]]$x) + (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2)
        if (ptchoice==1) {
          #message(paste0("Schritt ","i",": gewähltes Modell findet keine Nullstellen in Richtung PC ",jjj,". Messe an neuem Punkt außerhalb des bisherig betrachteten Bereichs."))
          if (min(abs(dat1d.PC12[[jjj]]$y[dat1d.PC12[[jjj]]$x==min(dat1d.PC12[[jjj]]$x)]-tgmean_norm)) < min(abs(dat1d.PC12[[jjj]]$y[dat1d.PC12[[jjj]]$x==max(dat1d.PC12[[jjj]]$x)]-tgmean_norm))) {
           newx1d.PC12[[jjj]] <- min(dat1d.PC12[[jjj]]$x) - (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2
          }
          else {
           newx1d.PC12[[jjj]] <- max(dat1d.PC12[[jjj]]$x) + (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2
          }
        } else if (ptchoice==2) {
          newx1d.PC12[[jjj]] <- c((head( sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))) ,-1) + diff(sort(c(0, unique(dat1d.PC12[[jjj]][[1]]))))*0.5)[which.max(diff(sort(c(0, unique(dat1d.PC12[[jjj]][[1]])))))],
                                  min(dat1d.PC12[[jjj]]$x) - (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2,
                                  max(dat1d.PC12[[jjj]]$x) + (max(dat1d.PC12[[jjj]]$x) - min(dat1d.PC12[[jjj]]$x))/2
                                  )
        } else if (ptchoice==3) {
          # nur 1 punkt für max abstand
          # wenn rand gegeben, dann wähle als abstand  2*(kleinstervorkommenderxwert-minpossiblexval) bzw
          #                                           2*(groestervorkommenderxwert-maxpossiblexval)
          # wenn keine grenzen dann wie ptchoice1
        }
      } else isprediction[[jjj]] <- TRUE
    }
  }
  # Vom 1D in den xdim-D: newx1d -> newx
  xcols <- grep("x.", names(xmeans), value = TRUE)
  newx.PC12 <- list()
  if (additive==TRUE) {
    if(FALSE){
    # # newx1d.add <- NULL #besser sowas? matrix(nrow=length(newx1d.PC12[[pcnr[1]]]),ncol=0)
    # # for (jjj in pcnr) {
    # #   newx1d.add <- cbind(newx1d.add,newx1d.PC12[[jjj]])
    # # } # oben evtl. besser als code unten, da es mit dimension 1 funktioniert?
    # newx1d.add <- newx1d.PC12[[pcnr[1]]]
    # for (jjj in pcnr[-1]) {
    #   newx1d.add <- cbind(newx1d.add,newx1d.PC12[[jjj]])
    # }
    # tpcnr <- c(1)
    #
    # #hier jetzt aber auch die nicht additiven punkte mitnehmen
    # tpcnr <- c(pcnr,max(pcnr)+1)
    # #der not additive teil:
    # for (jjj in pcnr) {
    #   newx.PC12[[jjj]] <- mkreg(matrix(newx1d.PC12[[jjj]])%*%solve(pls1$mod.wgs)[jjj,],xmeans[xcols],xsds[xcols])
    # }
    # # der additive teil
    # newx.PC12[[max(pcnr)+1]] <- mkreg(newx1d.add%*%solve(pls1$mod.wgs)[pcnr,],xmeans[xcols],xsds[xcols])
    }
  } else {
    tpcnr <- pcnr
    print("some cck")
    for (jjj in pcnr) {
      if (!is.null(newx1d.PC12[[jjj]])) { #if (length(newx1d.PC12)>=jjj) # passiert (nicht) wenn sequential
        newx.PC12[[jjj]] <- mkreg(matrix(newx1d.PC12[[jjj]])%*%solve(pls1$mod.wgs)[jjj,],xmeans[xcols],xsds[xcols])
      } else newx.PC12[[jjj]] <- list(NULL) # ok, macht das was kaputt?
    }
  }
  ###nur bestimmen, wenn neue Punkte um mehr als xeps verschieden
  ## PC12schleife
  dstncs.PC12 <- list()
  for (jjj in tpcnr) {
    ###nur bestimmen, wenn neue Punkte um mehr als xeps verschieden
    dstncs.PC12[[jjj]] <- matrix(nrow = dim(newx.PC12[[jjj]])[1],ncol=dim(dat$x)[1])
    #folgende Zeile nur wenn überhaupt punkte vorhanden in PCjjj:
    if (dim(newx.PC12[[jjj]])[1]>0)
      for (cnt in 1:dim(newx.PC12[[jjj]])[1]) {
        dstncs.PC12[[jjj]][cnt,] <- mvdistance(dat$x, newx.PC12[[jjj]][cnt,])
      }
  }
  for (jjj in tpcnr) {
    if (!is.null(newx.PC12[[jjj]])) # ?passt das?
    newx.PC12[[jjj]] <- newx.PC12[[jjj]][ suppressWarnings(apply(dstncs.PC12[[jjj]],1,min))>xeps ,,drop=F]
  }
  newx<-NULL
  for (jjj in tpcnr) {
    newx <- rbind(newx,newx.PC12[[jjj]] )
  }
  return(list(newx=newx, dat1d.PC12=dat1d.PC12, newx1d.PC12=newx1d.PC12))
}

#' newexp
#'
#' Creates new samples from model function. Used for simulations in the function autosolve as a helper function.
#'
#' @param n number of repeated measurements
#' @param xpos coordinates at which the measurement takes place
#' @param foo model for the real relationship
#'
#' @return a matrix containing the samples
#' @export
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


#' Autosolve / PLSolve.auto
#'
#' @param startx start values
#' @param tgmean target value vector
#' @param tgerr vector of acceptable error range
#' @param reps number of repeated measurements
#' @param maxit maximum number of iterations
#' @param reality real model
#' @param xeps smallest unit?
#' @param pplot do plots?
#' @param pcnr which principal directions will be considered
#' @param additive use soem additive approach (deprecated)
#'
#' @return
#' @export
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
#' #additive=TRUE
#' tmp<-autosolve(startx,tgmean,tgerr*0.125,reps=7,maxit=10,tfoo, xeps=0.01, F, pcnr=c(1,2), additive=T)
#' #nwexample
#' dstgmean <- tfoo(cbind(1.35,1.4))#,1.5))#tfoo(cbind(0.35,0.4,0.5))
#' tgerr <- c(0.2,0.2)#,0.5)
#' reps <- 5
#' maxit <-10
#' reality <- tfoo
#' xeps <- 0.01
#' pplot <- TRUE
#' additive=TRUE
#' startx <- expand.grid(x1=c(-1,3),x2=c(-1,3))
#' set.seed(123)
#' autosolve(startx,tgmean,tgerr*0.125,reps=7,maxit=10,tfoo, xeps, F, pcnr=c(1,2), additive=F)
#' # ? zurückgeben mit Welcher WK man im Ziel landet
autosolve <- function(startx, tgmean, tgerr, reps=25, maxit=10, reality=foo, xeps=0.01, pplot=TRUE, pcnr=c(1,2), additive=FALSE,maxarea=NULL, useweights=TRUE, mknormweights=F,gr2retlimit=T, mindeg=0,sequential=F,tgpcawg=1,yweights=F,datlim=NULL, knearest=NULL,tgdim=1,ylast=NULL,sto=T,mod.sd=NULL,...) {
  tgmean<-matrix(tgmean,nrow=1)
  tgerr<- matrix(tgerr,nrow=1)
  # Dimension?
  Xdim <- dim(startx)[2]
  Ydim <- dim(tgmean)[2]
  dat <- data.frame(x=numeric(),y=numeric(),nri=numeric()) # x und y eigentlich unsinnig definiert
  newx <- startx
  vartoobig <- -1;  solfound <- -1;  nomore <- -1
  modeg <- matrix(NA, ncol=2, nrow=0)
  for (i in 1:maxit){
    # wenn i>0 dann startx verändern: getroots(polymodel(data3,degByBIC(data3, mindeg=mindeg)),tgmean)
    if (i==1) {
      PI <- data.frame(pi.l=numeric(), pi.r=numeric(), pc=integer(), nr=integer(), nri=integer())
      ispred <- F
      is.pred <- NULL
    }
    if (i>1) {
      #print(ispred)
      #print(!ispred)
      flag <- TRUE
      tryCatch( { rs <- pred.solution(dat=dat,tgmean=tgmean,tgerr=tgerr,xeps=xeps,pcnr=pcnr,additive=additive,maxarea=maxarea, useweights=useweights, mknormweights=mknormweights,gr2retlimit=gr2retlimit,mindeg=mindeg,sequential=sequential,nptc=sum(!ispred),tgpcawg=tgpcawg,yweights=yweights,datlim=datlim,knearest=knearest,tgdim=1+((i-1-1) %% tgdim), ylast=ylast,sto=sto,...) }
                , error= function(e) { print("break out of for - pred.solution throws error:"); print(e); flag<<-FALSE})
      if (!flag) break
      #rs <- pred.solution(dat=dat,tgmean=tgmean,tgerr=tgerr,xeps=xeps,pcnr=pcnr,additive=additive,maxarea=maxarea, useweights=useweights, mknormweights=mknormweights,gr2retlimit=gr2retlimit,mindeg=mindeg,sequential=sequential,nptc=sum(!ispred),tgpcawg=tgpcawg,yweights=yweights,datlim=datlim,...)
      newx <- rs$newx
      # überprüfen ob NAs vorhanden
      if (anyNA(newx)){
        print(tail(dat))
        print(newx)
        print(tgmean)
        print(tgerr)
        print(rs)
        stop("datn hat NA")
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
      message(paste0("Schritt ",i,": keine neuen Punkte"))
      break
    }
  #  browser()
    datn <- newexp(n=reps,xpos=newx, foo=reality, sd=mod.sd)
    datn$nri <- i
    # man könnte noch speichern in welcher PC Richtung die Punkte hinzugefügt wurden

    # überprüfen ob NAs vorhanden
    if (anyNA(datn$y[,"y1",drop=F])){
      print(datn)
      stop("datn hat NA")
    }

    newdat <- rbind(dat,datn)
    #plot
    if (pplot==TRUE) {
      if (i>1) {
        #wo tauchen die neuen punkte im alten koordinatensystem auf?
        xcols <- grep("x.", names(xmeans), value = TRUE)
        tmpscores <- mknorm(newdat$x,xmeans[xcols],xsds[xcols])%*%pls1$mod.wgs
        ##### hier unten y falsch??, wir brauchen y nach alter pca##############################
        tmpscoresy <- mknorm(newdat$y,myypca$allobsmean,myypca$allobssd) %*% myypca$pca$rotation
        #funktioniert das für 1d-y?
        newdat1d.PC12 <- list()
        for (jjj in pcnr) {
          newdat1d.PC12[[jjj]] <- data.frame(x=tmpscores[,jjj], y=tmpscoresy[,1]) #newdat1d <-   data.frame(x=tmpscores[,1], y=newdat$y[,"y1"])
          print(plotexp(dat=newdat1d.PC12[[jjj]], target=myypca$pcatg[,"PC1"], myypca$pcatgerr[,"PC1"], prediction=NULL, model=polymodel(dat1d.PC12[[jjj]], degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)), limit0=TRUE) )
        }
      }
      # DIMENSIONSREDUKTION (erst für den folgenden Plot ermitteln)
      xmeans <- colMeans(newdat)
      xsds <- apply(newdat,2,sd)
      ### DIMENSIONSREDUKTION Y; pca für erklärte y
      if (Ydim>=2) {
        ypca <- tgpca(newdat$y,tgmean,tgerr,wg=tgpcawg)  #besser t(as.numeric(tgmean)) und code vereinfachen, vgl unten
        #browser()
        pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2]) #oder comps=2 genug hier?
      } else {
        if (anyNA(newdat$y[,"y1",drop=F])){
          print(newdat)
          stop("ZEIL HAT NA 1")
        }
       # browser()
        pls1 = plsdepot::plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2])
      }
      newdat1d.PC12<-list()
      for (jjj in pcnr) {
        if (Ydim>=2) {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=ypca$pca$x[,"PC1"])
        } else {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=newdat$y[,"y1"])
        }
      }
      #warning("affe1")
      if (Ydim>=2) {
        tgmean_norm <- ypca$pcatg[,"PC1"]
        tgerr_norm <- ypca$pcatgerr[,"PC1"]
      } else {
        tgmean_norm <- mknorm(tgmean, xmeans[c("y")], xsds[c("y")]) # evtl nochmal checken
        tgerr_norm <-  mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y")]) # evtl nochmal checken
      }
      #warning("affe2")



      #plot: auf Originalskala
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
      #supertoller plot mit plotly
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
      #Erfolg?
      # DIMENSIONSREDUKTION  wie oben
      xmeans <- colMeans(newdat)
      xsds <- apply(newdat,2,sd)
      if (Ydim>=2) {
        ypca <- tgpca(newdat$y,tgmean,tgerr,wg=tgpcawg)
        #ypca <- tgpca(newdat$y,tgmean,tgerr)  #besser t(as.numeric(tgmean)) und code vereinfachen, vgl unten
      #  browser()
        doCrosVal <- !(nrow(dat$x) < 10) # package internal check for plsreg1 is insufficient when comps != NULL - contact Gaston Sanchez for Bug report?
        pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2], crosval=doCrosVal) #oder comps=2? genug? #pls1 = plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2])
        #pls1 = plsdepot::plsreg1(newdat$x, ypca$pca$x[,"PC1",drop=F], comps = dim(newdat$x)[2])
      } else {
        if (anyNA(newdat$y[,"y1",drop=F])){
          print(tail(newdat,200))
          print(xmeans)
          print(xsds)
          stop("ZEIL HAT NA 2")
        }
     #   browser()
        pls1 = plsdepot::plsreg1(newdat$x, newdat$y[,"y1",drop=F], comps = dim(newdat$x)[2])
      }
            #fragwürdiges abbruchkriterium
      newdat1d.PC12 <- list()
      for (jjj in pcnr) {
        if (Ydim>=2) {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=ypca$pca$x[,"PC1"]) #data.frame(x=pls1$x.scores[,"t1"], y=newdat$y[,"y1"])
        } else {
          newdat1d.PC12[[jjj]] <- data.frame(x=pls1$x.scores[,jjj], y=newdat$y[,"y1"]) #data.frame(x=pls1$x.scores[,"t1"], y=newdat$y[,"y1"])
        }
      }
      #mv Fehler again
      #warning("affe1")
      if (Ydim>=2) {
        tgmean_norm <- ypca$pcatg # mknorm(tgmean, xmeans[c("y")], xsds[c("y")])#mknorm(tgmean, xmeans[c("y.y1","y.y2")], xsds[c("y.y1","y.y2")])
        tgerr_norm <- ypca$pcatgerr # mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y")])#mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y.y1","y.y2")])
      } else {
        tgmean_norm <- mknorm(tgmean, xmeans[c("y")], xsds[c("y")]) # evtl nochmal checken
        tgerr_norm <-  mknorm(tgerr, rep(0,length(tgerr)), xsds[c("y")]) # evtl nochmal checken
      }
      #warning("affe2")

      for (jjj in pcnr) { # nicht ganz klar wie gerechtfertigt das für additive == TRUE ist
        if (is.null(newx1d.PC12)) break # so ok???
       # print(paste("checkpoint3 - ",jjj))
        if (jjj>length(newx1d.PC12)) break # so ok???
      #  print(paste("checkpoint35 - ",jjj))
        if (is.null(newx1d.PC12[[jjj]])) next # so ok???
      #  print(paste("checkpoint4 - ",jjj))
        # müsste hier die summe der Fehler für die pcnr stehen?
        if (degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)==degByBIC(newdat1d.PC12[[jjj]], mindeg=mindeg)) {
          checkpred <- as.data.frame(predict.lm(polymodel(newdat1d.PC12[[jjj]],degByBIC(dat1d.PC12[[jjj]], mindeg=mindeg)), newdata = data.frame(x=newx1d.PC12[[jjj]]), interval="prediction") ) #,row.names=1:(dim(dat)[1]+dim(datn)[1])
          if (any(checkpred$lwr[checkpred$upr<tgmean_norm[1]+tgerr_norm[1]] > tgmean_norm[1]-tgerr_norm[1])) {
            if (solfound==-1) {solfound<-i}
            message(paste0("Punkt gefunden mit Modell für PC",jjj))

            dat <- newdat # oder im return newdat statt dat????
            break
          }
          rm(checkpred)
        }
      }
    }
    dat <- newdat
  }
  ## hinzufügen: Wahrscheinlichkeit unter wahrem Modell für Erfolg
  # achtung: erfolgswahrscheinlichkeit nicht für roots berechnen sondern für newx
  # roots <- getroots2(dat, degByBIC(dat, mindeg=mindeg), tgmean)
  erfolgwk0 <- suppressWarnings( max(pnorm(tgmean[[1]]+tgerr[[1]], mean=reality(newx)[,1], sd=1) - pnorm(tgmean[[1]]-tgerr[[1]], mean=reality(newx)[,1], sd=1)) )
  return(list(data=dat, data1d=newdat1d.PC12, #roots=roots,
              solfound=solfound, nomorepoints=nomore, variancetoobig=vartoobig, messreihen=max(dat$nri), messungen=dim(dat)[1], erfolgwk=erfolgwk0, PI=PI,  is.pred = is.pred, ispred=ispred, modeg=modeg#, pls1st=pls1st
  ))
}

setwd("C:/Users/maendle/OneDrive/workingr/tgopt")

library(jsonlite)
#document <- fromJSON(txt="all_data_metrics.json")
document <- fromJSON("http://test.sfb1232.de:85/all_metric_dl")
dat <- document$data

#library(tidyverse)
#document %>% select(data)

hts <- c("GKZ", "Q1",  "Q4" ,"Q3") #"QT1", "B3"
dat <- subset(dat, subset=dat$ht %in% hts)
alloys <- c("100Cr6")
dat <- subset(dat, subset=dat$alloys %in% alloys)

dat$alloy
dat$ht
dat$`micro-data` <- as.data.frame(lapply(dat$`micro-data`, function(x) x$global_avg))
dat$`macro-data` <- as.data.frame(lapply(dat$`macro-data`, function(x) x$global_avg))
dat$`mat-data`   <- as.data.frame(lapply(dat$`mat-data`, function(x) x$global_avg))

input <- "micro-data"
output <- "macro-data"
ipar <- c("D04.1.1.1.1", "U03.3.1.1.2", "U04.8.1.1.4")
opar <- c("D03.1.1.1.5", "U04.2.1.1.5", "U04.6.1.1.9")
objective <- c( -59.995950,  147.699394,    4.984113 )

idat <- dat[[input]]
odat <- dat[[output]]
idat <- idat[,ipar]
odat <- odat[,opar]

usedat <- data.frame(x=idat,y=odat[,1:2])

if (any(is.na(usedat))) stop(paste("NAs in coloumn", names(usedat)[apply(usedat,2,function(x) any(is.na(x)))],
                                  "- observations will be dropped"))
usedat <- na.omit(usedat)

mvTargetOpt::pred.solution(dat=usedat,#,nri=1),
                           tgmean=objective[1:2],
                           tgerr=t(rep(100,2)),
                           xeps=-1,#0.001,
                           pcnr=1:2,#length(opar),  # oder ipar???
                           additive=FALSE,
                           maxarea=NULL, ptchoice=3, useweights=TRUE, mknormweights=F, allpts=F, gr2retlimit=F,bpcenter=F,
                           mindeg=0,wfun=function(x) {(1/x)^2}, sequential=T, ptchng=F, nptc=0,tgpcawg=1,betterweights=3,yweights=F,
                           datlim=NULL,knearest=NULL,
                           tgdim=1,
                           ylast=5,sto=F,#...)
                           mdeg=1, messages=T)
mvTargetOpt::pred.solution(dat=data.frame(x=idat,y=odat[,1:3]),#,nri=1),
                           tgmean=objective[1:3],
                           tgerr=t(rep(100,3)),
                           xeps=-1,#0.001,
                           pcnr=1:2,#length(opar),  # oder ipar???
                           additive=FALSE,
                           maxarea=NULL, ptchoice=3, useweights=TRUE, mknormweights=F, allpts=F, gr2retlimit=F,bpcenter=F,
                           mindeg=0,wfun=function(x) {(1/x)^2}, sequential=T, ptchng=F, nptc=0,tgpcawg=1,betterweights=3,yweights=F,
                           datlim=NULL,knearest=NULL,
                           tgdim=1,
                           ylast=5,sto=F,#...)
                           mdeg=1, messages=T)
# startx=stval[[1]], dstgmean=dstgmean, tgerr=tgerr, mod=MODELL123, maxit=100,
# maxarea=cbind(rep(-5,2),rep(5,2)), betterweights=3,tgdim=1,
# ylast=5,sto=F,mod.sd=0, reps=2,
# mdeg=1, messages=F,
# xeps=-1, pplot=FALSE, pcnr=c(1,2), additive=F,
# useweights=T,  mknormweights=F, ptchoice=3, allpts=F, gr2retlimit=F, sequential=T, tgpcawg=1,
#yweights=F, bpcenter=F, datlim=NULL, knearest=NULL

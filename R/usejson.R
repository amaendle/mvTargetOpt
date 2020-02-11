#' usejson
#'
#' apply algorithm to some data from json files
#'
#' @return json file with suggested parameters and flag (is parameter set a prediction or simply a suggestion how to improve the model)
#' @export TRUE
#' @examples
#' setwd("C:/Users/maendle/OneDrive/workingr/tgopt")
#'
usejson <- function(processdat="processparameter.json", formdat="formdata_20200211.json", maindat="http://test.sfb1232.de:85/all_metric_dl", outjson="out.json", ignorehts=TRUE) {
  # get process parameter data
  formdat <- jsonlite::fromJSON(formdat)
  # filter for alloy and heat treatment
  alloys <- formdat$alloys    #  alloys <- "100Cr6"
  if (ignorehts==T)
    hts <- c("Q4",  "B3",  "GKZ", "Q3",  "Q1",  "QT1")
  else
    hts <- formdat$hts

  # get micro/macro descriptors & material data
  dat <- jsonlite::fromJSON(maindat)  #dat <- fromJSON(txt="all_data_metrics.json")
  dat <- dat$data
  dat <- subset(dat, subset =  (dat$alloy %in% alloys)&(dat$ht %in% hts))

  # get process parameter data
  procdat <- jsonlite::fromJSON(processdat)
  procdat <- procdat$data
  procdat <- subset(procdat, subset =  (procdat$alloy %in% alloys)&(procdat$ht %in% hts))
  procdat <- procdat$`process-param`
  procdat <- lapply(procdat,function(x) as.numeric(x$value))
  procdat <- as.data.frame(procdat)

  # simplify data.frame (use only global_avg for now)
  dat$`micro-data` <- as.data.frame(lapply(dat$`micro-data`, function(x) x$global_avg))
  dat$`macro-data` <- as.data.frame(lapply(dat$`macro-data`, function(x) x$global_avg))
  dat$`mat-data`   <- as.data.frame(lapply(dat$`mat-data`, function(x) x$global_avg))
  dat$`process-data` <- procdat
  #View(dat)

  # select coloumns/properties
  input <- "process-data"
  output <- "micro-data"
  ipar <- formdat$`process-params` #c("U03.2.P.3.1.1", "U03.2.P.3.1.5", "U03.2.P.3.1.7","U03.2.P.3.1.11") #formdat$`process-params` #c("D04.1.1.1.1", "U03.3.1.1.2", "U04.8.1.1.4")
  opar <- formdat$`char-values`$id #c("D01.1.1.1.1",  "D01.1.1.1.5",  "D01.5.1.1.1",  "D01.5.1.1.14")#formdat$`char-values`$id #c("D03.1.1.1.5", "U04.2.1.1.5", "U04.6.1.1.9") # opar <- selected$id


  if (sum(!(ipar %in% colnames(dat[[input]])))>0) {
    errmsg <- paste("Folgende Inputparameter wurden nicht gefunden:",paste(ipar[!(ipar %in% colnames(dat[[input]]))],collapse=", "))
    jsonlite::write_json(list(error=errmsg), outjson)
    stop(errmsg)
  }
  if (sum(!(opar %in% colnames(dat[[output]])))>0) {
    errmsg <- paste("Folgende Outputparameter wurden nicht gefunden:",paste(opar[!(opar %in% colnames(dat[[output]]))],collapse=", "))
    jsonlite::write_json(list(error=errmsg), outjson)
    stop(errmsg)
  }

  idat <- dat[[input]][,ipar]
  odat <- dat[[output]][,opar]

  dat.df <- data.frame(x=idat,y=odat)
  if (dim(dat.df)[1] < 4) {
    errmsg <- "Mindestens 4 Beobachtungen werden benötigt. Wähle mehr Wärmebehandlungen aus."
    jsonlite::write_json(list(error=errmsg), outjson)
    stop(errmsg)
  }
  # get objective
  objective <- formdat$`char-values`$value #c( 2000, 59,  16.5,    -0.5 ) #
  tgerr <- rep(100,length(objective))

  if (any(is.na(dat.df))) warning(paste("NAs in Spalte(n):",paste(names(dat.df)[apply(dat.df,2,function(x) any(is.na(x)))], collapse=", "),"- observations will be removed"))
  dat.df <- na.omit(dat.df)

  rs <- mvTargetOpt::pred.solution(dat=dat.df,
                            tgmean=objective,
                            tgerr=tgerr,
                            xeps=-1,#0.001,
                            pcnr=1:length(ipar),  # oder opar???
                            additive=FALSE,
                            maxarea=NULL, ptchoice=3, useweights=TRUE, mknormweights=F, allpts=F, gr2retlimit=F,bpcenter=F,
                            mindeg=0,wfun=function(x) {(1/x)^2}, sequential=T, ptchng=F, nptc=0,tgpcawg=1,betterweights=3,yweights=F,
                            datlim=NULL,knearest=NULL,
                            tgdim=1,
                            ylast=5,sto=F,#...)
                            mdeg=1, messages=T)
  colnames(rs$newx)<-NULL
  result <- list(newX = data.frame(id=ipar, value=rs$newx[1,]) ,isPrediction=rs$isprediction)
  jsonlite::write_json(result, outjson)
  return(result)
}

# startx=stval[[1]], dstgmean=dstgmean, tgerr=tgerr, mod=MODELL123, maxit=100,
# maxarea=cbind(rep(-5,2),rep(5,2)), betterweights=3,tgdim=1,
# ylast=5,sto=F,mod.sd=0, reps=2,
# mdeg=1, messages=F,
# xeps=-1, pplot=FALSE, pcnr=c(1,2), additive=F,
# useweights=T,  mknormweights=F, ptchoice=3, allpts=F, gr2retlimit=F, sequential=T, tgpcawg=1,
#yweights=F, bpcenter=F, datlim=NULL, knearest=NULL

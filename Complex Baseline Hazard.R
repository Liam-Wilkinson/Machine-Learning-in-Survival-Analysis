rm(list=ls())
####Simulations Code
library("foreach")
library("simsurv")
library("randomForestSRC")
library("pec")
library("survival")
library("survivalsvm")
library("caret")
library("plyr")
library("dplyr")
library("glmnet")
library("reshape2")
library("data.table")
library("Hmisc")
library("glmnetUtils")
library("rstpm2")
dgm <- 1:4
simdata <- function(i, dgm, 
                    fn,
                    age=0.013,tumour2=0.489,tumour3=0.867,menopause=-0.042,highgrade=0.33,
                    nodes=0.077,progesterone=-0.0005,eostrogen=-0.00008,
                    chemotherapy=0.073,hormon=-0.006)
{
  dataframe <- data.frame(id = 1:715,
                          age = rnorm(715,55, 12.95),
                          tumoursize = sample(c(1,2,3),size=715,replace=TRUE,prob=c(0.47,0.43,0.1)),
                          menopause = rbinom(715,1,0.56),
                          highgrade = rbinom(715,1,0.73),
                          nodes = rpois(715,1),
                          progesterone = rnorm(715,29.75,9.35),
                          eostrogen = rnorm(715,42.1,7.63),
                          chemotherapy = rbinom(715,1,0.19),
                          hormon = rbinom(715,1,0.11)
  )
  for (k in 1:nrow(dataframe)) {
    if (dataframe$tumoursize[k]== 1) {
      dataframe$tumour1[k]<- 1 
    }
    else {
      dataframe$tumour1[k]<-0
    }
    if (dataframe$tumoursize[k]== 2) {
      dataframe$tumour2[k]<- 1 
    }
    else {
      dataframe$tumour2[k]<-0
    }
    if (dataframe$tumoursize[k]== 3) {
      dataframe$tumour3[k]<- 1 
    }
    else {
      dataframe$tumour3[k]<-0
    }
  }
s <- simsurv(loghazard=fn,
betas = c(age=0.013,tumour2=0.489,tumour3=0.867,menopause=-0.042,highgrade=0.33,
                                       nodes=0.077,progesterone=-0.0005,eostrogen=-0.00008,
                                       chemotherapy=0.073,hormon=-0.006),
               x = dataframe,
               maxt = 10)
  dataframe <- merge(dataframe, s)
  dataframe <- na.omit(dataframe)
  dt <- createDataPartition(dataframe$status, p = .7,
                            list = FALSE,
                            times = 1,
  )
  traindata<-dataframe[dt,]
  testdata<-dataframe[-dt,]
  sj<-system.time({
    tune<- tryCatch(tune.rfsrc(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                                 progesterone+eostrogen+
                                 chemotherapy+hormon, data=traindata,doBest=TRUE),error = function(e) paste("error"))
  })
  sj<-as.matrix(sj)
  sj<-as.data.frame(sj)
  sjp<-sj$V1[[3]]
  tune<-tryCatch(tune.rfsrc(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                              progesterone+eostrogen+
                              chemotherapy+hormon, data=traindata,doBest=TRUE),error = function(e) paste("error"))
  tunebest<-tryCatch(as.data.frame(tune$optimal),error = function(e) paste("error"))
  optnode=tryCatch(tunebest$`tune$optimal`[1],error = function(e) paste("error"))
  optsplit=tryCatch(tunebest$`tune$optimal`[2],error = function(e) paste("error"))
  r1<-tryCatch(rfsrc(Surv(eventtime, status)~age+tumoursize+menopause+highgrade+nodes+
                       progesterone+eostrogen+
                       chemotherapy+hormon,mtry=optsplit,nodesize=optnode, data=traindata,ntree=500,
                     splitrule="logrank"),error = function(e) paste("error"))
  bcvCindex  <- tryCatch(pec::cindex(r1,formula=Surv(eventtime, status)~age+tumoursize+menopause+highgrade+nodes+
                                       progesterone+eostrogen+
                                       chemotherapy+hormon,
                                     data=testdata),error = function(e) paste("error"))
  rfsrc<-tryCatch(bcvCindex$AppCindex,error = function(e) paste("error"))
  
  trainmat=subset(traindata, select=c("age","tumoursize"
                                      ,"menopause","highgrade","nodes","progesterone",
                                      "eostrogen","chemotherapy","hormon"))
  trainmar=as.matrix(trainmat)
  testmat=subset(testdata, select=c("age","tumoursize"
                                    ,"menopause","highgrade","nodes","progesterone",
                                    "eostrogen","chemotherapy","hormon"))
  testmar=as.matrix(testmat)
  sdvp<-system.time({
    tryCatch(crossvalcox=cva.glmnet(trainmar, Surv(traindata$eventtime,traindata$status), family = "cox", maxit = 1000,alpha=seq(0, 1, by=0.1))
             ,error = function(e) paste("error"))})
  sdvp<-tryCatch(as.matrix(sdvp),error = function(e) paste("error"))
  sdvp<-tryCatch(as.data.frame(sdvp),error = function(e) paste("error"))
  sfs<-tryCatch(sdvp$V1[[3]],error = function(e) paste("error"))
  
  crossvalcox=tryCatch(cva.glmnet(trainmar, Surv(traindata$eventtime,traindata$status), family = "cox", maxit = 1000,alpha=seq(0, 1, by=0.1)),error = function(e) paste("error"))
  get_alpha <- function(crossvalcox) {
    alpha <- crossvalcox$alpha
    error <- sapply(crossvalcox$modlist, function(mod) {min(mod$cvm)})
    alpha[which.min(error)]
  }
  get_model_params <- function(crossvalcox) {
    alpha <- crossvalcox$alpha
    lambdaMin <- sapply(crossvalcox$modlist, `[[`, "lambda.min")
    lambdaSE <- sapply(crossvalcox$modlist, `[[`, "lambda.1se")
    error <- sapply(crossvalcox$modlist, function(mod) {min(mod$cvm)})
    best <- which.min(error)
    data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
               lambdaSE = lambdaSE[best], error = error[best])
  }
  gmp=tryCatch(get_model_params(crossvalcox),error = function(e) paste("error"))
  alphaopt=tryCatch(gmp$alpha,error = function(e) paste("error"))
  lambdaopt=tryCatch(gmp$lambdaMin,error = function(e) paste("error"))
  fit=tryCatch(glmnet(trainmar, Surv(traindata$eventtime,traindata$status), family = "cox",alpha=alphaopt,lambda=lambdaopt,nlambda=1,trace.it=1),error = function(e) paste("error"))
  preds=tryCatch(predict(fit, newx=testmar),error = function(e) paste("error"))
  response=subset(testdata, select=c("eventtime","status"))
  y = cbind(time=testdata$eventtime,status=testdata$status)
  
  ccox <- tryCatch(Cindex(preds,y),error = function(e) paste("error"))
  
  sjrk<-system.time({tryCatch(
    survsvm<- survivalsvm(Surv(eventtime,status)~age+
                            tumoursize+menopause+highgrade+nodes+
                            progesterone+eostrogen+
                            chemotherapy+hormon, data=traindata,type="regression",
                          gamma.mu=0.1,kernel="add_kernel"),error = function(e) paste("error"))
  })
  sjrk<-tryCatch(as.matrix(sjrk),error = function(e) paste("error"))
  sjrk<-tryCatch(as.data.frame(sjrk),error = function(e) paste("error"))
  sqx<-tryCatch(sjrk$V1[[3]],error = function(e) paste("error"))
  
  survsvm<- tryCatch(survivalsvm(Surv(eventtime,status)~age+
                                   tumoursize+menopause+highgrade+nodes+
                                   progesterone+eostrogen+
                                   chemotherapy+hormon, data=traindata,type="regression",
                                 gamma.mu=0.1,kernel="add_kernel"),error = function(e) paste("error"))
  
  p3<-tryCatch(predict(survsvm,testdata),error = function(e) paste("error"))
  ph<-tryCatch(Hmisc::rcorr.cens(p3$predicted, testdata$eventtime, testdata$status),error = function(e) paste("error"))
  ph<-tryCatch(as.data.frame(ph),error = function(e) paste("error"))
  svmc<-tryCatch(ph$ph[1],error = function(e) paste("error"))
  
  spq<-system.time({tryCatch(
    pfit<-pstpm2(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                   progesterone+eostrogen+
                   chemotherapy+hormon,data=traindata),error = function(e) paste("error"))
  })
  spq<-tryCatch(as.matrix(spq),error = function(e) paste("error"))
  spq<-tryCatch(as.data.frame(spq),error = function(e) paste("error"))
  sqfx<-tryCatch(spq$V1[[3]],error = function(e) paste("error"))
  
  pfit<-tryCatch(pstpm2(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                          progesterone+eostrogen+
                          chemotherapy+hormon,data=traindata),error = function(e) paste("error"))
  pfit1<-tryCatch(pstpm2(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                           progesterone+eostrogen+
                           chemotherapy+hormon,data=traindata,sp=pfit@sp),error = function(e) paste("error"))
  pnm<-tryCatch(predict(pfit1,testdata,type="cumhaz"),error = function(e) paste("error"))
  pnm<-tryCatch(as.data.frame(pnm),error = function(e) paste("error"))
  pnx<-tryCatch(cbind(pnm,testdata$eventtime,testdata$status),error = function(e) paste("error"))
  psh<-tryCatch(Hmisc::rcorr.cens(pnm$pnm, testdata$eventtime, testdata$status),error = function(e) paste("error"))
  psh<-tryCatch(as.data.frame(psh),error = function(e) paste("error"))
  pzsc<-tryCatch(psh$psh[[1]],error = function(e) paste("error"))
  
  results<- data.frame(
    i= i,
    dgm = dgm,
    rfsrc=rfsrc,
    sjp=sjp,
    ccox=ccox,
    sfs=sfs,
    svmc=svmc,
    sqx=sqx,
    pzsc=pzsc,
    sqfx=sqfx
  )
  return(results)
}
set.seed(811)
nsim <- 100
estimates <- data.frame(matrix(ncol = 10, nrow = (nsim*4)))
x <- c("i", "dgm","rfsrc","sjp","ccox","sfs","svmc","sqx","pzsc", "sqfx")
colnames(estimates) <- x
states <- matrix(ncol = 626 , nrow = (nsim*4))
for (r in 1:nsim) {
  # 1st data-generating mechanism
  states[r, ] <- .Random.seed
  estimates[r, ]<- simdata(i = r, dgm=1,fn <- function(t, x, betas, ...)
    (-2.7 + 0.2 * t + 0.03 * t ^ 2 -0.01 * t ^ 3))
  states[(nsim+r), ] <- .Random.seed
  estimates[(nsim+r), ] <- simdata(i = r, dgm=2,fn <- function(t, x, betas, ...)
    (-1.3 + 0.2 * t + 0.03 * t ^ 2 -0.01 * t ^ 3))
  states[((2*nsim)+r), ] <- .Random.seed
  estimates[((2*nsim)+r), ] <-simdata(i = r, dgm=3,fn <- function(t, x, betas, ...)
    (-3.5 + 0.2 * t + 0.03 * t ^ 2 -0.01 * t ^ 3))
  states[((3*nsim)+r), ] <- .Random.seed
  estimates[((3*nsim)+r), ] <- simdata(i = r, dgm=4,fn <- function(t, x, betas, ...)
    (-2.5 + 0.2 * t + 0.03 * t ^ 2 -0.01 * t ^ 3))
}
saveRDS(estimates,file="D:/Simulations/Simulated Data/simdata6.Rds")
data<-readRDS("D:/Simulations/Simulated Data/simdata6.Rds")

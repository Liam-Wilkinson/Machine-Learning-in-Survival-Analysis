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
lambdas <- c((4.6)^-10,10,(8.5)^-10,5)
gammas <- c(1.53,36,1.2,1)
simdata <- function(i, dgm, 
                    lambda = (4.6)^-10, gamma = 1.53,
                    age=-0.462,tumour2=0.495,tumour3=0.793,menopause=0.249,highgrade=0.215,
                    nodes=0.086,progesterone=-0.0004,eostrogen=-0.00005,
                    chemotherapy=-0.892,hormon=0.447,ageb=0.003,lage=8.36,
                    agechem=-0.004,canchem=0.409,
                    nodchem=-0.014,progmen=-3.27,agehorm=-0.008,ler=0.022,lpr=-0.12)
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
  dataframe$ageb<-dataframe$age*dataframe$age
  dataframe$lage<-log(dataframe$age)
  dataframe$lpr<-log(dataframe$progesterone)
  dataframe$ler<-log(dataframe$eostrogen)
  dataframe$agechem<-dataframe$age*dataframe$chemotherapy
  dataframe$canchem<-dataframe$highgrade*dataframe$chemotherapy
  dataframe$nodchem<-dataframe$nodes*dataframe$chemotherapy
  dataframe$progmen<-dataframe$progesterone*dataframe$menopause
  dataframe$agehorm<-dataframe$age*dataframe$hormon
  dataframe<-na.omit(dataframe)
  
  
  s <- simsurv(lambdas = lambda, gammas = gamma,
               betas = c(age=-0.462,tumour2=0.495,tumour3=0.793,menopause=0.249,highgrade=0.215,
                         nodes=0.086,progesterone=-0.0004,eostrogen=-0.00005,
                         chemotherapy=-0.892,hormon=0.447,ageb=0.003,lage=8.36,
                         agechem=-0.004,canchem=0.409,
                         nodchem=-0.014,progmen=-3.27,agehorm=-0.008,ler=0.022,lpr=-0.12
               ),
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
  agecoef<-tryCatch(as.numeric(pfit1@coef[2]),error = function(e) paste("error"))
  tumourcoef<-tryCatch(as.numeric(pfit1@coef[3]),error = function(e) paste("error"))
  menocoef<-tryCatch(as.numeric(pfit1@coef[4]),error = function(e) paste("error"))
  gradecoef<-tryCatch(as.numeric(pfit1@coef[5]),error = function(e) paste("error"))
  nodescoef<-tryCatch(as.numeric(pfit1@coef[6]),error = function(e) paste("error"))
  progcoef<-tryCatch(as.numeric(pfit1@coef[7]),error = function(e) paste("error"))
  ercoef<-tryCatch(as.numeric(pfit1@coef[8]),error = function(e) paste("error"))
  chemocoef<-tryCatch(as.numeric(pfit1@coef[9]),error = function(e) paste("error"))
  hormoncoef<-tryCatch(as.numeric(pfit1@coef[10]),error = function(e) paste("error"))
  testdata$linearpredictors<-tryCatch(testdata$age*agecoef+testdata$tumoursize*tumourcoef+testdata$menopause*menocoef+
    testdata$highgrade*gradecoef+testdata$nodes*nodescoef+testdata$progesterone*progcoef+
    testdata$eostrogen*ercoef+testdata$chemotherapy*chemocoef+testdata$hormon*hormoncoef,error = function(e) paste("error"))
  pzsc<-tryCatch(Cindex(testdata$linearpredictors,y),error = function(e) paste("error"))

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
  estimates[r, ]<- simdata(i = r, dgm=1,lambda=(4.6)^-10,gamma=1.53)
  states[(nsim+r), ] <- .Random.seed
  estimates[(nsim+r), ] <- simdata(i = r, dgm=2,lambda=10,gamma=36)
  states[((2*nsim)+r), ] <- .Random.seed
  estimates[((2*nsim)+r), ] <-simdata(i = r, dgm=3,lambda=(8.5)^-10,gamma=1.2)
  states[((3*nsim)+r), ] <- .Random.seed
  estimates[((3*nsim)+r), ] <- simdata(i = r, dgm=4,lambda=5,gamma=1)
}
saveRDS(estimates,file="D:/Simulations/Simulated Data/simdata4.Rds")
data<-readRDS("D:/Simulations/Simulated Data/simdata4.Rds")

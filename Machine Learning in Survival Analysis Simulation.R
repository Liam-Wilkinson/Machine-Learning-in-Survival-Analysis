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
lambdas <- c(0.0034, 0.025, 0.0015, 0.015)
gammas <- c(1.36,1.2,1.36,1)
simdata <- function(i, dgm, 
                    lambda = 0.0034, gamma = 1.36,
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
  s <- simsurv(lambdas = lambda, gammas = gamma,
               betas = c(age=0.013,tumour2=0.489,tumour3=0.867,menopause=-0.042,highgrade=0.33,
                         nodes=0.077,progesterone=-0.0005,eostrogen=-0.00008,
                         chemotherapy=0.073,hormon=-0.006),
               x = dataframe,
               maxt = 10)
  dataframe <- merge(dataframe, s)
  dt <- createDataPartition(dataframe$status, p = .7,
                            list = FALSE,
                            times = 1)
  traindata<-dataframe[dt,]
  testdata<-dataframe[-dt,]
  sj<-system.time({
  tune<-tune.rfsrc(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                     progesterone+eostrogen+
                     chemotherapy+hormon, data=traindata,doBest=TRUE)
  })
  sj<-as.matrix(sj)
  sj<-as.data.frame(sj)
  sjp<-sj$V1[[3]]
  tune<-tune.rfsrc(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+nodes+
                     progesterone+eostrogen+
                     chemotherapy+hormon, data=traindata,doBest=TRUE)
  tunebest<-as.data.frame(tune$optimal)
  optnode=tunebest$`tune$optimal`[1]
  optsplit=tunebest$`tune$optimal`[2]
  r1<-rfsrc(Surv(eventtime, status) ~ age+
              tumoursize+menopause+highgrade+nodes+
              progesterone+eostrogen+
              chemotherapy+hormon,mtry=optsplit,nodesize=optnode, data=traindata,ntree=500,
            splitrule="logrank")
  bcvCindex  <- pec::cindex(r1,formula=Surv(eventtime, status) ~ age+
                              tumoursize+menopause+highgrade+nodes+
                              progesterone+eostrogen+
                              chemotherapy+hormon,
                            data=testdata,
                            exact=TRUE)
  rfsrc<-bcvCindex$AppCindex
  
  trainmat=subset(traindata, select=c("age","tumoursize"
                                      ,"menopause","highgrade","nodes","progesterone",
                                      "eostrogen","chemotherapy","hormon"))
  trainmar=as.matrix(trainmat)
  testmat=subset(testdata, select=c("age","tumoursize"
                                    ,"menopause","highgrade","nodes","progesterone",
                                    "eostrogen","chemotherapy","hormon"))
  testmar=as.matrix(testmat)
  
  sdvp<-system.time({
  crossvalcox=cva.glmnet(trainmar, Surv(traindata$eventtime,traindata$status), family = "cox", maxit = 1000,alpha=seq(0, 1, by=0.1))
  })
  sdvp<-as.matrix(sdvp)
  sdvp<-as.data.frame(sdvp)
  sfs<-sdvp$V1[[3]]

  crossvalcox=cva.glmnet(trainmar, Surv(traindata$eventtime,traindata$status), family = "cox", maxit = 1000,alpha=seq(0, 1, by=0.1))
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
  gmp=get_model_params(crossvalcox)
  alphaopt=gmp$alpha
  lambdaopt=gmp$lambdaMin
  fit=glmnet(trainmar, Surv(traindata$eventtime,traindata$status), family = "cox",alpha=alphaopt,lambda=lambdaopt,nlambda=1,trace.it=1)
  preds=predict(fit, newx=testmar)
  response=subset(testdata, select=c("eventtime","status"))
  y = cbind(time=testdata$eventtime,status=testdata$status)
 
   ccox <- Cindex(preds,y)
   
  
  sjrk<-system.time({
  survsvm<- survivalsvm(Surv(eventtime,status)~age+
                           tumoursize+menopause+highgrade+nodes+
                           progesterone+eostrogen+
                           chemotherapy+hormon, data=traindata,type="regression",
                         gamma.mu=0.1,kernel="add_kernel")
  })
  sjrk<-as.matrix(sjrk)
  sjrk<-as.data.frame(sjrk)
  sqx<-sjrk$V1[[3]]
  
  survsvm<- survivalsvm(Surv(eventtime,status)~age+
                          tumoursize+menopause+highgrade+nodes+
                          progesterone+eostrogen+
                          chemotherapy+hormon, data=traindata,type="regression",
                        gamma.mu=0.1,kernel="add_kernel")

   p3<-predict(survsvm,testdata)
   ph<-Hmisc::rcorr.cens(p3$predicted, testdata$eventtime, testdata$status)
   ph<-as.data.frame(ph)
   svmc<-ph$ph[1]
   
   spq<-system.time({
   pfit<-pstpm2(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+
                  nodes+progesterone+eostrogen+chemotherapy+hormon,data=traindata)
   })
   spq<-as.matrix(spq)
   spq<-as.data.frame(spq)
   sqfx<-spq$V1[[3]]
   
   pfit<-pstpm2(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+
                  nodes+progesterone+eostrogen+chemotherapy+hormon,data=traindata)
   pfit1<-pstpm2(Surv(eventtime,status)~age+tumoursize+menopause+highgrade+
                   nodes+progesterone+eostrogen+chemotherapy+hormon,data=traindata,sp=pfit@sp)

   agecoef<-as.numeric(pfit1@coef[2])
   tumourcoef<-as.numeric(pfit1@coef[3])
   menocoef<-as.numeric(pfit1@coef[4])
   gradecoef<-as.numeric(pfit1@coef[5])
   nodescoef<-as.numeric(pfit1@coef[6])
   progcoef<-as.numeric(pfit1@coef[7])
   ercoef<-as.numeric(pfit1@coef[8])
   chemocoef<-as.numeric(pfit1@coef[9])
   hormoncoef<-as.numeric(pfit1@coef[10])
   testdata$linearpredictors<-testdata$age*agecoef+testdata$tumoursize*tumourcoef+testdata$menopause*menocoef+
     testdata$highgrade*gradecoef+testdata$nodes*nodescoef+testdata$progesterone*progcoef+
     testdata$eostrogen*ercoef+testdata$chemotherapy*chemocoef+testdata$hormon*hormoncoef
   y = cbind(time=testdata$eventtime,status=testdata$status)
   pzsc<-Cindex(testdata$linearpredictors,y)
   
   results<- data.frame(
    i = i,
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
set.seed(88888)
B<-100
dfs <- foreach(j = dgm) %do% {
 foreach(i = 1:B) %do% {
simdata(i = i, dgm = j, lambda = lambdas[j], gamma =
              gammas[j])
} 
}

dg<-as.data.frame(dfs)
tidy_data <- melt(setDT(dg), 
                  measure.vars = patterns('i','dgm','rfsrc','ccox','svmc','pzsc','sjp','sfs','sqx','sqfx'), 
                  value.name = c('i','dgm','rfsrc','ccox','svmc','pzsc','sjp','sfs','sqx','sqfx'))
tidy_data <- tidy_data[order(i),]
saveRDS(tidy_data,file="D:/Simulations/Simulated Data/simdata1.Rds")
data2<-readRDS("D:/Simulations/Simulated Data/simdata1.Rds")
   
  
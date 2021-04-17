rm(list=ls())
library("haven")
library("readxl")
library("randomForestSRC")
library("survival")

stset <- read_excel("D:/Simulations/Simulation Analysis/stset.xlsx")
x<-c("event","time")
colnames(stset)<-x
rotterdamdata<-readRDS("D:/Simulations/Simulations Code/rotterdamdata.Rds")
rotterdamdata<-cbind(rotterdamdata,stset)

####Kaplan-Meier plot
fit1 <- survfit(Surv(time,event) ~ 1, data = rotterdamdata)
ggsurvplot(fit1,data=rotterdamdata,title="Overall Survival",xlab="Years",ylab="S(t)",censor=FALSE,break.x.by=2,break.y.by=0.1,ylim=c(0.2,1),font.tickslab=35,font.main=50,font.x=40,font.y=40,font.legend=25,legend="top",ggtheme = theme_gray(),size=2)
attach(rotterdamdata)

###RSF
rf1<- tune.rfsrc(Surv(time,event)~age+meno+size+grade+nodes+pr+er+er+hormon+chemo,data=rotterdamdata,doBest=TRUE)
if (library("akima", logical.return = TRUE)) {
  
  ## Plot RSF results
  plot.tune <- function(rf1, linear = TRUE) {
    x <- rf1$results[,1]
    y <- rf1$results[,2]
    z <- rf1$results[,3]
    so <- interp(x=x, y=y, z=z, linear = linear)
    idx <- which.min(z)
    x0 <- x[idx]
    y0 <- y[idx]
    filled.contour(x = so$x,
                   y = so$y,
                   z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette =
                     colorRampPalette(c("yellow", "red")),
                   xlab = "Average number of unique survival times in the terminal nodes",
                   ylab = "Number of covariates used to split nodes",
                   main = "OOB error for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                     points(x,y,pch=16,cex=.25)})
  }
  
  ## plot the surface
  plot.tune(rf1)
  
}

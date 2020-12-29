rm(list=ls())
library("ggplot2")
library("reshape2")
library("gridExtra")
library("plyr")
library("Rmisc")
library("asbio")
simdata1<-readRDS("D:/Simulations/Simulated Data/simdata1correct.Rds")
simdata1$dgm<-as.factor(simdata1$dgm)

####Confidence Intervals for C-index
c1<-CI(simdata1$ccox[simdata1$dgm==1],ci=0.95)
c2<-CI(simdata1$ccox[simdata1$dgm==2],ci=0.95)
c3<-CI(simdata1$ccox[simdata1$dgm==3],ci=0.95)
c4<-CI(simdata1$ccox[simdata1$dgm==4],ci=0.95)
c5<-CI(simdata1$rfsrc[simdata1$dgm==1],ci=0.95)
c6<-CI(simdata1$rfsrc[simdata1$dgm==2],ci=0.95)
c7<-CI(simdata1$rfsrc[simdata1$dgm==3],ci=0.95)
c8<-CI(simdata1$rfsrc[simdata1$dgm==4],ci=0.95)
c9<-CI(simdata1$svmc[simdata1$dgm==1],ci=0.95)
c10<-CI(simdata1$svmc[simdata1$dgm==2],ci=0.95)
c11<-CI(simdata1$svmc[simdata1$dgm==3],ci=0.95)
c12<-CI(simdata1$svmc[simdata1$dgm==4],ci=0.95)
c13<-CI(simdata1$pzsc[simdata1$dgm==1],ci=0.95)
c14<-CI(simdata1$pzsc[simdata1$dgm==2],ci=0.95)
c15<-CI(simdata1$pzsc[simdata1$dgm==3],ci=0.95)
c16<-CI(simdata1$pzsc[simdata1$dgm==4],ci=0.95)

####CI's for Runtimes
ci.median(simdata1$sfs[simdata1$dgm==1])
ci.median(simdata1$sfs[simdata1$dgm==2])
ci.median(simdata1$sfs[simdata1$dgm==3])
ci.median(simdata1$sfs[simdata1$dgm==4])
ci.median(simdata1$sjp[simdata1$dgm==1])
ci.median(simdata1$sjp[simdata1$dgm==2])
ci.median(simdata1$sjp[simdata1$dgm==3])
ci.median(simdata1$sjp[simdata1$dgm==4])
ci.median(simdata1$sqx[simdata1$dgm==1])
ci.median(simdata1$sqx[simdata1$dgm==2])
ci.median(simdata1$sqx[simdata1$dgm==3])
ci.median(simdata1$sqx[simdata1$dgm==4])
ci.median(simdata1$sqfx[simdata1$dgm==1])
ci.median(simdata1$sqfx[simdata1$dgm==2])
ci.median(simdata1$sqfx[simdata1$dgm==3])
ci.median(simdata1$sqfx[simdata1$dgm==4])

####Analysis1-Cindex
dat.m <- melt(simdata1,id.vars='dgm', measure.vars=c('ccox','rfsrc','svmc','pzsc'))
data1<-dat.m[dat.m$dgm == 1, ]  
data2<-dat.m[dat.m$dgm == 2, ]  
data3<-dat.m[dat.m$dgm == 3, ]  
data4<-dat.m[dat.m$dgm == 4, ]  


attach(data1)
g1<-ggplot(data1,aes(x=variable,y=value,color=variable)) +
  geom_boxplot() +
scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="Rotterdam Baseline Hazard", y = "Concordance Index", x="Method") +
  scale_x_discrete(breaks=c("ccox", "rfsrc", "svmc" , "pzsc"),
                   labels=c("Cox-EN", "RSF", "SVM" , "Penalized Spline")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10)) +
  coord_cartesian(ylim = c(0.4,0.8)) +
  stat_summary(fun = mean, geom = "point", shape=17, size=2)



g2<-ggplot(data2,aes(x=variable,y=value,color=variable)) +
  geom_boxplot() +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="Low Censoring (~20%)", y = "Concordance Index", x="Method") +
  scale_x_discrete(breaks=c("ccox", "rfsrc", "svmc" , "pzsc"),
                   labels=c("Cox-EN", "RSF", "SVM" , "Penalized Spline")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10)) +
  coord_cartesian(ylim = c(0.4,0.8)) +
  stat_summary(fun = mean, geom = "point", shape=17, size=2)

g3<-ggplot(data3,aes(x=variable,y=value,color=variable)) +
  geom_boxplot() +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="High Censoring (~80%)", y = "Concordance Index", x="Method") +
  scale_x_discrete(breaks=c("ccox", "rfsrc", "svmc" , "pzsc"),
                   labels=c("Cox-EN", "RSF", "SVM" , "Penalized Spline")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10)) +
  coord_cartesian(ylim = c(0.4,0.8))  +
  stat_summary(fun = mean, geom = "point", shape=17, size=2)

g4<-ggplot(data4,aes(x=variable,y=value,color=variable)) +
  geom_boxplot() +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="Exponential", y = "Concordance Index", x="Method") +
  scale_x_discrete(breaks=c("ccox", "rfsrc", "svmc" , "pzsc"),
                   labels=c("Cox-EN", "RSF", "SVM" , "Penalized Spline")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10)) +
  coord_cartesian(ylim = c(0.4,0.8))   +
  stat_summary(fun = mean, geom = "point", shape=17, size=2)

grid.arrange(g1,g2,g3,g4, ncol=2, nrow=2)

####Analysis 1 - Runtimes

dat.t <- melt(simdata1,id.vars=c('i','dgm'), measure.vars=c('sfs','sjp','sqx','sqfx'))
datat1<-dat.t[dat.t$dgm == "1", ]  
datat2<-dat.t[dat.t$dgm == "2", ]  
datat3<-dat.t[dat.t$dgm == "3", ]  
datat4<-dat.t[dat.t$dgm == "4", ] 

g5<-ggplot(datat1) +
  geom_line(aes(x=i,y=value,color=variable)) +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="Rotterdam Baseline Hazard", y = "Run time (seconds)", x="Repitition") +
theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10))

g6<-ggplot(datat2) +
  geom_line(aes(x=i,y=value,color=variable)) +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="Low Censoring (~20%)", y = "Run time (seconds)", x="Repitition") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10))

g7<-ggplot(datat3) +
  geom_line(aes(x=i,y=value,color=variable)) +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="High Censoring (~80%)", y = "Run time (seconds)", x="Repitition") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10))

g8<-ggplot(datat4) +
  geom_line(aes(x=i,y=value,color=variable)) +
  scale_color_manual(values=c("Blue", "Darkgreen", "Darkred","Darkorange"),name="Method",labels=c("Cox-EN","RSF","SVM","PS")) +
  labs(title="Exponential", y = "Run time (seconds)", x="Repitition") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title = element_text(size = 10))

grid.arrange(g5,g6,g7,g8, ncol=2, nrow=2)



####Histogram-Analysis 1
rm(list=ls())
library("ggplot2")
library("reshape2")
library("gridExtra")
library("plyr")
simdata1<-readRDS("D:/Simulations/Simulated Data/simdata1correct.Rds")
simdata1$dgm<-as.factor(simdata1$dgm)
simdata1$variable<-as.factor(simdata1$variable)
simdata1$pzsc<-simdata1$V2
simdata1$sqfx<-simdata1$V3
levels(simdata1$dgm) <- c("Rotterdam","LowCens","HighCens","Exp")
dat.m <- melt(simdata1,id.vars='dgm', measure.vars=c('ccox','rfsrc','svmc','pzsc'))
levels(dat.m$variable) <-c("ccox"="Cox-EN", "rfsrc"="RSF", "svmc"="SVM","pzsc"="PS")

####Analysis 1 - Histogram
h1<-ggplot(dat.m,aes(x=value,fill=variable))+geom_histogram()
h1+facet_grid(variable~dgm) +
  scale_fill_manual(values=c("Blue","Darkgreen","Darkred","Darkorange")) +
  labs(title="Histograms of C-index by Method and Baseline Hazard Distribution", x = "Concordance Index", y="Count") +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(size = 0))
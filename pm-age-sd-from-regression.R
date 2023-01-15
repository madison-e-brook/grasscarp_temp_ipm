library(dplyr)
library(mgcv)

#linear regression + CI####

#we need the results of those K fold regressions
aamloc<-read.csv(file.choose()) #final complete data.csv from Brook et al, 2023
str(aamloc)

DATA<-aamloc[aamloc$Species!="Black",]
DATA$Species<-as.factor(DATA$Species)

#slope and intercept from regressions
dd.slope<--0.00016
dd.int<-2.3

sub<-DATA %>% group_by(Code) %>% sample_n(size=1)

#code for confidence intervals#

#dd
x.dd<-DATA$AnnualDD
y.dd<-log(DATA$AAM)
n<-nrow(sub)

pred.x.dd<-seq(min(DATA$AnnualDD), max(DATA$AnnualDD),10)
pred.y.dd<-dd.int+(dd.slope*pred.x.dd)
y.fitted.dd<-dd.int+(dd.slope*x.dd)

sse.dd <- sum((y.dd - y.fitted.dd)^2)
mse.dd <- sse.dd / (n - 2)
t.value<-qt(0.975, n - 2) 

mean.se.fit.dd <- (1 / n + (pred.x.dd - mean(x.dd))^2 / (sum((x.dd - mean(x.dd))^2)))
mean.conf.upper.dd <- pred.y.dd + t.value * sqrt(mse.dd * mean.se.fit.dd)
mean.conf.lower.dd <- pred.y.dd - t.value * sqrt(mse.dd * mean.se.fit.dd)
plot(x.dd,y.dd)
points(pred.x.dd,mean.conf.upper.dd,type="l",col="blue")
points(pred.x.dd,mean.conf.lower.dd,type="l",col="blue")
abline(dd.int,dd.slope)

#get degree day estimates for multiple ages/size#####
annual.dd<-(log(seq(2,12,0.1))-dd.int)/dd.slope

aam.dataframe<-data.frame(AnnualDD=annual.dd,LogAge=log(seq(2,12,0.1)))

plot(aam.dataframe$AnnualDD,aam.dataframe$LogAge)

#confidence intervals
pred.x.dd<-annual.dd
pred.y.dd<-dd.int+(dd.slope*pred.x.dd)

y.fitted.dd<-dd.int+(dd.slope*x.dd)

sse.dd <- sum((y.dd - y.fitted.dd)^2)
mse.dd <- sse.dd / (n - 2)
t.value<-qt(0.975, n - 2) 

mean.se.fit.dd <- (1 / n + (pred.x.dd - mean(x.dd))^2 / (sum((x.dd - mean(x.dd))^2)))
mean.conf.upper.dd <- pred.y.dd + t.value * sqrt(mse.dd * mean.se.fit.dd)
mean.conf.lower.dd <- pred.y.dd - t.value * sqrt(mse.dd * mean.se.fit.dd)


#get our standard deviation
aam.dataframe<-data.frame(AnnualDD=annual.dd,LogAge=seq(2,12,0.1),
                          Upper=exp(mean.conf.upper.dd),
                          Lower=exp(mean.conf.lower.dd))

plot(aam.dataframe$AnnualDD,aam.dataframe$LogAge,type="l",ylim=c(0,12))
points(aam.dataframe$AnnualDD,aam.dataframe$Upper,type="l")
points(aam.dataframe$AnnualDD,aam.dataframe$Lower,type="l")

range<-aam.dataframe$Upper-aam.dataframe$Lower

aam.dataframe<-data.frame(AnnualDD=annual.dd,LogAge=seq(2,12,0.1),
                          Upper=exp(mean.conf.upper.dd),
                          Lower=exp(mean.conf.lower.dd),
                          Range=range)


#in a true normal distribution, 68% of the data falls in abot 1 standard devation
#and 95% falls within 1.96 of the standard dev

standard.devs<-aam.dataframe$Range/(2*1.96) #these are our standard deviations

aam.dataframe<-data.frame(AnnualDD=annual.dd,LogAge=seq(2,12,0.1),
                          Upper=exp(mean.conf.upper.dd),
                          Lower=exp(mean.conf.lower.dd),
                          Range=range, Stand.Dev.age=standard.devs)

plot(aam.dataframe$AnnualDD,aam.dataframe$Stand.Dev)

write.csv(aam.dataframe,"parameter-relationship-with-dd.csv",row.names = F)


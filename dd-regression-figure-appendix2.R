library(colorBlindness)
library(sp)
library(raster)
library(maps)
library(maptools)
library(viridis)
library(RColorBrewer)
library(dplyr)

#final data from Brook et al., 2023
aamloc<-read.csv(file.choose())
str(aamloc)

DATA<-aamloc[aamloc$Species!="Black",]
DATA$Species<-as.factor(DATA$Species)

#the regression parameters
dd.slope<--0.00016
dd.int<-2.3

sub<-DATA %>% group_by(Code) %>% sample_n(size=1)


#dd confidence intervals
x.dd<-DATA$AnnualDD
y.dd<-log(DATA$AAM)
n<-nrow(sub)

pred.x.dd<-seq(range(DATA$AnnualDD)[1], range(DATA$AnnualDD)[2],10)
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

table(DATA$Code)[table(DATA$Code)>1]
uniques<-DATA[DATA$Code!="A" & DATA$Code!="AF"& DATA$Code!="AG"& DATA$Code!="AL"& DATA$Code!="AM" & DATA$Code!="B"& DATA$Code!="C" & DATA$Code!="D"& DATA$Code!="E"& DATA$Code!="F"& DATA$Code!="G"& DATA$Code!="H"& DATA$Code!="J" & DATA$Code!="M"& DATA$Code!="N"& DATA$Code!="O"& DATA$Code!="S"& DATA$Code!="Z",]
duplicates<-DATA[DATA$Code=="A" | DATA$Code=="AF"| DATA$Code=="AG"| DATA$Code=="AL"| DATA$Code=="AM" | DATA$Code=="B"| DATA$Code=="C"| DATA$Code=="D"|DATA$Code=="E"| DATA$Code=="F"| DATA$Code=="G"| DATA$Code=="H"| DATA$Code=="J"| DATA$Code=="M"| DATA$Code=="N"| DATA$Code=="O"| DATA$Code=="S"| DATA$Code=="Z",]

#figure#####
setwd("~/Documents/Thesis/good code")
png("appendix 3.2 figure 1.png", width= 2404, height= 1600, units="px", res = 300)
par(mar=c(6,4,1,3))
plot(uniques$AnnualDD,log(uniques$AAM),pch=c(0, 1, 2)[as.numeric(uniques$Species)],
     xlab="",
     ylab="",las=1,xlim=c(range(DATA$AnnualDD)[1],range(DATA$AnnualDD)[2]),
     ylim=c(range(log(DATA$AAM))[1],range(log(DATA$AAM))[2]),col=Blue2DarkOrange18Steps[18])
points(duplicates$AnnualDD,log(duplicates$AAM),pch=c(0, 1, 2)[as.numeric(duplicates$Species)],
       col=Blue2DarkOrange18Steps[2])
title(ylab=expression("ln Age at Maturity"),mgp=c(2.5,1,0))
title(xlab=expression("Annual Average Degree Days"),mgp=c(2.5,1,0))

points(pred.x.dd,mean.conf.upper.dd,type="l",lty=2)
points(pred.x.dd,mean.conf.lower.dd,type="l",lty=2)
curve(dd.int+dd.slope*x,range(DATA$AnnualDD)[1],range(DATA$AnnualDD)[2],add=T)
text(8000,2,expression(italic("y = 2.3 - 0.00016 * degree days")))


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",
       legend = c(levels(DATA$Species),"Unique","Subsampled"),
       col = c("black","black","black",Blue2DarkOrange18Steps[18],Blue2DarkOrange18Steps[2]),
       pch=c(0, 1, 2, 15, 15), xpd = TRUE, horiz = T, cex = 0.9, bty = 'n')

dev.off()

library(colorBlindness)
library(viridis)
library(RColorBrewer)

setwd("~/Documents/Thesis/good code")
surv.2330<-read.csv("surv range at 2330dd.csv")
surv.3560<-read.csv("surv range at 3560dd.csv")
surv.4280<-read.csv("surv range at 4280dd.csv")

#create labels for plot
x.label<-expression("Maximum Adult Survival ("*italic("max"["s"])*")")
y.label<-expression(lambda)

png("different max survival plot.png", width= 2404, height= 1600, units="px", res = 300)

plot(surv.2330$surv.max,surv.2330$lambdas,ylim=c(0.7,1.5),
     col=viridis(4)[1],type="l",
     xlab=x.label,ylab=y.label)
abline(h=1,lty=2)
points(surv.3560$surv.max,surv.3560$lambdas,
       col=viridis(4)[2],type="l")
points(surv.4280$surv.max,surv.4280$lambdas,
       col=viridis(4)[3],type="l")
legend("topleft",c("Superior (DD = 2330)",
                   "Erie (DD = 3560)",
                   "Mississippi (DD = 4280)"),
       lty=1,col=c(viridis(4)[1],
                    viridis(4)[2],
                    viridis(4)[3]),bty="n")
dev.off()


#test perccentage change at different degree days
one<-lm(surv.2330$lambdas~surv.2330$surv.max)
one

value1<-one$coef[1]+one$coef[2]*0.5
value2<-one$coef[1]+one$coef[2]*0.6

((value2-value1)/value1)*100

one<-lm(surv.3560$lambdas~surv.3560$surv.max)
one

value1<-one$coef[1]+one$coef[2]*0.5
value2<-one$coef[1]+one$coef[2]*0.6

((value2-value1)/value1)*100

one<-lm(surv.4280$lambdas~surv.4280$surv.max)
one
value1<-one$coef[1]+one$coef[2]*0.5
value2<-one$coef[1]+one$coef[2]*0.6

((value2-value1)/value1)*100



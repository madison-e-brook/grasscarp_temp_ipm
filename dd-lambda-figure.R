library(colorBlindness)
library(viridis)
library(RColorBrewer)

setwd("~/Documents/Thesis/good code")
#download csv of lambva values across temperature at two survivals
surv.09<-read.csv("lambda data survival 0.9.csv")
surv.09

surv.06<-read.csv("lambda data survival 0.6.csv")
surv.06

y.label<-expression(lambda)
png("aam lambda graph both survivals.png", width= 2404, height= 1600, units="px", res = 300)

plot(surv.09$DD,surv.09$default.lambda,pch=16,
     ylim=c(0.75,1.5),
     xlab="Annual Degree Days",
     ylab=y.label,
     col=Blue2DarkRed18Steps[1])
points(surv.06$DD,surv.06$default.lambda,pch=16,
       ylim=c(0.75,1.4),
       col=Blue2DarkRed18Steps[17])
points(surv.09$DD,surv.09$CIupper,type="l",
       col=Blue2DarkRed18Steps[1])
points(surv.09$DD,surv.09$CIlower,type="l",
       col=Blue2DarkRed18Steps[1])
points(surv.06$DD,surv.06$CIupper,type="l",
       col=Blue2DarkRed18Steps[17])
points(surv.06$DD,surv.06$CIlower,type="l",
       col=Blue2DarkRed18Steps[17])
text(2330,0.8,"Superior\nDD = 2330")
segments(2330,0.72,2330,0.75)
text(3560,0.8,"Erie\nDD = 3560")
segments(3560,0.72,3560,0.75)
text(4280,0.8,"Mississippi\nDD = 4280")
segments(4280,0.72,4280,0.75)
abline(h=1,lty=2)
legend("topleft",c("Max adult survival = 0.9",
                   "Max adult survival = 0.6"),
       col=c(Blue2DarkRed18Steps[1],Blue2DarkRed18Steps[18]),pch=16,
       bty="n")
dev.off()

#approximate percents####
#estimates of how much changing dd affects lambda
one<-lm(surv.09$default.lambda~surv.09$DD)
two<-lm(surv.06$default.lambda~surv.06$DD)

value1<-two$coef[1]+two$coef[2]*2000
value2<-two$coef[1]+two$coef[2]*3000

((value2-value1)/value1)*100

#changing the dd by 1000 at survival of 0.6 changes the lambda by about 2.8%
#changing the dd by 1000 at survival of 0.9 changes the lambda by about 2.8%



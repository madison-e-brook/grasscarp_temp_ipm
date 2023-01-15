eggs<-read.csv("fecundities with degree days.csv")
str(eggs)

png("fecundity boxplot.png", width= 2404, height= 1600, units="px", res = 300)
boxplot(eggs$Fecundity,ylab="Relative Fecundity (eggs/kg)")
dev.off()

#exclude 1st point for being an outlier
fecund.reg<-lm(Fecundity~AnnualDD,data=eggs[-1,])
fecund.int<-fecund.reg$coef[1]
fecund.slope<-fecund.reg$coef[2]
newdata.me<-data.frame(AnnualDD=seq(0,11000,1))
tests<-predict(fecund.reg,newdata=newdata.me,interval="confidence")
tests<-as.data.frame(tests)
summ<-summary(fecund.reg)

png("fecundity regression.png", width= 2404, height= 1600, units="px", res = 300)
plot(eggs$AnnualDD[-1],eggs$Fecundity[-1],
     xlab="Annual Average Degree Days",
     ylab="Relative Fecundity (eggs/kg)")

abline(fecund.reg)
points(seq(0,11000,1),tests$lwr,type="l",lty=3)
points(seq(0,11000,1),tests$upr,type="l",lty=3)

text(9000,130000,expression("Adj R"^2*" = 0.50"))
text(9000,120000,expression("p value = 0.031"))
text(9000,110000,expression("y = 130907 - 6.7x"))
dev.off()
#multiplied by degree days
#round to a whole number of eggs
degday<-5000
fecund<-round(fecund.int+fecund.slope*degday,0)


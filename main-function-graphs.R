library(colorBlindness)
library(viridis)
library(RColorBrewer)

displayAvailablePalette(color="white")
displayAllColors(Blue2DarkRed18Steps, color="white")

displayAllColors(c(Blue2DarkRed18Steps[1],
                   Blue2DarkRed18Steps[17]), color="white")



library(nlstools)

setwd("~/Documents/Population modelling/Data")

#growth and von bert####
fish.data<-read.csv("sullivan data raw.csv")
str(fish.data)

plot(fish.data$Age,fish.data$Length..mm.)

#anchor plot with t0 at 0
fish.data<-read.csv("sullivan data raw.csv")
str(fish.data)
plot(fish.data$Age,fish.data$Length..mm.)
von.bert.0<-nls(Length..mm. ~ a*(1-exp(-b*(Age-0))), 
                data = fish.data,
                start = list(a = 850, b = 0.36))
params.0<-von.bert.0$m$getPars()
curve(params.0[1]*(1-exp(-params.0[2]*(x-0))),0,15,
      ylim=c(200,1000),add=T,col="blue")
L_inf<-round(params.0[1],0)
k<-round(params.0[2],2)
#get growth increment
ages<-seq(0,15,0.5)
sizes<-L_inf*(1-exp(-k*(ages)))
plot(ages,sizes)
initial<-sizes[1:length(sizes)-1]
after<-sizes[-1]
change<-after-initial
#plot of growth increment
plot(initial,change,type="l")

curve((1-k)*x+k*L_inf,0,1500)

#length weight relationship####
#calculate from sullivan et al
fish.data<-read.csv("sullivan data raw.csv")
plot(log10(fish.data$Length..mm.),log10(fish.data$Weight..kg.*1000))
test<-lm(log10(fish.data$Weight..kg.*1000)~log10(fish.data$Length..mm.))
abline(test,col='red')

#other values taken from Wanner & Klumb, 2009
length.int<-round(mean(c(-4.59,-4.33,round(test$coefficients[1],2))),2)
length.slope<-round(mean(c(2.87,2.77,round(test$coefficients[2],2))),2)

#get weight of fish (in kg)
(10^(length.int+length.slope*log10(600)))/1000

curve((10^(length.int+length.slope*log10(x)))/1000,0,1500)


#fecundity####
eggs<-read.csv("fecundities with degree days.csv")
str(eggs)
boxplot(eggs$Fecundity)
#exclude 1st point for being an outlier
plot(eggs$AnnualDD[-1],eggs$Fecundity[-1])
fecund.reg<-lm(eggs$Fecundity[-1]~eggs$AnnualDD[-1])
fecund.int<-fecund.reg$coef[1]
fecund.slope<-fecund.reg$coef[2]
abline(fecund.reg)
shapiro.test(fecund.reg$residuals)

#multiplied by degree days
#round to a whole number of eggs
#(keep large number of sig dig because this is obtained by a relationship)
fecund.2330<-round(fecund.int+fecund.slope*2330,0)
fecund.3560<-round(fecund.int+fecund.slope*3560,0)
fecund.4280<-round(fecund.int+fecund.slope*4280,0)

curve((fecund.2330*10^(length.int+length.slope*log10(x)))/1000,0,1000)
curve((fecund.3560*10^(length.int+length.slope*log10(x)))/1000,0,1000,add=T)
curve((fecund.4280*10^(length.int+length.slope*log10(x)))/1000,0,1000,add=T)

#survival#####
surv.inflection<-1
surv.slope<--2
surv.min<-0.1
surv.max<-0.9

curve(((surv.max-surv.min)/(1+exp(surv.slope*(x-surv.inflection)))+surv.min),-1,15)


#probability of reproduction####
#relationship with parameters and dd
dd.para<-read.csv("degree day relationship with parameters.csv")


slope.A.lm<-lm(log(-dd.para$slope_beta)~ dd.para$degree_day)
slope.A<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*3560)

int.A.lm<-lm(log(dd.para$int_alpha)~ dd.para$degree_day)
int.A<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*3560)


curve(1/(1+exp(slope.A*(x-int.A))),0,20,col="blue")


#get actual parameters####

slope.2330<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*2330)
int.2330<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*2330)

slope.3560<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*3560)
int.3560<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*3560)

slope.4280<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*4280)
int.4280<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*4280)


curve(1/(1+exp(slope.2330*(x-int.2330))),-1,15,col="black")
curve(1/(1+exp(slope.3560*(x-int.3560))),-1,15,col="red",add=T)
curve(1/(1+exp(slope.4280*(x-int.4280))),-1,15,col="blue",add=T)

#make the graphs####

setwd("~/Documents/Thesis/good code")
png("model graphs.png", width= 2404, height= 1600, units="px", res = 300)
par(mfrow=c(2,2),mar=c(6,4,1,3))

curve(1/(1+exp(slope.2330*(x-int.2330))),-1,15,col=viridis(4)[1],
      xlab="Age (years)",ylab="Probability of maturing",
      lwd=1,ylim=c(0,1.1))
curve(1/(1+exp(slope.3560*(x-int.3560))),-1,15,col=viridis(4)[2],add=T,lwd=1)
curve(1/(1+exp(slope.4280*(x-int.4280))),-1,15,col=viridis(4)[3],add=T,lwd=1)
text(15,1.07,"a",cex=1)
legend("bottomright",c("Superior (DD = 2330)",
                       "Erie (DD = 3560)",
                       "Mississippi (DD = 4280)"),
       lty=1,col=c(viridis(4)[1],
                   viridis(4)[2],
                   viridis(4)[3]),bty="n",cex=0.7)


curve(((fecund.2330*10^(length.int+length.slope*log10(x)))/1000),0,1000,
      xlab="Length (mm)",ylab="Total number of eggs",col=viridis(4)[1],
      ylim=c(0,1300000))
curve((fecund.3560*10^(length.int+length.slope*log10(x)))/1000,0,1000,add=T,
      col=viridis(4)[2])
curve((fecund.4280*10^(length.int+length.slope*log10(x)))/1000,0,1000,add=T,
      col= viridis(4)[3])
legend("topleft",c("Superior (DD = 2330)",
                   "Erie (DD = 3560)",
                   "Mississippi (DD = 4280)"),
       lty=1,col=c(viridis(4)[1],
                   viridis(4)[2],
                   viridis(4)[3]),bty="n",cex=0.7)
text(1000,1280000,"b",cex=1)

curve(((surv.max-surv.min)/(1+exp(surv.slope*(x-surv.inflection)))+surv.min),-1,15,
      xlab="Age (years)",ylab="Probability of survival",
      ylim=c(0,1.1),lty=1)
curve(((0.6-surv.min)/(1+exp(surv.slope*(x-surv.inflection)))+surv.min),-1,15,add=T,
      lty=2)
legend("bottomright",c("Maximum Adult Survival = 0.9",
                       "Maximum Adult Survival = 0.6"),
       lty=c(1,2),bty="n",cex=0.8)
text(15,1.07,"c",cex=1)

me<-expression("Growth increment ("*italic("g"["inc"])*") (mm)")
plot(initial,change,type="l",xlab="Length (mm)",ylab=me,xlim=c(0,1000))
text(1000,159,"d",cex=1)
segments(813,0.1,1000,0.1)
dev.off()

expression(italics("max"["s"])*" = 0.9")

expression(italic("max"["s"])*" = 0.9")
expression(italic("max"["s"])*" = 0.6")

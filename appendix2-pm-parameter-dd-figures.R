#graph of age at maturity distribution#####
setwd("~/Documents/Population modelling/Data")
aam.dataframe<-read.csv("parameter relationship with dd.csv")
colnames(aam.dataframe)[2]<-"Age"
aam.dataframe<-aam.dataframe[aam.dataframe$AnnualDD>0,]

age.2<-which(aam.dataframe$Age==2)
age.5<-which(aam.dataframe$Age==5)
age.8<-which(aam.dataframe$Age==8)

#generate ages
age.crit.2<-rnorm(50000,aam.dataframe$Age[age.2],aam.dataframe$Stand.Dev.age[age.2])
age.crit.5<-rnorm(50000,aam.dataframe$Age[age.5],aam.dataframe$Stand.Dev.age[age.5])
age.crit.8<-rnorm(50000,aam.dataframe$Age[age.8],aam.dataframe$Stand.Dev.age[age.8])

#aam distribution graph####
setwd("~/Documents/Thesis/good code")
png("aam distribution graph.png", width= 2404, height= 1000, units="px", res = 300)

par(mfrow=c(1,3),cex.main=2)
hist(age.crit.2,
     xlab="Individual Age at Maturity",
     main="A")
hist(age.crit.5,
     xlab="Individual Age at Maturity",
     main="B")
hist(age.crit.8,
     xlab="Individual Age at Maturity",
     main="C")

dev.off()

#fish population simulation graph####
setwd("~/Documents/Thesis/good code")
png("fish pop distribution graph.png", width= 2404, height= 1600, units="px", res = 300)

fish<-rnorm(50000,mean=6,sd=4)
fish[fish<0]<-NA
hist(fish,
     xlab="Fish Age (years)",
     main="")
dev.off()

#which fish are older than their age at maturity
mature.2<-fish>age.crit.2
mature.5<-fish>age.crit.5
mature.8<-fish>age.crit.8


#fit functions####
#make a glm at each age
gm.2<-glm(mature.2~fish,family="binomial")
slope.2<-gm.2$coef[2]
int.2<-gm.2$coef[1]

gm.5<-glm(mature.5~fish,family="binomial")
slope.5<-gm.5$coef[2]
int.5<-gm.5$coef[1]

gm.8<-glm(mature.8~fish,family="binomial")
slope.8<-gm.8$coef[2]
int.8<-gm.8$coef[1]


setwd("~/Documents/Thesis/good code")
png("maturity distributions with fit functions.png", width= 2404, height= 1200, units="px", res = 300)
par(mfrow=c(1,3))
plot(fish,mature.2,
     xlab="Age of Fish (year)",
     ylab="Probability of Maturity")
text(20,0.1,"A",cex=2)
curve(((1)/(1+exp(-(int.2+slope.2*x)))),0,25,add=T,lwd=2,col="red")


plot(fish,mature.5,
     xlab="Age of Fish (year)",
     ylab="Probability of Maturity")
text(20,0.1,"B",cex=2)
curve(((1)/(1+exp(-(int.5+slope.5*x)))),0,25,add=T,lwd=2,col="red")

plot(fish,mature.8,
     xlab="Age of Fish (year)",
     ylab="Probability of Maturity")
text(20,0.1,"C",cex=2)
curve(((1)/(1+exp(-(int.8+slope.8*x)))),0,25,add=T,lwd=2,col="red")
dev.off()

#parameter relationship with degree days####
slopes.length.A<-NULL
intercepts.length.A<-NULL

for(i in 1:nrow(aam.dataframe)){
  age.crit<-rnorm(50000,aam.dataframe$Age[i],aam.dataframe$Stand.Dev.age[i])
  
  maturity<-data.frame(fish=fish,
                       age.crit=age.crit,
                       mature=fish>age.crit)
  m.model<-glm(mature~fish,data=maturity,family="binomial")
  new.slope<-m.model$coef[2]
  new.int<-m.model$coef[1]
  slopes.length.A<-c(slopes.length.A,new.slope)
  intercepts.length.A<-c(intercepts.length.A,new.int)
}

length(slopes.length.A)
#convert parameters to match equation
real.slopes<--slopes.length.A
real.ints<-intercepts.length.A/-slopes.length.A

#plots####
setwd("~/Documents/Thesis/good code")
png("parameter relationships with all degree days.png", width= 2404, height= 1200, units="px", res = 300)
par(mfrow=c(1,2))

plot(aam.dataframe$AnnualDD,real.ints,
     xlab="Annual Degree Days",
     ylab=expression(alpha))
text(9500,9.5,"A",cex=2)
plot(aam.dataframe$AnnualDD,real.slopes,
     xlab="Annual Degree Days",
     ylab=expression(beta))
text(9500,-1.7,"B",cex=2)

dev.off()


#only for less than 5000####

#clean the outputs (need dd<5000)
degree.clean<-aam.dataframe$AnnualDD[which.max(aam.dataframe$AnnualDD<=5000):length(aam.dataframe$AnnualDD)]
slopes.clean.A<-real.slopes[which.max(aam.dataframe$AnnualDD<=5000):length(aam.dataframe$AnnualDD)]
intercept.clean.A<-real.ints[which.max(aam.dataframe$AnnualDD<=5000):length(aam.dataframe$AnnualDD)]

#graphs####
setwd("~/Documents/Thesis/good code")
png("parameter relationships with less than 5000 dd.png", width= 2404, height= 1200, units="px", res = 300)
par(mfrow=c(1,2))
plot(degree.clean,log(intercept.clean.A),
     xlab="Annual Degree Days",
     ylab=expression(ln(alpha)))
abline(lm(log(intercept.clean.A)~degree.clean))
text(4800,1.9,"A",cex=2)

plot(degree.clean,log(-slopes.clean.A),
     xlab="Annual Degree Days",
     ylab=expression(ln(-beta)))
abline(lm(log(-slopes.clean.A)~degree.clean))
text(4800,1.1,"B",cex=2)
dev.off()




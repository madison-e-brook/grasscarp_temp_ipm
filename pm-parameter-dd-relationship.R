#simulate population####

fish<-rnorm(50000,mean=6,sd=4)
fish[fish<0]<-NA

#get the ages and standard devations

aam.dataframe<-read.csv("parameter-relationship-with-dd.csv")

#make empty vectors
slopes.length.A<-NULL
intercepts.length.A<-NULL

#for each age, generate the slope and intercept that describes the relationship
#of when the fish mature
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
#do equation conversion of parameters
real.slopes<--slopes.length.A
real.ints<-intercepts.length.A/-slopes.length.A

plot(aam.dataframe$AnnualDD,real.slopes)
plot(aam.dataframe$AnnualDD,real.ints)


#clean the outputs (need dd<5000)
degree.clean<-aam.dataframe$AnnualDD[which.max(aam.dataframe$AnnualDD<=5000):length(aam.dataframe$AnnualDD)]
slopes.clean.A<-real.slopes[which.max(aam.dataframe$AnnualDD<=5000):length(aam.dataframe$AnnualDD)]
intercept.clean.A<-real.ints[which.max(aam.dataframe$AnnualDD<=5000):length(aam.dataframe$AnnualDD)]

plot(degree.clean,log(-slopes.clean.A))

plot(degree.clean,log(intercept.clean.A))

#get the relationships
slope.A.lm<-lm(log(-slopes.clean.A)~degree.clean)
slope.A<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*3000)

int.A.lm<-lm(log(intercept.clean.A)~degree.clean)
int.A<-(exp(int.A.lm$coef[1]+int.A.lm$coef[2]*3000))


curve(1/(1+exp(real.slope*(x-real.int))),0,20,col="blue")


dd.parameter<-data.frame(degree_day=degree.clean,
                         slope_beta=slopes.clean.A,
                         int_alpha=intercept.clean.A)



setwd("~/Documents/Population modelling/Data")
write.csv(dd.parameter,"degree day relationship with parameters.csv",
          row.names = F)


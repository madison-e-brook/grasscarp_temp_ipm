#default model####
library(nlstools)
library(ipmr)
library(lhs)
setwd("~/Documents/Population modelling/Data")

#maximum survival of 0.9!!!!!#####
degday<-3000

#growth model####
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
plot(initial,change)
reg<-lm(change~initial)

max.growth.increment<-round(reg$coef[1],0)
max.size<-round(params.0[1],0)
curve(max.growth.increment+(-max.growth.increment/max.size)*x,col="red",
      add=T)
#growth increment function
growth.increment<-function(length_1,min.growth=0.1,max.growth.increment=161,max.size=813){
  incre<-max.growth.increment+(-max.growth.increment/max.size)*length_1
  for(i in 1:length(length_1)){
    if(incre[i]<=min.growth){
      incre[i]<-min.growth
    }
  }
  return(incre)
}
#standard deviation function
standard.dev<-function(length_1,min.std=5,max.growth.increment=161,max.size=813,division=4){
  incre<-max.growth.increment+(-max.growth.increment/max.size)*length_1
  std<-incre/division
  for(i in 1:length(length_1)){
    if(std[i]<=min.std){
      std[i]<-min.std
    }
  }
  return(std)
}
min.growth=0.1
min.std=5
division=4

#survival#####
surv.inflection<-1
surv.slope<--2
surv.min<-0.1
surv.max<-0.9

curve(((surv.max-surv.min)/(1+exp(surv.slope*(x-surv.inflection)))+surv.min),-3,20)

#length weight relationship####
#calculate from sullivan et al
fish.data<-read.csv("sullivan data raw.csv")
plot(log10(fish.data$Length..mm.),log10(fish.data$Weight..kg.*1000))
test<-lm(log10(fish.data$Weight..kg.*1000)~log10(fish.data$Length..mm.))
abline(test,col='red')

#other values taken from Wanner & Klumb, 2009
length.int<-round(mean(c(-4.59,-4.33,round(test$coefficients[1],2))),2)
length.slope<-round(mean(c(2.87,2.77,round(test$coefficients[2],2))),2)
abline(a=length.int,b=length.slope)
#get weight of fish (in kg)
(10^(length.int+length.slope*log10(600)))/1000
curve((10^(length.int+length.slope*log10(x)))/1000,0,1000)

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
#multiplied by degree days
#round to a whole number of eggs
#(keep large number of sig dig because this is obtained by a relationship)
fecund<-round(fecund.int+fecund.slope*degday,0)

de<-seq(0,5000,500)
fe<-signif(fecund.int+fecund.slope*de,2)
plot(de,fe)

#probability of reproducing####
#relationship with parameters and dd
dd.para<-read.csv("degree day relationship with parameters.csv")
plot(dd.para$degree_day,log(-dd.para$slope_beta))
plot(dd.para$degree_day,log(dd.para$int_alpha))

slope.A.lm<-lm(log(-dd.para$slope_beta)~dd.para$degree_day)
aam.slope<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*degday)

int.A.lm<-lm(log(dd.para$int_alpha)~dd.para$degree_day)
aam.int<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*degday)

curve(1/(1+exp(aam.slope*(x-aam.int))),0,20,col="blue")
abline(v=aam.int)

#percent eggs hatched####
percent.hact<-0.65

#Hatchling survival####
hatch.surv<-0.0002

#Hatchling size####
young.mean<-52
young.sd<-21

#model constant list####
#constant list#####
constant_list <- list(
  surv.inflection=surv.inflection,
  surv.slope=surv.slope,
  surv.min=surv.min,
  surv.max=surv.max,
  min.growth=min.growth,
  max.growth.increment=max.growth.increment,
  max.size=max.size,
  min.std=min.std,
  division=division,
  length.int=length.int,
  length.slope=length.slope,
  fecund=fecund,
  aam.int=aam.int,
  aam.slope=aam.slope,
  percent.hact= percent.hact,
  hatch.surv=hatch.surv,
  young.mean=young.mean,
  young.sd=young.sd
)


#model code#####
age_size_ipm <- init_ipm(sim_gen    = "general",
                         di_dd      = "di", 
                         det_stoch  = "det",
                         uses_age    = TRUE)

age_size_ipm<-define_kernel(
  proto_ipm = age_size_ipm,
  name          = "P_age",#define what IPM we're putting this kernel in
  family        = "CC",#the variable in this kernel starts and ends Continuous
  formula       = s_age * g * d_length, #the formula for the kernel
  #define formulas
  s_age         = ((surv.max-surv.min)/(1+exp(surv.slope*(age-surv.inflection)))+surv.min),
  grow_mean     = length_1+growth.increment(length_1,
                                            min.growth=min.growth,
                                            max.growth.increment=max.growth.increment,
                                            max.size=max.size),
  growth.sd     = standard.dev(length_1,min.std=min.std,
                               max.growth.increment=max.growth.increment,
                               max.size=max.size,
                               division=division),
  g             = dnorm(length_2, grow_mean, growth.sd),
  #the state variable is length
  states        = list(c("length")),
  data_list     = constant_list,
  #we don't have any parameters that vary by year
  uses_par_sets = FALSE,
  age_indices   = list(age = c(0:20)),
  evict_cor     = FALSE)

age_size_ipm<-define_kernel( 
  proto_ipm = age_size_ipm,
  name          = "F_age",
  family        = "CC",
  formula       =  (b/2) * pr_age * hp * hs * rs * d_length,
  b             = fecund*((10^(length.int+length.slope*log10(length_1)))/1000),
  pr_age        = 1/(1+exp(aam.slope*(age-aam.int))),
  hp            = percent.hact,
  hs            = hatch.surv,
  rs            = dnorm(length_2, young.mean, young.sd),
  states        = list(c("length")), #same variable
  data_list     = constant_list,
  uses_par_sets = FALSE,
  age_indices   = list(age = c(0:20)),
  evict_cor     = FALSE)

age_size_ipm<-define_impl(
  proto_ipm = age_size_ipm,
  make_impl_args_list(
    kernel_names = c("P_age", "F_age"), 
    int_rule     = rep("midpoint", 2),
    state_start    = c("length_age", "length_age"),
    state_end      = c("length_age", "length_0")
  )
)

age_size_ipm<-define_domains(
  proto_ipm = age_size_ipm,#define the upper and lower ranges, and how many midpoints
  length = c(0,
             2000,
             200))

age_size_ipm<-define_pop_state(proto_ipm = age_size_ipm,
                               n_length_age = rep(1/200,200))


age_size_ipm<-make_ipm(proto_ipm=age_size_ipm, iterations = 200,
                       usr_funs = list(growth.increment = growth.increment,
                                       standard.dev = standard.dev))


lam <- lambda(age_size_ipm)
lam

#Get errors around values#####
param.list<-c(surv.inflection,surv.slope,
              surv.min,surv.max,min.growth,
              max.growth.increment,max.size,
              min.std,division,length.int,
              length.slope,fecund,aam.int,aam.slope,
              percent.hact,hatch.surv,young.mean,young.sd)

param.var<-abs(param.list*0.2)

#max survival#####
#mean values of 0.6 and 0.9 from Sullivan et al

upper.06<-1-0.31
lower.06<-1-0.49

upper.09<-1-0.04
lower.09<-1-0.16

#get standard devs
sd.06<-(upper.06-lower.06)/(2*1.96)
sd.09<-(upper.09-lower.09)/(2*1.96)
#max size#####
setwd("~/Documents/Population modelling/Data")
fish.data<-read.csv("sullivan data raw.csv")
str(fish.data)
plot(fish.data$Age,fish.data$Length..mm.)
von.bert.0<-nls(Length..mm. ~ a*(1-exp(-b*(Age-0))), 
                data = fish.data,
                start = list(a = 850, b = 0.36))
params.0<-von.bert.0$m$getPars()

L_inf<-round(params.0[1],0)
k<-round(params.0[2],2)

summary(von.bert.0)

confint(von.bert.0)[1,2]

#sd for max size
sd.max.size<-(confint(von.bert.0)[1,2]-confint(von.bert.0)[1,1])/(2*1.96)
#length weight relationship####
fish.data<-read.csv("sullivan data raw.csv")
plot(log10(fish.data$Length..mm.),log10(fish.data$Weight..kg.*1000))
test<-lm(log10(fish.data$Weight..kg.*1000)~log10(fish.data$Length..mm.))
abline(test,col='red')

#other values taken from Wanner & Klumb, 2009
length.int<-round(mean(c(-4.59,-4.33,round(test$coefficients[1],2))),2)
length.slope<-round(mean(c(2.87,2.77,round(test$coefficients[2],2))),2)

sd.length.int<-sd(c(-4.59,-4.33,round(test$coefficients[1],2)))
sd.length.slope<-sd(c(2.87,2.77,round(test$coefficients[2],2)))
#fecundity####
eggs<-read.csv("fecundities with degree days.csv")
str(eggs)
fecund.reg<-lm(Fecundity~AnnualDD,data=eggs[-1,])
newdata.me<-data.frame(AnnualDD=seq(0,5000,1))
tests<-predict(fecund.reg,newdata=newdata.me,interval="confidence")
tests<-as.data.frame(tests)

standards<-(tests$upr-tests$lwr)/(2*1.96)
fecuncity.data<-data.frame(newdata.me$AnnualDD,tests,standards)
head(fecuncity.data)


#degree day sequence and params####
degree.days<-seq(2000,5000,500)

aam.slopes<-NULL
aam.ints<-NULL

for(i in 1:length(degree.days)){
  slope.A.lm<-lm(log(-dd.para$slope_beta)~dd.para$degree_day)
  temp.slope<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*degree.days[i])
  
  int.A.lm<-lm(log(dd.para$int_alpha)~dd.para$degree_day)
  temp.int<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*degree.days[i])
  
  aam.slopes<-c(aam.slopes,temp.slope)
  aam.ints<-c(aam.ints,temp.int)
}

curve(1/(1+exp(aam.slopes[1]*(x-aam.ints[1]))),0,20)
for(i in 2:length(aam.slopes)){
  curve(1/(1+exp(aam.slopes[i]*(x-aam.ints[i]))),0,20,add=T)
}

fecundities<-NULL

for(i in 1:length(degree.days)){
  fecundities<-c(fecundities,round(fecund.int+fecund.slope*degree.days[i],0))
}

curve(fecundities[1]*((10^(length.int+length.slope*log10(x)))/1000),0,1000)
for(i in 2:length(fecundities)){
  curve(fecundities[i]*((10^(length.int+length.slope*log10(x)))/1000),0,1000,add=T)
}

#get mean values for all degree days####
use_proto <- age_size_ipm$proto_ipm

lambdas<-NULL

for(i in 1:length(degree.days)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv.max,
    min.growth=min.growth,
    max.growth.increment=max.growth.increment,
    max.size=max.size,
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[i],
    aam.int=aam.ints[i],
    aam.slope=aam.slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas<-c(lambdas,lambda(test_ipm))
}

plot(degree.days,lambdas,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

lambda.data<-data.frame(DD=degree.days,default.lambda=lambdas)
lambda.data

#good version####
use_proto <- age_size_ipm$proto_ipm
trials<-1000
lambda.matrix<-matrix(NA,nrow=trials,ncol=length(degree.days))
colnames(lambda.matrix)<-degree.days
rownames(lambda.matrix)<-1:trials
head(lambda.matrix)

good.matrix<-matrix(NA,nrow=length(degree.days),ncol=3)
colnames(good.matrix)<-c("mean","CIupper","CIlower")
rownames(good.matrix)<-degree.days
good.matrix


for(i in 1:length(degree.days)){
  A<-randomLHS(trials,length(param.list))
  B<-matrix(nrow = nrow(A), ncol = ncol(A))
  B[,1]<-qnorm(A[,1],surv.inflection,sqrt(param.var[1]))
  B[,2]<-qnorm(A[,2],surv.slope,sqrt(param.var[2]))
  B[,3]<-qunif(A[,3],surv.min-surv.min*0.2,surv.min+surv.min*0.2)
  B[,4]<-qunif(A[,4],surv.max-sd.09^2,surv.max+sd.09^2)
  B[,5]<-qunif(A[,5],min.growth-min.growth*0.2,min.growth+min.growth*0.2)
  B[,6]<-qnorm(A[,6],max.growth.increment,sqrt(param.var[6]))
  B[,7]<-qnorm(A[,7],max.size,sd.max.size)
  B[,8]<-qnorm(A[,8],min.std,sqrt(param.var[8]))
  B[,9]<-qnorm(A[,9],division,sqrt(param.var[9]))
  B[,10]<-qnorm(A[,10],length.int,sd.length.int)
  B[,11]<-qnorm(A[,11],length.slope,sd.length.slope)
  B[,12]<-qnorm(A[,12],fecundities[i],fecuncity.data[fecuncity.data$newdata.me.AnnualDD==degree.days[i],5])
  B[,13]<-qnorm(A[,13],aam.ints[i],sqrt(aam.ints[i]*0.2))
  B[,14]<-qnorm(A[,14],aam.slopes[i],sqrt(abs(aam.slopes[i]*0.2)))
  B[,15]<-qunif(A[,15],percent.hact-percent.hact*0.2,percent.hact+percent.hact*0.2)
  B[,16]<-qunif(A[,16],hatch.surv-hatch.surv*0.2,hatch.surv+hatch.surv*0.2)
  B[,17]<-qnorm(A[,17],young.mean,sqrt(param.var[17]))
  B[,18]<-qnorm(A[,18],young.sd,sqrt(param.var[18]))
  for(j in 1:nrow(B)){
    constant_list_new <- list(
      surv.inflection=B[j,1],
      surv.slope=B[j,2],
      surv.min=B[j,3],
      surv.max=B[j,4],
      min.growth=B[j,5],
      max.growth.increment=B[j,6],
      max.size=B[j,7],
      min.std=B[j,8],
      division=B[j,9],
      length.int=B[j,10],
      length.slope=B[j,11],
      fecund=B[j,12],
      aam.int=B[j,13],
      aam.slope=B[j,14],
      percent.hact= B[j,15],
      hatch.surv=B[j,16],
      young.mean=B[j,17],
      young.sd=B[j,18]
    )
    parameters(use_proto)<-constant_list_new
    test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
    lambda.matrix[j,i]<-lambda(test_ipm)
  }
  good.matrix[i,1]<-mean(lambda.matrix[,i])
  good.matrix[i,2]<-mean(lambda.matrix[,i])+1.96*(sd(lambda.matrix[,i])/sqrt(trials))
  good.matrix[i,3]<-mean(lambda.matrix[,i])-1.96*(sd(lambda.matrix[,i])/sqrt(trials))
}

setwd("~/Documents/Thesis/good code")
write.csv(lambda.matrix,"raw lambda values 0.9.csv",row.names = F)
write.csv(good.matrix,"clean lambda and CI values 0.9.csv",row.names = F)


#0.6 version#####
use_proto <- age_size_ipm$proto_ipm
trials<-1000
lambda.matrix.6<-matrix(NA,nrow=trials,ncol=length(degree.days))
colnames(lambda.matrix.6)<-degree.days
rownames(lambda.matrix.6)<-1:trials
head(lambda.matrix.6)

good.matrix.6<-matrix(NA,nrow=length(degree.days),ncol=3)
colnames(good.matrix.6)<-c("mean","CIupper","CIlower")
rownames(good.matrix.6)<-degree.days
good.matrix.6


for(i in 1:length(degree.days)){
  A<-randomLHS(trials,length(param.list))
  B<-matrix(nrow = nrow(A), ncol = ncol(A))
  B[,1]<-qnorm(A[,1],surv.inflection,sqrt(param.var[1]))
  B[,2]<-qnorm(A[,2],surv.slope,sqrt(param.var[2]))
  B[,3]<-qunif(A[,3],surv.min-surv.min*0.2,surv.min+surv.min*0.2)
  B[,4]<-qunif(A[,4],0.6-sd.06^2,0.6+sd.06^2)
  B[,5]<-qunif(A[,5],min.growth-min.growth*0.2,min.growth+min.growth*0.2)
  B[,6]<-qnorm(A[,6],max.growth.increment,sqrt(param.var[6]))
  B[,7]<-qnorm(A[,7],max.size,sd.max.size)
  B[,8]<-qnorm(A[,8],min.std,sqrt(param.var[8]))
  B[,9]<-qnorm(A[,9],division,sqrt(param.var[9]))
  B[,10]<-qnorm(A[,10],length.int,sd.length.int)
  B[,11]<-qnorm(A[,11],length.slope,sd.length.slope)
  B[,12]<-qnorm(A[,12],fecundities[i],fecuncity.data[fecuncity.data$newdata.me.AnnualDD==degree.days[i],5])
  B[,13]<-qnorm(A[,13],aam.ints[i],sqrt(aam.ints[i]*0.2))
  B[,14]<-qnorm(A[,14],aam.slopes[i],sqrt(abs(aam.slopes[i]*0.2)))
  B[,15]<-qunif(A[,15],percent.hact-percent.hact*0.2,percent.hact+percent.hact*0.2)
  B[,16]<-qunif(A[,16],hatch.surv-hatch.surv*0.2,hatch.surv+hatch.surv*0.2)
  B[,17]<-qnorm(A[,17],young.mean,sqrt(param.var[17]))
  B[,18]<-qnorm(A[,18],young.sd,sqrt(param.var[18]))
  for(j in 1:nrow(B)){
    constant_list_new <- list(
      surv.inflection=B[j,1],
      surv.slope=B[j,2],
      surv.min=B[j,3],
      surv.max=B[j,4],
      min.growth=B[j,5],
      max.growth.increment=B[j,6],
      max.size=B[j,7],
      min.std=B[j,8],
      division=B[j,9],
      length.int=B[j,10],
      length.slope=B[j,11],
      fecund=B[j,12],
      aam.int=B[j,13],
      aam.slope=B[j,14],
      percent.hact= B[j,15],
      hatch.surv=B[j,16],
      young.mean=B[j,17],
      young.sd=B[j,18]
    )
    parameters(use_proto)<-constant_list_new
    test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
    lambda.matrix.6[j,i]<-lambda(test_ipm)
  }
  good.matrix.6[i,1]<-mean(lambda.matrix.6[,i])
  good.matrix.6[i,2]<-mean(lambda.matrix.6[,i])+1.96*(sd(lambda.matrix.6[,i])/sqrt(trials))
  good.matrix.6[i,3]<-mean(lambda.matrix.6[,i])-1.96*(sd(lambda.matrix.6[,i])/sqrt(trials))
}

setwd("~/Documents/Thesis/good code")
write.csv(lambda.matrix.6,"raw lambda values 0.6.csv",row.names = F)
write.csv(good.matrix.6,"clean lambda and CI values 0.6.csv",row.names = F)

#mean lambda with dd####
use_proto <- age_size_ipm$proto_ipm

lambdas.6<-NULL

for(i in 1:length(degree.days)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=0.6,
    min.growth=min.growth,
    max.growth.increment=max.growth.increment,
    max.size=max.size,
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[i],
    aam.int=aam.ints[i],
    aam.slope=aam.slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.6<-c(lambdas.6,lambda(test_ipm))
}

plot(degree.days,lambdas.6,ylim=c(0.85,1.01),
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

#plot all together now####
plot(degree.days,lambdas,ylim=c(0.85,1.4),
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=18)
points(degree.days,lambdas.6)
abline(h=1)
points(degree.days,good.matrix.6[,2],type="l")
points(degree.days,good.matrix.6[,3],type="l")

points(degree.days,good.matrix[,2],type="l")
points(degree.days,good.matrix[,3],type="l")

all.data.9<-data.frame(DD=degree.days,default.lambda=lambdas,
                       CIupper=good.matrix[,2],
                       CIlower=good.matrix[,3])

all.data.9

all.data.6<-data.frame(DD=degree.days,default.lambda=lambdas.6,
                       CIupper=good.matrix.6[,2],
                       CIlower=good.matrix.6[,3])

all.data.6

write.csv(all.data.9,"lambda data survival 0.9.csv",row.names = F)
write.csv(all.data.6,"lambda data survival 0.6.csv",row.names = F)



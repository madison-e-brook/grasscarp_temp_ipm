install.packages("nlstools")
install.packages("ipmr")
library(nlstools)
library(ipmr)
library(RColorBrewer)
library(scales)
setwd("~/Documents/Population modelling/Data")


#growth and von bert####
fish.data<-read.csv("sullivan data raw.csv")
str(fish.data)

plot(fish.data$Age,fish.data$Length..mm.)

#anchor plot with t0 at 0
von.bert.0<-nls(Length..mm. ~ a*(1-exp(-b*(Age-0))), 
                data = fish.data,
                start = list(a = 850, b = 0.36))

summary(von.bert.0)

params.0<-von.bert.0$m$getPars()
curve(params.0[1]*(1-exp(-params.0[2]*(x-0))),0,15,
      ylim=c(200,1000),add=T,col="blue")

growth.sd<-round(sd(fish.data$Length..mm.),0)

L_inf<-round(params.0[1],0)

k<-round(params.0[2],2)

temperatures<-c(2500,5000)
k.values<-c(0.17,0.615)


#make relationship of temp and k####
plot(temperatures,k.values)
k.lm<-lm(k.values~temperatures)
abline(k.lm)


setwd("~/Documents/Thesis/good code")
png("appendix 3 temp depen k plot.png", width= 2404, height= 1600, units="px", res = 300)
plot(temperatures,k.values,pch=16,
     xlab="Annual Degree Days",
     ylab="k value")
abline(k.lm)

text(3500,0.55,expression(italic("k = 0.00018(degree days) - 0.28")))
dev.off()
setwd("~/Documents/Population modelling/Data")

#relationship with L_inf
L.values<-c(810,1297)
temp.Lvalue<-c(5000,2500)

#make relationship of temp and k####
plot(temp.Lvalue,L.values)
L.lm<-lm(L.values~temp.Lvalue)
abline(L.lm)
expression(italic("L"["inf"]* "= 1784 - 0.19(degree days)"))
#plot relationship of temp and L#####
setwd("~/Documents/Thesis/good code")
png("appendix 3 temp depen Linf plot.png", width= 2404, height= 1600, units="px", res = 300)
plot(temp.Lvalue,L.values,pch=16,
     xlab="Annual Degree Days",
     ylab=expression("L"["inf"]*" value"))
abline(L.lm)
text(4000,1250,expression(italic("L"["inf"]* "= 1784 - 0.19(degree days)")))
dev.off()

#survival and temp####
survs.0.6<-c(0.6+0.09,
             0.6-0.09)

survs.0.9<-c(0.9+0.09,
             0.9-0.09)
temperatures<-c(2500,5000)


plot(temperatures,survs.0.6)
plot(temperatures,survs.0.9)

s.lm.06<-lm(survs.0.6~temperatures)
abline(s.lm.06)

s.lm.09<-lm(survs.0.9~temperatures)
abline(s.lm.09)

setwd("~/Documents/Thesis/good code")
png("appendix 3 temp depen survival plot.png", width= 2404, height= 1600, units="px", res = 300)
plot(temperatures,survs.0.6,pch=16,ylim=c(0.5,1),
     xlab="Annual Degree Days",
     ylab="Maximum Adult Survival")
points(temperatures,survs.0.9)
abline(s.lm.06)
abline(s.lm.09,lty=2)
text(3300,0.51,expression(italic("max"["s 0.6"]* "= 0.87 - 7.2x10"^-5*"(degree days)")))
text(4250,0.99,expression(italic("max"["s 0.9"]* "= 1.2 - 7.2x10"^-5*"(degree days)")))
dev.off()
setwd("~/Documents/Population modelling/Data")
#look at temperatures between 2000 and 5000
temp<-seq(2000,5000,100)

k.values<-k.lm$coef[2]*temp + k.lm$coef[1]
L.inf.values<-L.lm$coef[2]*temp + L.lm$coef[1]
surv.values.06<-s.lm.06$coef[2]*temp + s.lm.06$coef[1]
surv.values.09<-s.lm.09$coef[2]*temp + s.lm.09$coef[1]
max.growths<-numeric(length(temp))

#just k changes first####
for(i in 1:length(temp)){
  ages<-seq(0,15,0.5)
  sizes<-L_inf*(1-exp(-k.values[i]*(ages)))
  initial<-sizes[1:length(sizes)-1]
  after<-sizes[-1]
  change<-after-initial
  reg<-lm(change~initial)
  max.growths[i]<-round(reg$coef[1],0)
}

length(max.growths)

#if both change####

max.growths.both<-numeric(length(temp))

for(i in 1:length(temp)){
  ages<-seq(0,15,0.5)
  sizes<-L.inf.values[i]*(1-exp(-k.values[i]*(ages)))
  initial<-sizes[1:length(sizes)-1]
  after<-sizes[-1]
  change<-after-initial
  reg<-lm(change~initial)
  max.growths.both[i]<-round(reg$coef[1],0)
}

length(max.growths.both)


#run model with just degree days changing####
setwd("~/Documents/Population modelling/Data")

#THIS IS FOR A SURVIVAL MAX OF 0.9!!!!!##
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
#increase upper limit for larger max sizes
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
             1500,
             150))

age_size_ipm<-define_pop_state(proto_ipm = age_size_ipm,
                               n_length_age = rep(1/150,150))


age_size_ipm<-make_ipm(proto_ipm=age_size_ipm, iterations = 200,
                       usr_funs = list(growth.increment = growth.increment,
                                       standard.dev = standard.dev))


lam <- lambda(age_size_ipm)
lam



#just with the normal dd changing####
length(temp)
aam.slopes<-NULL
aam.ints<-NULL

for(i in 1:length(temp)){
  slope.A.lm<-lm(log(-dd.para$slope_beta)~dd.para$degree_day)
  temp.slope<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*temp[i])
  
  int.A.lm<-lm(log(dd.para$int_alpha)~dd.para$degree_day)
  temp.int<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*temp[i])
  
  aam.slopes<-c(aam.slopes,temp.slope)
  aam.ints<-c(aam.ints,temp.int)
}

curve(1/(1+exp(aam.slopes[1]*(x-aam.ints[1]))),0,20)
for(i in 2:length(aam.slopes)){
  curve(1/(1+exp(aam.slopes[i]*(x-aam.ints[i]))),0,20,add=T)
}

fecundities<-NULL

for(i in 1:length(temp)){
  fecundities<-c(fecundities,round(fecund.int+fecund.slope*temp[i],0))
}

curve(fecundities[1]*((10^(length.int+length.slope*log10(x)))/1000),0,1000)
for(i in 2:length(fecundities)){
  curve(fecundities[i]*((10^(length.int+length.slope*log10(x)))/1000),0,1000,add=T)
}

#0.6 version
#model with no temp dependence######
surv_good.06<-surv.values.06
use_proto <- age_size_ipm$proto_ipm

lambdas.no.temp<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.no.temp<-c(lambdas.no.temp,lambda(test_ipm))
}

plot(temp,lambdas.no.temp,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only survival######
use_proto <- age_size_ipm$proto_ipm


lambdas.just.surv<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[i],
    min.growth=min.growth,
    max.growth.increment=max.growth.increment,
    max.size=max.size,
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.surv<-c(lambdas.just.surv,lambda(test_ipm))
}

plot(temp,lambdas.just.surv,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

#only growth######
use_proto <- age_size_ipm$proto_ipm


lambdas.just.growth<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.growth<-c(lambdas.just.growth,lambda(test_ipm))
}

plot(temp,lambdas.just.growth,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only linf#####
use_proto <- age_size_ipm$proto_ipm

lambdas.just.Linf<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[i],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.Linf<-c(lambdas.just.Linf,lambda(test_ipm))
}

plot(temp,lambdas.just.Linf,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

#only growth and linf######
use_proto <- age_size_ipm$proto_ipm

lambdas.growth.Linf<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[i],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.growth.Linf<-c(lambdas.growth.Linf,lambda(test_ipm))
}

plot(temp,lambdas.growth.Linf,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only aam######
use_proto <- age_size_ipm$proto_ipm

lambdas.just.aam<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[i],
    aam.slope=aam.slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.aam<-c(lambdas.just.aam,lambda(test_ipm))
}

plot(temp,lambdas.just.aam,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only fecund######
use_proto <- age_size_ipm$proto_ipm

lambdas.just.fecund<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[i],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.fecund<-c(lambdas.just.fecund,lambda(test_ipm))
}

plot(temp,lambdas.just.fecund,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)



#all positive factors#####
use_proto <- age_size_ipm$proto_ipm

lambdas.just.pos<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[i],
    aam.slope=aam.slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.pos<-c(lambdas.just.pos,lambda(test_ipm))
}

plot(temp,lambdas.just.pos,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)



#all negative factors#####
use_proto <- age_size_ipm$proto_ipm

lambdas.just.neg<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[i],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[i],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[i],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.neg<-c(lambdas.just.neg,lambda(test_ipm))
}

plot(temp,lambdas.just.neg,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)




#model with all factors#####
use_proto <- age_size_ipm$proto_ipm

lambdas.all.factor<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.06[i],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[i],
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
  lambdas.all.factor<-c(lambdas.all.factor,lambda(test_ipm))
}

plot(temp,lambdas.all.factor,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#~~~~~~~~repeat with 0.9########
surv_good.09<-surv.values.09
use_proto <- age_size_ipm$proto_ipm

lambdas.no.temp.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.no.temp.09<-c(lambdas.no.temp.09,lambda(test_ipm))
}

plot(temp,lambdas.no.temp.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only survival######
use_proto <- age_size_ipm$proto_ipm


lambdas.just.surv.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[i],
    min.growth=min.growth,
    max.growth.increment=max.growth.increment,
    max.size=max.size,
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.surv.09<-c(lambdas.just.surv.09,lambda(test_ipm))
}

plot(temp,lambdas.just.surv.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

#only growth######
use_proto <- age_size_ipm$proto_ipm


lambdas.just.growth.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.growth.09<-c(lambdas.just.growth.09,lambda(test_ipm))
}

plot(temp,lambdas.just.growth.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only linf#####
use_proto <- age_size_ipm$proto_ipm

lambdas.just.Linf.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[i],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.Linf.09<-c(lambdas.just.Linf.09,lambda(test_ipm))
}

plot(temp,lambdas.just.Linf.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

#only growth and linf######
use_proto <- age_size_ipm$proto_ipm

lambdas.growth.Linf.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[i],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.growth.Linf.09<-c(lambdas.growth.Linf.09,lambda(test_ipm))
}

plot(temp,lambdas.growth.Linf.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only aam######
use_proto <- age_size_ipm$proto_ipm

lambdas.just.aam.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[i],
    aam.slope=aam.slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.aam.09<-c(lambdas.just.aam.09,lambda(test_ipm))
}

plot(temp,lambdas.just.aam.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)


#only fecund######
use_proto <- age_size_ipm$proto_ipm

lambdas.just.fecund.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[i],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.fecund.09<-c(lambdas.just.fecund.09,lambda(test_ipm))
}

plot(temp,lambdas.just.fecund.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)



#all positive factors#####
use_proto <- age_size_ipm$proto_ipm

lambdas.just.pos.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[16],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[16],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[16],
    aam.int=aam.ints[i],
    aam.slope=aam.slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.pos.09<-c(lambdas.just.pos.09,lambda(test_ipm))
}

plot(temp,lambdas.just.pos.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)



#all negative factors#####
use_proto <- age_size_ipm$proto_ipm

lambdas.just.neg.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[i],
    min.growth=min.growth,
    max.growth.increment=max.growths[16],
    max.size=L.inf.values[i],
    min.std=min.std,
    division=division,
    length.int=length.int,
    length.slope=length.slope,
    fecund=fecundities[i],
    aam.int=aam.ints[16],
    aam.slope=aam.slopes[16],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.just.neg.09<-c(lambdas.just.neg.09,lambda(test_ipm))
}

plot(temp,lambdas.just.neg.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)




#model with all factors#####
use_proto <- age_size_ipm$proto_ipm

lambdas.all.factor.09<-NULL

for(i in 1:length(temp)){
  constant_list_new <- list(
    surv.inflection=surv.inflection,
    surv.slope=surv.slope,
    surv.min=surv.min,
    surv.max=surv_good.09[i],
    min.growth=min.growth,
    max.growth.increment=max.growths[i],
    max.size=L.inf.values[i],
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
  lambdas.all.factor.09<-c(lambdas.all.factor.09,lambda(test_ipm))
}

plot(temp,lambdas.all.factor.09,
     xlab="Annual Degree days",
     ylab="Predicted Population Growth Rate",pch=16)

#plot#####

png("simulation output figure.png", width= 2400, height= 2400, units="px", res = 300)
layout(matrix(c(seq(1,9),10,10,10), ncol=3, byrow=TRUE),heights = c(1.5,1.5,1.5,1))

par(oma=c(5,4,0,0)+0.1,mar=c(0,0,0.5,0.5)+0.1,cex=1.2)
plot(temp,lambdas.just.growth,ylim=c(0.65,1.6),
     xlab="",yaxt="n",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.growth.09,type="l",col=cols[1],
       lwd=2)
box()
#axis(1,at=seq(2250,4750,length.out=7),labels=F,cex=1.5)
axis(2,at=seq(0.85,1.45,by=0.15),labels=round(seq(0.85,1.45,by=0.15),1),cex=1.5)
text(4900,1.55,"a)")
text(2100,1.55,expression(paste(italic(k))))


plot(temp,lambdas.just.Linf,ylim=c(0.65,1.6),
     xlab="  ",yaxt="n",
     ylab="",type="l",lwd=2.5,col=cols[2], axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.Linf.09,type="l",col=cols[1],
       lwd=2)

#axis(1,at=seq(2250,4750,length.out=7),labels=F,cex=1.5)
#axis(2,at=seq(0.85,1.45,by=0.15),labels=F,cex=1.5)
box()
text(4900,1.55,"b)")
text(2250,1.55,expression(paste(italic(L[inf]))))


plot(temp,lambdas.growth.Linf,ylim=c(0.65,1.6),
     xlab=" ",yaxt="nn",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.growth.Linf.09,type="l",col=cols[1],
       lwd=2)
#axis(1,at=seq(2250,4750,length.out=7),labels=F,cex=1.5)
#axis(2,at=seq(0.85,1.45,by=0.15),labels=F,cex=1.5)
box()
text(4900,1.55,"c)")
text(2750,1.55,expression(paste(italic(k)," and ",italic(L[inf]),"")))


plot(temp,lambdas.just.aam,ylim=c(0.65,1.6),
     xlab="  ",yaxt="nn",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.aam.09,type="l",col=cols[1],
       lwd=2)
#axis(1,at=seq(2250,4750,length.out=7),labels=F,cex=1.5)
axis(2,at=seq(0.85,1.45,by=0.15),labels=round(seq(0.85,1.45,by=0.15),1),cex=1.5)
box()

text(4900,1.55,"d)")
text(2200,1.55,expression(paste(italic(alpha[m]))))



plot(temp,lambdas.just.fecund,ylim=c(0.65,1.6),
     xlab="",yaxt="nn",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.fecund.09,type="l",col=cols[1],
       lwd=2)
#axis(1,at=seq(2250,4750,length.out=7),labels=F,cex=1.5)
#axis(2,at=seq(0.85,1.45,by=0.15),labels=F,cex=1.5)
box()

text(4900,1.55,"e)")
text(2050,1.55,expression(paste(italic(f))))


plot(temp,lambdas.just.surv,ylim=c(0.65,1.6),
     xlab="",yaxt="n",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.surv.09,type="l",col=cols[1],
       lwd=2)

#axis(1,at=seq(2250,4750,length.out=7),labels=F,cex=1.5)
#axis(2,at=seq(0.85,1.45,by=0.15),labels=F,cex=1.5)
box()


text(4900,1.55,"f)")
text(2400,1.55,expression(paste(italic(max[s]))))



par(mar=c(0,0,0.5,0.5)+0.1)
plot(temp,lambdas.just.pos,ylim=c(0.65,1.6),
     xlab="",yaxt="n",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.pos.09,type="l",col=cols[1],
       lwd=2)
axis(1,at=seq(2250,4750,length.out=7),labels=seq(2250,4750,length.out=7),cex=1.5)
axis(2,at=seq(0.85,1.45,by=0.15),labels=round(seq(0.85,1.45,by=0.15),1),cex=1.5)
box()

text(4900,1.55,"g)")
text(3150,1.55,expression(paste(italic(Only)," ",italic(Positive))))



plot(temp,lambdas.just.neg,ylim=c(0.65,1.6),
     xlab="",yaxt="n",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.just.neg.09,type="l",col=cols[1],
       lwd=2)
axis(1,at=seq(2250,4750,length.out=7),labels=seq(2250,4750,length.out=7),cex=1.5)
#axis(2,at=seq(0.85,1.45,by=0.15),labels=F,cex=1.5)
box()

text(4900,1.55,"h)")
text(3200,1.55,expression(paste(italic(Only)," ",italic(Negative))))


plot(temp,lambdas.all.factor,ylim=c(0.65,1.6),
     xlab="",yaxt="n",
     ylab="",type="l",lwd=2.5,col=cols[2],axes=F,
     panel.first = c(points(temp,lambdas.no.temp,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[2],0.5))
                     ,abline(h=1,lty=3,col="grey",lwd=2.5),
                     points(temp,lambdas.no.temp.09,type="l",lty=2,lwd=2.5,
                            col=alpha(cols[1],0.5))))
points(temp,lambdas.all.factor.09,type="l",col=cols[1],
       lwd=2)
axis(1,at=seq(2250,4750,length.out=7),labels=seq(2250,4750,length.out=7),cex=1.5)
#axis(2,at=seq(0.85,1.45,by=0.15),labels=F,cex=1.5)
box()

text(4900,1.55,"i)")
text(3300,1.55,expression(paste(italic(All)," ",
                                italic(Interactions))))

title(ylab=expression(lambda),
      outer=T,line=3)
title(xlab="Annual Degree Days",
      outer=T,line=-2.5)

#plot.new()
par(fig = c(0, 1, 0, 1),oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),new=T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="bottom",inset=0,
       legend = c("Model Prediction","Prediction without Temperature Dependence",expression(lambda*" = 1"),
                  expression(paste("Default ",italic(max[s])," = 0.6")),
                  expression(paste("Default ",italic(max[s])," = 0.9"))),
       col = c("black","darkgrey","grey",cols[2],cols[1]),
       lty=c(1, 2,3,0,0),pch=c(NA, NA,NA,15,15), xpd = TRUE, horiz = F,cex=1,bty="n",
       lwd=2.5)
dev.off()



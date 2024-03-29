library(raster)
library(maps)
library(viridis)


setwd("~/Documents/Thesis/good code")
#import raster map with only degree days less than 5000
NA.map.dd<-raster("degree day map only over 5000.grd")
plot(NA.map.dd)

library(nlstools)
library(ipmr)
setwd("~/Documents/Population modelling/Data")

#maximum survival of 0.9#####
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
             1100,
             150))

age_size_ipm<-define_pop_state(proto_ipm = age_size_ipm,
                               n_length_age = rep(1/150,150))


age_size_ipm<-make_ipm(proto_ipm=age_size_ipm, iterations = 200,
                       usr_funs = list(growth.increment = growth.increment,
                                       standard.dev = standard.dev))


lam <- lambda(age_size_ipm)
lam


#getting slopes and intercepts for all degree days####

#get slope and int for all degree days from 0 to 5000
temps<-seq(0,5000,1)


slopes<--1*exp(slope.A.lm$coef[1]+slope.A.lm$coef[2]*temps)
ints<-exp(int.A.lm$coef[1]+int.A.lm$coef[2]*temps)
fecunds<-round(fecund.int+fecund.slope*temps,0)

#get the many lambdas#####
use_proto <- age_size_ipm$proto_ipm

lambdas<-NULL


#change slope and int for each location
for(i in 1:length(temps)){
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
    fecund=fecunds[i],
    aam.int=ints[i],
    aam.slope=slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas<-c(lambdas,lambda(test_ipm))
}

storage.results<-cbind(temps,slopes,ints,fecunds,lambdas)
storage.results<-as.data.frame(storage.results)
str(storage.results)
plot(temps,lambdas)
#2000 was used as a minimum degree days
storage.results<-storage.results[storage.results$temps>=2000,]

setwd("~/Documents/Thesis/good code")
write.csv(storage.results,"0.9 survival lambda for all degree days.csv",row.names = F)

storage.results<-storage.results[storage.results$temps>=2000,]


map.09<-NA.map.dd
plot(map.09)
map.09[map.09<2000]<-NA
plot(map.09)
map.09

for(i in 1:nrow(storage.results)){
  #each cell greater than or equal to lower value but less than that value+1 gets
  #the lambda associated with lower value
  map.09[map.09>=storage.results[i,1] & map.09<storage.results[i+1,1]]<-storage.results[i,5]
}

writeRaster(map.09,"0.9 survival North america lambda map.grd",overwrite=T)


#repeat for lower survival####
lambdas.6<-NULL


#change slope and int for each location
for(i in 1:length(temps)){
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
    fecund=fecunds[i],
    aam.int=ints[i],
    aam.slope=slopes[i],
    percent.hact= percent.hact,
    hatch.surv=hatch.surv,
    young.mean=young.mean,
    young.sd=young.sd
  )
  parameters(use_proto)<-constant_list_new
  test_ipm<-make_ipm(proto_ipm=use_proto, iterations = 100)
  lambdas.6<-c(lambdas.6,lambda(test_ipm))
}

storage.results.6<-cbind(temps,slopes,ints,fecunds,lambdas.6)
storage.results.6<-as.data.frame(storage.results.6)
str(storage.results)


setwd("~/Documents/Thesis/good code")
write.csv(storage.results.6,"0.6 survival lambda for all degree days.csv",row.names = F)
head(storage.results.6)

storage.results.6<-storage.results.6[storage.results.6$temps>=2000,]


map.06<-NA.map.dd

map.06<-NA.map.dd
plot(map.06)
map.06[map.06<2000]<-NA
plot(map.06)


for(i in 1:nrow(storage.results.6)){
  #each cell greater than or equal to lower value but less than that value+1 gets
  #the lambda associated with lower value
  map.06[map.06>=storage.results.6[i,1] & map.06<storage.results.6[i+1,1]]<-storage.results.6[i,5]
}

writeRaster(map.06,"0.6 survival North america lambda map.grd",overwrite=T)


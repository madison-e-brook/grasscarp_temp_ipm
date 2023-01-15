setwd("~/Documents/Population modelling/Data")

#download Sullivan et al., 2020 data
fish.data<-read.csv("sullivan data raw.csv")
str(fish.data)


fish.data$Dam<-as.factor(fish.data$Dam)

#the relationship for the change in growth increment####
plot(fish.data$Age,fish.data$Length..mm.)
von.bert.0<-nls(Length..mm. ~ a*(1-exp(-b*(Age-0))), 
                data = fish.data,
                start = list(a = 850, b = 0.36))

summary(von.bert.0)$coef[1,2]


params.0<-von.bert.0$m$getPars()
curve(params.0[1]*(1-exp(-params.0[2]*(x-0))),0,15,
      ylim=c(200,1000),add=T,col="blue")

#von bert plot####
setwd("~/Documents/Thesis/good code")
png("von bert raw.png", width= 2404, height= 1600, units="px", res = 300)
plot(fish.data$Age,fish.data$Length..mm.,
     xlab="Age (year)",ylab="Length (mm)")
curve(params.0[1]*(1-exp(-params.0[2]*(x-0))),0,15,
      ylim=c(200,1000),add=T,col="blue")
text(10,550,expression(italic("L"[t]*" = 813 "*"(1 - e"^"-0.44t"*")")))
dev.off()



#growth increment####
ages<-seq(0,15,0.5)
sizes<-params.0[1]*(1-exp(-params.0[2]*(ages)))
points(ages,sizes)
initial<-sizes[1:length(sizes)-1]
after<-sizes[-1]
change<-after-initial

#growth rate
lm(change~initial)


png("growth increment for appendix1.png", width= 2404, height= 1600, units="px", res = 300)
plot(initial,change,type="l",xlab="Length (mm)",
     ylab="Expected growth increment (mm)",xlim=c(0,1000))
segments(813,0.1,1000,0.1)
text(700,145,expression(italic("Length < 813; y = 161 – "*frac(161,813)*" * Length")))
text(605,120,expression(italic("Length "*"\u2265"*" 813; y = 0.1")))
dev.off()

#distribution graph####
growth.increment<-function(length_1,min.growth=0.1,max.growth.increment=161,max.size=813){
  incre<-max.growth.increment+(-max.growth.increment/max.size)*length_1
  for(i in 1:length(length_1)){
    if(incre[i]<=min.growth){
      incre[i]<-min.growth
    }
  }
  return(incre)
}


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
sizes<-seq(100,900,200)
growth.means<-growth.increment(sizes)
growth.sd<-standard.dev(sizes)


x<-seq(0,950,1)

plot(x,dnorm(x,growth.means[1]+sizes[1],growth.sd[1]),type="l",ylim=c(0,0.08),
     col=2,xlab="Length (mm)",ylab="Density")
for(i in 2:length(sizes)){
  points(x,dnorm(x,growth.means[i]+sizes[i],growth.sd[i]),type="l",col=i+1)
}

abline(v=sizes,col=2:6)
abline(h=0)
legend("topleft",paste(as.character(sizes),"mm"),lty=1,
       col=2:6)

png("growth distributions appendix1.png", width= 2404, height= 1600, units="px", res = 300)
plot(x,dnorm(x,growth.means[1]+sizes[1],growth.sd[1]),type="l",ylim=c(0,0.08),
     col=2,xlab=expression("Length"["t+1"]),ylab="Density",
     lwd=1.5)
for(i in 2:length(sizes)){
  points(x,dnorm(x,growth.means[i]+sizes[i],growth.sd[i]),type="l",col=i+1,lwd=1.5)
}

abline(v=sizes,col=2:6)
abline(h=0)
legend("topleft",paste(as.character(sizes),"mm"),lty=1,lwd=1.5,
       col=2:6,title=expression("Length"["t"]))
dev.off()





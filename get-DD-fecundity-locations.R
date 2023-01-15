library(raster)
library(ncdf4)

#import our location data
setwd("~/Documents/Thesis/good code")
loc<-read.csv("locations of relative fecundities.csv")
str(loc)



#set wd to where CPC data was downloaded to
setwd("/Volumes/My Passport/CPC temperature data/")

#download one file to get the spatial points
tmin.1979 <- brick("tmin.1979.nc", varname = "tmin")
tmin.1979<-rotate(tmin.1979)
plot(tmin.1979$X1979.01.01)

#spatial points####
lat<-loc$Lat
lon<-loc$Lon
#turn them to spatial points
coords<-data.frame(longitude=lon, latitude=lat)
locations<-SpatialPoints(coords, proj4string = tmin.1979@crs)

points(locations)

#we are interested in 41 years (1979-2019)
degree_day<-matrix(NA,length(loc$Loc),42)
degree_day[,1]<-loc$Loc
#we want 1979-2019
years<-seq(1979,2019,1)

#for loop####
for(i in 1:41){
  #download tmax
  tmax.name<-paste0("tmax.",years[i],".nc")
  tmax<- brick(tmax.name, varname = "tmax")
  tmax<-rotate(tmax)
  #extract tmax
  max<-extract(tmax,locations,layer=1)
  #download tmin
  tmin.name<-paste0("tmin.",years[i],".nc")
  tmin<- brick(tmin.name, varname = "tmin")
  tmin<-rotate(tmin)
  #extract tmin
  min<-extract(tmin,locations,layer=1)
  
  #calculate the degree days
  degree<-(max+min)/2
  degree[degree<0]<-0
  calculated.degrees<-rowSums(degree)
  degree_day[,i+1]<-calculated.degrees
  removeTmpFiles(h=0.001)
}


degree_day<-as.data.frame(degree_day)
str(degree_day)
for(i in 1:41){
  degree_day[,i+1]<-as.numeric(degree_day[,i+1])
}

annual.average.dd<-rowMeans(degree_day[,2:42],na.rm = TRUE)
annual.average.dd.sd<-apply(degree_day[,2:42],1,sd,na.rm=T)


loc$AnnualDD<-annual.average.dd

#save file####
setwd("~/Documents/Thesis/good code")
write.csv(loc,"fecundities with degree days.csv",row.names = F)

boxplot(loc$Fecundity)

plot(loc$AnnualDD,loc$Fecundity)
abline(lm(loc$Fecundity~loc$AnnualDD))
summary(lm(loc$Fecundity~loc$AnnualDD))

abline(lm(loc$Fecundity[-1]~loc$AnnualDD[-1]))
summary(lm(loc$Fecundity[-1]~loc$AnnualDD[-1]))






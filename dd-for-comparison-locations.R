library(raster)
setwd("~/Documents/Thesis/good code")
#upload degree day base 0 raster map
original.map<-raster("GDD0_mean.tif")
plot(original.map)

names<-c("Mississippi","Erie","Superior")
#coordinates match locations of interest
lat<-c(40.379319,42.2,47.7)
lon<-c(-91.419208,-81.2,-87.5)

coords<-data.frame(x=lon,y=lat)
points(loc)

#convert to spatial points
loc<-SpatialPoints(coords,proj4string =original.map@crs)
points(loc)
#extract temperature data
extract(original.map,loc)

round(extract(original.map,loc),-1)
plot(1:3,round(extract(original.map,loc),-1))

location.data<-data.frame(loc=names,lat=lat,lon=lon,DD=round(extract(original.map,loc),-1))
location.data

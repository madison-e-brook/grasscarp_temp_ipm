library(ncdf4)
library(sp)
library(raster)

setwd("~/Desktop/CPC Global Daily Temperature")

#downloading the datafiles into the computer (approx 7GB)####

#make the outline of the url
start<-"ftp://ftp.cdc.noaa.gov/Datasets/cpc_global_temp/"
end<-".nc"

#the data contains two variables (tmax and tmin)
#for years 1979-2020 (but I'm just downloading until the end of 2019)
years<-seq(1979,2019,1)
var<-c("tmax.","tmin.") #periods added to make url easier

#download all of the data files (took me 2.5-3 hours)####
for(i in 1:length(var)){
  for(j in 1:length(years)){
    myurl<-paste0(start,var[i],years[j],end)
    download.file(myurl,
                  destfile = paste0("~/Desktop/CPC Global Daily Temperature/",
                                    var[i],years[j],end),mode="wb")
  }
}
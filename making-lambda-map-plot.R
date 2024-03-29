library(colorBlindness)
library(viridis)
library(RColorBrewer)

library(raster)
library(maps)
library(maptools)
setwd("~/Documents/Thesis/good code")

original.map<-raster("GDD0_mean.tif")

over.2000<-original.map
over.2000[over.2000>2000]<-NA
plot(over.2000)
over.2000<-reclassify(over.2000,matrix(c(0,2000,10000),1,3))
plot(over.2000)


over.5000<-original.map
over.5000[over.5000<5000]<-NA
plot(over.5000)
over.5000<-reclassify(over.5000,matrix(c(0,100000,10000),1,3))
plot(over.5000)

#import lambda graphs####
map.9<-raster("0.9 survival North america lambda map.grd")
map.6<-raster("0.6 survival North america lambda map.grd")

plot(map.9)
plot(map.6)
map('world', regions=c('usa','canada','mexico'),add=T)


#adding layers

#good plot with both#####
setwd("~/Documents/Thesis/good code")
png("both maps on one figure.png", width= 2404, height= 1000, units="px", res = 300)

test<-colorRampPalette(c("darkred","orange","wheat"))
par(fig = c(0, 0.45, 0, 1),bty = 'n',mar=c(6,4,1,2))
plot(map.9,col=rev(test(100)),zlim=c(summary(map.9)[1],summary(map.9)[5]),xlab="Longitude",ylab="Latitude")
par(fig = c(0, 0.45, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2))
plot(over.5000,col="grey",legend=F,axes=F)
par(fig = c(0, 0.45, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2))
plot(over.2000,col="black",legend=F,axes=F)
par(fig = c(0, 0.45, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2))
map('world', regions=c('usa','canada','mexico'),add=T,xlim=c(-167,-60),
    ylim=c(25,100))
map('lakes',add=T,xlim=c(-180,-60))
par(fig = c(0, 0.45, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2),xpd=NA)
text(-60,75,"a",cex=1)

test<-colorRampPalette(c("darkblue","cornflowerblue","azure"))
par(fig = c(0.5, 0.95, 0, 1),new=T,mar=c(6,4,1,2))
plot(map.6,col=rev(test(100)),zlim=c(c(summary(map.6)[1],summary(map.6)[5])),xlab="Longitude",ylab="Latitude")
par(fig = c(0.5, 0.95, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2))
plot(over.5000,col="grey",legend=F,axes=F)
par(fig = c(0.5, 0.95, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2))
plot(over.2000,col="black",legend=F,axes=F)
par(fig = c(0.5, 0.95, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2))
map('world', regions=c('usa','canada','mexico'),add=T,xlim=c(-167,-60),
    ylim=c(25,85))
map('lakes',add=T,xlim=c(-180,-60))
par(fig = c(0.5, 0.95, 0, 1),new=T,bty = 'n',mar=c(6,4,1,2),xpd=NA)
text(-60,75,"b",cex=1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",
       legend = c("Unlikely to survive (DD<2000)","Outside of model parameters"),
       col = c("black","grey"),
       pch=c( 15, 15), xpd = TRUE, horiz = T, cex = 0.9, bty = 'n')
dev.off()



#at high survival
setwd("~/Documents/Thesis/good code")
test<-colorRampPalette(c("darkred","orange","wheat"))
png("lambda at 0.9 surv map.png", width= 2404, height= 1600, units="px", res = 300)
par(mfrow=c(1,1),mar=c(4,2,2,2))
plot(map.9,col=rev(test(100)))
par(mfrow=c(1,1),new=T)
plot(over.5000,col="grey",legend=F)
par(mfrow=c(1,1),new=T)
plot(over.2000,col="black",legend=F)
par(mfrow=c(1,1),new=T)
map('world', regions=c('usa','canada','mexico'),add=T)
map('lakes',add=T)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",
       legend = c("Unlikely to survive (DD<2000)","Outside of model parameters"),
       col = c("black","grey"),
       pch=c( 15, 15), xpd = TRUE, horiz = T, cex = 0.9, bty = 'n')
dev.off()

#at 0.6 survival#####
test<-colorRampPalette(c("darkblue","cornflowerblue","azure"))
png("lambda at 0.6 surv map.png", width= 2404, height= 1600, units="px", res = 300)
par(mfrow=c(1,1),mar=c(4,2,2,2))
plot(map.6,col=rev(test(100)))
par(mfrow=c(1,1),new=T)
plot(over.5000,col="grey",legend=F)
par(mfrow=c(1,1),new=T)
plot(over.2000,col="black",legend=F)
par(mfrow=c(1,1),new=T)
map('world', regions=c('usa','canada','mexico'),add=T)
map('lakes',add=T)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",
       legend = c("Unlikely to survive (DD<2000)","Outside of model parameters"),
       col = c("black","grey"),
       pch=c( 15, 15), xpd = TRUE, horiz = T, cex = 0.9, bty = 'n')
dev.off()

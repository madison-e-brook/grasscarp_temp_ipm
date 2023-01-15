setwd("~/Documents/Thesis/good code")


e.2330.06.10<-read.csv("elas.2330.06 minus 10.csv")
e.2330.06<-read.csv("elas.2330.06 plus 10.csv")

e.2330.09.10<-read.csv("elas.2330.09 minus 10.csv")
e.2330.09<-read.csv("elas.2330.09 plus 10.csv")

e.3560.06.10<-read.csv("elas.3560.06 minus 10.csv")
e.3560.06<-read.csv("elas.3560.06 plus 10.csv")

e.3560.09.10<-read.csv("elas.3560.09 minus 10.csv")
e.3560.09<-read.csv("elas.3560.09 plus 10.csv")


e.4280.06.10<-read.csv("elas.4280.06 minus 10.csv")
e.4280.06<-read.csv("elas.4280.06 plus 10.csv")

e.4280.09.10<-read.csv("elas.4280.09 minus 10.csv")
e.4280.09<-read.csv("elas.4280.09 plus 10.csv")


#minus length and weight variables
#and the two other parameters that had very small elasticities

e.2330.06.10<-e.2330.06.10[-c(5,8,10,11),]
e.2330.06<-e.2330.06[-c(5,8,10,11),]

e.3560.06.10<-e.3560.06.10[-c(5,8,10,11),]
e.3560.06<-e.3560.06[-c(5,8,10,11),]

e.4280.06.10<-e.4280.06.10[-c(5,8,10,11),]
e.4280.06<-e.4280.06[-c(5,8,10,11),]

e.2330.09.10<-e.2330.09.10[-c(5,8,10,11),]
e.2330.09<-e.2330.09[-c(5,8,10,11),]

e.3560.09.10<-e.3560.09.10[-c(5,8,10,11),]
e.3560.09<-e.3560.09[-c(5,8,10,11),]

e.4280.09.10<-e.4280.09.10[-c(5,8,10,11),]
e.4280.09<-e.4280.09[-c(5,8,10,11),]

#combiine data
dd.2330<-rowMeans(cbind(e.2330.06.10$elas,e.2330.06$elas,e.2330.09.10$elas,e.2330.09$elas))
dd.3560<-rowMeans(cbind(e.3560.06.10$elas,e.3560.06$elas,e.3560.09.10$elas,e.3560.09$elas))
dd.4280<-rowMeans(cbind(e.4280.06.10$elas,e.4280.06$elas,e.4280.09.10$elas,e.4280.09$elas))

all.data<-cbind(dd.2330,dd.3560,dd.4280)

plot(rowMeans(cbind(e.2330.06.10$elas,e.2330.06$elas,e.2330.09.10$elas,e.2330.09$elas)),col="black",pch=1,
     ylim=c(-0.3,0.9))
points(rowMeans(cbind(e.3560.06.10$elas,e.3560.06$elas,e.3560.09.10$elas,e.3560.09$elas)),col="red",pch=1)
points(rowMeans(cbind(e.4280.06.10$elas,e.4280.06$elas,e.4280.09.10$elas,e.4280.09$elas)),col="blue",pch=1)
abline(h=0)

#which parameters have postive or negative elasticities
neg<-c(1,2,7,9,10,14)
pos<-c(3,4,5,6,8,11,12,13)


plot(rowMeans(cbind(e.2330.06.10$elas,e.2330.06$elas,e.2330.09.10$elas,e.2330.09$elas)),col="black",pch=1,
     ylim=c(-0.3,0.9))
points(rowMeans(cbind(e.3560.06.10$elas,e.3560.06$elas,e.3560.09.10$elas,e.3560.09$elas)),col="red",pch=1)
points(rowMeans(cbind(e.4280.06.10$elas,e.4280.06$elas,e.4280.09.10$elas,e.4280.09$elas)),col="blue",pch=1)
abline(h=0)

for(i in neg){
  text(i,min(all.data[i,])-0.05,signif(rowMeans(all.data)[i],2))
}

for(i in pos){
  text(i,max(all.data[i,])+0.05,signif(rowMeans(all.data)[i],2))
}

#good graph####
paras<-c(expression(alpha[s]),
         expression(beta[s]),
         expression(italic(min[s])),
         expression(italic(max[s])),
         expression(italic(max[g])),
         expression(italic(L[inf])),
         expression(italic(d[g])),
         expression(italic(f)),
         expression(alpha[m]),
         expression(beta[m]),
         expression(italic(hp)),
         expression(italic(hs)),
         expression(mu[r]),
         expression(sigma[r]))

setwd("~/Documents/Thesis/good code")
png("elasticities.png", width= 3000, height= 1600, units="px", res = 300)
par(xpd=T,mar = c(6, 4, 4, 2) + 0.1)
plot(rowMeans(cbind(e.2330.06.10$elas,e.2330.06$elas,e.2330.09.10$elas,e.2330.09$elas)),pch=1,
     ylim=c(-0.3,0.9),xaxt = "n",
     xlab="",
     ylab="Elasticity",
     col=viridis(4)[1])
axis(1, at = 1:14,
     labels = paras)
text(4,-0.7,"Growth \n and \n Survival")
segments(3,-0.7,1,-0.7)
segments(1,-0.7,1,-0.6)
segments(5,-0.7,7,-0.7)
segments(7,-0.7,7,-0.6)

text(11,-0.7,"Reproduction")
segments(9.5,-0.7,8,-0.7)
segments(8,-0.7,8,-0.6)
segments(12.5,-0.7,14,-0.7)
segments(14,-0.7,14,-0.6)

points(rowMeans(cbind(e.3560.06.10$elas,e.3560.06$elas,e.3560.09.10$elas,e.3560.09$elas)),col=viridis(4)[2],pch=1)
points(rowMeans(cbind(e.4280.06.10$elas,e.4280.06$elas,e.4280.09.10$elas,e.4280.09$elas)),col=viridis(4)[3],pch=1)
segments(0.5,0,14.5,0)



for(i in neg){
  text(i,min(all.data[i,])-0.05,signif(rowMeans(all.data)[i],2))
}

for(i in pos){
  text(i,max(all.data[i,])+0.05,signif(rowMeans(all.data)[i],2))
}
legend("topright",c("Superior (DD = 2330)",
                    "Erie (DD = 3560)",
                    "Mississippi (DD = 4280)"),
       pch=1,col=c(viridis(4)[1],
                   viridis(4)[2],
                   viridis(4)[3]),bty="n")
dev.off()


#means of the variables not in the figure#####
e.2330.06.10<-read.csv("elas.2330.06 minus 10.csv")
e.2330.06<-read.csv("elas.2330.06 plus 10.csv")

e.2330.09.10<-read.csv("elas.2330.09 minus 10.csv")
e.2330.09<-read.csv("elas.2330.09 plus 10.csv")

e.3560.06.10<-read.csv("elas.3560.06 minus 10.csv")
e.3560.06<-read.csv("elas.3560.06 plus 10.csv")

e.3560.09.10<-read.csv("elas.3560.09 minus 10.csv")
e.3560.09<-read.csv("elas.3560.09 plus 10.csv")


e.4280.06.10<-read.csv("elas.4280.06 minus 10.csv")
e.4280.06<-read.csv("elas.4280.06 plus 10.csv")

e.4280.09.10<-read.csv("elas.4280.09 minus 10.csv")
e.4280.09<-read.csv("elas.4280.09 plus 10.csv")

e.2330.06.10<-e.2330.06.10[c(5,8,10,11),]
e.2330.06<-e.2330.06[c(5,8,10,11),]

e.3560.06.10<-e.3560.06.10[c(5,8,10,11),]
e.3560.06<-e.3560.06[c(5,8,10,11),]

e.4280.06.10<-e.4280.06.10[c(5,8,10,11),]
e.4280.06<-e.4280.06[c(5,8,10,11),]

e.2330.09.10<-e.2330.09.10[c(5,8,10,11),]
e.2330.09<-e.2330.09[c(5,8,10,11),]

e.3560.09.10<-e.3560.09.10[c(5,8,10,11),]
e.3560.09<-e.3560.09[c(5,8,10,11),]

e.4280.09.10<-e.4280.09.10[c(5,8,10,11),]
e.4280.09<-e.4280.09[c(5,8,10,11),]

e.2330<-rowMeans(cbind(e.2330.06.10$elas,e.2330.06$elas,e.2330.09.10$elas,e.2330.09$elas))
e.3560<-rowMeans(cbind(e.3560.06.10$elas,e.3560.06$elas,e.3560.09.10$elas,e.3560.09$elas))
e.4280<-rowMeans(cbind(e.4280.06.10$elas,e.4280.06$elas,e.4280.09.10$elas,e.4280.09$elas))

rowMeans(cbind(e.2330,e.3560,e.4280))



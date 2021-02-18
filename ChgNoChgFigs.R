#
library(RStoolbox)
library(raster)
library(rgdal)
library(aRn)
library(fitdistrplus)
library(ggplot2)

# MAC
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Orinoquia"
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Pucallpa"
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Mexico"

# WINDOWS
path=#("X:/VictorShare/s3dFiles/Pucallpa")
#("X:/VictorShare/s3dFiles/MontesTest")
("X:/VictorShare/s3dFiles/Orinoquia")
("X:/VictorShare/s3dFiles/Mexico")#
setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

# Extract data from files
## path
## shape, rate, kd over time


# Using masks
## ECDF before and after iterations (gamma for s3d and normal for bands)
## 

# Retrieve infomrmation from R objects
rdata=list.files('.', pattern='.RData')

#organize as Gamma, Gamma CCA, Chisq CCA
rdata=rdata[c(2,3,1)] 
names=c("Gamma", "Gamma CCA", "Chisquare CCA")#, "Chisquare2 CCA")

for(i in 1:length(rdata)){
  load(rdata[i])
  data=s3dmod$paramstats[[2]]
  data=data[,-c(3, 4)]
  #data$location=rep(substr(path, 25, nchar(path)), nrow(data))
  data$names=rep(names[i], nrow(data))
  if(i==1){dataset=data} else{dataset=rbind(dataset,data)}
  # 
  # print(nrow(s3dmod$data[[6]]))
  # print(max(data$iter))
  # plot(data$ksD~data$iter, xlab="iter", ylab="ksD")
  # ksDstd=(max(data$ksD)-data$ksD)/max(data$ksD)
  # lag=ksDstd[2:length(ksDstd)]
  # lag-ksDstd[1:length(ksDstd)-1]
  # ksDlag=data$ksD[2:length(data$ksD)]
  # diff=ksDlag-data$ksD[1:length(data$ksD)-1]
  # diff2=diff[2:length(diff)]
  # print(diff2/diff[1:length(diff2)])
  # print(data$ksD)
  # print(data$ksD[nrow(data)]-data$ksD[nrow(data)-1])
  # print(data$ksD[13]-data$ksD[12])
  #plot(data$shape~data$iter, xlab="iter", ylab="shape")
  #plot(data$rate~data$iter, xlab="iter", ylab="rate")
# 
#  lmparamtot=s3dmod$paramstats[[1]]
#  par(mar = c(4, 4, 1, 1) + 0.1)
#  par(mfrow=c(round(sqrt(6))+1,round(sqrt(6))))
#  for(b in 1:6){
#    plot(lmparamtot$iter[which(lmparamtot$band==b)],
#         lmparamtot$intercept[which(lmparamtot$band==b)], xlab="band", ylab="intercept")
#  }
 # par(mfrow=c(round(sqrt(6))+1,round(sqrt(6))))
 # for(b in 1:6){
 #   plot(lmparamtot$iter[which(lmparamtot$band==b)],
 #        lmparamtot$slope[which(lmparamtot$band==b)], xlab="band", ylab="slope")
 # }
}
pdf(file="Convergence.pdf",width=4,height=4,paper='special')
  ggplot(dataset, aes(x=dataset$iter, y=dataset$ksD, col=names)) + 
   geom_line() + ylim(0, 0.5)
 #ggplot(dataset, aes(x=dataset$iter, y=dataset$rate, col=names)) +
 #  geom_line() + ylim(0, 0.4)
 #ggplot(dataset, aes(x=dataset$iter, y=dataset$shape, col=names)) + 
 #  geom_line() + ylim(0, 2.5)
dev.off()




s3dmod$paramstats[[2]]


calp=s3dmod
par(mfrow=c(2,3))
pdf("reg.pdf")
for (b in 1:6){
  plot(calp$data[[b]][,1]~calp$data[[b]][,2], xlab="target", ylab="reference")
  abline(unlist(calp$parameters[[b]][1])[1], unlist(calp$parameters[[b]][1])[2])
}
dev.off()
lmparamtot=s3dmod$paramstats[[1]]
strips=stacks
par(mar = c(4, 4, 1, 1) + 0.1)
par(mfrow=c(2,3))
for(b in 1:6){
  plot(lmparamtot$iter[which(lmparamtot$band==b)], 
       lmparamtot$intercept[which(lmparamtot$band==b)], xlab="band", ylab="intercept")
}
par(mfrow=c(2,3))
for(b in 1: 6){
  plot(lmparamtot$iter[which(lmparamtot$band==b)], 
       lmparamtot$slope[which(lmparamtot$band==b)], xlab="band", ylab="slope")
}

#
library(RStoolbox)
library(raster)
library(rgdal)
library(aRn)
library(fitdistrplus)
library(ggplot2)

# MAC
#path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ImistManuscript/Orinoquia"
# WINDOWS
path=#("X:/VictorShare/s3dFiles/Pucallpa")
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
names=c("05057CCA", "06057CCA", "05057", "06057")  # ORINOQUIA
names=c("06066", "07066") # PUCALLPA
names=c("026046", "026047") # MEXICO
for(i in 1:length(rdata)){
  load(rdata[i])
  data=s3dmod$paramstats[[2]]
  data$location=rep(substr(path, 25, nchar(path)), nrow(data))
  data$names=rep(names[i], nrow(data))
  if(i==1){dataset=data} else{dataset=rbind(dataset,data)}
}
pdf(file="gammaResults.pdf",width=4,height=4,paper='special') 
  ggplot(dataset, aes(x=dataset$iter, y=dataset$ksD, col=names)) + 
    geom_line() + ylim(0, 0.3)
  ggplot(dataset, aes(x=dataset$iter, y=dataset$rate, col=names)) +
    geom_line() + ylim(0, 0.4)
  ggplot(dataset, aes(x=dataset$iter, y=dataset$shape, col=names)) + 
    geom_line() + ylim(0, 2.5)
dev.off()



s3dmod$paramstats[[2]]


calp=s3dmod
par(mfrow=c(2,3))
for (b in 1:nlayers(strips[[1]])){
  plot(calp$data[[b]][,1]~calp$data[[b]][,2], xlab="target", ylab="reference")
  abline(unlist(calp$parameters[[b]][1])[1], unlist(calp$parameters[[b]][1])[2])
}

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

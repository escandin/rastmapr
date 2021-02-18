# This code processes a workflow for producing time series of land cover classifications
# and perform accuracy assessment using landsat 8 collection 1 level 2 data

# WORKFLOW
# 1 Downloading (pending)
# 2 Stacking
# 3 Cloud filtering
# 4 Topographic correction (if needed)
# 5 Normalization
# 6 Mosaicking
# 7 LC calibration
# 8 LC validation
# 9 LC extrapolation

####### ENVIRONMENT SETTINGS
library(RStoolbox)
library(raster)
library(rgdal)
library(aRn)
library(randomForest)
library(caret)
library(maptools)
library(stringr)
library(ggplot2)

path=#("X:/VictorShare/aRnFiles/Pucallpa")
("X:/VictorShare/aRnFiles/Orinoquia")
#("X:/VictorShare/aRnFiles/Mexico")#

savepath=(paste(path, "outputsBandsAll", sep="/"))
savepath1=(paste(path, "outputsBands3_7", sep="/"))
setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

##### 1. DATA INPUTS
tares <- list.files('.', pattern='tar.gz')

calibdata=readOGR(".","calibdata_2016_NEW_2.0")# Orinoquia
readOGR(".", "Training_DF_2018_5.0")               # Mexico
readOGR(".", "Training_Pucallpa_2015_5.0")    #Pucallpa

validdata= readOGR(".", "Valdata_2013_NEW_2.0")  #Orinoquia
readOGR(".", "Training_DF_2010_5.0") # for Mexico    # Mexico
readOGR(".", "Training_Pucallpa_2011_5.0_")           #Pucallpa

# Check if levels are the same
levels(validdata$CLASS_NAME)==levels(calibdata$CLASS_NAME)

# FOR ORINOQUIA
  levels(validdata$CLASS_NAME)=levels(calibdata$CLASS_NAME)

# FOR MEXICO # correct mislabeling of sand class
  #calibdata@polygons=calibdata@polygons[which(calibdata@data$class!="sand")]
  calibdata@data$class[which(calibdata@data$class=="Baresoil")]="baresoil"
  calibdata@data$class=droplevels(calibdata@data$class, "Baresoil")
  
  calibdata@data$class[which(calibdata@data$class=="Burn Scar (new)")]="Burn Scar (New)"
  calibdata@data$class=droplevels(calibdata@data$class, "Burn Scar (new)")
  
  calibdata@polygons=calibdata@polygons[which(calibdata@data$class!="Airport_1")]
  calibdata@data=calibdata@data[which(calibdata@data$class!="Airport_1"),]
  calibdata@data$class=droplevels(calibdata@data$class, "Airport_1")
  
  calibdata@polygons=calibdata@polygons[which(calibdata@data$class!="Airport_2")]
  calibdata@data=calibdata@data[which(calibdata@data$class!="Airport_2"),]
  calibdata@data$class=droplevels(calibdata@data$class, "Airport_2")

  #validdata@polygons=validdata@polygons[which(validdata@data$class!="snow")]
  #validdata@data=validdata@data[which(validdata@data$class!="snow"),]
  #validdata@data$class=droplevels(validdata@data$class, "snow")
  levels(validdata$class)==levels(calibdata$class)
  
#stack with ALOS PALSAR
HHref=raster( "palsar_HH_2018_CDMX_Spk.tif" )
#("palsar_HH_2015_Peru_Spk.tif")
HVref=raster("palsar_HV_2018_CDMX_Spk.tif")
#("palsar_HV_2015_Peru_Spk.tif")
HHtar=raster("palsar_HH_2010_CDMX_Spk.tif")
#("palsar_HH_2010_Peru_Spk.tif")
HVtar=raster("palsar_HV_2010_CDMX_Spk.tif")
#("palsar_HV_2010_Peru_Spk.tif")

##### 2 create original stacks ###########################################

HHHVref=stack(HHref, HVref)
HHHVtar=stack(HHtar, HVtar)
bnames=c("HH", "HV")
names(HHHVref)=bnames
names(HHHVtar)=bnames

writeRaster(HHHVref, "HHHVref.tif")
writeRaster(HHHVtar, "HHHVtar.tif")
rm(HHref, HHtar, HVref, HVtar, HHHVref, HHHVtar)#, HHHVref,HHHVtar)

##### 3-6 Cloud filtering, normalizing, mosaicking ####
# create mosaics for landsat

# FOR OLI
stacks<- EEstackWithoutMeta(tares, sat.nm="LO08")
normethods=c("imist", "ed", "none")
# This is to use only the data that is common between all layers for normalization
#e=commonExtent(c(list(calibdata, validdata),stacks))
for (i in 1:length(normethods)){
  start=Sys.time()
  mosaick<- smg(inlist=stacks, method= normethods[i], refitem=2, mosaicitems=c(1,3),
                pval.pif = 1e-02, pval.chg=0.99, QAbandname="pixel_qa", indem=NA,#demfile,  # DEM file only for Mexico
                normbands=seq(3, 7), verbose=TRUE, sensor="OLI", norm.ext=NULL, cloudbuff=NA)
  print(Sys.time()-start)  # Processing time: 1.39 hrs aRn, 24.88 min cor, 10.4 min none
  writeRaster(mosaick, filename=paste(paste("mosaicTAR", normethods[i], sep=""), "tif", sep="."))
}
# FOR TM
stacks<- EEstackWithoutMeta(tares, sat.nm="LT05")
substacks=list(stacks[[1]][[3:nlayers(stacks[[1]])]], stacks[[2]][[3:nlayers(stacks[[2]])]],
               stacks[[3]][[2:nlayers(stacks[[3]])]], stacks[[4]][[2:nlayers(stacks[[4]])]])
names(substacks)=names(stacks)
for (i in 1: length(normethods)){
  # remove band 1 from OLI stacks. THIS IS A WORKAROUND THAT NEEDS TO BE FIXED IN THE FUNCTION
  start=Sys.time()
  mosaick<- smg(inlist=substacks, method= normethods[i], refitem=1, mosaicitems=c(3,4), 
                pval.pif = 1e-02, pval.chg=0.99, QAbandname="pixel_qa", indem= NA, # demfile # only for mexico, 
                  normbands=seq(1, 5), verbose=TRUE, sensor="TM")
  print(Sys.time()-start)  # Processing time: 1.39 h imist, 1.39 hrs aRn, 24.88 min cor, 10.4 min none
  writeRaster(mosaick, filename=paste(paste("mosaicTAR", normethods[i], sep=""), "tif", sep="."))
}
rm(mosaick, stacks)
plotRGB(stack("mosaicTARnone.tif"), r=3,g=2,b=1, stretch="lin")

# For radar (optional)
for (i in 1:1){
  start=Sys.time()
  mosaick<- smg(inlist=HHHVstacks, method= normethods[i], refitem=2, mosaicitems=1, 
               normbands=seq(1, 2), pval.pif=0.05, verbose=TRUE)
  print(Sys.time()-start)  # Processing time: 1.39 hrs aRn, 24.88 min cor, 10.4 min none
  writeRaster(mosaick, filename=paste(paste("mosaicHHHVTAR", normethods[i], sep=""), "tif", sep="."))
}
rm(mosaick)

# FOR MEXICO: Both landsat scenes are from the same date in both REF and TAR years
# Therefore I will mosaic them first with no normalization
  demfile=raster("SRTM_CDMXint.tif") # For Mexico only
  stacks<- EEstackWithoutMeta(tares, sat.nm="LT05")
  #stacks=list(stack('mosaicREFnone.tif'), stack('mosaicTARnone.tif') )
  mosaick<- smg(inlist=stacks, method="none", refitem=1, mosaicitems=2, 
                QAbandname="pixel_qa", indem=NA,#demfile, 
                normbands=seq(1, 5), # (1,7), #
                verbose=TRUE, sensor="TM")
  writeRaster(mosaick, filename=paste(paste("mosaicTAR", "none", sep=""), "tif", sep="."), datatype="INT2S")
  
  # Then I will normalize the TAR mosaic based on REF mosaic
for (i in 1:2){
  start=Sys.time()
  mosaick<- smg(inlist=instacks, method= normethods[i], refitem=1, mosaicitems=2, 
                QAbandname=NA, indem=NA, normbands=seq(1, 6), verbose=TRUE, sensor="TM")
  print(Sys.time()-start)  # Processing time: 1.39 hrs aRn, 24.88 min cor, 10.4 min none
  writeRaster(mosaick, filename=paste(paste("mosaicTAR", normethods[i], sep=""), "tif", sep="."))
}

# crop all the layers to the same extent
mosaiclist=list(stack("mosaicREFnone.tif"), stack("mosaicTARaRn.tif"), 
                stack("HHHVref.tif"), stack("HHHVtar.tif"))
e=commonExtent(mosaiclist)
rm (mosaiclist)

MosaicREFnoneAllbands=stack(crop(stack("mosaicREFnone.tif"), e), crop(stack("HHHVref.tif"), e))
MosaicREFedAllbands  =stack(crop(stack("mosaicREFed.tif"), e),   crop(stack("HHHVref.tif"), e))
MosaicREFaRnAllbands =stack(crop(stack("mosaicREFaRn.tif"), e),  crop(stack("HHHVref.tif") ,e))
MosaicTARnoneAllbands=stack(crop(stack("mosaicTARnone.tif"), e), crop(stack("HHHVtar.tif"), e))
MosaicTARedAllbands  =stack(crop(stack("mosaicTARed.tif"), e),   crop(stack("HHHVtar.tif"), e))
MosaicTARaRnAllbands =stack(crop(stack("mosaicTARaRn.tif"), e),  crop(stack("HHHVtar.tif"), e))
writeRaster(MosaicREFnoneAllbands, "MosaicREFnoneAllbands.tif")
writeRaster(MosaicREFedAllbands, "MosaicREFedAllbands.tif")
writeRaster(MosaicREFaRnAllbands, "MosaicREFaRnAllbands.tif")
writeRaster(MosaicTARnoneAllbands, "MosaicTARnoneAllbands.tif")
writeRaster(MosaicTARedAllbands, "MosaicTARedAllbands.tif")
writeRaster(MosaicTARaRnAllbands, "MosaicTARaRnAllbands.tif")
plotRGB(stack("MosaicTARaRnAllbands.tif"), r=8,g=7,b=6, stretch="lin")
rm(MosaicREFnoneAllbands, MosaicREFedAllbands, MosaicREFaRnAllbands, MosaicTARnoneAllbands,
   MosaicTARedAllbands, MosaicTARaRnAllbands)


########################################################
####### Classification and validation

#### For Reference year
#mosaicREF=stack("MosaicREFaRnAllbands.tif")
selbands=c(1:5)#(3:7) # this was the best band combination for OLI
classcolname= "CLASS_NAME" 
#"class"#
normethods=c("imist", "ed", "none")
calibvalid=crossvalidsp(dataset=calibdata, classcolname=classcolname, validprop=0.3)

for (i in 1:length(normethods)){
  start=Sys.time()
  #name=paste(paste("mosaicREF", normethods[i], sep=""), "Allbands.tif", sep="")
  name=paste(paste("mosaicREF", normethods[i], sep=""), ".tif", sep="")
  mosaicREF=stack(name)[[selbands]]
  tm=tmg(inrast=mosaicREF, calibObject=calibvalid[[1]], ntrees=500, 
         classcolname=classcolname, 
         validObject=calibvalid[[2]],  maptype="class", savefiles=TRUE, 
         outfile=paste("classREF", normethods[i], sep=""), 
         plot=TRUE, tmproduce=TRUE)
  Sys.time()-start  #Time difference of 57.60239 mins (excluding stack2df())
  # Time difference of 1.509214 hours (including stack2df())
  #Time difference of 30.14512 mins (excluding tmproduce)
}
#plotRGB(mosaicREF, r=9,g=8,b=4, stretch="lin")

##### For validation year
iter=10
selbandsval= c(1:5) # for TM or ETM only
  #c(3:7)   #for OLI only
  
#setwd(savepath)
#classcolname="CLASS_NAME"
  #"class"# 
for (j in 1:length(normethods)){
  rfname=paste(paste("classREF", normethods[j], sep=""),"RF.RData", sep="_")
  #rfname=paste(paste("mosaicREF", normethods[3], sep=""),"RF.RData", sep="_") # FOR MEXICO
  
  load(rfname)
  #attributes(class.RF$terms)$term.labels
  #levels(validdata[[classcolname]])=rownames(class.RF$confusion)# This is to match the names of the classes in calibdata and validdata
  #setwd(path)
  #mosname=paste(paste("MosaicTAR", normethods[j], sep=""), "Allbands.tif", sep="")
  mosname=paste(paste("mosaicTAR", normethods[j], sep=""), ".tif", sep="")
  mosaicTAR=stack(mosname)[[selbandsval]]
  
  #setwd(savepath)
  for (i in 1:iter){
    outname=paste(paste(paste("classTAR", normethods[j], sep="_"), "iter", sep="_"), i, sep="_")
    calibvalid1=crossvalidsp(dataset=validdata, classcolname=classcolname, validprop=0.3)
    start=Sys.time()
    tm=tmg(inrast=mosaicTAR, calibObject=class.RF, classcolname=classcolname, 
          validObject=calibvalid1[[1]],  maptype="class", savefiles=TRUE, 
           outfile=outname, plot=TRUE, tmproduce=FALSE)
    print(Sys.time()-start)
  }
}


### Calculate statistics from all iteration
#setwd(paste(getwd(), "outputsBandsAllV2", sep="/"))
imiststats <-AccuStats(prefix="classTAR_imist_iter", suffix="accumatrix.RData", iter=iter)
edStats <-AccuStats(prefix="classTAR_ed_iter", suffix="accumatrix.RData", iter=iter)
noneStats <-AccuStats(prefix="classTAR_none_iter", suffix="accumatrix.RData", iter=iter)
overalldata <- data.frame(
  name=normethods,
  value=c((1-mean(imiststats[[3]])), (1-mean(edStats[[3]])), (1-mean(noneStats[[3]]))),
  sd=c(sd(imiststats[[3]]), sd(edStats[[3]]), sd(noneStats[[3]]))
  )

crossusers=data.frame(cbind(imiststats[[4]][,c(1,2)],edStats[[4]][,2],noneStats[[4]][2]),
                      stringsAsFactors=FALSE)
names(crossusers)=c("classnames", "imist", "ed", "none")
crossproducers=cbind(imiststats[[4]][,c(1,4)],edStats[[4]][,4],noneStats[[4]][4])
names(crossproducers)=c("classnames", "imist", "ed", "none")

pdf("ProducersUsersPuc.pdf", width=4,height=4,paper='special') 
par(mar = c(7, 4, 2, 2) + 0.2)
barplot(rbind(crossproducers[,2], crossproducers[,3],crossproducers[,4]),
        col=c("aquamarine3","coral", "cornflowerblue"), 
        names.arg=crossproducers$classnames,  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Producers accuracy")
legend("top", legend=c("imist", "ed", "none"), box.lty=0, bg= "transparent",
       col=c("aquamarine3","coral", "cornflowerblue"), lty=1:2, cex=0.8)
barplot(rbind(crossusers[,2], crossusers[,3],crossusers[,4]),col=c("aquamarine3","coral", "cornflowerblue"), 
        names.arg=crossusers$classnames,  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Users accuracy")
legend("top", legend=c("imist", "ed", "none"), box.lty=0, bg= "transparent",
       col=c("aquamarine3","coral", "cornflowerblue"), lty=1:2, cex=0.8)
ggplot(overalldata) +
  geom_bar(aes(x=name, y=value), stat="identity", fill=c("aquamarine3","coral", "cornflowerblue"), alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="gray28", alpha=0.9, size=1.3)+
  coord_cartesian(ylim = c(0, .35))
dev.off()
#############################################################
# POST CLASSIFICATION ANALYSIS

######################################################
# Post classification class merging
path=#("X:/VictorShare/aRnFiles/Orinoquia")
  #("X:/VictorShare/aRnFiles/Pucallpa")
  ("X:/VictorShare/aRnFiles/Mexico")#
setwd(path)

fromtoclass=c(4,6, 5,10, 9,10, 8,7, 11,10, 12,6, 17,10)# ORINOQUIA
fromtoclass=c(1,6, 2,5, 4,9, 15,10)        # PUCALLPA
fromtoclass=c(9,1)                                     # MEXICO
#fromtoclass=c(5,2, 7,6)                                     # MEXICO
fromtomat=matrix(fromtoclass, ncol=2, byrow=T)

normethods=c("imist", "ed", "none")
iter=10
for (j in 1:length(normethods)){
  for(i in 1:iter){
    load(paste((paste(paste(paste("classTAR", normethods[j], sep="_"), 
                            "iter", sep="_"), i, sep="_")), "accumatrix.RData", sep="_"))
    test=conmatrixmerge(reclassx=conmatrix$table, rcl=fromtomat)
    if(i==1){
      overall=test$overall
      users=test[[2]][,1]
      producers=test[[2]][,2]
    } else {
      overall=c(overall, test$overall)
      users=cbind(users, test[[2]][,1])
      producers=cbind(producers, test[[2]][,2])
  }
  }
  if (j==1){
    alloverallmean=mean(overall)
    alloverallsd=sd(overall)
    allusersmean=apply(users, 1, "mean")
    alluserssd=apply(users, 1, "sd")
    allproducersmean=apply(producers, 1, "mean")
    allproducerssd=apply(producers, 1, "sd")
  }else{
    alloverallmean=c(alloverallmean, mean(overall))
    alloverallsd=c(alloverallsd, sd(overall))
    allusersmean=cbind(allusersmean,apply(users, 1, "mean"))
    alluserssd=cbind(alluserssd, apply(users, 1, "sd"))
    allproducersmean=cbind(allproducersmean, apply(producers, 1, "mean"))
    allproducerssd=cbind(allproducerssd, apply(producers, 1, "sd"))
  }
}

# compile data
overalldata <- data.frame(
  name=normethods,
  value=(1-alloverallmean)*100,
  sd=alloverallsd*100
)

producersdata <- data.frame(
  name=c(rownames(allproducersmean), rownames(allproducersmean),rownames(allproducersmean)),
  type=c(rep(normethods[1], nrow(allproducersmean)), 
         rep(normethods[2], nrow(allproducersmean)), 
         rep(normethods[3], nrow(allproducersmean))),
  value=c((1-allproducersmean[,1])*100, (1-allproducersmean[,2])*100, (1-allproducersmean[,3])*100),
  sd=(c(allproducerssd[,1], allproducerssd[,2], allproducerssd[,3]))*100
)

usersdata <- data.frame(
  name=c(rownames(allusersmean), rownames(allusersmean),rownames(allusersmean)),
  type=c(rep(normethods[1], nrow(allusersmean)), 
         rep(normethods[2], nrow(allusersmean)), 
         rep(normethods[3], nrow(allusersmean))),
  value=(1-c(allusersmean[,1], allusersmean[,2], allusersmean[,3]))*100,
  sd=c(alluserssd[,1], alluserssd[,2], alluserssd[,3])*100
)

pdf("ProducersUsersMergeMex.pdf", width=4,height=4,paper='special') 
par(mar = c(7, 4, 2, 2) + 0.2)
ggplot(producersdata, aes(x = name, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha=0.7) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),                            
                width = 0.2,
                position = position_dodge(0.9))+ 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_fill_manual(values=c("aquamarine3","coral", "cornflowerblue"))

ggplot(usersdata, aes(x = name, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha=0.7) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),                            
                width = 0.2,
                position = position_dodge(0.9))+ 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_fill_manual(values=c("aquamarine3","coral", "cornflowerblue"))
ggplot(overalldata) +
  geom_bar(aes(x=name, y=value), stat="identity", fill=c("aquamarine3","coral", "cornflowerblue"), alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="gray28", alpha=0.9, size=1.3)+
  coord_cartesian(ylim = c(0, 20))
dev.off()






# ASSESS OPTIMAL RESAMPLING ITERATIONS
optstats=aRnStats
pref="classTAR_aRn_iter"
for (i in 2: iter){
  optstats <-AccuStats(prefix=pref, suffix="accumatrix.RData", iter=i)
  test=cbind(mean(optstats[[3]]), sd(optstats[[3]]), sd(optstats[[3]])/mean(optstats[[3]]))
  if (i==2){metrics=test} else {metrics=rbind(metrics, test)}
}
metrics=data.frame(cbind(seq(2,iter), metrics))
names(metrics)=c("iter", "mean", "sd", "cv")
plot(metrics$iter, metrics$mean, type="l")
plot(metrics$iter, metrics$sd, type="l")
plot(metrics$iter, metrics$cv, type="l")

#### This is to produce a correlation matrix between the accuracy of different categories
# Transpose data frames
allusersT=as.data.frame(t(as.matrix(allusers[,-1])))
names(allusersT)=paste(as.character(allusers[,1]), "U", sep="_")

allproducersT=as.data.frame(t(as.matrix(allproducers[,-1])))
names(allproducersT)=paste(as.character(allproducers[,1]), "P", sep="_")

allaccu=cbind(allusersT, allproducersT)

# Correlation matrix
coran=cor(allaccu, use="complete.obs", method="pearson")
write.csv(coran, file='coran.csv')

plot(allaccu$River_U, allaccu$Lake_P)

### THIS IS TO ASSESS HOW MANY ITERATIONS ARE NEEDED FOR STABILIZING ACCURACY
for (iter in 3:10){
  aRnStats <-AccuStats(prefix="classTAR_aRn_iter", suffix="accumatrix.RData", iter=iter)
  stats=cbind(mean(aRnStats[[3]]), sd(aRnStats[[3]]))
  if (iter==3){cumstats=stats} 
  else {
    cumstats=rbind(cumstats, stats)}
}

attributes(class.RF$terms)$term.labels



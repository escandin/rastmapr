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
#("X:/VictorShare/aRnFiles/Orinoquia")
#("X:/VictorShare/aRnFiles/Mexico")#
("X:/VictorShare/s3dFiles/MontesTest")#

savepath=(paste(path, "outputsAllBands", sep="/"))
savepath1=(paste(path, "outputsBands3_7", sep="/"))
setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

##### 1. DATA INPUTS

calibdata=#readOGR(".","Training_Orinoquia_2016_3.0")# Orinoquia
#readOGR(".", "Training_Pucallpa_2015_6.0")    #Pucallpa
readOGR(".", "Training_DF_2018_6.0")               # Mexico

validdata=#readOGR(".", "Training_Orinoquia_2013_3.0")  #Orinoquia
#readOGR(".", "Training_Pucallpa_2011_6.6")           #Pucallpa
readOGR(".", "Training_DF_2010_6.0") # for Mexico    # Mexico

# Check if levels are the same
levels(validdata$CLASS_NAME)==levels(calibdata$CLASS_NAME)
levels(validdata$class)==levels(calibdata$class)

# FOR ORINOQUIA
  levels(validdata$CLASS_NAME)=levels(calibdata$CLASS_NAME)
  validdata@data$CLASS_NAME[which(validdata@data$CLASS_NAME=="B urn Scar (Old)")]="Burn Scar (Old)"
  validdata@data$CLASS_NAME=droplevels(validdata@data$CLASS_NAME, "B urn Scar (Old)")
  validdata@data$CLASS_NAME[which(validdata@data$CLASS_NAME=="Burn Scar (old)")]="Burn Scar (Old)"
  validdata@data$CLASS_NAME=droplevels(validdata@data$CLASS_NAME, "Burn Scar (old)")
  

# FOR PUCALLPA
  validdata=spTransform(validdata, crs(calibdata))
  
# FOR MEXICO # correct mislabeling of sand class
  #calibdata@polygons=calibdata@polygons[which(calibdata@data$class!="sand")]
  #calibdata@data$class[which(calibdata@data$class=="Baresoil")]="baresoil"
  #calibdata@data$class=droplevels(calibdata@data$class, "Baresoil")
  
  #calibdata@data$class[which(calibdata@data$class=="Burn Scar (new)")]="Burn Scar (New)"
  #calibdata@data$class=droplevels(calibdata@data$class, "Burn Scar (new)")
  
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
  
########################################################
#### CLASS REFERENCE YEAR

#settngs for landsat only
selbands=c(1:5)#
suffix=".tif"
preffix="classREF"

#settings for landsat and ALOS
#selbands=c(1:7)
#suffix="Allbands.tif"
#preffix="classREFAllbands"

classcolname=# "CLASS_NAME" 
"class"#
normethods=c("imist", "ed", "cor", "none")
calibvalid=crossvalidsp(dataset=calibdata, classcolname=classcolname, validprop=0.3)

for (i in 1:length(normethods)){
  start=Sys.time()
  name=paste(paste("mosaicREF", normethods[i], sep=""), suffix, sep="")
  mosaicREF=stack(name)[[selbands]]
  if(i==1){mapping=TRUE}else{mapping=FALSE}
  tm=tmg(inrast=mosaicREF, calibObject=calibvalid[[1]], ntrees=500, 
         classcolname=classcolname, 
         validObject=calibvalid[[2]],  maptype="class", savefiles=TRUE, 
         outfile=paste(preffix, normethods[i], sep=""), 
         plot=TRUE, tmproduce=mapping)
  Sys.time()-start  #Time difference of 57.60239 mins (excluding stack2df())
  # Time difference of 1.509214 hours (including stack2df())
  #Time difference of 30.14512 mins (excluding tmproduce)
}
#plotRGB(mosaicREF, r=9,g=8,b=4, stretch="lin")

########################################################
##### CLASS FOR VALIDATION YEAR
iter=3

#settings for landsat only
preffixtar="classTAR"
#settings for landsat and ALOS
#preffixtar="classTARAllbands"

for (i in 1:length(normethods)){
  rfname=paste(paste(preffix, normethods[i], sep=""),"RF.RData", sep="_")
  load(rfname)
  #attributes(class.RF$terms)$term.labels
  #levels(validdata[[classcolname]])=rownames(class.RF$confusion)# This is to match the names of the classes in calibdata and validdata
  #setwd(path)
  mosname=paste(paste("mosaicTAR", normethods[i], sep=""), suffix, sep="")
  mosaicTAR=stack(mosname)[[selbands]]
  
  #setwd(savepath)
  for (j in 1:iter){
    outname=paste(paste(paste(preffixtar, normethods[i], sep="_"), "iter", sep="_"), j, sep="_")
    calibvalid1=crossvalidsp(dataset=validdata, classcolname=classcolname, validprop=0.3)
    start=Sys.time()
    if(i==1 & j==1){mapping=TRUE}else{mapping=FALSE}
    tm=tmg(inrast=mosaicTAR, calibObject=class.RF, classcolname=classcolname, 
          validObject=calibvalid1[[1]],  maptype="class", savefiles=TRUE, 
           outfile=outname, plot=TRUE, tmproduce=mapping)
    print(Sys.time()-start)
  }
}

########################################################
##### UPDATED CLASS VALIDATION YEAR
#ORINOQUIA
mask=merge(raster(paste(paste("LC08_L1TP_005057_20140128_20170426_01_T1_removeChgmsk", 3, sep=""),
                        "tif", sep=".")),
           raster(paste(paste("LC08_L1TP_006057_20140220_20170425_01_T1V2_removeChgmsk", 10, sep=""),
                        "tif", sep=".")))

# MEXICO
mask=merge(raster(paste(paste("LT05_L1TP_026046_20100205_20161016_01_T1_removeChgmsk", 15, sep=""),
                        "tif", sep=".")),
           raster(paste(paste("_removeChgmsK", 14, sep=""),
                        "tif", sep=".")))

tmref=raster("classREFimist.tif")
tmtar=raster(paste(paste("classTAR_imist_iter", 1, sep="_"),"tif", sep="."))
inlist=(list(tmref, tmtar, mask))
e=commonExtent(inlist)
tmref=crop(tmref, e)
mask=crop(mask, e)

for (i in 2:iter){
  tmtar=raster(paste(paste("classTAR_imist_iter", i, sep="_"),"tif", sep="."))
  tmtar=crop(tmtar, e)
  # Update target year with the no change mask
  classTARupdated=lcupdate(refclass=tmref, tarclass=tmtar, nochmsk=mask)
  writeRaster(classTARupdated, file=paste(paste("classTARupdated", i, sep="_"), "tif", sep="."))
  validdf=stack2df(inrast=classTARupdated, invec=validdata, classcolname=classcolname)
  validdf$class_ID=as.factor(validdf$class_ID)
  validdf$layer=as.factor(validdf$layer)
  levels(validdf$layer)=levels(validdf$class_name)
  
  conmatrix=confusionMatrix(validdf$class_name ,validdf$layer)
  save(conmatrix, file=paste(paste(paste("classTAR_updated", i, sep="_"), "accumatrix", sep="_"),
                        "RData", sep="."))
  print(paste(i, "iterations processed", sep=" "))
}




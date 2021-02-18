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
("X:/VictorShare/aRnFiles/Mexico")#

savepath=(paste(path, "outputsBandsAll", sep="/"))
savepath1=(paste(path, "outputsBands3_7", sep="/"))
setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

##### 1. DATA INPUTS
tares <- list.files('.', pattern='tar.gz')
#stack with ALOS PALSAR
HHref=raster("palsar_HH_2016_Orinoquia_Spk.tif")
#raster("palsar_HH_2015_Peru_Spk.tif")
( "palsar_HH_2018_CDMX_Spk.tif" )

HVref=raster("palsar_HV_2016_Orinoquia_Spk.tif")
#raster("palsar_HV_2015_Peru_Spk.tif")
("palsar_HV_2018_CDMX_Spk.tif")

HHtar=raster("palsar_HH_2015_Orinoquia_Spk.tif")
#raster("palsar_HH_2010_Peru_Spk.tif")
("palsar_HH_2010_CDMX_Spk.tif")

HVtar=raster("palsar_HV_2015_Orinoquia_Spk.tif")
#raster("palsar_HV_2010_Peru_Spk.tif")
("palsar_HV_2010_CDMX_Spk.tif")


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
normethods=c("imist", "ed", "cor", "none")
normethods=c("imist", "cancor", "pca")
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

# # For radar (optional)
# for (i in 1:1){
#   start=Sys.time()
#   mosaick<- smg(inlist=HHHVstacks, method= normethods[i], refitem=2, mosaicitems=1, 
#                normbands=seq(1, 2), pval.pif=0.05, verbose=TRUE)
#   print(Sys.time()-start)  # Processing time: 1.39 hrs aRn, 24.88 min cor, 10.4 min none
#   writeRaster(mosaick, filename=paste(paste("mosaicHHHVTAR", normethods[i], sep=""), "tif", sep="."))
# }
# rm(mosaick)

# Create no change masks for the non-intersecting area in adjacent image
normitem=#3 #Orinoquia
4 # Pucallpa, Mexico
  
#substacks=list(stack("mosaicREFimist.tif"), stacks[[normitem]][[3:7]]) # For OLI
substacks=list(stack("mosaicREFimist.tif"), stacks[[normitem]][[2:6]]) # For TM
names(substacks)=c("mosaicREF", paste(names(stacks[normitem]),'V2', sep=""))
names(substacks[[1]])=names(substacks[[2]])
rm(stacks)
i=1 # selects normalization method for mosaic creation
# remove band 1 from OLI stacks. THIS IS A WORKAROUND THAT NEEDS TO BE FIXED IN THE FUNCTION
mosaick<- smg(inlist=substacks, method= normethods[i], refitem=1, mosaicitems=2, 
              pval.pif = 1e-02, pval.chg=0.99, QAbandname=NA, indem= NA, # demfile # only for mexico, 
              normbands=seq(1, 5), verbose=TRUE, sensor="OLI")
writeRaster(mosaick, filename=paste(paste(names(substacks)[[2]], normethods[i], sep=""), "tif", sep="."))


#######################################################################
# Produce combined mosaics using Landsat + ALOS-PALSAR
# crop all the layers to the same extent
mosaiclist=list(stack("mosaicREFimist.tif"), stack("mosaicTARimist.tif"), 
                stack("HHHVref.tif"), stack("HHHVtar.tif"))
e=commonExtent(mosaiclist, method="intersection")
rm (mosaiclist)

#MosaicREFnoneAllbands=stack(crop(stack("mosaicREFnone.tif"), e), crop(stack("HHHVref.tif"), e))
#MosaicREFedAllbands  =stack(crop(stack("mosaicREFed.tif"), e),   crop(stack("HHHVref.tif"), e))
MosaicREFimistAllbands =stack(crop(stack("mosaicREFimist.tif"), e),  crop(stack("HHHVref.tif") ,e))
#MosaicTARnoneAllbands=stack(crop(stack("mosaicTARnone.tif"), e), crop(stack("HHHVtar.tif"), e))
#MosaicTARedAllbands  =stack(crop(stack("mosaicTARed.tif"), e),   crop(stack("HHHVtar.tif"), e))
MosaicTARimistAllbands =stack(crop(stack("mosaicTARimist.tif"), e),  crop(stack("HHHVtar.tif"), e))
#writeRaster(MosaicREFnoneAllbands, "MosaicREFnoneAllbands.tif")
#writeRaster(MosaicREFedAllbands, "MosaicREFedAllbands.tif")
writeRaster(MosaicREFimistAllbands, "MosaicREFimistAllbands.tif")
#writeRaster(MosaicTARnoneAllbands, "MosaicTARnoneAllbands.tif")
#writeRaster(MosaicTARedAllbands, "MosaicTARedAllbands.tif")
writeRaster(MosaicTARimistAllbands, "MosaicTARimistAllbands.tif")
plotRGB(stack("MosaicTARimistAllbands.tif"), r=8,g=7,b=6, stretch="lin")
rm(MosaicREFnoneAllbands, MosaicREFedAllbands, MosaicREFaRnAllbands, MosaicTARnoneAllbands,
   MosaicTARedAllbands, MosaicTARaRnAllbands)




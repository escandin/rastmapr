# This code processes a workflow for producing time series of land cover classifications
# and perform accuracy assessment using landsat 8 collection 1 level 2 data

####### ENVIRONMENT SETTINGS
library(RStoolbox)
library(raster)
library(rgdal)
library(randomForest)
library(caret)
library(maptools)
library(stringr)
library(ggplot2)
library(rgeos)
library(RColorBrewer)
library(R.utils)
sourcedir <- '/Users/tug61163/Documents/Repositories/rastmapr/rasterMapR/R'
sourceDirectory(sourcedir, modifiedOnly=FALSE)

# MAC
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Orinoquia"
#path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Pucallpa"
#path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Mexico"

setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

##### 1. DATA INPUTS
shapes <- str_remove(list.files('.', pattern='.shp'),".shp")
dir()
calibdata=readOGR(".", shapes[2])# Orinoquia
#readOGR(".", "Training_DF_2018_6.0")               # Mexico
#readOGR(".", "Training_Pucallpa_2015_5.0")    #Pucallpa

validdata= readOGR(".", shapes[1]) #Orinoquia
#readOGR(".", "Training_DF_2010_6.0") # for Mexico    # Mexico
#readOGR(".", "Training_Pucallpa_2011_5.0_")           #Pucallpa

# Check if levels are the same. 
sort(unique(validdata$CLASS_NAME))==sort(unique(calibdata$CLASS_NAME))

# FOR ORINOQUIA
  validdata@data$CLASS_NAME[which(validdata@data$CLASS_NAME=="B urn Scar (Old)" )]="Burn Scar (Old)"
  validdata@data$CLASS_NAME[which(validdata@data$CLASS_NAME=="Burn Scar (old)" )]="Burn Scar (Old)"
  
# FOR MEXICO # check labeling
  sort(unique(calibdata@data$class)) == sort(unique(validdata@data$class))
  
  calibdata@polygons=calibdata@polygons[which(calibdata@data$class!="Airport_1")]
  calibdata@data=calibdata@data[which(calibdata@data$class!="Airport_1"),]
  #calibdata@data$class=droplevels(calibdata@data$class, "Airport_1")
  
  calibdata@polygons=calibdata@polygons[which(calibdata@data$class!="Airport_2")]
  calibdata@data=calibdata@data[which(calibdata@data$class!="Airport_2"),]
 # calibdata@data$class=droplevels(calibdata@data$class, "Airport_2")
########################################################
####### Classification and validation

#### For Reference year
grds=  list.files('.', pattern='.grd')
mosaicREF=stack(grds[2]) # ORINOQUIA
#plotRGB(mosaicREF, r=5, g=4, b=3, stretch='lin')

#Select only polygons within the area of influence of the image: 
#extent(calibdata)
#extent(mosaicREF)

classcolname= "CLASS_NAME" #ORINOQUIA
#"class"# Mexico
e=commonExtent(inlist=list(calibdata, mosaicREF))
p <- as(e, 'SpatialPolygons') 
crs(p)=crs(calibdata)
selitems=gIntersects(calibdata,p, byid=TRUE)
calibcrop <- calibdata[as.vector(selitems),]
length(calibcrop[[classcolname]])
length(calibdata[[classcolname]])

#selbands=c(1:5)#(3:7) # this was the best band combination for OLI

normethods="s3d"#c("imist", "ed", "none")
calibvalid=crossvalidsp(dataset=calibcrop, classcolname=classcolname, validprop=0.3)

for (i in 1:length(normethods)){
  start=Sys.time()
  #name=paste(paste("mosaicREF", normethods[i], sep=""), "Allbands.tif", sep="")
  #name=paste(paste("mosaicREF", normethods[i], sep=""), ".tif", sep="")
  #mosaicREF=stack(name)[[selbands]]
  tm=tmg(inrast=mosaicREF, calibObject=calibvalid[[1]], ntrees=500, 
         classcolname=classcolname, 
         validObject=calibvalid[[2]],  maptype="class", savefiles=TRUE, 
         outfile=paste("classREF", normethods[i], sep=""), 
         plot=TRUE, tmproduce=TRUE)
  Sys.time()-start  #Time difference of 57.60239 mins (excluding stack2df())
  # Time difference of 1.509214 hours (including stack2df())
  #Time difference of 30.14512 mins (excluding tmproduce)
}
# Generate random colors
n=length(unique(calibdata[[classcolname]]))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
colors=sample(col_vector, length(unique(calibdata[[classcolname]])))
plot(tm, col=colors)
#plotRGB(mosaicREF, r=9,g=8,b=4, stretch="lin")

# Reclassify
recmat= c(4,8, 6,14, 7,14, 9,14, 10,5, 12,2, 15,14)
recmat=matrix(recmat, ncol=2, byrow=TRUE)
conmatrixrec=conmatrixmerge(conmatrix$table, recmat)
tmrec=reclassify(tm, recmat)
plot(tmrec,colors[1:8])
writeRaster(tmrec, "classREFreclass.tif", datatype="INT1U")

###### RADNORM AND CHANGE PER CLASSES
# Open and stack images
grds=  list.files('.', pattern='_01_T1_mskd_.grd')
filenames=c("LC08_L1TP_005057_20140128_20170426_01_T1_mskd_.grd", 
            "LC08_L1TP_005057_20160118_20180130_01_T1_mskd_.grd")
instacks=list(stack(filenames[1]), 
              stack(filenames[2]))
e=commonExtent(instacks)
instacks=Map(function(x)(crop(x, e)), instacks)
names(instacks)=str_remove(filenames, "_.grd")

#Open classified image
classimg=raster("classREFs3d.tif")#("classREFreclass.tif")
classimg=crop(classimg, e)
plot(classimg)
classval=unique(classimg)
i=1 #1 New burn scars
    #2 old burn scars
    #5 forest
    #8 grassland
    #10 Riparian
#for (i in 1:  length(classval)){
  v <- raster::getValues(classimg)
  v[v!=classval[which(classval==i)]]=NA
  v[v==classval[which(classval==i)]]=1
  msk <- raster::setValues(classimg, v)
  writeRaster(msk, paste(paste("MaskClass_", i),"tif", sep="."), datatype="INT1U")
  instackmskd=Map(function(x)(mask(x, msk)), instacks)
  
  s3dmod=s3d(strips=instackmskd, thres=1e-2, distype="gamma",
             pval.pif=1e-3,  pval.chg=0.99, cca=FALSE, 
             prefix=paste(names(instackmskd)[1], paste("class", i, sep=""), sep="_"))
  save(s3dmod, file=paste(paste(paste(str_remove(names(instackmskd)[2], "_.grd"), i, sep="_class"),'s3d_gamma', sep="_"),  
                          "RData", sep="."))
#}
# Compare normalizations
s3dlist=list()
j=1
for (i in c(1,2,5,8,10)){
  load(paste(paste("LC08_L1TP_005057_20160118_20180130_01_T1_mskd_class", 
                    i, sep=""),"s3d_gamma.RData", sep="_"))
  s3dlist[[j]]=s3dmod
  names(s3dlist[[j]])=names(s3dmod)
  j=j+1
  }  

scatterclass=4
abclass=3
plot(s3dlist[[scatterclass]]$data[[4]]$y~s3dlist[[scatterclass]]$data[[4]]$x)
  abline(a=s3dlist[[abclass]]$paramstats$lmparam$intercept[length(s3dlist[[abclass]]$paramstats$lmparam$intercept)],
         b=s3dlist[[abclass]]$paramstats$lmparam$slope[length(s3dlist[[abclass]]$paramstats$lmparam$slope)])

#Join all radnorm datapoints in one list
scatterlist=s3dlist[[1]]$data
for (j in 1:6){ # iterate over bands#
  for (i in 1:5){  #iterate over the rest of the classes
    #scatterlist[[j]]$class=rep(1, nrow(scatterlist[[j]]))
    s3dlist[[i]]$data[[j]]$class=rep(i, nrow(s3dlist[[i]]$data[[j]]))
    if (i==1){scatterlist[[j]]=s3dlist[[i]]$data[[j]]} else {
    scatterlist[[j]]=rbind(scatterlist[[j]],s3dlist[[i]]$data[[j]])}
  }
}
par(mfrow=c(2,3))
for(band in 1:6){
with(scatterlist[[band]], plot(x, y, col = as.factor(class)))
}
colores=c('black', "red", "green", "cyan")
load("LC08_L1TP_005057_20160118_20180130_01_T1_mskd_class8_s3d_gamma.RData")
grass=s3dmod$paramstats$lmparam[which(s3dmod$paramstats$lmparam$iter==5),]
for (abclass in c(1,2,3,5)){
  abline(a=s3dlist[[abclass]]$paramstats$lmparam$intercept[length(s3dlist[[abclass]]$paramstats$lmparam$intercept)],
       b=s3dlist[[abclass]]$paramstats$lmparam$slope[length(s3dlist[[abclass]]$paramstats$lmparam$slope)], col=colores[abclass])
  abline(a=grass$intercept[which(grass$band==abclass)], b=0.9903301, col="blue")
}
plotRGB(instackmskd[[2]], r=4,g=3,b=2, stretch="lin")

for(band in 1:6){
  with(scatterlist[[band]], plot(x, y, col = as.factor(class)))
  dif=scatterlist[[band]]$x-scatterlist[[band]]$y
  hist(dif)
  difsel=dif[which(dif<0)]
  difsel=difsel[which(difsel>0)]
  hist(difsel)
}

ndvistack=Map(function(x)(spectralIndices(x, 
                                          blue=1, green=2, red=3, nir=4, indices="NDVI")),instackmsk)
ndvidiff=ndvistack[[1]]-ndvistack[[2]]
hist(ndvistack[[1]])
hist(ndvistack[[1]])
hist(ndvidiff[[1]])

par(mfrow=c(3,2))
for (i in 1:6){
#s3dmod$data[[i]]=rbind(s3dmod$data[[i]],s3dmodClass8$data[[i]])
  plot(rbind(s3dmod$data[[i]],s3dmodClass8$data[[i]]))
}

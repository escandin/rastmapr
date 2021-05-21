smg <- structure(function #Seamless mosaic generator
### This function produces a cloud masked and radiometrically normalized 
  # raster mosaic from landsat images.
                       ##details<< 
(
    inlist, ##<<\code{list}. List of input rasters
    cca="FALSE", ##<<\code{logical}. 
    method="none", ##<<\code{character}. "none" for no normalization,
                   ##"aRn" for std diff normalization or "cor" for
                   ##pearson correlation (pending)
    refitem=NA, ##<<\code{}. Item in the list corresponding to the
                ##raster to use as reference
    mosaicitems=c(1,3), ##<<\code{numeric}. Items to be normalized and
                        ##mosaiced. If mosaicitems=1, no mosaics are
                        ##built.
    QAbandname="pixel_qa", ##<<\code{character}
    cloudbuff=NA, ##<<\code{}.
    savelcloudmsk=TRUE, ##<<\code{logical}.
    indem=NA, ##<<\code{}.
    normbands=seq(1,7), ##<<\code{numeric}.
    pval.pif=1e-02, ##<<\code{numeric}.
    pval.chg=0.99, ##<<\code{numeric}.
    sensor="OLI", ##<<\code{character}
    verbose=TRUE, ##<<\code{logical}.
    norm.ext=NULL ##<<\code{}.
){
  # smg (seamless mosaic generator) produces a cloud masked and radiometrically normalized 
  # raster mosaic from landsat images.
  # inlist: list of input rasters
  # method: "none" for no normalization, "aRn" for std diff normalization or "cor" for pearson correlation (pending)
  # refitem= item in the list corresponding to the raster to use as reference 
  # for the radiometric normalization. This image should overlap at least partially
  # with all mosaicitems in order to work.
  # mosaicitems: items to be normalized and mosaiced. If mosaicitems=1, no mosaics are built
  # QAband (inherited): position of the landsat QA band in each rasterstack contained in 
  #         inlist. If QA=NA, no cloud masking is performed
  # sensor (inherited): Sensor to encode. Options: c("OLI", "TIRS", "ETM+", "TM", "MSS")
  # transform: whether to apply a canonical transformation to the data before normalizing
  # required functions and packages
  require(rgdal, raster, RStoobox, aRn)
  is.integer0 <- function(x){
    is.integer(x) && length(x) == 0L
  }
  validmethods=c("s3d", "ed", "cor", "sam", "none")
  if (is.integer0(which(validmethods==method))){
    stop("please provide a valid method: s3d, ed, cor, sam, none")
  }
  
  stacknames=names(inlist)
  maskitems<-unique(c(refitem, mosaicitems))
  maskitems<-maskitems[!is.na(maskitems)]
 
  
  # select items to mask from the list
  #bnames=names(inlist[[refitem]])
  
  #Mask clouds and subset stacks to include only normbands:
  for(i in maskitems){
    layernames=names(inlist[[i]])
    # perform cloud masking
    if(!is.na(QAbandname)){
      print('masking clouds')
      msk=RStoolbox::classifyQA(inlist[[i]][[QAbandname]], sensor=sensor) #do Sys.time to know whether it's faster with get/setValues
      if(is.numeric(cloudbuff)){
         msk=buffer(msk, width=cloudbuff)}
      msk=namaskplus(msk, 1, NA)
      #msk=namask(msk) # see function above
      if(savelcloudmsk==TRUE){
      writeRaster(msk, filename=paste(paste("CloudMsk", stacknames[i], sep="_"), "tif", sep="."))
      }
    }
    inlist[[i]]=inlist[[i]][[normbands]] # only select bands to normalize
    names(inlist[[i]])=layernames[normbands]
    if(!is.na(QAbandname)){
      inlist[[i]]=raster::mask(inlist[[i]], msk)
      rm(msk)
      print(paste(paste("item", i, sep=" "), "masked", sep=" "))
    }
  }
  
  # perform topographic correction
  if(class(indem)[1]=="RasterLayer"){
    print ("topographic correction started")
    for (j in maskitems){ 
      # ensure that the dem is aligned with the raster stack
      if(as.character(crs(indem))!=as.character(crs(inlist[[j]]))){
        print("reprojecting raster dem")
        indem<-raster::projectRaster(indem, inlist[[j]], 
                                     datatype=dataType(inlist[[j]][[normbands[1]]]))
      } else if(as.character(extent(indem))!=as.character(extent(inlist[[j]]))){
        print("resampling raster dem") 
        # OJO: dataType is not working. It only allows 
        subdem<-raster::resample(indem, inlist[[j]], 
                                 datatype=dataType(inlist[[j]][[normbands[1]]]))
      } else {
        subdem<-indem
      }
      # This is a workaround to change the data type. 
      # The format of the file changes to FLT4S after resampling regardless of
      # whether I add the datatype argument or change the global rasterOptions 
      # writeRaster(subdem, paste(tmpDir(), "subdem.tif", sep=""), 
      #            datatype=dataType(inlist[[j]][[normbands[1]]]))
      #subdem=raster(paste(tempdir, "subdem.tif", sep="/"))
      print("producing common mask between dem and raster stack")
      msk<-stackmask(stack(inlist[[j]], subdem), maskrast=TRUE)
      inlist[[j]]<-msk[[1]][[1:nlayers(stack(inlist[[j]]))]]
      subdem<-msk[[1]][[nlayers(stack(inlist[[j]]))+1]]
      rm(msk)
      #print("extracting topographic variables")
      #slope<-raster::terrain(dem, opt='slope')
      #aspect<-raster::terrain(dem, opt='aspect')
      #terrainvar<-raster::stack(slope,aspect)
      #names(terrainvar)<-c('slope', 'aspect')
      #rm(slope,aspect,dem)
      print ("performing topographic correction")
      meta<-RStoolbox::readMeta(paste(names(inlist)[[j]], "xml", sep="."))
      inlist[[j]]<-RStoolbox::topCor(inlist[[j]], dem=subdem, metaData=meta, method="C")
      rm(subdem)
      # THIS PART OF THE WORKAROUND: remove temporary dem file so that a new one can be saved
      #file.remove(paste(tempdir, "subdem.tif", sep="/"))
      print(paste(paste("item", j, sep=" "), "topo-corrected", sep=" "))
    }
  }

  # normalize
  normitems=mosaicitems[which(mosaicitems!=refitem)]
  if(method=="s3d"){
    print("normalizing")

    for (i in normitems){
      PIFlist=c(inlist[[i]],inlist[[refitem]])
      names(PIFlist)=names(inlist)[c(i, refitem)]
      models <- s3d(strips=PIFlist, cca=cca)
      # make sure that the items in the list have a name
      save(models, file=paste(paste("s3dModels", stacknames[i], sep="_"), "RData", sep=".")) 
      #save(models, file=paste("modelsItem", i, sep="_"))
      if(verbose == TRUE) {print(summary(models))}
      names(models$parameters)=names(PIFlist[[1]])
      inlist[[i]]<-aRn::CalibrateRaster(models$parameters, PIFlist)
      # if(verbose == TRUE) {plot(models)} # I stoped this to reduce processing time
      print(paste(paste("item", i, sep=" "), "normalized", sep=" "))
    }
  } else if(is.integer0(which(validmethods[2:4]==method))==FALSE){
    print("normalizing")
    for (i in normitems){
      normobj<-RStoolbox::pifMatch(inlist[[i]], inlist[[refitem]], method=method, 
                                   quantile=1-pval.pif, returnModels=TRUE)
      save(normobj, file=paste(paste(paste(method, "Models", sep=""), stacknames[i], sep="_"), "RData", sep=".")) 
      inlist[[i]]<-normobj[[1]]
      #if(verbose == TRUE) {plot(models)}
      print(paste(paste("item", i, sep=" "), "normalized", sep=" "))
    }
  }
  
  # mosaic TRY THE FUNCTION gdalUtils::mosaic_rasters IT MIGHT BE FASTER!!
  mosaicYr=inlist[[mosaicitems[1]]]
  if (length(mosaicitems)>1){
    print("mosaicking layers")
    mosaicYr=raster::merge(inlist[[mosaicitems[[1]]]], inlist[[mosaicitems[[2]]]])
    print("2 layers mosaicked")
  }
  
  if(length(mosaicitems)>2){
    for (i in 3:length(mosaicitems)){
      mosciacYr=raster::merge(mosaicYr, inlist[[mosaicitems[[i]]]])
      print(paste(i, "layers mosaicked", sep=" "))
    }
  }
  names(mosaicYr)=names(inlist[[refitem]])[normbands]
  return(mosaicYr)
### \code{}... 
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1:4))
    strips <- RasterIntersection(stack)
    ## thrs <- thresraster(strips[[2L]], strips[[1L]])
    ## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
    ## model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    ## calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})

stack2df <- structure(function #Extracts into data frame 
### This function extracts into a data frame the pixel values for
### all bands from different classes defined in a spatial dataframe
                       ##details<< This function ...
(
    inrast, ##<<\code{}. 
    invec, ##<<\code{}. 
    classcolname="class" ##<<\code{}.
){
  # extracts into a data frame the pixel values for all bands from different classes
  # defined in a spatial dataframe
  # inrast: the raster dataset containing pixel values to extract in [[1]]
  # invec: spatial dataframe object defining the locations where the data
  # should be extracted from 
  # classcolname: the column in the spatial dataframe containing the names 
  # of the attributes associated to those areas
  # value: a data frame with columns representing the pixel values in each band for
  # the areas labeled as defined by classcolname
  if (is.null(raster::intersect(extent(invec), extent(inrast)))){
    stop("the extents of inrast and invec do not overlap")
  }
  if(as.character(crs(inrast))!= as.character(crs(invec))){
    stop("inrast and invec should have the same projection")
  }
  # required function
  extractval=function(inraster=inrast, msk=msk){
    outvector=raster::mask(inraster, msk) 
    outvector=na.omit(raster::getValues(outvector))
    return(outvector)
  }
  
  # assign class ID to each class
  invec$class_ID=rep(NA, nrow(invec@data))
  for (i in 1:length(invec[[classcolname]])){
    invec$class_ID[which(invec[[classcolname]]==unique(invec[[classcolname]])[i])]=i
  }
  
  # mask the input raster including  pixels with valid values in all bands only
  inrast=stackmask(inrast)
  
  # create a raster of class_ids. TRY gdalUtils::gdal_rasterize. It might be faster!!!
  calibrast=raster::rasterize(invec, inrast[[2]], field=invec$class_ID)
  calibmsk<-maskfun(calibrast, 0, 1, NA)
  calibmsk=raster::mask(calibmsk, inrast[[2]])
  
  # Extract pixel values into a dataframe
  class_ID=(extractval(calibrast, calibmsk))
  dataset=data.frame(matrix(data=NA, nrow=length(class_ID), ncol=nlayers(inrast[[1]])))
  
  # add a column with a class name
  dbclassname=rep(NA, length(class_ID))
  for (i in 1:length(unique(invec[[classcolname]]))){
    dbclassname[which(class_ID==i, arr.ind=TRUE)] = unique(invec[[classcolname]])[i]
  }
  commonclasses= match(sort(unique(dbclassname)), sort(levels(invec[[classcolname]])))
  if(length(commonclasses)< length(levels(invec[[classcolname]]))){
    missing=sort(levels(invec[[classcolname]]))[-commonclasses]
    warning(paste(paste("the class", missing, sep= " "), 
                  "has no valid pixels in input raster", sep=" "))
    print(paste(paste("Warning: the class", missing, sep= " "), 
                "has no valid pixel values in the input raster", sep=" "))
  }
  
  dataset=cbind(class_ID, dbclassname, dataset)
  rm(class_ID, dbclassname)
  
  #I HAVE TO REMOVE DATA THAT HAS BEEN MASKED FROM THE RASTER
  # I MIGHT HAVE TO BRING BACK THE STACKMASK FUNCTION
  # ranges=c(-Inf,1,NA, 1,99999,1)
  # calibmsk=reclassify(calibrast,ranges)
  # calibmsk=calibmsk*inrastmskd[[2]]
  
  for (i in 1:nlayers(inrast[[1]])){
    dataset[,i+2]=extractval(inrast[[1]][[i]], calibmsk)
    print(paste(i, "layers extracted", sep=" "))
  }
  names(dataset)=c("class_ID", "class_name", names(inrast[[1]]))
  return(dataset)
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


thresraster2Updated <- structure(function #Threshold raster
### This function produces a threshold raster by computing
### standardized differences between two overlaying \code{Raster*}
### objects: reference and target, according to a range of valid pixel
### values.
                         ##details<<                         
(
    refstack, ##<<\code{RasterBrick}. Reference raster stack
    tarstack, ##<<\code{RasterBrick}. Target raster stack
    cca = FALSE, ##<<\code{logical}. Whether to perform a canonical
                 ##transformation to input rasters and then use it to
                 ##produce SSSD.
    propsamp=1 ##<<\code{numeric}. proportion of pixels to sample for
               ##fitting the cca.
){
  print("calculating standardized differences")
 
  # Apply canonical transformation to the data
  if(cca==TRUE){
    cct=rasterCCA(refstack, tarstack, propsamp=propsamp)
    dif=cct[[2]]-cct[[3]]} else {
    dif=refstack-tarstack}

  standardizediff=dif

  # Check if it is better to do this directly on the raster
  for (k in 1:nlyr(dif)) {
    difvec = terra::values(dif[[k]])
    difvec = difvec[!is.na(difvec)]
    meandif = mean(difvec, na.rm = TRUE)
    sddif = sd(difvec, na.rm = TRUE)
    standardizediff[[k]] = (meandif - dif[[k]])/sddif
    print(paste("standardized differences calculated for band", k, sep = " "))
  }
  sumstandardizediff = terra::app(standardizediff^2, sum, cores=4)#, na.rm = TRUE)
  return(sumstandardizediff)
  return(sumstandardizediff)
### \code{RasterLayer}. 
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    thrs <- thresraster2(strips[[2L]], strips[[1L]])
})

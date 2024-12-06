stackmaskUpdated <- structure(function #Stack mask
### This function ...
                       ##details<< This function ...
(
    inrast, ##<<\code{RasterBrick}. Reference raster.
    maskrast=TRUE ##<<\code{logical}. Apply the mask to the raster
                  ##stackoutput: If maskrast=true a list with the
                  ##stack masked and the mask produced. If
                  ##maskrat=FALSE, the raster mask.
){
  # THIS FUNCTION HAS TO BE OPTIMIZED. IT TAKES TOO LONG
  # masks out any pixels that have NAs in at least 1 band 
  # in a raster stack
  # maskrast: TRUE applies the mask to the raster stack
  # output: If maskrast=true a list with the stack masked and
  # the mask produced. If maskrat=FALSE, the raster masck)
  msk<-max(inrast)
  msk<-maskfunUpdated(msk, 0,1,NA)
  #inrast<-raster::mask(inrast, msk)
  if (maskrast==TRUE){
    inrast<-terra::mask(inrast, msk) # this is more than a minute faster than using the raster::mask function
    out=list(inrast,msk)} else {
      out=msk}
  return(out)
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
    model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})

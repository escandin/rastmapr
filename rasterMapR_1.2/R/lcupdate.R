lcupdate <- structure(function #Update target land cover map
### This function produces an update of a target land cover map
### (tarclass) by replacing the classification values of pixels
### located in areas considered as no change as denoted in a no change
### mask (nochmsk). See s3d for methods to produce no change masks).
                       ##details<< 
(
    refclass=NA, ##<<\code{}. Reference classification
    tarclass=NA, ##<<\code{}. Target classification
    nochmsk=NA ##<<\code{}. mask representing areas labeled as no
                ##change between refclass and tarclass.
){
      # lcupdate produces an update of a target land cover map
      # (tarclass) by replacing the classification values of pixels
      # located in areas considered as no change as denoted in a no
      # change mask (nochmsk). See s3d for methods to produce no
      # change masks).
  # refclass: reference classification
  # tarclass: target classification
  # nochmsk: #mask representing areas labeled as no change between refclass and tarclass
  chmsk=namaskplus(msk, 1, NA)
  #chmsk=namask(nochmsk) 
  refclass=raster::mask(refclass, nochmsk)
  tarclass=raster::mask(tarclass, chmsk)
  tarclass=raster::merge(tarclass,refclass)
  return(tarclass)
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

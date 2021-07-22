maskfun <- structure(function #Stack mask
### This function sets a threshold value (thresh). Any pixels above
### thresh are converted to aboveval. Any pixels bellow or equal to
### thresh are converted to belowval
                       ##details<< This function ...
(
    x, ##<<\code{}...
    thresh, ##<<\code{}...
    aboveval, ##<<\code{}...
    belowval ##<<\code{}...
){
  # sets a threshold value (thresh). Any pixels above 
  # thresh are converted to aboveval
  # Any pixels bellow or equal to thresh are converted to belowval
  require(raster)
  v <- raster::getValues(x)
  v[v>thresh]=aboveval
  v[v<=thresh]=belowval
  x <- raster::setValues(x, v)
  return(x)
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

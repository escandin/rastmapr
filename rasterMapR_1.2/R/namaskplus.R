namaskplus <- structure(function  # convert NAs 
### This function convertes NAs to to a value indicated by toval and
### any other values to a new value indicated by otherval

                        ##details<< This function ...
(
    x, ##<<\code{}. 
    toval, ##<<\code{}.
    ##\code{PIFmodel2}.
    otherval ##<<\code{}. 
){
  v <- raster::getValues(x)
  isna=is.na(v)
  isnotna=!is.na(v)
  v[isna]=toval
  v[isnotna]=otherval
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

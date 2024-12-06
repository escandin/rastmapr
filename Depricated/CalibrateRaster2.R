pifs=test$parameters
mlayer=tifstacks[[1]]
CalibrateRaster2 <- structure(function #Calibrate raster
### This function implements parameters from \code{\link{PIFmodel}} to
### calibrate multi- or hiper-spectral layers.
(
    pifs, ##<<\code{List}. Set of linear models such as these produced
          ##by \code{\link{PIFmodel}}.
    mlayer, ##<<\code{RasterLayer}. Mult-layer object such as that
            ##returned by \code{\link{EEstackWithoutMeta}}. If missing
            ##then the set used to compute the \code{pifs} is
            ##recycled.
    round. = 0 ##<<\code{numeric} or \code{NULL}. Integer indicating
             ##the number of decimal places. If \code{NULL} then this
             ##argument is ignored.
){
 coefs <- Map("coefficients", pifs)
  intercepts <- mapply(function(x) x[1L], coefs)
  names(pifs)=names(mlayer)
  slopes <- mapply(function(x) x[2L], coefs)
  if (missing(mlayer)) 
    mlayer <- as.list(attr(pifs, "env"))[["strips"]]
  #radiostack <- mlayer[[1L]][[names(pifs)]]
  radiostack <- mlayer[[names(pifs)]]
  b = 1
  normedstack = intercepts[b] + (slopes[b] * radiostack[[b]])
  print(paste(b, "band processed", sep = " "))
  for (b in 2:nlayers(radiostack)) {
    normedstack = stack(normedstack, intercepts[b] + (slopes[b] * 
                                                        radiostack[[b]]))
    print(paste(b, "bands processed", sep = " "))
  }
  if (!is.null(round.) | round. >= 0) 
    normedstack = round(normedstack, digits = round.)
  names(normedstack) <- names(radiostack)
  nms <- !names(mlayer[[1L]]) %in% names(pifs)
  if (any(nms)) {
    normedstack <- stack(normedstack, raster::subset(mlayer[[1L]], 
                                                     names(mlayer[[1L]])[nms]))
  }
  return(normedstack)
## RasterBrick.
} , ex=function(){
    ## \donttest{
    ##     tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
    ##                   'LT050060661988072201T1-SC20181031160603.tar.gz')
    ##     tarPaths <- system.file(tarFiles, package = 'aRn')
    ##     stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    ##     ## model <- PIFmodel(stack, pvalue = 3E-1, brackets = 7)
    ##     ## or
    ##     mlayer <- RasterIntersection(stack)
    ##     model <- PIFmodel(mlayer, pvalue = 3E-1, brackets = 7)
    ##     ## then
    ##     calib <- CalibrateRaster(model, stack)
    ## }
 
})

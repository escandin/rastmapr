PredictFromPIFmodel <- structure(function #DEPRECATED PIF-model prediction
### Use \code{\link{CalibrateRaster}} instead. The function predicts
### hyper-spectral data using PIF model
(
    stacks, ##<<\code{list}. Set of \code{RasterBrick} data such as that
            ##produced by \code{\link{EEstackWithoutMeta}}.
    pif.mod ##<<\code{list}. PIF models such as these
                        ##produced by
                        ##\code{\link{PIFmodel}}.
){
        if(is.null(stacks))
        return(NULL)
refs. <- lapply(stacks, raster::as.list)[[1L]]
stkLs <- Map(function(x,y)
    raster::calc(x, function(z)
        y$'coefficients'[1L] +
        y$'coefficients'[2L] * z),
    refs., pif.mod)
names(stkLs) <- names(stacks[[1L]])
predicts <- stack(stkLs)
### \code{RasterBrick}. Predicted features.
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    model <- PIFmodel(strips, pvalue = 3E-1, brackets = 7)
    plot(model)
    calib <- CalibrateRaster(model, strips)
})

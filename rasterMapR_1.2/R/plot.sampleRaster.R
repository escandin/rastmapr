plot.sampleRaster <- structure(function #sampleRaster plots
### This function plots \code{\link{sampleRaster}} objects.
(
    x, ##<<Raster* object.
    y = 1, ##<<If x is a RasterStack or RasterBrick: integer,
           ##character (layer name(s)), see \code{plot} method in the
           ##\code{raster} package.
    ... ##<< Additional arguments in \code{\link{plot}}. 
){
    raw. <- as.list(attr(x, 'env'))
    raster::plot(raw.[['rs']], y, ...)
    ## raster::plot(raw.[['pols']], border = 'gray80', add = TRUE)
    raster::text(raw.[['smp']], labels = raw.[['code']], col = 'gray20', cex = 0.7, adj = c(0,0), offset = 0.0)
    raster::plot(raw.[['smp']], border = 'black', lwd = 2, add = TRUE)

### \code{plot}. Raster with sample units.
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    set.seed(2)
    sr <- sampleRaster(stack[[1L]], side = 200, n = 2)
    ## plot(sr)
})

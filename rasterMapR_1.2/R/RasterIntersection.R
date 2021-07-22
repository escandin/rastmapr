RasterIntersection <- structure(function #Raster Intersection
### This function extracts intersections between two Earth-Explorer
### Landsat data sets.
                                ##details<< The function implements
                                ##\code{\link{intersect}} to derive a
                                ##list of two raster objects sharing
                                ##the same region.
(
    rasterImp ##<< Either \code{list} of two rasters, or
              ##\code{caracter} vector of paths to the data.
){
      if(is.null(rasterImp))
          return(NULL)

      nms <- Map('names', rasterImp)
      intr <- Reduce('intersect', nms)
      rasterImp <- lapply(rasterImp,function(x) x[[intr]])

      
    ## ri <- Reduce(as.list, rasterImp)
    ## ri <- Reduce(raster::as.list, rasterImp)
    inters <- Reduce(raster::intersect, rasterImp)
    newCrop <- Map(function(x)
        raster::crop(x, inters[[1L]]),
        rasterImp)
      m <- rclMatrix(0, FALSE)
    msk <- Map(function(x)
        raster::reclassify(x, m), newCrop)
    msk1 <- Reduce('*',msk)[[1L]]
    strips <- Map(function(x)
        raster::mask(x, msk1),
        newCrop)
      print(paste(nlayers(strips[[1L]]),
                  "intersections extracted", sep= " "))
      ## class(strips) <- c('RasterIntersection',
      ##                    class(strips))
    return(strips)
### \code{list} of rasters.
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    ## model <- PIFmodel(strips, pvalue = 3E-1, brackets = 7)
    ## plot(model)
    ## calib <- CalibrateRaster(model, strips)
})

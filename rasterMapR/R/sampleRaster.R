sampleRaster <- structure(function #Sample Raster
### This function processes Raster* objects to sample \code{n}
### square polygons. Sampling without replacement is
### implemented. Polygons having \link{NA} pixels are excluded.
(
    rl, ##<<\code{Raster*}. Raster Object with \code{UTM} Coordinate
        ##Reference System.
    side = sqrt(10) * 1E3, ##<<\code{numeric}. Side of the
                           ##square. Default produces squares of 10
                           ##km2.
    n = 30 ##<<\code{numeric}. Number of samples.
){
    oldw <- getOption("warn")
    options(warn = -1)
    pr <- projectRaster(rl[[1L]], crs = crs(rl),
                        res = c(side,side), method = 'ngb')
    r2pol <- rasterToPolygons(pr)
    sm <- sample(1:length(r2pol), n)
    smp <- r2pol[sm,]
    cropRaster <- function(rst, br){
        crp <- crop(rst, br)
        msk <- rasterize(br, crp, mask = TRUE)
        return(msk)}
    mrs <- Map(function(x)
        cropRaster(rl, smp[x,]), 1:length(smp))
    options(warn = oldw)
    names(mrs) <- sm
    mrs[['polygons']] <- smp
    ## raw <- list2env(list(rs = rl[[1L]],
    raw <- list2env(list(rs = rl,
                         pols = r2pol,
                         smp = smp,
                         code = sm))
    attributes(mrs) <- c(attributes(mrs),list(env=raw))
    class(mrs) <- c('sampleRaster',class(mrs))
    return(mrs)
### \code{list}. Set of \code{n} fixed-area rasters and the polygon geometry.
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    set.seed(2)
    sr <- sampleRaster(stack[[1L]], side = 200, n = 2)
})

nochg <- structure(function #Change and no change mask
### This function computes a binary change/no-change mask processing a
### threshold raster, see \code{\link{thresraster}}
                       ##details<< This function is implemented by
                         ##\code{\link{PIFmodel}}. The mask is
                         ##computed comparing values in the threshold
                         ##raster with the chi-square statistic, see
                         ##\code{\link{qchisq}}.
(
    thresraster, ##<<\code{RasterLayer}, such as that produced by
                 ##\code{\link{thresraster}}.
    degfree, ##<<\code{numeric}. Degrees of freedom.
    pvalue = 1E-4 ##<<\code{numeric}. Probability threshold
){
    print("calculating change/no-change mask")
    threshold=qchisq(pvalue,degfree)
    m <- c(-Inf,threshold,1,  threshold,999999,NA)
    m <- matrix(m, ncol=3, byrow=TRUE)
    nochgmsk = reclassify(thresraster, m)
    return(nochgmsk)
### \code{RasterLayer}. 
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    thrs <- thresraster(strips[[2L]], strips[[1L]])
    noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## plot(noch)
})

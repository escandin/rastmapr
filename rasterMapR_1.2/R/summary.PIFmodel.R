summary.PIFmodel <- structure(function #summarize PIFmodel object
### A set of summaries of a \code{\link{PIFmodel}} object is obtained.
(
    object, ##<< an object of class \code{"lm"}, usually, a result of
          ##a call to \code{\link{lm}}.
    correlation = FALSE, ##<< \code{logical}; if \code{TRUE}, the
                         ##correlation matrix of the estimated
                         ##parameters is returned and printed
    symbolic.cor = FALSE, ##<< \code{logical}. If \code{TRUE}, print
                          ##the correlations in a symbolic form (see
                          ##symnum) rather than as numbers.
    ... ##<< further arguments passed to or from other methods.
    
) {
    sm <- Map(function(object, correlation, symbolic.cor, ...)
        summary(object, correlation, symbolic.cor, ...), object, correlation, symbolic.cor,...)
    return(sm)
    ## \code{list}. Set of  summaries.
} , ex=function() {
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,4:6))
    strips <- RasterIntersection(stack)
    ## thrs <- thresraster(strips[[2L]], strips[[1L]])
    ## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
    model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    summary(model)
})

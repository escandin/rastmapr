PIFmodel <- structure(function #PIF modeling
### Modeling Pseudo-Invariant Features (\code{PIF}) using two
### hyper- or multi-spectral  rasters which overlay each other.
                      ##details<< Pseudo-invariant features
                          ##(\code{PIF}) are defined here as those
                          ##pixels sharing a common region of two
                          ##adjacent hyper-spectral data sets which
                          ##are similar according to a threshold in a
                          ##chi-square distribution, see
                          ##\code{\link{qchisq}}. The (\code{PIF}) are
                          ##computed implementing three in-package
                          ##functions: \code{\link{thresraster}},
                          ##\code{\link{nochg}}, and
                          ##\code{\link{calibrationParameters}},
                          ##establishing thresholds between rasters
                          ##(standardized differences), creating
                          ##change and no-change masks, and extracting
                          ##calibration parameters. The standardized
                          ##differences help to model the (\code{PIF})
                          ##into classes of reflectance. The
                          ##(\code{PIF}) classes are modeled fitting
                          ##linear models: \code{y ~ x}, where
                          ##\code{y} are the \code{PIF} classes in the
                          ##target hyper-spectral bands (first raster
                          ##in \code{strips}) and \code{x} are
                          ##corresponding \code{PIF} values in the
                          ##reference hyper-spectral band (second
                          ##raster in \code{strips}). The internal
                          ##method for modeling \code{PIF} classes
                          ##could produce blank intervals having empty
                          ##brackets (\code{NaN}), affecting the
                          ##number of observations in the regression
                          ##model. Finally, only bands having
                          ##alphanumeric names (e.g.,
                          ##\code{'sr_band1'}, ..., \code{sr_band10})
                          ##are included in the modeling procedure
                          ##while other bands such as these for
                          ##quality assessmend (e.g.,
                          ##\code{'pixel_qa', 'radsat_qa',
                          ##'radsat_qa','sr_aerosol'}) are excluded.
(
    strips, ##<<\code{list}. Set of two \code{Raster*} objects:
            ##target and reference, sharing a common region such as
            ##these produced by \code{\link{EEstackWithoutMeta}} or
            ##\code{\link{RasterIntersection}}.
    minvalid = 0, ##<<\code{numeric}. Minimum valid value in the
                  ##\code{PIF} classes.
    maxvalid = Inf, ##<<\code{numeric}. Maximum valid value in the
                  ##\code{PIF} classes.
    degfree = nlayers(strips[[2L]]) - 1, ##<<\code{numeric}. Degrees
                                         ##of freedom.  
    pvalue = 1E-4, ##<<\code{numeric}. Probability threshold affecting
                   ##the number of extracted \code{PIF}
    brackets = 30 ##<<\code{numeric}. Approximate number of
                  ##brackets. Blank intervals produce empty brackets
                  ##(\code{NaN}) affecting the number of observations
                  ##in the regression model.
){
    if(is.null(strips))
        return(NULL)
    if(!Reduce('==', Map(raster::extent, strips)))
        strips <- RasterIntersection(strips)

    tos. <- names(strips[[1L]])[!is.na(names(strips[[1L]]))]
    ## tos. <- names(strips[[1L]])[grep('[[:digit:]]', names(strips[[1L]]))]
    strips <- Map(function(x)
        raster::subset(x,tos.), strips)

    
    thrs <- thresraster(strips[[2L]], strips[[1L]],
                        minvalid, maxvalid)
    noch <- nochg(thrs, degfree, pvalue)
    calp <- calibrationParameters(strips[[2L]], strips[[1L]],
                                  noch, brackets)

    psimap <- Map(function(x)
        mask(x, noch), strips)
    psi <- Map(function(x)
        as.data.frame(getValues(x)),
        psimap) # <- long run; don't check it out!
    dats <- lapply(
        seq_len(ncol(psi[[1L]])),function(i)
            data.frame(y = psi[[2L]][,i],
                       x = psi[[1L]][,i]))
    nmst. <- names(strips[[1L]])
    names(dats) <- nmst.
    dats <- Map(function(x)
        na.omit(x), dats)

    raw <- list2env(list(raw = dats[names(calp)],
                         nm = names(strips),
                         strips = strips))
    attributes(calp) <- c(attributes(calp),list(env=raw))

    
    class(calp) <- c('PIFmodel',class(calp))
    return(calp)
### \code{list}. Set of \code{\link{lm}} models
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    ## thrs <- thresraster(strips[[2L]], strips[[1L]])
    ## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
    model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    ## plot(model)
})

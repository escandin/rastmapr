rasterCCA <- structure(function #Canonical correlation transformation
### This function applies a canonical correlation transformation to
### two raster stacks. THIS FUNCTION NEEDS TO BE MORE EFFICIENT. I
### have to figure out how to do CCA directly in the raster and ignore
### NAs. SELECTING A PROPSAMP DIFFERENT THAN 1 PRODUCES SUBSTANTIALLY
### DIFFERNT RESULTS #
                       ##details<< 
(
    xstack, ##<<\code{}. 
    ystack, ##<<\code{}. 
    propsamp=1
){
                                        # applies a canonical correlation transformation to two raster stacks
                                        # THIS FUNCTION NEEDS TO BE MORE EFFICIENT. I have to figure out how
                                        # to do CCA directly in the raster and ignore NAs.
                                        # SELECTING A PROPSAMP DIFFERENT THAN 1 PRODUCES SUBSTANTIALLY DIFFERNT RESULTS
    if(nlayers(xstack)!=nlayers(ystack)){
        stop("the xstack and ystack stacks should have the same number of bands")
    }
    msk=mask(stackmask(xstack, maskrast=FALSE),
             stackmask(ystack, maskrast=FALSE))
    xstack=mask(xstack, msk)
    ystack=mask(ystack, msk)
    
    xm=getValues(xstack)
    ym=getValues(ystack)
    if(is.numeric(propsamp) & propsamp>=0 & propsamp<1){
        xtrain=raster::sampleRandom(xstack, size = nrow(xstack)*ncol(xstack)*propsamp, na.rm = TRUE)
        ytrain=raster::sampleRandom(ystack, size = nrow(ystack)*ncol(ystack)*propsamp, na.rm = TRUE)
    } else if (propsamp==1){
        xtrain=na.omit(xm)
        ytrain=na.omit(ym)
    } else {
        stop("propsamp should be a number between 0 and 1")
    }
    CCA=cancor(xtrain,ytrain)
    xCCA=xm %*% CCA$xcoef
    yCCA=ym %*% CCA$xcoef
    xstack=setValues(xstack, xCCA)
    ystack=setValues(ystack, yCCA)
    return(list(CCA, xstack, ystack))
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
    ## model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    ## calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})

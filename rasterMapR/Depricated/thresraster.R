thresraster <- structure(function #Threshold raster
### This function produces a threshold raster by computing
### standardized differences between two overlaying \code{Raster*}
### objects: reference and target, according to a range of valid pixel
### values.
                         ##details<< This function is implemented by
                         ##\code{\link{PIFmodel}}. The standardized
                         ##differences in the threshold raster will be
                         ##contrasted with the chi-square statistic to
                         ##identify significant differences between
                         ##the overlaying rasters, see
                         ##\code{\link{nochg}}.
                         
(
    ref, ##<<\code{RasterBrick}. Reference raster
    target, ##<<\code{RasterBrick}. Target raster
    minvalid = 0, ##<<\code{numeric}. Minimum valid value.
    maxvalid = Inf ##<<\code{numeric}. Maximum valid value.
){
    ## ref <- strips[[2L]]
    ## target <- strips[[1L]]
    print("calculating standardized differences")
    ##drops not valid and negative values
    m <- c(-Inf,minvalid,NA, minvalid, maxvalid, 1, maxvalid,Inf,NA )
    m <- matrix(m, ncol=3, byrow=TRUE)
    msk1 <- reclassify(ref, m)
    msk2 <- reclassify(target, m)
    msk=msk1*msk2
    ## msk=msk[[1]]*msk[[2]]*msk[[3]]
    msk <- Reduce('*', raster::as.list(msk))
    rm(msk1, msk2)
    refmskd=mask(ref,msk)
    tarmskd=mask(target, msk)
    rm(ref,target)

    ## calculates standardized differences or the valid pixels
    dif=refmskd-tarmskd
    dif=mask(dif, msk)
    standardizediff=dif
    sampix=Which(msk == 1, cells=TRUE, na.rm=TRUE)
    for (k in 1:nlayers(dif)) {
        difvec=getValues(dif[[k]])
        difvec=difvec[!is.na(difvec)]
        meandif=mean(difvec, na.rm=TRUE)
        sddif=sd(difvec, na.rm=TRUE)
        standardizediff[[k]]= ((meandif-dif[[k]])/sddif)^2
        ## print(paste(k, "bands processed", sep=" "))

        ##the standardized difference is the value that will be
        ##contrasted with the chi-square statistic to identify
        ##significant differences\
    }
    print(paste(nlayers(dif), "bands processed", sep=" "))
    standardizediff=mask(standardizediff,msk)
    sumstandardizediff=mask(calc(standardizediff, sum, na.rm=TRUE),msk)
    return(sumstandardizediff)
### \code{RasterLayer}. 
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    thrs <- thresraster(strips[[2L]], strips[[1L]])
})

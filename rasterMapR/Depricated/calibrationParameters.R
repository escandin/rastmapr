calibrationParameters <- structure(function #Calibration parameters
### This function generates calibration parameters comparing bands in
### two overlaying \code{Raster*} objects: reference and target, while
### using a change/no-change mask, see \code{\link{nochg}}.
                       ##details<< This function is implemented by
                         ##\code{\link{PIFmodel}}. The function
                         ##develops four procedures. First, values in
                         ##the \code{Raster*} objects are extracted
                         ##using the binary mask. Second, the selected
                         ##values are reduced to classes according to
                         ##a number of brackets in
                         ##\code{nbrackets}. Third, average classes of
                         ##values are modeled generating intercepts
                         ##and slopes for each of the bands. Fourth,
                         ##non-significant models are subtracted.
(
    ref, ##<<\code{RasterBrick}. Reference raster.
    target, ##<<\code{RasterBrick}. Target raster.
    threshMask, ##<<\code{Raster*}. Change/no-change mask such as that
                ##produced by \code{\link{nochg}}.
    nbrackets = 30 ##<<\code{numeric}. Number of brackets, see Section
                   ##of \code{Details}.
){

    bracksplit=function(xval, yval, nbrackets){
    ## bracksplit=function(xval=seq(3,102), yval=xval+rnorm(length(xval), mean=0, sd=10), nbrackets){
        ## reduce the number of data points to the averages of different data ranges
        indata=cbind(xval, yval)
        bracketlength=(max(xval)-min(xval))/nbrackets
        bracketbins=seq(1,nbrackets-1)
        bracketranges=c(min(xval), bracketbins*bracketlength+min(xval), max(xval))
        
        labels=cbind(rep(NA, length(bracketbins)),rep(NA, length(bracketbins)))
        for(i in 1:length(bracketbins)) {
            datasubset=subset(indata, indata[,1]>=bracketranges[i] & indata[,1]<bracketranges[i+1])
            labels[i,1]=mean(datasubset[,1], na.rm=TRUE)
            labels[i,2]=mean(datasubset[,2], na.rm=TRUE)
        }
        return(labels)
    }
    
    ##builds mask. Not necessary to create externally anymore
    ranges =c(-Inf, 0, NA,0,99999999999,1)
    ranges= matrix(ranges, nrow=2, byrow=TRUE)
    mask= reclassify(target,ranges)
    mask<-mask[[1L]]
    
    weighted=1
    
    ##Masks the NA pixels
    threshMask[is.na(threshMask[])] <- 0
    
    ranges <-c(-Inf,0.0001,NA, 0.0001,Inf,1)
    plotranges<- matrix(ranges, ncol=3, byrow=TRUE)
    maskNA<- reclassify(mask, ranges)
    inmskd=mask(target,maskNA)
    refmskd=mask(ref,maskNA)
    rm(mask)
    ## This converts the mask with the pseudo invariant features into a matrix
    allPIF=raster::as.matrix(threshMask)
    allPIF1=allPIF
    allPIF1[which(allPIF==0)]=NaN
    
    if (nlayers(inmskd)==nlayers(refmskd)) {nb = (nlayers(inmskd))}
    else {print("different number of bands in  input files, aborting")}
    
    intercepts = rep(0, nb)
    slopes = rep(0, nb)
    
    
    ##Finds parameters for relative radiometric normalization
    
#############################################################
    insampnormw <- list()
    brackdata <- list()

    ## jpeg(paste(substr(inname, 1, 22) ,"scatterplot.jpg", sep="_"))
    nb <- nlayers(ref)

    ## par(mfrow = n2mfrow(nb))
    for (b in 1:nb){
        inb=raster::as.matrix(inmskd[[b]])
        inb=inb[which(allPIF1==1)]
        inb=as.vector(inb)
        refb=raster::as.matrix(refmskd[[b]])
        refb=refb[which(allPIF1==1)]
        refb=as.vector(refb)
        brackdata[[b]]=bracksplit(xval=inb, yval=refb, nbrackets)

        ## brackdata=bracksplit(xval=inb, yval=refb, nbrackets=30)

        ## insampnormw=lm(brackdata[,2]~brackdata[,1],
        ##                na.action=na.omit)

        
        ## intercepts[b]=as.numeric(insampnormw$coefficients[1])
        ## slopes[b]=as.numeric(insampnormw$coefficients[2])

        
                                        #save scatterplot
        ## plot(inb, refb, col='dimgray')
        ## lines(inb, intercepts[b] + slopes[b]*inb)
                                        #rm(insampnormw)
    }

    brackdata <- lapply(brackdata,function(x)
        data.frame(y = x[,2L], x = x[,1L]))
    names(brackdata) <- names(ref)
    print('fitting models')
    funcReg <- function(z){
        lm(y ~ x, data = z)}
    lms <- Map(function(z)
        tryCatch(
            funcReg(z),
            error = function(e) NULL),
        brackdata)
    nmst. <- names(target)
    names(lms) <- nmst.
    lms <- Filter(Negate(is.null), lms)

    nm <- lapply(lms, function(x)
        any(is.na(x$'coefficients')))
    nm. <- names(nm)[unlist(nm)]
    
    nnm <- intersect(nmst., names(lms))
    cnm <- setdiff(nmst., nnm)
    cnm <- c(cnm, nm.)

    lms <- lms[!names(lms)%in%cnm]
    
    if(length(cnm) >= 1){
        cnm. <- paste(cnm, collapse = ', ')
        rms <- paste(
            'excluding non-significant models: (',
            cnm., ')', sep = "")
        print(rms)
    }

    
    return(lms)


    
    ## return(brackdata)
    ## dev.off()
    param=data.frame(cbind(intercepts, slopes))
    ##save parameters
    ## write.table(param, (paste(substr(radioname, 1, 22),"parameters.csv", sep="")),col.names=T,row.names=T, sep= ",")
    
    return(param)
### \code{RasterLayer}. 
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1:4))
    strips <- RasterIntersection(stack)
    ## thrs <- thresraster(strips[[2L]], strips[[1L]])
    ## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
    model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})

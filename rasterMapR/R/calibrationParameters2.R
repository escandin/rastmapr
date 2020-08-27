calibrationParameters2 <- structure(function #Calibration parameters
### This function ...
                       ##details<< This function ...
(
    ref, ##<<\code{RasterBrick}. Reference raster.
    tar, ##<<\code{RasterBrick}. Target raster.
    threshMask ##<<\code{Raster*}. Change/no-change mask such as that
                ##produced by \code{\link{nochg}}.
){
      print('extracting reflectance values from no change pixels')
  refmsk=stackmask(ref)[[2]]
  tarmsk=stackmask(tar)[[2]]
  msk=refmsk*tarmsk*threshMask
  refmskd=ref*msk
  tarmskd=tar*msk
  rm(ref, tar, msk)
  gc()

  if (nlayers(tarmskd) == nlayers(refmskd)) {
    nb = (nlayers(tarmskd))
  } else {
    print("different number of bands in  input files, aborting")
  }
  intercepts = rep(0, nb)
  slopes = rep(0, nb)
  insampnormw <- list()
  indata<-list()
  nb <- nlayers(refmskd)
  
  extractvalid <- function(x){
    v=getValues(x)
    v=v[!is.na(v)]
    return(v)
  }
  
  inrefmat=lapply(raster::as.list(refmskd), function(x)
    extractvalid(x))
  intarmat=lapply(raster::as.list(tarmskd), function(x)
    extractvalid(x))
  
  indata <- Map(function(x,y)
    data.frame(y, x), y=inrefmat, x=intarmat)
  
  print("fitting models")
  funcReg <- function(z) 
    lm(y ~ x, data = z)#}
  
  # This function demands a large memory allocation. It might have to be
  # changed to process each band separately
  lms <- Map(function(z) tryCatch(funcReg(z), error = function(e) NULL), 
             indata)
  
  nmst. <- names(tarmskd)
  names(lms) <- nmst.
  lms <- Filter(Negate(is.null), lms)
  nm <- lapply(lms, function(x) any(is.na(x$coefficients)))
  nm. <- names(nm)[unlist(nm)]
  nnm <- intersect(nmst., names(lms))
  cnm <- setdiff(nmst., nnm)
  cnm <- c(cnm, nm.)
  lms <- lms[!names(lms) %in% cnm]
  if (length(cnm) >= 1) {
    cnm. <- paste(cnm, collapse = ", ")
    rms <- paste("excluding non-significant models: (", 
                 cnm., ")", sep = "")
    print(rms)
  }
  out=list(indata,lms)
  names(out)=c("data", "parameters")
  return(out)

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
    model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})

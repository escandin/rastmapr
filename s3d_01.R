EEstackWithoutMeta2 <- function (tarPath, bands = NULL, sat.nm = "LO08") {
  fimp <- function(tarPath, bands) {
    if (is.null(tarPath)) 
      return(NULL)
    tarNm <- basename(tarPath)
    tmp. <- tempdir()
    temp. <- file.path(tmp., "EEstack")
    lstar <- untar(tarPath, exdir = temp., list = TRUE)
    if (!all(lstar %in% list.files(temp.))) 
      untar(tarPath, exdir = temp.)
    cmmPth <- Reduce(PTXQC::LCS, lstar)
    tif <- lstar[grepl(".tif", lstar, ignore.case = TRUE)]
    ftiford <- function(tif) {
      tifa <- tif
      cmmPth <- Reduce(PTXQC::LCS, tif)
      tif <- gsub(cmmPth, "", tif)
      names(tif) <- 1:length(tif)
      fOrd <- function(x) x[order(nchar(x), x)]
      tif. <- fOrd(tif[grepl("[0-9]", tif)])
      tif.. <- fOrd(tif[!tif %in% tif.])
      tif.. <- as.numeric(names(c(tif., tif..)))
      return(tifa[tif..])
    }
    tif <- ftiford(tif)
    pathtifs <- file.path(temp., tif)
    rasterImp <- Map(function(x) raster::raster(x), pathtifs)
    exts <- Map(extent, rasterImp)
    isext <- sapply(exts, function(x) extent(rasterImp[[1L]]) == 
                      x)
    rasterImp <- rasterImp[isext]
    tifNm. <- gsub(cmmPth, "", tif)
    tifNm <- gsub("^_|.tif", "", tifNm., ignore.case = TRUE)
    names(rasterImp) <- tifNm[isext]
    stkk <- stack(rasterImp)
    attributes(stkk) <- c(attributes(stkk), list(nam. = cmmPth))
    return(stkk)
  }
  mta <- Map(function(x, y) fimp(x, y), tarPath, MoreArgs = list(y = bands))
  mtt <- lapply(mta, function(x) attr(x, "nam."))
  names(mta) <- unlist(mtt)
  nm. <- names(mta)
  if (!is.null(sat.nm)) 
    mta <- equateEEnames(mta, sat.nm)
  if (!is.null(bands)) 
    mta <- Map(function(x) x[[bands]], mta)
  class(mta) <- c("EEstack", class(mta))
  return(mta)
}
thresraster2<-function (refstack, tarstack, cca=FALSE, propsamp=1){
  # refstack: reference raster stack
  # tarstack: target raster stack
  # cca: whether to perform a canonical transformation to input rasters and then use it to produce SSSD
  # ccasamp: when cca=TRUE, the proportion of pixels to sample for cca transformation
  # propsamp: proportion of pixels to sample for fitting the cca.
  print("calculating standardized differences")
 
  # Apply canonical transformation to the data
  if(cca==TRUE){
    cct=rasterCCA(refstack, tarstack, propsamp=propsamp)
    dif=cct[[2]]-cct[[3]]} else {
    dif=refstack-tarstack}

  standardizediff=dif
  # Check if it is better to do this directly on the raster
  for (k in 1:nlayers(dif)) {
    difvec = raster::getValues(dif[[k]])
    difvec = difvec[!is.na(difvec)]
    meandif = mean(difvec, na.rm = TRUE)
    sddif = sd(difvec, na.rm = TRUE)
    standardizediff[[k]] = (meandif - dif[[k]])/sddif
    print(paste("standardized differences calculated for band", k, sep = " "))
  }
  
  sumstandardizediff = calc(standardizediff^2, sum)#, na.rm = TRUE)
  return(sumstandardizediff)
}
nochg2<-function (thresraster, pvalue = 1e-02, distype="gamma",
                  fitmethod="mle", propsamp=0.01, degfree=NULL){#degfree, filename="gammadist") {
  # thresraster: a file containing the sum of the square of
  #       the standardized differences (SSSD) of changes between target and reference images.
  # pvalue: pvalue set to define the threshold for selecting the PIFs
  # distype: the type of satistical distribution assumed for the sum of the square of
  #       the standardized differences (SSSD) of changes between target and reference images.
  #       It is very likely gamma but also accepts "chisquare". I SHOULD GENERALIZE IT TO 
  #       ACCEPT ALSO negative exponential and weibull
  # fitmethod: the method to be used to fit the statistical distribution. "moments" 
  #       might also work but it has not been tested
  # propsamp: proportion of pixels to sample for fitting the statistical distribution.
  # degfree:  only applies for chisquare. It sould be equal to the number of layers in the 
  #      input raster stacks
  # distype: "gamma", "chisq", ("weibull", "nomal" "negative exponential") 

  print("calculating change/no-change mask")
  if (!is.numeric(propsamp) & propsamp<0 & propsamp>1){
    stop("propsamp should be a number between 0 and 1")
  }
  samp=raster::sampleRandom(thresraster, length(thresraster)*propsamp)
  if(distype=="chisq"){
    threshold = qchisq(pvalue, degfree)
    simulated=rchisq(length(samp), degfree)
    ks<-ks.test(samp, simulated)
    fit.stats=list(samp, simulated, ks)} else {
      require(fitdistrplus)
    #if(cca==TRUE) {thresraster=thresraster[[1]]}
      fit.stats <- fitdistrplus::fitdist(samp, dist =distype)
      threshold= qgamma(pvalue, shape=fit.stats$estimate[1],
                        rate=fit.stats$estimate[2])
  }
  if(threshold < minValue(thresraster)){
    warning("the pvalue selected produces a threshold that is larger that the minimum 
        value in thresraster producing an NA. Select a higher pvalue to produce a valid threshold mask")
  nochgmsk=NA} else {nochgmsk=maskfun(thresraster, threshold, NA, 1)}
  nochgmsk=list(nochgmsk, fit.stats)
  return(nochgmsk)
}
calibrationParameters2<-function (ref, tar, threshMask){#, 
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
}
s3d <- function(strips, cca=FALSE, distype="gamma", distsamp=0.01, pval.pif = 1e-02, pval.chg=0.99,
                minPIF=5, thres=1e-2, maxiter=20, prefix="", 
                norm.ext=NULL,fitline=TRUE, writemasks=TRUE){ #degfree = nlayers(strips[[2L]]),
  # s3d: Iterated Sum of the Squared Standaradized Differences iterates the selection of
  # PIFs defined by a pval.tif threshold by excluding in each iteration the pixels 
  # with the largest change based on the level of probability specified by the pval.chg parameter.
  # strips:  strips: a list containing the target and reference raster stacks respectively
  # cca: same argument as for PIFmodel2
  # strips: a list containing the target and reference raster stacks respectively
  # distype: same as noch2
  # distsamp: same as noch2
  # cca: same as thresraster2
  # ccasamp: same as propsamp in thresraster2
  # pvalue: same as noch2
  # degfree: same as noch2
  # minPIF: minimum number of pixels that should be selected to produce the
  #   normalization. The function will abort if the number of pixels selected is lower
  # fitline: if true, it will label pseudoinvariant targets and use them to fit a regression
  #   line in each iteration to be used later for radiometric normalization
  # writemasks: If true, it writes in the hard drive, 
      # the no change masks produced in each iteration
  
  # Crop rasters
  if (is.null(strips)) 
    return(NULL)
  if (!Reduce("==", Map(raster::extent, strips))) 
    if (is.null(norm.ext)){
      if (as.character(raster::extent(strips[[1]]))!=as.character(raster::extent(strips[[2]]))){  
        print("cropping rasters to a common extent")
        e=commonExtent(strips)
        strips=Map(function(x)
          crop(x, e), raster::as.list(strips))
      }
    } else {
      strips=Map(function(x)
        crop(x, norm.ext), raster::as.list(strips))
    }
 
  # obtain a first no change mask
  print("calculating initial parameters")
  thrs <- thresraster2(strips[[2L]], strips[[1L]], cca=cca)#, degfree = nlayers(strips[[2L]])-1,
  if (distype=="chisq"){
    noch <- nochg2(thrs, pvalue= pval.pif, distype=distype, 
                   propsamp=distsamp, degfree=nlayers(strips[[1]]))
  } else {
    noch <- nochg2(thrs, pvalue= pval.pif, distype=distype, propsamp=distsamp)
  }
  # iterate until the intercept values are stabilized
  i=1
  if(fitline==TRUE){
    lmparam=data.frame(matrix(nrow=nlayers(strips[[1]]),ncol=4))
    names(lmparam)=c("iter","band", "intercept", "slope")
    if (class(noch[[1]])== "RasterLayer"){
      nPIFs=cellStats(noch[[1]], stat="sum", na.rm=TRUE)
      print(paste("the number of PIFs selected for normalization is equal to", 
                  nPIFs, sep=" "))
      if(nPIFs<minPIF){
        warning("the number of PIFs to produce the regression is lower than 
         the number specified in minPIF, producing NAs for regression parameters. 
         Either decrease minPIF or  increase the pvalue used to select PIFs")
        lmparam$iter <- rep(i,nlayers(strips[[1]]))
        lmparam$band <- seq(1:nlayers(strips[[1]]))
        lmparam$intercept <- rep(NA, nlayers(strips[[1]])) # extract coefficients
        lmparam$slope <- rep(NA, nlayers(strips[[1]])) # extract slope
      } else {calp <- calibrationParameters2(strips[[2L]], strips[[1L]], noch[[1]])
      lmparam$iter <- rep(i,nlayers(strips[[1]]))
      lmparam$band <- seq(1:nlayers(strips[[1]]))
      lmparam$intercept <- unlist(lapply(lapply(calp$parameters, '[[', 1), '[[', 1)) # extract coefficients
      lmparam$slope <- unlist(lapply(lapply(calp$parameters, '[[', 1), '[[', 2)) # extract slope
      }} else {
      warning("thresMask is not a valid raster layer mask. NAs produced for regression parameters")
      lmparam$iter <- rep(i,nlayers(strips[[1]]))
      lmparam$band <- seq(1:nlayers(strips[[1]]))
      lmparam$intercept <- rep(NA, nlayers(strips[[1]])) # extract coefficients
      lmparam$slope <- rep(NA, nlayers(strips[[1]])) # extract slope
    }
    lmparamtot=lmparam
  }
  
  if(distype=="chisq"){
    distparam=data.frame(cbind(i, noch[[2]][[3]][[1]], noch[[2]][[3]][[2]]))
    names(distparam)=c('iter', 'ksD', "pval")
    } else {
      distparam=data.frame(matrix(nrow=1,ncol=4))
      names(distparam)=c("iter", "ksD", "shape", "rate")
      diststat=gofstat(noch[[2]])
      distparam$iter <- i
      distparam$ksD=diststat$ks
      distparam$shape=noch[[2]]$estimate[1]
      distparam$rate=noch[[2]]$estimate[2]
  }
  distparamtot=distparam
 
  print('iterating PIF extraction until the change in the distribution parameters are below thres values')
  repeat{
    i <- i+1
    if (distype=="chisq"){
      noch <- nochg2(thrs, pvalue= pval.chg, distype=distype, 
                     propsamp=distsamp, degfree=nlayers(strips[[1]]))
    } else { noch <- nochg2(thrs, pvalue= pval.chg, distype=distype, propsamp=distsamp)}
 
    if (writemasks==TRUE){
      raster::writeRaster(noch[[1]], overwrite=TRUE,
                  paste(prefix, paste(paste("removeChgMsk", i-1, sep=""), "tif", sep="."), sep="_"))
    }
      #plot(noch)
    stripsmskd=list(strips[[1]]*noch[[1]], strips[[2]]*noch[[1]])
    # release memory
    gc()
    
    thrs <- thresraster2(stripsmskd[[2L]], stripsmskd[[1L]], cca=cca)#, degfree = nlayers(strips[[2L]])-1,
    if (distype=="chisq"){
      noch <- nochg2(thrs, pvalue= pval.pif, distype=distype, 
                     propsamp=distsamp, degfree=nlayers(strips[[1]]))
    } else {
      noch <- nochg2(thrs, pvalue= pval.pif, distype=distype, propsamp=distsamp)
    }
    
    if(fitline==TRUE){
      if (class(noch[[1]])== "RasterLayer"){      #is.na(noch[[1]])==FALSE){
        nPIFs=cellStats(noch[[1]], stat="sum", na.rm=TRUE)
        print(paste("the number of PIFs selected for normalization is equal to", 
                    nPIFs, sep=" "))
        if(nPIFs<minPIF){
          warning("the number of PIFs to produce the regression is lower than 
          the number specified in minPIF, producing NAs for regression parameters. 
          Either decrease minPIF or  increase the pvalue used to select PIFs")
          lmparam$iter <- rep(i,nlayers(strips[[1]]))
          lmparam$band <- seq(1:nlayers(strips[[1]]))
          lmparam$intercept <- rep(NA, nlayers(strips[[1]])) # extract coefficients
          lmparam$slope <- rep(NA, nlayers(strips[[1]])) # extract slope
        } else {calp <- calibrationParameters2(strips[[2L]], strips[[1L]], noch[[1]])
          lmparam$iter <- rep(i,nlayers(strips[[1]]))
          lmparam$band <- seq(1:nlayers(strips[[1]]))
          lmparam$intercept <- unlist(lapply(lapply(calp$parameters, '[[', 1), '[[', 1)) # extract coefficients
          lmparam$slope <- unlist(lapply(lapply(calp$parameters, '[[', 1), '[[', 2)) # extract slope
          }
      } else {
        warning("thresMask is not a valid raster layer mask. NAs produced for regression parameters")
        lmparam$iter <- rep(i,nlayers(strips[[1]]))
        lmparam$band <- seq(1:nlayers(strips[[1]]))
        lmparam$intercept <- rep(NA, nlayers(strips[[1]])) # extract coefficients
        lmparam$slope <- rep(NA, nlayers(strips[[1]])) # extract slope
      }
      lmparamtot=rbind(lmparamtot, lmparam)
    }
    if(distype=="chisq"){
      distparam=data.frame(i, cbind(noch[[2]][[3]][[1]], noch[[2]][[3]][[2]]))
      names(distparam)=c('iter', 'ksD', "pval")
    } else {
      diststat=gofstat(noch[[2]])
      distparam$iter <- i
      distparam$ksD=diststat$ks
      distparam$shape=noch[[2]]$estimate[1]
      distparam$rate=noch[[2]]$estimate[2]
    }
    distparamtot=rbind(distparamtot, distparam)
    print(paste(i, "iterations processed", sep=" "))
    
    # if (distype=="chisqs"){
    # chg.intp=abs(lmparamtot$intercept[which(lmparamtot$iter==i)]-
    #                lmparamtot$intercept[which(lmparamtot$iter==(i-1))])
    # chg.slope=abs(lmparamtot$slope[which(lmparamtot$iter==i)]-
    #                 lmparamtot$slope[which(lmparamtot$iter==(i-1))])
    # if (max(chg.intp)<thres.intp & max(chg.slope)<thres.slope) {
    #   break}
    # }
    # else {
    #   chg.shape=abs(distparamtot$shape[which(distparamtot$iter==i)]-
    #                 distparamtot$shape[which(distparamtot$iter==i-1)])
    #   chg.rate=abs(distparamtot$rate[which(distparamtot$iter==i)]-
    #                 distparamtot$rate[which(distparamtot$iter==i-1)])
    # if (max(chg.shape)<thres.shape & max(chg.rate)<thres.rate) 
    
    # Setup threshold for convergence to stop the iterations
    ksDstd=(max(distparamtot$ksD)-distparamtot$ksD)/max(distparamtot$ksD)
    delta=ksDstd[2:length(ksDstd)]-ksDstd[1:length(ksDstd)-1]
    if (delta[length(delta)]<thres){
      break
    }
    if (i==maxiter){
      break
    }
  }
  par(mar = c(4, 4, 1, 1) + 0.1)
  par(mfrow=c(round(sqrt(nlayers(strips[[1]])))+1,round(sqrt(nlayers(strips[[1]])))))
  for(b in 1:nlayers(strips[[1]])){
    plot(lmparamtot$iter[which(lmparamtot$band==b)], 
         lmparamtot$intercept[which(lmparamtot$band==b)], xlab="band", ylab="intercept")
  }
  par(mfrow=c(round(sqrt(nlayers(strips[[1]])))+1,round(sqrt(nlayers(strips[[1]])))))
  for(b in 1:nlayers(strips[[1]])){
    plot(lmparamtot$iter[which(lmparamtot$band==b)], 
         lmparamtot$slope[which(lmparamtot$band==b)], xlab="band", ylab="slope")
  }
  
  par(mfrow=c(round(sqrt(nlayers(strips[[1]])))+1,round(sqrt(nlayers(strips[[1]])))))
  for (b in 1:nlayers(strips[[1]])){
    plot(calp$data[[b]][,1]~calp$data[[b]][,2], xlab="target", ylab="reference")
    abline(unlist(calp$parameters[[b]][1])[1], unlist(calp$parameters[[b]][1])[2])
  }
  if (distype=="chisq"){
    par(mfrow=c(1,1))
    plot(distparamtot$ksD~distparamtot$iter, xlab="iter", ylab="ksD")
  } else {
    par(mfrow=c(2,2))
    plot(distparamtot$ksD~distparamtot$iter, xlab="iter", ylab="ksD")
    plot(distparamtot$shape~distparamtot$iter, xlab="iter", ylab="shape")
    plot(distparamtot$rate~distparamtot$iter, xlab="iter", ylab="rate")
    plot(distparamtot$rate~distparamtot$shape, xlab="shape", ylab="rate")
  } #else {paramstats=lmparamtot}
  if(fitline==TRUE){
    paramstats=list(lmparamtot, distparamtot)
    names(paramstats)=c('lmparam', "distparam")} else{
      paramstats=distparamtot}
  
  out=list(thrs[[1]], noch, calp[[1]], calp[[2]], paramstats)
  names(out)=c("sumstandardizediff", "noch", "data", "parameters", "paramstats")
  return(out)
}
CalibrateRaster2 <- function (pifs, mlayer, round. = 0) {
  coefs <- Map("coefficients", pifs)
  intercepts <- mapply(function(x) x[1L], coefs)
  slopes <- mapply(function(x) x[2L], coefs)
  if (missing(mlayer)) 
    mlayer <- as.list(attr(pifs, "env"))[["strips"]]
  radiostack <- mlayer[[1L]][[names(pifs)]]
  b = 1
  normedstack = intercepts[b] + (slopes[b] * radiostack[[b]])
  print(paste(b, "band processed", sep = " "))
  for (b in 2:nlayers(radiostack)) {
    normedstack = stack(normedstack, intercepts[b] + (slopes[b] * 
                                                        radiostack[[b]]))
    print(paste(b, "bands processed", sep = " "))
  }
  if (!is.null(round.) | round. >= 0) 
    normedstack = round(normedstack, digits = round.)
  names(normedstack) <- names(radiostack)
  nms <- !names(mlayer[[1L]]) %in% names(pifs)
  if (any(nms)) {
    normedstack <- stack(normedstack, raster::subset(mlayer[[1L]], 
                                                     names(mlayer[[1L]])[nms]))
  }
  return(normedstack)
}


# PIFmodel2<-function (strips, distype= "gamma", distsamp=0.01, cca=FALSE,#minvalid = -Inf, maxvalid = Inf, 
#                      ccasamp=1, pvalue = 2e-02, degfree=NULL, minPIF=2){# degfree = nlayers(strips[[2L]])\
#   # strips: a list containing the target and reference raster stacks respectively
#   # distype: same as noch2
#   # distsamp: same as noch2
#   # cca: same as thresraster2
#   # ccasamp: same as propsamp in thresraster2
#   # pvalue: same as noch2
#   # degfree: same as noch2
#   # minPIF: minimum number of pixels that should be selected to produce the
#   # normalization. The function will abort if the number of pixels selected is lower
#   
#   thrs <- thresraster2(strips[[2L]], strips[[1L]], cca=cca, propsamp=ccasamp)#, quant.rm=quant.rm)
#   
#   if(cca==TRUE){thrs=thrs[[1]]}
#   noch <- nochg2(thrs, pvalue, distype=distype, propsamp=distsamp)
#   nPIFs=cellStats(noch[[1]], stat="sum", na.rm=TRUE)
#   print(paste("the number of PIFs selected for normalization is equal to", 
#               nPIFs, sep=" "))
#   if(nPIFs<minPIF){
#     stop("the number of PIFs to produce the regression is lower than 
#          the number specified in minPIF. Either decrease minPIF or 
#          increase the pvalue used to select PIFs")
#   }
#   calp <- calibrationParameters2(strips[[2L]], strips[[1L]], noch[[1]])#, brackets, weightVal=NULL)
#  
#   out=list(thrs, noch, calp[[1]], calp[[2]])
#   names(out)=c("sumstandardizediff", "noch", "data", "parameters")
#   return(out)
# }
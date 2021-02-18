s3d <- structure(function #Iterated Sum of the Squared Standaradized
                          #Differences
### This function iterates a selection of \code{PIFs} defined by a
### \code{pval.tif} threshold by excluding in each iteration the
### pixels with the largest change based on the level of probability
### specified by a \code{pval.chg} parameter strips
                       ##details<< This function ...
(
    strips, ##<<\code{list}. Target and reference raster stacks
            ##respectively.
    cca = FALSE, ##<<\code{logical}. Same argument as for
                 ##\code{PIFmodel2}.
    distype = "gamma", ##<<\code{character}. Same as noch2
    distsamp = 0.01, ##<<\code{numeric}. Same as noch2
    pval.pif = 1e-02, ##<<\code{numeric}.
    pval.chg=0.99, ##<<\code{numeric}.
    minPIF = 5, ##<<\code{numeric}. Minimum number of pixels that
                ##should be selected to produce the normalization. The
                ##function will abort if the number of pixels selected
                ##is lower.
    thres=1e-2, ##<<\code{numeric}.
    maxiter=20, ##<<\code{numeric}.
    prefix="", ##<<\code{character}.
    norm.ext=NULL, ##<<\code{} or \code{NULL}.
    fitline=TRUE, ##<<\code{logical}. Label pseudoinvariant targets
                  ##and use them to fit a regression.
    writemasks=TRUE ##<<\code{logical}. Write in the hard drive the no
                    ##change masks produced in each iteration
){
 #degfree = nlayers(strips[[2L]]),
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
    distparam=data.frame(cbind(i, noch[[2]][[1]], noch[[2]][[2]]))
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
      distparam=data.frame(i, cbind(noch[[2]][[1]], noch[[2]][[2]]))
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
### \code{}... 
} , ex=function(){
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'rasterMapR')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
sssd <- s3d(strips, distsamp = 0.5,pval.pif = 0.4, pval.chg = 0.5, maxiter = 1)
})

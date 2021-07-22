AccuStats <- structure(function #Extract accuracy results 
### This function extracts accuracy results from confusion
### matrices. Add this to accuplot!
                       ##details<< This function ...
(
    prefix="class2013_AP_raw_iter", ##<<\code{}. 
    suffix="accumatrix.RData", ##<<\code{}. 
    iter=30 ##<<\code{}.
){
  for (i in 1:iter){
    outname=paste(paste(prefix, i, sep="_"), suffix, sep="_") 
    load(outname)
    accuresults=data.frame(conmatrix$byClass, stringsAsFactors=FALSE)
    users=accuresults$Precision
    producers=accuresults$Recall
    overall=conmatrix$overall[[1]]
    classnames=rownames(accuresults)
    classnames=str_remove(classnames, "Class: ")
    if (i==1){
      allusers=users
      allproducers=producers
      alloverall=overall} else {
        allusers=data.frame(cbind(allusers, users), stringsAsFactors=FALSE)
        allproducers=data.frame(cbind(allproducers, producers), stringsAsFactors=FALSE)
        alloverall=rbind(alloverall, overall)
      }
  }
  allusers=data.frame(cbind(classnames, allusers), stringsAsFactors=FALSE)
  allproducers=data.frame(cbind(classnames, allproducers), stringsAsFactors=FALSE)
  alloverall=alloverall
  mean(alloverall)
  sd(alloverall, na.rm=TRUE)
  allusers=data.frame(allusers, stringsAsFactors=FALSE)
  allproducers=data.frame(allproducers,stringsAsFactors=FALSE)
  
  # summarize user's and producer's mean and sd accuracy in a table
  finalstats=data.frame(cbind(round(apply(allusers[,2:ncol(allusers)], 1,  mean), digits=2), 
                              round(apply(allusers[,2:ncol(allusers)], 1,  sd), digits=2), 
                              round(apply(allproducers[,2:ncol(allproducers)], 1,  mean), digits=2),
                              round(apply(allproducers[,2:ncol(allproducers)], 1,  sd), digits=2)), stringsAsFactors=FALSE)
  # I have to do this twice otherwise data.frame converts numbers into factors.
  finalstats=cbind(classnames,finalstats)
  names(finalstats)=c("classnames", "users mean", "users sd", "producers mean", "producers sd")
  
  return(list(allusers, allproducers, alloverall, finalstats))
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

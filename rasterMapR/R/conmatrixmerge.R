conmatrixmerge <- structure(function #Merges pairs of classes in an
                                     #accuracy matrix
### This function merges pairs of classes in an accuracy matrix and
### recalculates accuracy metrics.
                       ##details<< This function ...
(
    reclassx, ##<<\code{}. The input accuracy matrix as a table. 
    rcl ##<<\code{}. A nx2 reclassification matrix with a number
        ##representing the target class object in the first column and
        ##the class to be merged to in the second column.
){
  # This function merges pairs of classes in an accuracy matrix and
  # recalculates accuracy metrics.  reclassx= the input accuracy
  # matrix as a table rcl= a nx2 reclassification matrix with a number
  # representing the target class object in the first column and the
  # class to be merged to in the second column
  for (i in 1:nrow(rcl)){
    reclassx[rcl[i,2],]=reclassx[rcl[i,2],]+reclassx[rcl[i,1],]
    reclassx[,rcl[i,2]]=reclassx[,rcl[i,2]]+reclassx[,rcl[i,1]]
  }
  reclassx=reclassx[,-rcl[,1]]
  reclassx=reclassx[-rcl[,1],]
  overall=sum(diag(reclassx))/sum(reclassx)
  producers=diag(reclassx)/apply(reclassx, 2,sum)
  users=diag(reclassx)/apply(reclassx, 1,sum)
  accustats=cbind(users,producers)
  out=list(overall, accustats, reclassx)
  names(out)=c("overall","users/producers", "matrix")
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
    ## model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    ## calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})

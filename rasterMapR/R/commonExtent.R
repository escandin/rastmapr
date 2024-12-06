commonExtent <- structure(function #Common extent between raster objects
### This function calculates the common extent between all the raster
### objects within a list. The objects should all overlap at least
### partially and have the same projection, spatial resolution and be
### aligned. This can be written in a function. It can calculate
### either the largest extent that includes all the listed elements
### (union) other minimum common extent (intersection).
                       ##details<< This function ...
(
    inlist, ##<<\code{}. 
    method="intersection" ##<<\code{}. 
){
  # calculates the common extent between all the raster objects within
  # a list.  The objects should all overlap at least partially and
  # have the same projection, spatial resolution and be aligned This
  # can be written in a function. It can calculate either the largest
  # extent that includes all the listed elements (union) orthe minimum
  # common extent (intersection).
  for (i in 1: length(inlist)){
    ext=as.matrix(ext(inlist[[i]]))
    if(i==1){
      extall=ext} else {extall=rbind(extall, ext)}
  }
  
  extall=data.frame(extall)
  if (method=="intersection"){
  minx=max(extall$min[c(seq(1,nrow(extall), by=2))])
  maxx=min(extall$max[c(seq(1,nrow(extall), by=2))])
  miny=max(extall$min[c(seq(2,nrow(extall), by=2))])
  maxy=min(extall$max[c(seq(2,nrow(extall), by=2))]) 
  } else if (method=="union"){
    minx=min(extall$min[c(seq(1,nrow(extall), by=2))])
    maxx=max(extall$max[c(seq(1,nrow(extall), by=2))])
    miny=min(extall$min[c(seq(2,nrow(extall), by=2))])
    maxy=max(extall$max[c(seq(2,nrow(extall), by=2))])
  } else {stop("unvalid method name")
    }
  e=ext(minx, maxx, miny, maxy)
  return(e)
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

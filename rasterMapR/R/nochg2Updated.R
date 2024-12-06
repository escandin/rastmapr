
nochg2Updated <- structure(function #Change and no change mask
### This function computes a binary change/no-change mask processing a
### threshold raster, see \code{\link{thresraster}}
                       ##details<< This function is implemented by
                         ##\code{\link{PIFmodel}}. The mask is
                         ##computed comparing values in the threshold
                         ##raster with the chi-square statistic, see
                         ##\code{\link{qchisq}}.
(
    thresraster, ##<<\code{RasterLayer}, such as that produced by
                 ##\code{\link{thresraster}}.
    pvalue = 1e-02,
    distype="gamma",
    fitmethod="mle",
    propsamp=0.01,
    degfree=NULL ##<<\code{numeric}. Degrees of freedom.
){
    ## filename="gammadist") {
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
      require(fitdistrplus)
      samp=terra::spatSample(thresraster, ncell(thresraster)*propsamp, na.rm=TRUE)
      if(distype=="chisq"){
    threshold = qchisq(pvalue, degfree)
    #simulated=rchisq(length(samp), degfree)
    fit.stats<-ks.test(samp,  pchisq, df=6)} else {
    #fit.stats=list(samp, simulated, ks)} else {
                                           ## require(fitdistrplus)
    #if(cca==TRUE) {thresraster=thresraster[[1]]}
      fit.stats <- fitdistrplus::fitdist(samp$sum, dist =distype)
      threshold= qgamma(pvalue, shape=fit.stats$estimate[1],
                        rate=fit.stats$estimate[2])
                                       }
  if(threshold < terra::global(thresraster, fun="min", na.rm=TRUE)){
    warning("the pvalue selected produces a threshold that is larger that the minimum 
        value in thresraster producing an NA. Select a higher pvalue to produce a valid threshold mask")
  nochgmsk=NA} else {nochgmsk=maskfunUpdated(thresraster, threshold, NA, 1)}
  nochgmsk=list(nochgmsk, fit.stats)
  return(nochgmsk)
### \code{}. 
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    thrs <- thresraster2(strips[[2L]], strips[[1L]])
    noch <- nochg2Updated(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1, propsamp=1)
    ## plot(noch)
})

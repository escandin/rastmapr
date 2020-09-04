ecdf.comp<-function (thresraster, distype="gamma",
                  fitmethod="mle", propsamp=0.01, degfree=NULL){#degfree, filename="gammadist") {
  # ecdf.comp: a file containing the empirical vs theoretical distribution values
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
  
  #print("calculating change/no-change mask")
  if (!is.numeric(propsamp) & propsamp<0 & propsamp>1){
    stop("propsamp should be a number between 0 and 1")
  }
  samp=raster::sampleRandom(thresraster, length(thresraster)*propsamp)
  if(distype=="chisq"){
    #threshold = qchisq(pvalue, degfree)
    simulated=rchisq(length(samp), degfree)
    ks<-ks.test(samp, simulated)
    fit.stats=list(samp, simulated, ks)} else {
      require(fitdistrplus)
      #if(cca==TRUE) {thresraster=thresraster[[1]]}
      fit.stats <- fitdistrplus::fitdist(samp, dist =distype)
      simulated=rgamma(length(samp), shape=fit.stats$estimate[1], rate = fit.stats$estimate[2])
      fit.stats <- list(samp, simulated, fit.stats)
      #threshold= qgamma(pvalue, shape=fit.stats$estimate[1],
       #                 rate=fit.stats$estimate[2])
    }
  # if(threshold < minValue(thresraster)){
  #   warning("the pvalue selected produces a threshold that is larger that the minimum 
  #       value in thresraster producing an NA. Select a higher pvalue to produce a valid threshold mask")
  #   nochgmsk=NA} else {nochgmsk=maskfun(thresraster, threshold, NA, 1)}
  # nochgmsk=list(nochgmsk, fit.stats)
  return(fit.stats)
}

library(raster)
library(rgdal)
dir="X:/VictorShare/s3dFiles/Orinoquia"
setwd(dir)
tar.ref=list(stack("LC08_L1TP_005057_20140128_20170426_01_T1_mskd_.grd"),
             stack("LC08_L1TP_005057_20160118_20180130_01_T1_mskd_.grd"))
common=commonExtent(tar.ref)
tar.ref=list(crop(tar.ref[[1]], common), crop(tar.ref[[2]], common))

# Before masking
sssd=thresraster2(tar.ref[[2]], tar.ref[[1]], propsamp=0.01)
compare=ecdf.comp(sssd)

# After masking
#msk=raster(paste(dir, "Outputsgamma/LC08_L1TP_005057_20140128_20170426_removeChgMsk14.tif", sep="/"))

tar.ref.mskd=list(mask(tar.ref[[1]], msk), mask(tar.ref[[2]], msk))
sssdmskd=thresraster2(tar.ref.mskd[[2]], tar.ref.mskd[[1]], propsamp=0.01)
comparemskd=ecdf.comp(sssdmskd)

# Plot before vs after
xmax.=40
pdf("ECDFbeforeAfterGamma.pdf", width=4,height=4,paper='special')
  plot(ecdf(compare[[1]]), xlim=c(0,xmax.), do.points=FALSE)
  lines(ecdf(compare[[2]]), col="red")
  plot(ecdf(comparemskd[[1]]), xlim=c(0,xmax.), do.points=FALSE)
  lines(ecdf(comparemskd[[2]]), col="red")
dev.off()
             


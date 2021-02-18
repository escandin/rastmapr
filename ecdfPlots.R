library(raster)
library(rgdal)
library(fitdistrplus)
ecdf.comp<-function (thresraster, distype="gamma",
                  fitmethod="mle", propsamp=0.01, degfree=NULL){#degfree, filename="gammadist") {
  # thresraster: produces a comparison of the empirical cummulative density functions
  # between an empirical vs a theoretical distribution specified by distype. The memprical distribution
  # is calculated using a sample of pixels randomly extraced from a raster layer specified by thresraster 
  # with an intensity specified by propsamp.
  # This function is a modification of the function noch2
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

#SETTINGS
# MAC
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Orinoquia"
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Pucallpa"
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Mexico"
# WINDOWS
#path=#("X:/VictorShare/s3dFiles/Pucallpa")
#("X:/VictorShare/s3dFiles/Orinoquia")
#("X:/VictorShare/s3dFiles/Mexico")#

tarRast=#"LC08_L1TP_005057_20140128_20170426_01_T1_mskd_.grd" #Orinoquia
#"LT05_L1TP_006066_20110807_20161007_01_T1_mskd_.grd" #Pucallpa
"LT05_L1TP_026046_20100205_20161016_01_T1_mskd_.grd" #Mexico
  
refRast= #"LC08_L1TP_005057_20160118_20180130_01_T1_mskd_.grd" #Orinoquia
#"LC08_L1TP_006066_20150818_20170406_01_T1_mskd_.grd" #Pucallpa
"LC08_L1TP_026046_20180110_20180119_01_T1_mskd_.grd" #Mexico

distr="chisq"
canoco=TRUE

#ORINOQUIA
mask.=  #"OutputsGamma/LC08_L1TP_005057_20140128_20170426_removeChgMsk14.tif" # Gamma
#"OutputsGammaCCA/LC08_L1TP_005057_20140128_20170426_removeChgMsk8.tif" # GammaCCA
  "OutputsChsqCCA/LC08_L1TP_005057_20140128_20170426_removeChgMsk12.tif" #Chsq

#PUCALLPA#
mask.=#"OutputsGamma/LT05_L1TP_006066_20110807_20161007_removeChgMsk8.tif" #Gamma
  #"OutputsGammaCCA/LT05_L1TP_006066_20110807_20161007_removeChgMsk8.tif" #GammaCCA
  "OutputsChsqCCA/LT05_L1TP_006066_20110807_20161007_removeChgMsk6.tif" #Chsq
  
# MEXICO
mask.=#"OutputsGamma/LT05_L1TP_026046_20100205_20161016_removeChgMsk8.tif" #Gamma
  #"OutputsGammaCCA/LT05_L1TP_026046_20100205_20161016_removeChgMsk5.tif" #GammaCCA
  "OutputsChsqCCA/LT05_L1TP_026046_20100205_20161016_removeChgMsk6.tif" #Chsq

pdfFilename="ECDFbeforeAfterChsqCCA.pdf" #ChsqCCA
outnameRast="sssdChsqCCAfirst.tif"#"sssdChsqCCAfirst.tif"
outnameRastmskd="sssdChsqCCAskd5.tif"#"sssdChsqCCAMskd8.tif"

setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

tar.ref=list(stack(tarRast),
             stack(refRast))

common=commonExtent(tar.ref)
tar.ref=list(crop(tar.ref[[1]], common), crop(tar.ref[[2]], common))

# Before masking
sssd=thresraster2(tar.ref[[2]], tar.ref[[1]], cca=canoco)
compare=ecdf.comp(sssd, distype=distr, degfree=6)

# After masking
msk=raster(paste(path, mask., sep="/"))

tar.ref.mskd=list(mask(tar.ref[[1]], msk), mask(tar.ref[[2]], msk))
sssdmskd=thresraster2(tar.ref.mskd[[2]], tar.ref.mskd[[1]],  cca=canoco)
comparemskd=ecdf.comp(sssdmskd, distype=distr, degfree=6)

# Plot before vs after
xmax.=40
pdf(pdfFilename, width=4,height=4,paper='special')
  plot(ecdf(compare[[1]]), xlim=c(0,xmax.), do.points=FALSE)
  lines(ecdf(compare[[2]]), col="red")
  plot(ecdf(comparemskd[[1]]), xlim=c(0,xmax.), do.points=FALSE)
  lines(ecdf(comparemskd[[2]]), col="red")
dev.off()
writeRaster(sssd, outnameRast)             
writeRaster(sssdmskd, outnameRastmskd)   

gofstat(compare[[3]])
gofstat(comparemskd[[3]])
s3dmod$paramstats[[2]]


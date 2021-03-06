library(RStoolbox)
library(raster)
library(rgdal)
#library(rasterMapR)
library(fitdistrplus)
library(R.utils)
library(profvis)
sourcedir <- '/Users/tug61163/Documents/Repositories/rastmapr/rasterMapR/R'
sourceDirectory(sourcedir, modifiedOnly=FALSE)

# MAC
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Orinoquia"
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Pucallpa"
path="/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/Mexico"
#paths=c(path1, path2, path3)
# WINDOWS
path=#("X:/VictorShare/s3dFiles/Pucallpa")
  #("X:/VictorShare/s3dFiles/Orinoquia")
#("X:/VictorShare/s3dFiles/Mexico")#

setwd(path)
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

tares <- list.files('.', pattern='tar.gz')
#stacks<- EEstackWithoutMeta(tares, sat.nm="LO08")
#stacks<- EEstackWithoutMeta(tares, sat.nm="LT05")

stacknames=names(stacks)
layernames=names(stacks[[1]])

#PRODUCE MASKED STACKS
# SPECIFY DIFFERENT SETTINGS FOR OLI, TM
sat="TM" 
normbands=seq(1,6) # for TM
#sat="OLI" 
#normbands=seq(2,7) # for OLI

#cloud masking
for (i in 3:4){ # MAKE SURE THE selected elements correspond to the sat.nm in EEstackWithoutMeta
  stacks[[i]]=smg(inlist=stacks, method= "none", mosaicitems=i, 
                  normbands=normbands, sensor=sat)
  names(stacks[[i]])=layernames[normbands]
  writeRaster(stacks[[i]], filename=paste(names(stacks)[i], "mskd_", sep="_"), 
              format="raster",  datatype="INT2S") # the raster format is to preserve the band names
}

####### Retrieve datasets as a stack
stacks=list()
normbands=seq(1,6)
setwd(path)
stacknames<- list.files('.', pattern='mskd_.grd')

for (i in 1:length(stacknames)){
    stacks[[i]]=raster::stack(stacknames[i], bands=normbands)
    #names(stacks[[i]])=layernames[normbands]
    #writeRaster(stacks[[i]], filename=paste(names(stacks)[i], "mskd_", sep="_"), 
    #            format="raster",  datatype="INT2S")
 }

stacknames= substr(stacknames, 1, nchar(stacknames)-16)
names(stacks)=stacknames
plotRGB(stacks[[1]], r=4, g=3, b=2, stretch="lin")

ref=c(2,4) #Orinoquia
c(1,2) 
tar=c(1,3) #Orinoquia
c(3,4) 
i=1

for(i in 1:length(ref)){
  instacks=list(stacks[[tar[[i]]]], stacks[[ref[[i]]]])
  names(instacks)=c(stacknames[[tar[[i]]]], stacknames[[ref[[i]]]])
  #a=Sys.time()
test=  profvis({
  s3dmod=s3d(strips=instacks, thres=1e-2, distype="gamma",
             pval.pif=1e-3,  pval.chg=0.99, cca=FALSE, fitline=FALSE,
             prefix=names(instacks)[1], writemasks=FALSE)
  }, height = "400px")
  #b=Sys.time()-a
  save(test, file=paste(path, (paste(paste(names(instacks)[1],
                                           's3d_GammaPerform', sep="_"),  
                                     "RData", sep=".")), sep="/"))  
  #save(s3dmod, file=paste(path, (paste(paste(names(instacks)[1],
                                           #  's3d_GammaItertime', sep="_"),  
                                       #"#RData", sep=".")), sep="/"))
}


#######################################################
##### DEBUGGING TESTS
## For test data

###### Arguments from s3d
plotRGB(stacks[[2]], r=4, g=3, b=2, stretch="lin")
instacks=list(stacks[[tar[[i]]]], stacks[[ref[[i]]]])
e=drawExtent()
instacks=Map(function(x)
  crop(x, e), raster::as.list(instacks))
plotRGB(instacks[[2]], r=4, g=3, b=2, stretch="lin")
### Adjusted
cca=TRUE
pval.pif = 1e-02 
distype="chisq"
thres=0.01
### Default
strips=instacks
distsamp=0.01
pval.chg=0.99
minPIF=5
maxiter=20 
prefix="test"
norm.ext=NULL
fitline=TRUE
writemasks=TRUE
#pvalpif=1e-3 # Pucallpa
#refvec=seq(reftar[1], reftar[length(reftar)-1], 2)

#Arguments for thresraster2
refstack=instacks[[2]]
tarstack=instacks[[1]]
propsamp=1

#Arguments for nochg
thresrater=sumstandardizediff
pvalue=pval.chg
fitmethod="mle"
propsamp=0.01
degfree=6


s3dmod=s3d(strips, thres=0.01, pval.pif=0.005,  pval.chg=0.99, cca=TRUE, 
           prefix="test", distype="chisq")


# FIGURES
# Distribution before correction

thresrasterTemp<-function (refstack, tarstack, cca=FALSE, propsamp=0.01){
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
    # stdif1[[k]]=stadif(dif, k)
    print(paste("standardized differences calculated for band", k, sep = " "))
  }
  sumstandardizediff = calc(standardizediff^2, sum)#, na.rm = TRUE)
  out=list(sumstandardizediff, standardizediff)
  if(cca==TRUE){
    out=list(out, cct)
    names(out)=c("sumstandardizediff", "CCT")} 
  return(out)
}
setwd(paste(path, "Outputs06057", sep="/"))

# BEFORE
test=thresrasterTemp(instacks[[2]], instacks[[1]])
for(i in 1:nlayers(instacks[[2]])){
  samp=sampleRandom(test[[2]][[i]], 
                    ncell(test[[2]][[i]])*0.01, na.rm=TRUE)
  fit.norm1 <- fitdistrplus::fitdist(samp,distr = "norm", method = "mle")
  summary(fit.norm1)
  simulated=rnorm(length(samp), mean=fit.norm1$estimate[1], 
                  sd = fit.norm1$estimate[2])
  png(file=paste(paste("NormBefore", i, sep="_"),
        "png", sep=".")) 
    plot(fit.norm1)
    plot(ecdf(samp), xlim=c(-3,3), do.points=FALSE)
    lines(ecdf(simulated), col="darkgreen")
  dev.off()
}
samp=sampleRandom(test[[1]], ncell(test[[1]])*0.01, na.rm=TRUE)
fit.gamma <- fitdist(samp, distr = "gamma", method = "mle")
gammastats=gofstat(fit.gamma)
refgamma=rgamma(length(samp), shape=fit.gamma$estimate[1], rate = fit.gamma$estimate[2])
refchisq=rchisq(length(samp), df=nlayers(instacks[[2]]))
png(file="GammaBefore.png") 
  plot(fit.gamma)
  plot(ecdf(samp), xlim=c(0,20), do.points=FALSE)
  lines(ecdf(refgamma), col="darkgreen")
  lines(ecdf(refchisq), col="red")
dev.off()

#  AFTER
masks <- list.files('.', pattern='_removeChgmsk')
msk=raster(masks[length(masks)])
tarmskd=mask(instacks[[1]], msk)  
refmskd=mask(instacks[[2]], msk)
test=thresrasterTemp(refmskd, tarmskd)

for(i in 1:nlayers(instacks[[2]])){
  samp=sampleRandom(test[[2]][[i]], 
                    ncell(test[[2]][[i]])*0.01, na.rm=TRUE)
  fit.norm1 <- fitdistrplus::fitdist(samp,distr = "norm", method = "mle")
  summary(fit.norm1)
  simulated=rnorm(length(samp), mean=fit.norm1$estimate[1], 
                  sd = fit.norm1$estimate[2])
  png(file=paste(paste("NormAfter", i, sep="_"),
                 "png", sep="."))
    plot(fit.norm1)
    plot(ecdf(samp), xlim=c(-3,3), do.points=FALSE)
    lines(ecdf(simulated), col="darkgreen")
  dev.off()
}
samp=sampleRandom(test[[1]], ncell(test[[1]])*0.01, na.rm=TRUE)
fit.gamma <- fitdist(samp, distr = "gamma", method = "mle")
gammastats=gofstat(fit.gamma)
refgamma=rgamma(length(samp), shape=fit.gamma$estimate[1], rate = fit.gamma$estimate[2])
refchisq=rchisq(length(samp), df=nlayers(instacks[[2]]))
png(file="GammaAfter.png") 
  plot(fit.gamma)
  plot(ecdf(samp), xlim=c(0,20), do.points=FALSE)
  lines(ecdf(refgamma), col="darkgreen")
  lines(ecdf(refchisq), col="red")
dev.off()



test=PIFmodel2(strips=instacks, pvalue=3e-2,  pval.chg=0.99, cca=FALSE, prefix=names(instacks)[1])
library(rasterMapR)
setwd("/Users/tug61163/Documents/MyResearch/Collaborations")
dir.create('tempfiledir')
tempdir=paste(getwd(),'tempfiledir', sep="/")
rasterOptions(tmpdir=tempdir)

tares <- list.files('.', pattern='tar.gz')
stacks<- EEstackWithoutMeta(tares, sat.nm="LO08")

stacks[[1]][stacks[[1]]==0]<-NA
stacks[[2]][stacks[[2]]==0]<-NA

mosaick<- smg(inlist=stacks, method= "none", refitem=2, mosaicitems=c(2,1),
              pval.pif = 1e-02, pval.chg=0.99, QAbandname=NA, indem=NA,#demfile,  # DEM file only for Mexico
              normbands=seq(1, 7), verbose=TRUE, sensor="OLI", norm.ext=NULL, cloudbuff=NA)
plotRGB(mosaick, r=5, g=4, b=3, stretch="lin")

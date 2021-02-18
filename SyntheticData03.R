
# FROM https://rdrr.io/github/debruine/faux/src/R/rnorm_multi.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(faux)
library(fitdistrplus)
library(Rmisc)

setwd("/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/ModelEvaluation")

#Produce a dataset with random variables sorted with three levels of correlation. 
# The first column denotes the level of correlation between variables (1= low, 2= middle or 3= high).
# The correlation levels are between 0 and 1/3 for 1, 1/3 to 2/3 for 2 and 2/3 to 0.95 for 3

nbands=6
# in gg plot format
for (k in 1:20){
  for (i in 1:100){
    for(j in 1:(nbands/2)){ # It will produce pairs of correlated bands
      #rval=runif(1, min=rmin[selgroup[j]], max=rmax[selgroup[j]])
      rval=0#(k-1)/20
      out <- mvrnorm(100, mu = c(0,0), Sigma = matrix(c(1,rval,rval,1), ncol = 2),
                     empirical = TRUE)
      if(j==1){syntdata=out} else {syntdata=cbind(syntdata, out)}
    }
    
    s3d=apply(syntdata^2,1, sum) 
    chisqdata=rchisq(length(s3d), 6)
    #kschisq<-ks.test(s3d, chisqdata)
    kschisq1<-ks.test(s3d, pchisq, df=6)
    
    fit.stats <- fitdistrplus::fitdist(s3d, dist ="gamma")
    gammadata=rgamma(length(s3d), shape=fit.stats$estimate[1], rate = fit.stats$estimate[2])
    ksgamma=gofstat(fit.stats) 
    dval=cbind(kschisq1$statistic, ksgamma$ks)
    if (i ==1){dvaltot=dval} else {dvaltot=rbind(dvaltot,dval)}
  }
  means=apply(dvaltot, 2, mean)
  sds=apply(dvaltot, 2, sd)
  outstats=data.frame(cbind(rbind(rval, rval),
                            rbind(as.numeric(means[1]),as.numeric(means[2])), 
                            rbind(as.numeric(sds[1]),as.numeric(sds[2]))),stringsAsFactors=FALSE)
  #outstats=data.frame(cbind(rval, means[1],means[2], sds[1],sds[2]))
  names(outstats)=c("rval","meanval", "sdval")
  if (k==1){totstats=outstats} else {totstats=rbind(totstats,outstats)}
}  
# I have to add this a-posteriori because otherwise it converts all columnns to character
totstats$dist=as.factor(rep(c("chisq","gamma"), nrow(totstats)/2))
save(totstats, file="syntdataGGplot.RData") 

pdfFilename="syntGammaVsXhisq.pdf"
pdf(pdfFilename, width=4,height=4,paper='special')
# Based on http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(totstats, aes(x=rval, y=meanval, colour=dist)) + 
  geom_errorbar(aes(ymin=meanval-sdval, ymax=meanval+sdval), width=.1) +
  geom_line() +
  geom_point()
dev.off()

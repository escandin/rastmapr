
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
for (k in 1:20){
  for (i in 1:100){
    for(j in 1:(nbands/2)){ # It will produce pairs of correlated bands
      #rval=runif(1, min=rmin[selgroup[j]], max=rmax[selgroup[j]])
      rval=(k-1)/20
      out <- mvrnorm(100, mu = c(0,0), Sigma = matrix(c(1,rval,rval,1), ncol = 2),
                     empirical = TRUE)
      if(j==1){syntdata=out} else {syntdata=cbind(syntdata, out)}
    }
    
    s3d=apply(syntdata^2,1, sum) 
    #chisqdata=rchisq(length(s3d), 6)
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

############################################################
### with correlated bands
for (k in 1:20){
  for (i in 1:100){
    for(j in 1:(nbands/2)){ # It will produce pairs of correlated bands
      #rval=runif(1, min=rmin[selgroup[j]], max=rmax[selgroup[j]])
      rval=(k-1)/20
      out <- mvrnorm(100, mu = c(0,0), Sigma = matrix(c(1,rval,rval,1), ncol = 2),
                     empirical = TRUE)
      #cor(out, method="pearson")
      # Rescale the second term to avoid negative values, apply log and then re-scale to zero
      #summary(out[,2])
      out[,2]= out[,2]-min(out[,2])
      out[,2]=log(out[,2])
      out[,2][which(out[,2]==-Inf)]=0
      out[,2]= (out[,2]-mean(out[,2]))/sd(out[,2])
      rval=cor(out, method="pearson")[2,1]
      if(j==1){syntdata=out} else {syntdata=cbind(syntdata, out)}
    }
    
    s3d=apply(syntdata^2,1, sum) 
    #chisqdata=rchisq(length(s3d), 6)
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
  if (k==1){totstats1=outstats} else {totstats1=rbind(totstats1,outstats)}
}  
# I have to add this a-posteriori because otherwise it converts all columnns to character
totstats1$dist=as.factor(rep(c("chisq","gamma"), nrow(totstats1)/2))

totstats1$distypes=rep("norm_lognorm", nrow(totstats1))
totstats$distypes=rep("norm", nrow(totstats))

allstats=rbind(totstats, totstats1)
names(allstats)=c("rval", "meanval", "sdval",  "dist", "distypes")
pdf(file='synthdata.pdf',
    width = 16, height=9, paper='special')
ggplot(allstats, aes(x=rval,y=meanval, color=dist))+#, fill=algorithm))+
  geom_errorbar(aes(ymin=meanval-sdval, ymax=meanval+sdval), width=.1) +
  facet_grid(cols=vars(distypes))+
  #theme(text=element_text(size=20))+
  #scale_y_continuous(breaks = seq(0, 1, by = .2))+
  theme_bw(base_size=24)+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Correlation (r)") + ylab("KS D statistic")+
  geom_line() + 
  scale_colour_manual(values = c("#FB9A99", "#8DA0CB"))+
  geom_point() #+ theme_bw()
dev.off()

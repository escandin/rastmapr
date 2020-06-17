library(fitdistrplus)
# Synthetic data and functions to test:
# 1. the effect of different emprical distributions and interdependence 
# of variables on the parameters of a gamma distribution
# 2. the effectiveness of different methods for mitigating the effect of 
# change on deviations from normality in the empirical distribution 
# of variables on a gamma distribution

z.score=function(x){
  x.z=(x-mean(x))/sd(x)
  return(x.z)
}
test.gamma=function(nvar=4, meanvar=0, sdvar=1, ndraws=1000,
                    vardep=c(1,2,3,4), sdvartimes=c(2,1)){
  # nvar: number of input variables to simulate
  # meanvar, sdvar: parameters of the normal distribution for each input variables
  # ndraws: number of entries per variable
  # vardep: index for pairs of interdependent variables. The first term of
  # each pair is the independent and the second is the dependent
  # corlevel: the level of correlation between dependent variables.
  # the length of this vector should be equal to ahlf the size of vardep
  # chgvar: variable to represent change
  for(i in 1:nvar){
    varvec=rnorm(ndraws,meanvar,sdvar)
  if(i==1){syndata=varvec} else {syndata=cbind(syndata,varvec)}
  syndata=data.frame(syndata)
  # assign correlation level between variables and renormalize them
  }
  if (class(sdvartimes)=="numeric"){
    if(length(sdvartimes)!=length(vardep)/2){
      return("sdvartimes should have a length equal to the length of vardep/2")}
    for(i in seq(2, length(vardep), by=2)){
      syndata[,i]=syndata[,i-1]+rnorm(ndraws, meanvar, (sdvar*sdvartimes[i/2]))
      syndata[,i]=z.score(syndata[,i])
    }
  }
  sumsyndatasq=apply(syndata^2,1, sum)
  gammadist=fitdist(sumsyndatasq, dist="gamma", method="mle")
  #fitdist(sumsyndatasq[[2]], dist="gamma", method="mle")
  return(list(syndata, sumsyndatasq, gammadist))
}

# IMPLEMENTATION
sdvec=c(0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8)
iter=(length(sdvec))
for(i in seq(1, iter, by=1)){
gammafun=test.gamma(nvar=4, meanvar=0, sdvar=1, ndraws=1000, 
                    vardep=c(1,2,3,4), 
                    sdvartimes=c(sdvec[i], (12.8-sdvec[i])))
hist(gammafun[[2]], freq=FALSE, 
     breaks=seq(0,ceiling(max(gammafun[[2]])), by=1), 
     xlim=c(1,20), ylim=c(0, 0.6))
out=cbind(cor(gammafun[[1]][,2], gammafun[[1]][,1], method="pearson"),
      gammafun[[3]]$estimate[1],gammafun[[3]]$estimate[2])
if(i==1){outot=out}else{outot=rbind(outot,out)}
}
outot=data.frame(outot)
names(outot)=c("cor", "shape", "rate")

plot(outot$cor, outot$shape)
plot(outot$cor, outot$rate)
plot(outot$shape, outot$rate)






a=rnorm(1000,0,1)
b=rnorm(1000,0,1)
c=rnorm(1000,0,1)
d=rnorm(1000,0,1)
abcd=data.frame(cbind(a,b,c,d))
sumabcdsq=apply(abcd^2, 1, sum)
test=fitdist(sumabcdsq, dist="gamma", method="mle")

gammafun=test.gamma(4, 0, 1, 1000, sdvartimes=NA)
test=fitdist(gammafun[[2]], dist="gamma", method="mle")
summary(test)
plot(test)

cor(syndata[,3], syndata[,1], method="pearson")
test=fitdist(sumsyndatasq, dist="gamma", method="mle")

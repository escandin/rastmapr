
calp=s3dmod
par(mfrow=c(2,3))
for (b in 1:nlayers(strips[[1]])){
  plot(calp$data[[b]][,1]~calp$data[[b]][,2], xlab="target", ylab="reference")
  abline(unlist(calp$parameters[[b]][1])[1], unlist(calp$parameters[[b]][1])[2])
}

lmparamtot=s3dmod$paramstats[[1]]
strips=stacks
par(mar = c(4, 4, 1, 1) + 0.1)
par(mfrow=c(2,3))
for(b in 1:6){
  plot(lmparamtot$iter[which(lmparamtot$band==b)], 
       lmparamtot$intercept[which(lmparamtot$band==b)], xlab="band", ylab="intercept")
}
par(mfrow=c(2,3))
for(b in 1: 6){
  plot(lmparamtot$iter[which(lmparamtot$band==b)], 
       lmparamtot$slope[which(lmparamtot$band==b)], xlab="band", ylab="slope")
}

####### ENVIRONMENT SETTINGS
library(stringr)

path=#("X:/VictorShare/aRnFiles/Pucallpa")
#  ("X:/VictorShare/aRnFiles/Orinoquia")
("X:/VictorShare/aRnFiles/Mexico")#

savepath=(paste(path, "outputsBandsAll", sep="/"))
savepath1=(paste(path, "outputsBands3_7", sep="/"))
setwd(savepath1)

# Compare users' producer's overall accuracy for REF and TAR years
normnames=c("classnames", "aRn", "ed", "none")

# For REF year
load("mosaicREFaRn_accumatrix.RData")
aRnStats<-AccuPlot(conmatrix)
load("mosaicREFed_accumatrix.RData")
edStats <-AccuPlot(conmatrix)
load("mosaicREFnone_accumatrix.RData")
noneStats <-AccuPlot(conmatrix)

#load overall accuracy
c(aRnStats[1], edStats[1], noneStats[1])

# Graph users, producers for all the methods
crossusers=data.frame(cbind(aRnStats[[2]][,c(1,2)],edStats[[2]][,2],noneStats[[2]][2]),
                      stringsAsFactors=FALSE)
names(crossusers)=normnames
crossproducers=cbind(aRnStats[[2]][,c(1,3)],edStats[[2]][,3],noneStats[[2]][,3])
names(crossproducers)=normnames

pdf("UsersProducersREFBands3_7DF.pdf", width=6,height=4,paper='special') 
par(mar = c(7, 4, 2, 2) + 0.2)
barplot(rbind(crossusers[,2], crossusers[,3],crossusers[,4]),col=c("aquamarine3","coral", "cornflowerblue"), 
        names.arg=crossusers[,1],  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Users accuracy")
legend("top", legend=c("aRn", "ed", "none"), box.lty=0, bg= "transparent",
       col=c("aquamarine3","coral", "cornflowerblue"), lty=1:2, cex=0.8)

barplot(rbind(crossproducers[,2], crossproducers[,3],crossproducers[,4]),col=c("aquamarine3","coral", "cornflowerblue"), 
        names.arg=crossproducers[,1],  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Producers accuracy")
legend("top", legend=c("aRn", "ed", "none"), box.lty=0, bg= "transparent",
       col=c("aquamarine3","coral", "cornflowerblue"), lty=1:2, cex=0.8)
dev.off()

# For TAR year
### Calculate statistics from all iteration
setwd(savepath1)
iter=10
aRnStats <-AccuStats(prefix="classTAR_aRn_iter", suffix="accumatrix.RData", iter=iter)
edStats <-AccuStats(prefix="classTAR_ed_iter", suffix="accumatrix.RData", iter=iter)
noneStats <-AccuStats(prefix="classTAR_none_iter", suffix="accumatrix.RData", iter=iter)
mean(aRnStats[[3]])
sd(aRnStats[[3]])
mean(edStats[[3]])
sd(edStats[[3]])
mean(noneStats[[3]])
sd(noneStats[[3]])

crossusers=cbind(aRnStats[[4]][,c(1,2)],edStats[[4]][,2],noneStats[[4]][2])
names(crossusers)=normnames
crossproducers=cbind(aRnStats[[4]][,c(1,4)],edStats[[4]][,4],noneStats[[4]][4])
names(crossproducers)=normnames

pdf("UsersProducersTARBands3_7Orinoquia.pdf", width=6,height=4,paper='special') 

par(mar = c(7, 4, 2, 2) + 0.2)
barplot(rbind(crossusers[,2], crossusers[,3],crossusers[,4]),col=c("aquamarine3","coral", "cornflowerblue"), 
        names.arg=crossusers[,1],  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Users accuracy")
legend("top", legend=c("aRn", "ed", "none"), box.lty=0, bg= "transparent",
       col=c("aquamarine3","coral", "cornflowerblue"), lty=1:2, cex=0.8)

barplot(rbind(crossproducers[,2], crossproducers[,3],crossproducers[,4]),col=c("aquamarine3","coral", "cornflowerblue"), 
        names.arg=crossproducers[,1],  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Producers accuracy")
legend("top", legend=c("aRn", "ed", "none"), box.lty=0, bg= "transparent",
       col=c("aquamarine3","coral", "cornflowerblue"), lty=1:2, cex=0.8)
dev.off()

#############################################################
# POST CLASSIFICATION ANALYSIS

setwd(savepath1)

# ASSESS OPTIMAL RESAMPLING ITERATIONS
optstats=aRnStats
pref="classTAR_aRn_iter"
for (i in 2: iter){
  optstats <-AccuStats(prefix=pref, suffix="accumatrix.RData", iter=i)
  test=cbind(mean(optstats[[3]]), sd(optstats[[3]]), sd(optstats[[3]])/mean(optstats[[3]]))
  if (i==2){metrics=test} else {metrics=rbind(metrics, test)}
}
metrics=data.frame(cbind(seq(2,iter), metrics))
names(metrics)=c("iter", "mean", "sd", "cv")
plot(metrics$iter, metrics$mean, type="l")
plot(metrics$iter, metrics$sd, type="l")
plot(metrics$iter, metrics$cv, type="l")

######################################################
############### sampling size testing
# Post classification class merging
#fromtoclass=c(2,1, 4,6, 5,10, 9,10, 8,7, 11,10, 12,6, 17,10) # ORINOQUIA
fromtoclass=c(2,1, 4,3, 5,10, 9,10, 8,7, 11,10, 12,6, 17,10) # ORINOQUIA
fromtomat=matrix(fromtoclass, ncol=2, byrow=T)

load("classTAR_ed_iter_2_accumatrix.RData")
test=conmatrixmerge(reclassx=conmatrix$table, rcl=fromtomat)
test$overall
test[[2]]

#### This is to produce a correlation matrix between the accuracy of different categories
# Transpose data frames
allusersT=as.data.frame(t(as.matrix(allusers[,-1])))
names(allusersT)=paste(as.character(allusers[,1]), "U", sep="_")

allproducersT=as.data.frame(t(as.matrix(allproducers[,-1])))
names(allproducersT)=paste(as.character(allproducers[,1]), "P", sep="_")

allaccu=cbind(allusersT, allproducersT)

# Correlation matrix
coran=cor(allaccu, use="complete.obs", method="pearson")
write.csv(coran, file='coran.csv')

plot(allaccu$River_U, allaccu$Lake_P)

### THIS IS TO ASSESS HOW MANY ITERATIONS ARE NEEDED FOR STABILIZING ACCURACY
for (iter in 3:10){
  aRnStats <-AccuStats(prefix="classTAR_aRn_iter", suffix="accumatrix.RData", iter=iter)
  stats=cbind(mean(aRnStats[[3]]), sd(aRnStats[[3]]))
  if (iter==3){cumstats=stats} 
  else {
    cumstats=rbind(cumstats, stats)}
}

attributes(class.RF$terms)$term.labels


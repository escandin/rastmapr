
### Calculate statistics from all iterations
setwd("/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/aRnManuscript/Data/Mexico")
#setwd(paste(getwd(), "outputsBandsAllV2", sep="/"))

normethods=c("imist", "imist2","ed", "cor", "none")
iter=3
imistStats <-AccuStats(prefix="classTAR_imist_iter", suffix="accumatrix.RData", iter=iter)
imist2Stats <-AccuStats(prefix="classTAR_AllBands_imist_iter", suffix="accumatrix.RData", iter=iter)
edStats <-AccuStats(prefix="classTAR_ed_iter", suffix="accumatrix.RData", iter=iter)
corStats <-AccuStats(prefix="classTAR_cor_iter", suffix="accumatrix.RData", iter=iter)
noneStats <-AccuStats(prefix="classTAR_none_iter", suffix="accumatrix.RData", iter=iter)

overalldata <- data.frame(
  name=normethods,
  value=c((1-mean(imistStats[[3]])), (1-mean(imist2Stats[[3]])), (1-mean(edStats[[3]])), (1-mean(corStats[[3]])),
          (1-mean(noneStats[[3]]))),
  sd=c(sd(imistStats[[3]]), sd(imist2Stats[[3]]), sd(edStats[[3]]), sd(corStats[[3]]), sd(noneStats[[3]]))
)

crossusers=data.frame(cbind(imistStats[[4]][,c(1,2)],imist2Stats[[4]][,2], edStats[[4]][,2], corStats[[4]][,2], noneStats[[4]][2]),
                      stringsAsFactors=FALSE)
names(crossusers)=c("classnames", "imist", "imist2", "ed", "cor", "none")
crossproducers=cbind(imistStats[[4]][,c(1,4)], imist2Stats[[4]][,4], edStats[[4]][,4], corStats[[4]][,4], noneStats[[4]][4])
names(crossproducers)=names(crossusers)

graphcol=col=c("aquamarine3","coral", "cornflowerblue", "darkgoldenrod3", "gray")
pdf("ProducersUsersMex.pdf", width=4,height=4,paper='special') 
par(mar = c(7, 4, 2, 2) + 0.2)
barplot(rbind(crossproducers[,2], crossproducers[,3],crossproducers[,4], crossproducers[,5], crossproducers[,6]),
        col=graphcol, 
        names.arg=crossproducers$classnames,  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Producers accuracy")
legend("top", legend=c("imist", "imist2", "ed", "cor", "none"), box.lty=0, bg= "transparent",
       col=graphcol, lty=1:2, cex=0.8)
barplot(rbind(crossusers[,2], crossusers[,3],crossusers[,4], crossusers[,5], crossusers[,6]),
        col=graphcol, 
        names.arg=crossusers$classnames,  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7, 
        main="Users accuracy")
legend("top", legend=c("imist", "imist2", "ed", "cor", "none"), box.lty=0, bg= "transparent",
       col=graphcol, lty=1:2, cex=0.8)
ggplot(overalldata) +
  geom_bar(aes(x=name, y=value), stat="identity", fill=graphcol, alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="gray28", alpha=0.9, size=1.3)+
  #coord_cartesian(ylim = c(0, .7)) #ORINOQUIA
  coord_cartesian(ylim = c(0, .2)) #PUCALLPA, MEXICO
dev.off()

#############################################################
# POST CLASSIFICATION ANALYSIS
######################################################

fromtoclass=c(4,3, 5,10, 9,10, 8,7, 11,3, 12,6, 17,10)# ORINOQUIA
fromtoclass=c(4,3, 5,10, 9,10, 8,7, 11,3, 12,6, 13,14, 17,10)# ORINOQUIA V2
fromtoclass=c(1,6, 2,5, 4,9, 15,10)        # PUCALLPA
fromtoclass=c(1,6, 2,5, 4,9, 12,14, 13,9, 15,10)        # PUCALLPA V2
fromtoclass=c(9,8)                                     # MEXICO
fromtoclass=c(3,2, 7,4, 8,1, 9,1, 11,10)                                     # MEXICO V2
#fromtoclass=c(5,2, 7,6)                                     # MEXICO
fromtomat=matrix(fromtoclass, ncol=2, byrow=T)

normethods=c("imist", "AllBands_imist", "ed", "cor", "none")
iter=3
for (j in 1:length(normethods)){
  for(i in 1:iter){
    load(paste((paste(paste(paste("classTAR", normethods[j], sep="_"), 
                            "iter", sep="_"), i, sep="_")), "accumatrix.RData", sep="_"))
    test=conmatrixmerge(reclassx=conmatrix$table, rcl=fromtomat)
    if(i==1){
      overall=test$overall
      users=test[[2]][,1]
      producers=test[[2]][,2]
    } else {
      overall=c(overall, test$overall)
      users=cbind(users, test[[2]][,1])
      producers=cbind(producers, test[[2]][,2])
    }
  }
  if (j==1){
    alloverallmean=mean(overall)
    alloverallsd=sd(overall)
    allusersmean=apply(users, 1, "mean")
    alluserssd=apply(users, 1, "sd")
    allproducersmean=apply(producers, 1, "mean")
    allproducerssd=apply(producers, 1, "sd")
  }else{
    alloverallmean=c(alloverallmean, mean(overall))
    alloverallsd=c(alloverallsd, sd(overall))
    allusersmean=cbind(allusersmean,apply(users, 1, "mean"))
    alluserssd=cbind(alluserssd, apply(users, 1, "sd"))
    allproducersmean=cbind(allproducersmean, apply(producers, 1, "mean"))
    allproducerssd=cbind(allproducerssd, apply(producers, 1, "sd"))
  }
}

# compile data
overalldata <- data.frame(
  name=normethods,
  value=(1-alloverallmean)*100,
  sd=alloverallsd*100
)

producersdata <- data.frame(
  name=rep(rownames(allproducersmean), 5),
  type=c(rep(normethods[1], nrow(allproducersmean)), 
         rep(normethods[2], nrow(allproducersmean)), 
         rep(normethods[3], nrow(allproducersmean)),
         rep(normethods[4], nrow(allproducersmean)),
         rep(normethods[5], nrow(allproducersmean))),
  value=c((1-allproducersmean[,1])*100, (1-allproducersmean[,2])*100, 
          (1-allproducersmean[,3])*100, (1-allproducersmean[,4])*100, (1-allproducersmean[,5])*100),
  sd=(c(allproducerssd[,1], allproducerssd[,2], allproducerssd[,3], allproducerssd[,4], allproducerssd[,5]))*100
)

usersdata <- data.frame(
  name=rep(rownames(allusersmean), 5),
  type=c(rep(normethods[1], nrow(allusersmean)), 
         rep(normethods[2], nrow(allusersmean)), 
         rep(normethods[3], nrow(allusersmean)),
         rep(normethods[4], nrow(allusersmean)),
         rep(normethods[5], nrow(allusersmean))),
  value=(1-c(allusersmean[,1], allusersmean[,2], allusersmean[,3], allusersmean[,4], allusersmean[,5]))*100,
  sd=c(alluserssd[,1], alluserssd[,2], alluserssd[,3], alluserssd[,4], alluserssd[,5])*100
)

pdf("ProducersUsersMergeMex2.pdf", width=4,height=4,paper='special') 
par(mar = c(7, 4, 2, 2) + 0.2)
ggplot(producersdata, aes(x = name, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha=0.7) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),                            
                width = 0.2,
                position = position_dodge(0.9))+ 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_fill_manual(values=graphcol)

ggplot(usersdata, aes(x = name, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha=0.7) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),                            
                width = 0.2,
                position = position_dodge(0.9))+ 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_fill_manual(values=graphcol)
ggplot(overalldata) +
  geom_bar(aes(x=name, y=value), stat="identity", fill=graphcol, alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="gray28", alpha=0.9, size=1.3)+
  #coord_cartesian(ylim = c(0, 65)) # ORINOQUIA
  #coord_cartesian(ylim = c(0, 10)) # PUCALLPA
  coord_cartesian(ylim = c(0, 20)) # MEXICO
dev.off()






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



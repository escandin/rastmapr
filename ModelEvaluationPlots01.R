#SETTINGS
# MAC
setwd("/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/ModelEvaluation")
dir()
rdata=list.files('.', pattern='.rdata')
load(rdata[[2]])
ori=binded_or
rm(binded_or)
load(rdata[[3]])
puc=binded_pu
rm(binded_pu)
load(rdata[[1]])
mex=binded
rm(binded)

site=puc
#site=site[site$samplesID!=14,] #remove outlier from Mexico only
attach(site)
pdfFilename="scatterPucV4.pdf"
pdf(pdfFilename, width=4,height=4,paper='special')
for (i in c(5,4,6)){ #gamma, gammaCCA, chisq
  plot(c(0:100), c(0:100), type="l", lty=2, xlab=names(site[i]), ylab="Observed", 
       xlim=c(0,100), ylim=c(0,100))
  # points(User1, site[,i], col="darkgreen")
  # points(User2,site[,i], col="red")
  # abline(lm(site[,i]~User1), col="darkgreen")
  # abline(lm(site[,i]~User2), col="red")
  points(site[,i], User1, col="black")
  points(site[,i], User2, col="grey")
  abline(lm(User1~site[,i]), col="black")
  abline(lm(User2~site[,i]), col="grey")
}
dev.off()

# Calculate fit and bias
# Fit (RMSE)
fit=cbind(rep(NA, 3),rep(NA, 3))
bias=fit
j=1
for (i in c(5,4,6)){ #gamma, gammaCCA, chisq
  # predicted1=predict(lm(site[,i]~User1))
  # predicted2=predict(lm(site[,i]~User2))
  predicted1=predict(lm(User1~site[,i]))
  predicted2=predict(lm(User2~site[,i]))
  fit[j,1]=sqrt(sum((predicted1-User1)^2)/length(predicted1))
  fit[j,2]=sqrt(sum((predicted2-User2[!is.na(User2)])^2)/length(predicted2))
  # bias[i-3,1]=sum(predicted1-User1)/length(predicted1)
  # bias[i-3,2]=sum(predicted2-User2[!is.na(User2)])/length(predicted2)
  # bias[j,1]=sum(predicted1-site[,i])/length(predicted1)
  # bias[j,2]=sum(predicted2-site[,i][!is.na(User2)])/length(predicted2)
  bias[j,1]=sum(abs(User1-site[,i]))/length(User1)
  bias[j,2]=sum(abs(User2-site[,i]))/length(User2)
  j=j+1
}
detach(site)

pdfFilename="FitBiasPucV4.pdf"
pdf(pdfFilename, width=4,height=4,paper='special')
  plot(fit[,1], col="black",ylim=c(0,40), xlab="", ylab="RMSE")
  points(fit[,2], col="grey")
  plot(bias[,1], col="black",ylim=c(0,40), xlab="", ylab="bias")
  points(bias[,2], col="grey")
dev.off()


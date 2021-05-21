library("ggplot2")
library("tidyverse")
library("gridExtra")
library("grid")
setwd("/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/ModelEvaluation")
load("processingTimes.RData")

# convert elements that are in hours to minutes
hr2min=function(x){
  x[x<3]=x[x<3]*60
  x=c(x[1], diff(x))
  return(x)
}

proctimes2=Map(hr2min, proctimes)
# means=unlist(Map('mean', proctimes2))
# sds=unlist(Map(sd, proctimes2))
# iter=unlist(Map(length, proctimes2))
# total=unlist(Map(sum, proctimes2))
#stats=as.data.frame(cbind(location, dist, means, sds, iter, total))

vars  = c('mean','sd','length','sum')
stats=Map(function(x)
  unlist(Map(x, proctimes2)), x= vars)

stats=do.call('cbind',stats)
#stats=Reduce('cbind',tmp)
stats=as.data.frame(stats)

location=factor(c("Orinoquia", "Orinoquia", "Pucallpa", "Pucallpa", "Mexico", "Mexico"),
                levels=c("Orinoquia", "Pucallpa", "Mexico"))
dist=rep(c("Chi-square", "Gamma"), 3)
stats=cbind(location, dist, stats)
names(stats)=c("Location", "Distribution", "Time", "sd", "length", "Total" )

pdf(file='Performance.pdf',
    width = 18, height=9, paper='special')
p<- ggplot(stats, aes(x=Location, y=Time, fill=Distribution)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())+#, show.legend = FALSE) +
  geom_errorbar(aes(ymin=Time-sd, ymax=Time+sd), width=.2,
                position=position_dodge(.9)) +
  labs(y = "Minutes per iteration")+
  labs(x = "")+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme_classic(base_size=24)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(title=NULL))+
  #theme_bw(base_size=18)+
  scale_fill_manual(values=c("#FB9A99", "#8DA0CB"))

q<- ggplot(stats, aes(x=Location, y=Total, fill=Distribution)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  labs(y = "Minutes")+
  labs(x = "")+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme_classic(base_size=24)+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = c(0.8, .9))+
  #theme_bw(base_size=18)+
  scale_fill_manual(values=c("#FB9A99", "#8DA0CB"))
#grid.arrange(p, q, nrow = 1)
#print(p)
#print(q)
grid.arrange(p,q, nrow = 1)
dev.off()
#grid.newpage()
#grid.draw(rbind(ggplotGrob(p), ggplotGrob(q), size = "last"))


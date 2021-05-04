setwd("/Users/tug61163/Documents/PROJECTS/NASAGeo/Manuscripts/ChgNoChgManuscript/ModelEvaluation")
library(tidyverse)
library(plyr)
library(RColorBrewer)

n=8
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
samp=sample(col_vector, n)
pie(rep(1,n), col=samp)

load('acc_change.RData')
load('acc_no_change.RData')

names(accu_no_ch)=c("Accuracy", "Interpreter", "Type", "Parameterization", "Location", "Proportion", "Window")
names(accu_ch)=c("Accuracy", "Interpreter", "Type", "Parameterization", "Location", "Proportion", "Window")
#mapvalues(accu_no_ch$Parameterization, from = c("cc", "ch_sq", "gamma"), to = c("Gamma CCA", "Chisq", "Gammma"))
levels(accu_no_ch$Parameterization)=c("Gamma CCA", "Chisq", "Gamma")
levels(accu_ch$Parameterization)=c("Gamma CCA", "Chisq", "Gamma")
##############################
#Boxplot no-change V2
pdf(file='NoChgAccuracyV2.pdf',
    width = 16, height=9, paper='special')
ggplot(subset(accu_no_ch, Parameterization!="Gamma CCA"), 
       aes(x=Interpreter,y=Accuracy, color=Parameterization))+#, fill=algorithm))+
  geom_boxplot()+
  facet_grid(vars(Location),vars(Type))+
  theme_bw(base_size=24)+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = seq(0, 1, by = .25))+
  scale_colour_manual(values = c("#FB9A99", "#8DA0CB"))+
  ggtitle('')
dev.off()

#Boxplot change
pdf(file='ChgAccuracyV2.pdf',
    width = 16, height=9, paper='special')
ggplot(subset(accu_ch, Parameterization!="Gamma CCA"), aes(x=Interpreter,y=Accuracy, color=Parameterization))+#, fill=algorithm))+
  geom_boxplot()+
  facet_grid(vars(Location),vars(Type))+
  #theme(text=element_text(size=24))+ 
  theme_bw(base_size=24)+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = seq(0, 1, by = .2))+
  scale_colour_manual(values = c("#FB9A99", "#8DA0CB"))+
  ggtitle('')
dev.off()

##############################
#Acc/Area no change
pdf(file='accuracy_area_noch.pdf',
    width = 16, height=9)
ggplot(subset(accu_no_ch, Parameterization!="Gamma CCA"), aes(x=Proportion, y=Accuracy, color=interaction(Interpreter, Parameterization), 
                       shape=Interpreter, group=interaction(Type,Interpreter,Parameterization)))+
  geom_point(alpha = 0.3)+
  facet_grid(vars(Location),vars(Type))+
  geom_smooth( method='gam', alpha=0.1, size=0.5, aes(linetype=Interpreter))+
  scale_colour_manual(values = c("#FB9A99", "#FB9A99", "#8DA0CB", "#8DA0CB"))+
  ylim(0,1)+
  #scale_x_continuous(labels = scales::percent_format(scale = 100))+
  scale_x_continuous(breaks = seq(0, 1, by = .2))+
  scale_y_continuous(breaks = seq(0, 1, by = .2))+
  theme_bw(base_size=20)+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x='Proportion no change')
dev.off()

#Acc/Area  change
pdf(file='accuracy_area_chg.pdf',
    width = 16, height=9)
ggplot(subset(accu_ch, Parameterization!="Gamma CCA"), aes(x=Proportion, y=Accuracy, color=interaction(Interpreter, Parameterization), 
                    shape=Interpreter, group=interaction(Type,Interpreter,Parameterization)))+
  geom_point(alpha = 0.3)+
  facet_grid(vars(Location),vars(Type))+
  geom_smooth( method='gam', alpha=0.1, size=0.5, aes(linetype=Interpreter))+
  scale_colour_manual(values = c("#FB9A99", "#FB9A99", "#8DA0CB", "#8DA0CB"))+
  ylim(0,1)+
  #scale_x_continuous(labels = scales::percent_format(scale = 100))+
  scale_x_continuous(breaks = seq(0, 1, by = .2))+
  scale_y_continuous(breaks = seq(0, 1, by = .2))+
  theme_bw(base_size=20)+
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x='Proportion change')
dev.off()







#Acc/Area no change V2
pdf(file='accuracy_areaV2.pdf',
    width = 16, height=9)
ggplot(accu_no_ch, aes(x=Proportion, y=Accuracy, color=interaction(Interpreter, Parameterization), 
                       shape=Interpreter, group=interaction(Type,Interpreter,Parameterization)))+
  geom_point(alpha = 0.3)+
  facet_grid(vars(Location),vars(Type))+
  geom_smooth( method='gam', alpha=0.1, size=0.5, aes(linetype=Interpreter))+
  ylim(0,1)+
  scale_x_continuous(labels = scales::percent_format(scale = 100))+
  labs(x='No change')

ggplot(accu_ch, aes(x=pixels, y=accuracy, color=algorithm, shape=user, group=interaction(type,user,algorithm)))+
  geom_point(alpha=0.3)+
  #geom_line(aes(linetype=user))+
  facet_grid(vars(location),vars(type))+
  geom_smooth( method='gam', alpha=0.1, size=0.5, aes(linetype=user))+
  ylim(0,1)+
  scale_x_continuous(labels = scales::percent_format(scale = 100))+
  labs(x='Change')
dev.off()









###################################
#Boxplot no-change
pdf(file='accuracy.pdf',
    width = 16, height=9)
ggplot(accu_no_ch, aes(x=algorithm,y=accuracy, color=user))+#, fill=algorithm))+
  geom_boxplot()+
  facet_grid(vars(location),vars(type))+
  ggtitle('Accuracies no-change')

#Boxplot change
ggplot(accu_ch, aes(x=algorithm,y=accuracy, color=user))+#, fill=algorithm))+
  geom_boxplot()+
  facet_grid(vars(location),vars(type))+
  ggtitle('Accuracies change')
dev.off()

#Acc/Area no change
pdf(file='accuracy_area.pdf',
    width = 16, height=9)
ggplot(accu_no_ch, aes(x=pixels, y=accuracy, color=interaction(user, algorithm), shape=user, group=interaction(type,user,algorithm)))+
  geom_point()+
  facet_grid(vars(location),vars(type))+
  geom_smooth( method='gam', alpha=0.1, size=0.5)+
  ylim(0,1)+
  scale_x_continuous(labels = scales::percent_format(scale = 100))+
  labs(x='no-Change (share)')

#Acc/Area  change
ggplot(accu_ch, aes(x=pixels, y=accuracy, color=interaction(user, algorithm), shape=user, group=interaction(type,user,algorithm)))+
  geom_point()+
  facet_grid(vars(location),vars(type))+
  geom_smooth( method='gam', alpha=0.1, size=0.5)+
  ylim(0,1)+
  scale_x_continuous(labels = scales::percent_format(scale = 100))+
  labs(x='Change (share)')
dev.off()


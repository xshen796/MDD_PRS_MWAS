setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')
library(dplyr)
library(ggplot2)
library(ggpubr)
options(bitmapType='cairo')

load('result/MDDprs_assoc_mPrediction/MRS_MDDprs_assoc_noassoc_lifestyle_MDD.RData')
targetresult=result.lifestyle
targetresult$factor=gsub('_predictor','',targetresult$factor)
targetresult$factor=paste0('MRS: ',targetresult$factor)
targetresult$factor=gsub('_',' ',targetresult$factor)
targetresult$ord=1:nrow(targetresult)

# MDD bar graph -------------------------------------------------------------------------------------------

fig.dat=filter(targetresult,dependent=='MDD_status')
fig.dat$ord=1:nrow(fig.dat)

fig.beta<-ggplot(fig.dat, aes(y=beta, x=factor, fill=factor)) +
  geom_bar(stat="identity",position=position_dodge(),width=0.5)+
  geom_errorbar(aes(ymin=beta-std, ymax=beta+std), width=.1,
                 position=position_dodge()) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position = "none")+
  coord_flip()+
  ylab('Log OR')+
  xlab('Independent variable')
  
fig.p<-ggplot(fig.dat, aes(y=-log10(p.value), x=factor, fill=factor)) +
  geom_bar(stat="identity",position=position_dodge(),width=0.5)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  coord_flip()+
  ylab('-log10(p)')+
  xlab('Independent variable')  
  
fig.r2<-ggplot(fig.dat, aes(y=R2*100, x=factor, fill=factor)) +
  geom_bar(stat="identity",position=position_dodge(),width=0.5)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  coord_flip()+
  ylab('R2(%)')+
  xlab('Independent variable')  

fig.total=ggarrange(fig.beta,fig.p,fig.r2,ncol=3,widths=c(1.5,1,1))


ggsave(plot = fig.total,filename = 'Figs/compare_MRS_MDD.png',
       width = 32, height = 8,
       units = 'cm', dpi = 200)
  

# lifestyle bar graph -------------------------------------------------------------------------------------------
fig.dat=filter(targetresult,dependent!='MDD_status')
fig.dat=filter(fig.dat,dependent!='bmi')
fig.dat$factor=gsub('_',' ',fig.dat$factor)
fig.dat$ord=1:nrow(fig.dat)

fig.beta<-ggplot(fig.dat, aes(y=beta, x=dependent, fill=factor)) +
  geom_bar(stat="identity",position=position_dodge(.8),width=.8)+
  geom_errorbar(aes(ymin=beta-std, ymax=beta+std), width=.1,colour='grey',
                 position=position_dodge(.7)) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #theme(legend.position = "none")+
  guides(fill=guide_legend(title="MRS type"))+
  coord_flip()+
  ylab('Log OR')+
  xlab('Independent variable')
  
fig.p<-ggplot(fig.dat, aes(y=-log10(p.value), x=dependent, fill=factor)) +
  geom_bar(stat="identity",position=position_dodge(.8),width=.8)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(#legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
        guides(fill=guide_legend(title="MRS type"))+
  coord_flip()+
  ylab('-log10(p)')+
  xlab('Independent variable')  
  
fig.r2<-ggplot(fig.dat, aes(y=R2*100, x=dependent, fill=factor)) +
  geom_bar(stat="identity",position=position_dodge(.8),width=.8)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(#legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
        guides(fill=guide_legend(title="MRS type"))+
  coord_flip()+
  ylab('R2(%)')+
  xlab('Independent variable')  

fig.total=ggarrange(fig.beta,fig.p,fig.r2,ncol=3,widths=c(1.5,1,1),common.legend=T,legend='bottom')


ggsave(plot = fig.total,filename = 'Figs/compare_MRS_lifestyle.png',
       width = 33, height = 12,
       units = 'cm', dpi = 200)


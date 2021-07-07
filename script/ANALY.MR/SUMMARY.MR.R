setwd('/sdata/images/projects/UKBIOBANK/users/Shen/ii.PGRS/IDP_4thrls_20k/pheWAS')
library(scales)
library(dplyr)
### correct for p values for MD to IM MR results:

# Correct p values MR -----------------------------------------------------


res.MD_to_IM=read.table('result//ii.MR//MDD_to_IM//outcomes.MRsummary',
                        sep='\t',header=T,stringsAsFactors=F)
colnames(res.MD_to_IM)[11:12]=c('Egger.intercept_se','Egger.intercept_pval')
res.MD_to_IM$MRpresso_Pval=gsub('<','',res.MD_to_IM$MRpresso_Pval)
res.MD_to_IM$MRpresso_Pval=as.numeric(res.MD_to_IM$MRpresso_Pval)

table.tocorrect=res.MD_to_IM
table.tocorrect=table.tocorrect[order(table.tocorrect$method),]

ls.MRmethod=as.character(unique(table.tocorrect$method))
ls.dmri=c('^gFA_','^gMD_','^gICVF_','^tFA_','^tMD_','tICVF_','^rsamp')
ls.total=data.frame(rep(ls.MRmethod,each=length(ls.dmri)),ls.dmri)

for (m in 1:nrow(ls.total)){
      MRmethod=ls.total[m,1]
      dmri=ls.total[m,2]
      tmp.dat=table.tocorrect[grep(MRmethod,table.tocorrect$method),]
      tmp.dat=tmp.dat[grep(dmri,tmp.dat$outcome),]
      if (m==1){
            p.new=p.adjust(tmp.dat$pval,method='fdr')
            p.intercept=p.adjust(tmp.dat$Egger.intercept_pval,method='fdr')
            p.Q=p.adjust(tmp.dat$Q_pval,method='fdr')
            p.presso=p.adjust(tmp.dat$MRpresso_Pval,method='fdr')
      }else{
            p.new=c(p.new,p.adjust(tmp.dat$pval,method='fdr'))
            p.intercept=c(p.intercept,p.adjust(tmp.dat$Egger.intercept_pval,method='fdr'))
            p.Q=c(p.Q,p.adjust(tmp.dat$Q_pval,method='fdr')) 
            p.presso=c(p.presso,p.adjust(tmp.dat$MRpresso_Pval,method='fdr'))
      }
}
table.tocorrect$p.corrected=p.new
table.tocorrect$p.intercept=p.intercept
table.tocorrect$p.Q=p.Q
table.tocorrect$p.presso=p.presso

table.tocorrect=table.tocorrect[order(table.tocorrect$outcome),]

new.table.MD_to_IM=table.tocorrect[,c('exposure','outcome','method','nsnp','b','se','pval','p.corrected',
                                      'egger_intercept','Egger.intercept_se','Egger.intercept_pval','p.intercept',
                                      'Q','Q_df','Q_pval','p.Q','MRpresso_RSSobs','MRpresso_Pval','p.presso')]
colnames(new.table.MD_to_IM)=c('Exposure','Outcome','MR method','Nsnp','Beta','SE','MR.pval','MR.pFDR',
                               'Egger_Intercept','Egger_Intercept.SE','Egger_Intercept.pval','Egger_intercept.pFDR',
                               'Q','Q.df','Q.pval','Q_pFDR','MRpresso.RSS','MRpresso.pval','MRpresso.pFDR')
saveRDS(new.table.MD_to_IM,file='result/ii.MR//MR_summary.MD_to_IM.rds')

rownames(new.table.MD_to_IM)=NULL
new.table.MD_to_IM$Exposure='Depression'
new.table.MD_to_IM[,c(5,6,10,13,17)]=round(new.table.MD_to_IM[,c(5,6,10,13,17)],digits = 3)
new.table.MD_to_IM[,grep('pval|pFDR|Egger_Intercept',colnames(new.table.MD_to_IM))]=
      apply(new.table.MD_to_IM[,grep('pval|pFDR|Egger_Intercept',colnames(new.table.MD_to_IM))],FUN=scientific,digits=3,MARGIN = 2)
new.table.MD_to_IM$Outcome=gsub('^amp','Amplitude in ',new.table.MD_to_IM$Outcome)
new.table.MD_to_IM$Outcome=gsub('^FA','FA in ',new.table.MD_to_IM$Outcome)
new.table.MD_to_IM$Outcome=gsub('^gFA','gFA-',new.table.MD_to_IM$Outcome)
new.table.MD_to_IM$Outcome=gsub('^MD','MD in ',new.table.MD_to_IM$Outcome)
new.table.MD_to_IM$Outcome=gsub('^gMD','gMD-',new.table.MD_to_IM$Outcome)

save(new.table.MD_to_IM,file='result/ii.MR/MR.MD_to_IM.RData')



### correct for p values for IM to MD MR results:

res.IM_to_MD=read.table('result//ii.MR//IM_to_MDD//outcomes.MRsummary',
                        sep='\t',header=T,stringsAsFactors=F)
res.IM_to_MD=filter(res.IM_to_MD,exposure!='tFA_IFOF')
#res.IM_to_MD=filter(res.IM_to_MD,method!='MR-Presso_outlier-corrected')
colnames(res.IM_to_MD)[11:12]=c('Egger.intercept_se','Egger.intercept_pval')
res.IM_to_MD$MRpresso_Pval=gsub('<','',res.IM_to_MD$MRpresso_Pval)
res.IM_to_MD$MRpresso_Pval=as.numeric(res.IM_to_MD$MRpresso_Pval)

table.tocorrect=res.IM_to_MD
table.tocorrect=table.tocorrect[order(table.tocorrect$method),]

ls.MRmethod=as.character(unique(table.tocorrect$method))
ls.dmri=c('^gFA_','^gMD_','^gICVF_','^tFA_','^tMD_','tICVF_','^rsamp')
ls.total=data.frame(rep(ls.MRmethod,each=length(ls.dmri)),ls.dmri)

for (m in 1:nrow(ls.total)){
      MRmethod=ls.total[m,1]
      dmri=ls.total[m,2]
      tmp.dat=table.tocorrect[grep(MRmethod,table.tocorrect$method),]
      tmp.dat=tmp.dat[grep(dmri,tmp.dat$exposure),]
      if (m==1){
            p.new=p.adjust(tmp.dat$pval,method='fdr')
            p.intercept=p.adjust(tmp.dat$Egger.intercept_pval,method='fdr')
            p.Q=p.adjust(tmp.dat$Q_pval,method='fdr')
            p.presso=p.adjust(tmp.dat$MRpresso_Pval,method='fdr')
      }else{
            p.new=c(p.new,p.adjust(tmp.dat$pval,method='fdr'))
            p.intercept=c(p.intercept,p.adjust(tmp.dat$Egger.intercept_pval,method='fdr'))
            p.Q=c(p.Q,p.adjust(tmp.dat$Q_pval,method='fdr')) 
            p.presso=c(p.presso,p.adjust(tmp.dat$MRpresso_Pval,method='fdr'))
      }
}
table.tocorrect$p.corrected=p.new
table.tocorrect$p.intercept=p.intercept
table.tocorrect$p.Q=p.Q
table.tocorrect$p.presso=p.presso

table.tocorrect=table.tocorrect[order(table.tocorrect$exposure),]

new.table.IM_to_MD=table.tocorrect[,c('exposure','outcome','method','nsnp','b','se','pval','p.corrected',
                                      'egger_intercept','Egger.intercept_se','Egger.intercept_pval','p.intercept',
                                      'Q','Q_df','Q_pval','p.Q','MRpresso_RSSobs','MRpresso_Pval','p.presso')]
colnames(new.table.IM_to_MD)=c('Exposure','Outcome','MR method','Nsnp','Beta','SE','MR.pval','MR.pFDR',
                               'Egger_Intercept','Egger_Intercept.SE','Egger_Intercept.pval','Egger_intercept.pFDR',
                               'Q','Q.df','Q.pval','Q_pFDR','MRpresso.RSS','MRpresso.pval','MRpresso.pFDR')
rownames(new.table.IM_to_MD)=NULL
new.table.IM_to_MD$Outcome='Depression'
saveRDS(new.table.IM_to_MD,file='result/ii.MR//MR_summary.IM_to_MD.rds')
new.table.IM_to_MD[,c(5,6,10,13,17)]=round(new.table.IM_to_MD[,c(5,6,10,13,17)],digits = 3)
new.table.IM_to_MD[,grep('pval|pFDR',colnames(new.table.IM_to_MD))]=
      apply(new.table.IM_to_MD[,grep('pval|pFDR',colnames(new.table.IM_to_MD))],FUN=scientific,digits=3,MARGIN = 2)

new.table.IM_to_MD$Exposure=gsub('^amp','Amplitude in ',new.table.IM_to_MD$Exposure)
new.table.IM_to_MD$Exposure=gsub('^FA','FA in ',new.table.IM_to_MD$Exposure)
new.table.IM_to_MD$Exposure=gsub('^gFA','gFA-',new.table.IM_to_MD$Exposure)
new.table.IM_to_MD$Exposure=gsub('^MD','MD in ',new.table.IM_to_MD$Exposure)
new.table.IM_to_MD$Exposure=gsub('^gMD','gMD-',new.table.IM_to_MD$Exposure)


save(new.table.IM_to_MD,file='result/ii.MR/MR.IM_to_MD.RData')


### make a table for a graph
table.forgraph=cbind(new.table.IM_to_MD[new.table.IM_to_MD$`MR method`!='MR-Presso_outlier-corrected',
                                        c('Exposure','MR method','MR.pval','MR.pFDR')],
                     new.table.MD_to_IM[new.table.MD_to_IM$`MR method`!='MR-Presso_outlier-corrected',
                                        c('MR.pval','MR.pFDR')])
rownames(table.forgraph)=NULL
colnames(table.forgraph)=c('IM.pheno','method','IMexposure_pval','IMexposure_pcorr','IMoutcome_pval','IMoutcome_pcorr')

# select the ones with at least sig corr p in two methods 
ls.pheno=table(table.forgraph[table.forgraph$IMoutcome_pcorr<0.05,'IM.pheno'])
ls.pheno=names(ls.pheno[ls.pheno>=2])
ls.pheno=c(ls.pheno,'MD_ATR')

table.p_plot=table.forgraph[grep(paste0(ls.pheno,collapse = '|'),table.forgraph$IM.pheno),]

write.table(table.p_plot,file='result/ii.MR/table.p_plot.txt',sep='\t',quote=F,row.names = F,col.names = T)





### make tables that report the bi-directional results (supplementary materials)
Table.report.IM_to_MDD=new.table.IM_to_MD[,c('Exposure','Outcome','MR method','Nsnp',
                                             'Beta','SE','MR.pval','MR.pFDR',
                                             'Egger_Intercept','Egger_Intercept.SE','Egger_intercept.pFDR',
                                             'Q','Q.df','Q_pFDR','MRpresso.pFDR')]
colnames(Table.report.IM_to_MDD)=c('Exposure','Outcome','MR method','Nsnp','Beta','SE','P','Pcorr',
                                   'Egger intercept','SE Egger intercept','Pcorr Egger intercept',
                                   'Q','Q df','Pcorr Q','Pcorr MRpresso')
save(Table.report.IM_to_MDD,file='result/ii.MR/table_report.IM_to_MDD.RData')



Table.report.MDD_to_IM=new.table.MD_to_IM[,c('Exposure','Outcome','MR method','Nsnp',
                                             'Beta','SE','MR.pval','MR.pFDR',
                                             'Egger_Intercept','Egger_Intercept.SE','Egger_intercept.pFDR',
                                             'Q','Q.df','Q_pFDR','MRpresso.pFDR')]
colnames(Table.report.MDD_to_IM)=c('Exposure','Outcome','MR method','Nsnp','Beta','SE','P','Pcorr',
                                   'Egger intercept','SE Egger intercept','Pcorr Egger intercept',
                                   'Q','Q df','Pcorr Q','Pcorr MRpresso')
save(Table.report.MDD_to_IM,file='result/ii.MR//table_report.MDD_to_IM.RData')


### make a table reporting the bi-directional sig results

# table.tmp1=new.table.IM_to_MD[,c('Exposure','MR method','MR.pval','MR.pFDR')]
# table.tmp2= new.table.MD_to_IM[,c('Outcome','MR.pval','MR.pFDR')]
# 
# table.forgraph=merge(table.tmp1,table.tmp2,by.x='Exposure',by.y='Outcome')
# rownames(table.forgraph)=NULL
# colnames(table.forgraph)=c('IM.pheno','method','IMexposure_pval','IMexposure_pcorr','IMoutcome_pval','IMoutcome_pcorr')
# # select the ones with at least sig corr p in two methods 
# ls.pheno=table(table.forgraph[table.forgraph$IMoutcome_pcorr<0.05,'IM.pheno'])
# ls.pheno=names(ls.pheno[ls.pheno>=2])
# ls.pheno=c(ls.pheno,'MD_ATR')
# 
# Table.report.IM_to_MDD=Table.report.IM_to_MDD[grep(paste0(ls.pheno,collapse = '|'),Table.report.IM_to_MDD$Exposure),]
# Table.report.MDD_to_IM=Table.report.MDD_to_IM[grep(paste0(ls.pheno,collapse = '|'),Table.report.MDD_to_IM$Outcome),]
#   
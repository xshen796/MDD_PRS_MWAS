setwd('/sdata/images/projects/UKBIOBANK/users/Shen/ii.PGRS/IDP_4thrls_20k/pheWAS')
library(scales)
library(dplyr)
### correct for p values for MD to IM MR results:

# Correct p values MR -----------------------------------------------------


res.neuroSNPremoved=read.table('result//ii.MR//sensitivity_analysis/remove_neuroticism_SNPs/outcomes.MRsummary',
                        sep='\t',header=T,stringsAsFactors=F)
colnames(res.neuroSNPremoved)[11:12]=c('Egger.intercept_se','Egger.intercept_pval')
res.neuroSNPremoved$MRpresso_Pval=gsub('<','',res.neuroSNPremoved$MRpresso_Pval)
res.neuroSNPremoved$MRpresso_Pval=as.numeric(res.neuroSNPremoved$MRpresso_Pval)

table.tocorrect=res.neuroSNPremoved
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

save(new.table.MD_to_IM,file='result/ii.MR/MR.sensitivity_neuroticism_SNP_removed.RData')


Table.report.MDD_to_IM=new.table.MD_to_IM[,c('Exposure','Outcome','MR method','Nsnp',
                                             'Beta','SE','MR.pval','MR.pFDR',
                                             'Egger_Intercept','Egger_Intercept.SE','Egger_intercept.pFDR',
                                             'Q','Q.df','Q_pFDR','MRpresso.pFDR')]
colnames(Table.report.MDD_to_IM)=c('Exposure','Outcome','MR method','Nsnp','Beta','SE','P','Pcorr',
                                   'Egger intercept','SE Egger intercept','Pcorr Egger intercept',
                                   'Q','Q df','Pcorr Q','Pcorr MRpresso')
save(Table.report.MDD_to_IM,file='result/ii.MR//table_report.sensitivity_neuroticism_SNP_removed.RData')


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
library(dplyr)
library(data.table)
library(ggplot2)
setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD')

mdd.ewas=fread('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/EWAS_MDDprs_fromKathryn/mdd_waveincluded_covs.toptable.txt',stringsAsFactors=T,fill=T)
mdd.ewas=mdd.ewas[,c('ID','CHR','P.Value','adj.P.Val','beta','se','N')]
colnames(mdd.ewas)[1]='CpG'

mddprs.ewas=fread('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/EWAS_MDDprs_fromKathryn/meta_wave1_wave3/mdd-pgrs-snp_ch/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.metal.out1.tbl',stringsAsFactors=T,fill=T)
mddprs.ewas=mddprs.ewas[,c('MarkerName','Effect','StdErr','P-value','Direction')]
colnames(mddprs.ewas)=c('CpG','beta.meta_prs','se.meta_prs','P.Value.meta_prs','Direction.meta_prs')


summstats2=merge(mddprs.ewas,mdd.ewas,by='CpG',all.x=T)
r.beta=cor(summstats2[,c('beta','beta.meta_prs')],use="complete.obs")

fig.dat=summstats2[,c('beta','beta.meta_prs')]
fig.dat=fig.dat[complete.cases(fig.dat),]

fig.beta=ggplot(fig.dat, aes(x=beta,y=beta.meta_prs)) +
  geom_point(alpha=0.7,size=1)+
  theme(legend.title = element_blank())+
  annotate(geom="text",label=paste0('r = ',r.beta), x=1, y=3)+
  xlab('MDD EWAS')+
  ylab('MDD PRS EWAS')
  geom_smooth(method=lm,size=0.3)


tiff("Figs/MDDprsEWAS_MDDewas_beta_correlation.tiff", width = 12, height = 12, units = 'in', res = 200)
fig.beta # Make plot
dev.off()


p.samedir=(sum(mddprs.ewas$Direction.meta_prs=='++',na.rm=T)+sum(mddprs.ewas$Direction.meta_prs=='--',na.rm=T))/nrow(mddprs.ewas)
mddprs.ewas.sig=filter(mddprs.ewas,P.Value.meta_prs<5*10^-8)
p.samedir=(sum(mddprs.ewas.sig$Direction.meta_prs=='++',na.rm=T)+sum(mddprs.ewas.sig$Direction.meta_prs=='--',na.rm=T))/nrow(mddprs.ewas.sig)
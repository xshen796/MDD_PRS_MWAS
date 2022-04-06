setwd("/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/")
library(dplyr)
library(pbapply)
library(data.table)


# Create phenotype file  ---------------------------------------------------------------------------------
# p+T method:

ls.f=list.files('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/',pattern = '.profile',full.names = T)

prs = as.list(ls.f) %>%
  lapply(.,fread,stringsAsFactors=F,header=T) %>%
  lapply(.,select,IID,starts_with('SCORE')) %>%
  Reduce(function(x,y) left_join(x,y,by='IID'),.) %>%
  mutate_at(vars(matches('SCORE')), ~scale(.)) %>% 
  as.data.frame %>% 
  rename(id=IID)

colnames(prs) = colnames(prs) %>% gsub('SCORE','pT',.)

write.table(prs,file='data/GS_phenotype/MDDprs_MWAS_input.txt',col.names=T,row.names=F,sep=',',quote=F)  

ls.prs = colnames(prs) %>% .[!. %in% 'id']
write.table(ls.prs,file='data/GS_phenotype/ls_MDDprs.txt',col.names=F,row.names=F,sep='\n',quote=F)  

# unrelated, p+T method:
w1.id = covs=read.table('result/EWAS_MDDprs_fromKathryn/wave1-pgrs-mdd/mdd-pgrs-wave1-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',') %>% 
  .$Sample_Name
kinship = read_tsv('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/GSrelated.kin_forSHEN') %>% 
  .[(.$ID1 %in% w1.id) & (.$ID2 %in% w1.id),]
id.unrelated.kinship =  kinship %>% filter(Kinship < 0.05) %>% c(.$ID1,.$ID2) %>% unique 
id.unrelated = w1.id %>% .[!. %in% c(kinship$ID1,kinship$ID2)] %>% c(.,id.unrelated.kinship) %>% unique

prs.unrelated.w1 = prs %>% .[.$id %in% id.unrelated,]

write.table(prs.unrelated.w1,file='data/GS_phenotype/MDDprs_MWAS_input_unrelated.txt',col.names=T,row.names=F,sep=',',quote=F)  


# Create covariate file  ---------------------------------------------------------------------------------

covs=read.table('result/EWAS_MDDprs_fromKathryn/wave1-pgrs-mdd/mdd-pgrs-wave1-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',')
covs=covs[,c('Sample_Name','sex','age','ever_smoke','pack_years')]
covs$Batch=as.factor(covs$Batch)
covs$ever_smoke=as.factor(covs$ever_smoke)
colnames(covs)[1]='id'

saveRDS(covs,file='data/EWAS_MDDprs_Shen/covs_RosieData_wave1.rds',version=2)

covs=read.table('data/EWAS_MDDprs_fromKathryn/wave3/mdd-pgrs-wave3-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',')
covs=covs[,c('Sample_Name','sex','age','ever_smoke','pack_years')]
covs$Batch=as.factor(covs$Batch)
covs$ever_smoke=as.factor(covs$ever_smoke)
colnames(covs)[1]='id'

saveRDS(covs,file='data/EWAS_MDDprs_Shen/covs_RosieData_wave3.rds',version=2)


# phenotype file for PRS without MHC region -------------------------------

prs.nomhc = fread('/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC.all_score',sep = ' ')

prs.nomhc = prs.nomhc %>% 
  select(id=IID,Pt_5e.08=`Pt_5e-08`)

write.table(prs.nomhc,file = '/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC_forEWAS.txt',
            col.names=T,row.names=F,sep=',',quote=F)


# phenotype file for PRS created using the top indenpendent hits ----------

prs.gwsigSNP = fread('/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_5e_08.all_score',sep = ' ')

prs.gwsigSNP = prs.gwsigSNP %>% 
  select(id=IID,Pt_5e.08=`Pt_5e-08`)

write.table(prs.gwsigSNP,file = '/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_gwsigSNP_forEWAS.txt',
            col.names=T,row.names=F,sep=',',quote=F)

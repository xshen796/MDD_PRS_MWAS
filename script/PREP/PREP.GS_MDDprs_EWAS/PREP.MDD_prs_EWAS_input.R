setwd("/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/")
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


# Create covariate file  ---------------------------------------------------------------------------------
#PC=read.table('data/HM3mds.mds',header=T,stringsAsFactors=F)
#PC=PC[!duplicated(PC$IID),]
#PC.tomerge=PC[,c('IID',paste0('C',1:10))]
#colnames(PC.tomerge)[1]='id'

covs=read.table('data/EWAS_MDDprs_fromKathryn/wave1/mdd-pgrs-wave1-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',')
covs=covs[,c('Sample_Name','sex','age','ever_smoke','pack_years','Batch')]
covs$Batch=as.factor(covs$Batch)
covs$ever_smoke=as.factor(covs$ever_smoke)
colnames(covs)[1]='id'

saveRDS(covs,file='data/EWAS_MDDprs_Shen/covs_RosieData_wave1.rds')

covs=read.table('data/EWAS_MDDprs_fromKathryn/wave3/mdd-pgrs-wave3-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',')
covs=covs[,c('Sample_Name','sex','age','ever_smoke','pack_years','Batch')]
covs$Batch=as.factor(covs$Batch)
covs$ever_smoke=as.factor(covs$ever_smoke)
colnames(covs)[1]='id'

saveRDS(covs,file='data/EWAS_MDDprs_Shen/covs_RosieData_wave3.rds')


# meta ewas for all SNPs 102  -----------------------------------------------------------------------------

ls.SNP=dir(path = "/gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/EWAS_MDDprs_fromKathryn/meta_wave1_wave3/mdd-pgrs-snp_ch/", pattern = 'mdd_SNPCH_rs', all.files = FALSE,full.names = F)

ls.folder=dir(path = "/gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/EWAS_MDDprs_fromKathryn/meta_wave1_wave3/mdd-pgrs-snp_ch/", pattern = 'mdd_SNPCH_rs', all.files = FALSE,full.names = T)

ls.file=paste0(ls.folder,'/',ls.SNP,'.metal.out1.tbl')

comm_metal=paste0('PROCESS ',ls.file)

write.table(comm_metal,sep='\n',file='data/EWAS_MDDprs_Shen/MDDprs_ewas_102SNP/files.to.process.metal.txt',quote=F,row.names=F,col.names=F)
# Manually copy the file list to the file: /gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/EWAS_MDDprs_Shen/MDDprs_ewas_102SNP/metal_input_ewas_SNPs.txt
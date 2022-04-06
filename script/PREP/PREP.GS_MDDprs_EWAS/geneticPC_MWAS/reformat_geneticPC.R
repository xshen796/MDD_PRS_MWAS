library(dplyr)
library(readr)

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS')

# load data ---------------------------------------------------------------

pc.wave1 = read_delim('data/GS_genetic_meth/Genetic_PC_w1.eigenvec',col_names=F) %>% 
  rename(FID=X1,id=X2)
colnames(pc.wave1)[3:ncol(pc.wave1)]=paste0('C',1:20) 


pc.wave3 = read_delim('data/GS_genetic_meth/Genetic_PC_w3.eigenvec',col_names=F) %>% 
  rename(FID=X1,id=X2)
colnames(pc.wave3)[3:ncol(pc.wave1)]=paste0('C',1:20)

PC.tomerge = rbind(pc.wave1,pc.wave3) %>% 
  select(id,paste0('C',1:10))


# Additional covariates with genetic PCs ----------------------------------


covs = readRDS('data/GS_phenotype/covs_RosieData_wave1.rds') %>% 
  select(-Batch) %>% 
  left_join(.,PC.tomerge,by='id')

saveRDS(covs,'data/GS_phenotype/covs_geneticPC_RosieData_wave1.rds',version=2)

covs = readRDS('data/GS_phenotype/covs_RosieData_wave3.rds') %>% 
  select(-Batch) %>% 
  left_join(.,PC.tomerge,by='id')

saveRDS(covs,'data/GS_phenotype/covs_geneticPC_RosieData_wave3.rds',version = 2)
  

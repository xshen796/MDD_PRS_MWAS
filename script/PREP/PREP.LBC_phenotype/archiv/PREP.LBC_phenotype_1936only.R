
# Basic settings ----------------------------------------------------------

setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

library(dplyr)
library(data.table)


# Load data ---------------------------------------------------------------

# Phenotype
LBC1921.wave1.pheno=
  read.csv("/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/LBC_phenotype/LBC1921_XueyiShen_EWAS_DepressionDisorderPGR_AM_24AUG2020.csv",
           stringsAsFactors=F,header=T)
LBC1936.wave1.pheno=
  read.csv("/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/LBC_phenotype/LBC1936_XueyiShen_EWAS_DepressionDisorderPGR_AM_24AUG2020.csv",
           stringsAsFactors=F,header=T)

# PRS
LBC1921.PRS=read.table('data/LBC_phenotype/MDDprs/MDDprs.LBC1921.imputed.all_score',header=T,stringsAsFactors = F)
LBC1936.PRS=read.table('data/LBC_phenotype/MDDprs/MDDprs.LBC1936.imputed.all_score',header=T,stringsAsFactors = F)
LBC1921.PRS.tophits=read.table('data/LBC_phenotype/MDDprs/MDDprs.LBC1921.imputed.topSNP.all_score',header=T,stringsAsFactors = F)
LBC1936.PRS.tophits=read.table('data/LBC_phenotype/MDDprs/MDDprs.LBC1936.imputed.topSNP.all_score',header=T,stringsAsFactors = F)

# tech-covars, methylation data
LBC_w1_cellcount=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC_w1_estWBCs_for_Shen_26Aug2020.rds') %>% mutate_if(is.factor, as.character)
LBC21.tech=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_w1_DNAm_target_for_Shen_26Aug2020.rds') %>% mutate_if(is.factor, as.character)
LBC36.tech=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC36_w1_DNAm_target_for_Shen_26Aug2020.rds') %>% mutate_if(is.factor, as.character)



# Genetic PC1-20
PC=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC_both/LBC21_n_36_PCs.eigenvec',
              header=F,stringsAsFactors = F)

# Create a dataset for covariates -----------------------------------------

#  Phenotype: age, sex, smoking category
colnames(LBC1921.wave1.pheno)[1]='ID'
colnames(LBC1936.wave1.pheno)[1]='ID'

LBC21.pheno.tomerge = LBC1921.wave1.pheno[,c('ID','smoking')]
LBC36.pheno.tomerge = LBC1936.wave1.pheno[,c('ID','smokcat_w1')]

colnames(LBC36.pheno.tomerge) <- colnames(LBC21.pheno.tomerge)

LBC.pheno.covar = rbind(LBC21.pheno.tomerge,LBC36.pheno.tomerge)

# Tech covars
LBC.meth.tech.covar = rbind(LBC21.tech,LBC36.tech) %>%
  select(Basename,ID,age,sex,plate) %>%
  merge(.,LBC_w1_cellcount[,!colnames(LBC_w1_cellcount) %in% colnames(.)[2:ncol(.)]],by='Basename')

# PCs
PC.tokeep = 
  PC[(PC$V1 %in% LBC.pheno.covar$ID) |(PC$V1 %in% LBC.meth.tech.covar$ID),] %>%
  .[,c(1,3:12)] 
colnames(PC.tokeep)=c('ID',paste0('PC',1:10))


# Harmonise IDs (keep methylation study ID)

input.covar = merge(LBC.meth.tech.covar,LBC.pheno.covar,by='ID') %>%
  select( -matches("^Basename$|WAVE|cohort")) %>%
  mutate(smoking=as.factor(smoking),ID=as.character(ID),sex=as.factor(sex))
colnames(input.covar)[1]='id'

input.PC = merge(PC.tokeep,LBC_w1_cellcount,by='ID') %>%
  select(Basename, matches('^PC'))
input.PC$Basename=paste0('X',input.PC$Basename)
colnames(input.PC)[1]='Sample_Sentrix_ID'

saveRDS(input.covar,file='data/LBC_phenotype/ewas_covar.rds',version = 2)
saveRDS(input.PC,file='data/LBC_phenotype/ewas_gPC1_20.rds',version = 2)


# Create target file for phenotype of interest ----------------------------

# pT PRS
LBC.PRS.pt = rbind(LBC1921.PRS,LBC1936.PRS) %>%
  select(IID,matches('Pt'))
LBC.PRS.tophits = rbind(LBC1921.PRS.tophits,LBC1936.PRS.tophits) %>%
  select(IID,matches('Pt'))
colnames(LBC.PRS.tophits)[2]='Pt_gwsig'

LBC.PRS=merge(LBC.PRS.pt,LBC.PRS.tophits,by='IID') %>%
  merge(.,LBC_w1_cellcount,by.x='IID',by.y='ID') %>%
  select(IID, matches('^Pt'))
colnames(LBC.PRS)[1]='id'

LBC.PRS[,2:ncol(LBC.PRS)]=scale(LBC.PRS[,2:ncol(LBC.PRS)]) 

LBC.PRS = LBC.PRS[LBC.PRS$id %in% LBC1936.wave1.pheno$ID,]

write.csv(LBC.PRS,file=paste0('data/LBC_phenotype/LBC_MDDprs.csv'),
          quote = F,row.names = F)



# Linkage file ------------------------------------------------------------

linkage.file=rbind(LBC21.tech,LBC36.tech) %>%
  select(ID,array,pos,Basename)
colnames(linkage.file)=c('id','Sentrix_ID','Sentrix_Position','Sample_Sentrix_ID')
linkage.file$Sample_Sentrix_ID=paste0('X',linkage.file$Sample_Sentrix_ID)

saveRDS(linkage.file,file='data/LBC_phenotype/sentrix.rds',version = 2)


# Prep ewas command -------------------------------------------------------
ls.pheno = colnames(LBC.PRS)[2:ncol(LBC.PRS)]

cmmd=paste0('lbc_ewas --pdata /gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/LBC_MDDprs.csv --pheno ',
            ls.pheno,
            ' --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_LBC/MDDprs_', ls.pheno,'_ewas')

write.table(cmmd,file='script/ANALY.MDDprs_EWAS/LBC/ewas_cmmd.sh',
            sep='\n\n',row.names=F,quote = F,col.names = F)

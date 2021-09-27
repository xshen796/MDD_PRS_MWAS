args <- commandArgs(trailingOnly = T)

library(dplyr)
library(pbapply)


# Combine two waves and create M-values

#setwd('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/')
# LBC21_w1=readRDS('LBC21_w1_DNAm_for_Shen_26Aug2020.rds')
# LBC36_w1=readRDS('LBC36_w1_DNAm_for_Shen_26Aug2020.rds')
# 
# LBC21_LBC36_w1_DNAm=data.frame(LBC21_w1,LBC36_w1)
# 
# saveRDS(LBC21_LBC36_w1_DNAm,file='LBC21_LBC36_w1_DNAm.rds')

# mval_LBC21_LBC36_w1_DNAm=log2(LBC21_LBC36_w1_DNAm/(1-LBC21_LBC36_w1_DNAm))
# saveRDS(mval_LBC21_LBC36_w1_DNAm,file='mval_LBC21_LBC36_w1_DNAm.rds')


# Calculate PCs in scratch space
# Copy data to scratch space
# system('cp -r /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr $myscratch')

# Basic settings
#setwd('/exports/eddie/scratch/xshen33/LBC21_n_36_w1_DNAm_bychr')

# Load data
sentrix.id=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/LBC_phenotype/sentrix.rds') %>%
  .[,c(1,4)]
covariates=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/LBC_phenotype/ewas_covar.rds') %>%
  .[,1:(ncol(.)-1)]

# Define functions
residualise_dat <- function(col.input,dat) {
  coln=as.character(col.input[1])
  cov=as.character(col.input[2])
  colnames(dat)[1]='ID'
  
  eq=paste0(coln,'~',cov)
  fit=glm(as.formula(eq),dat=dat)
  new.dat=resid(fit)
  ls.cov=strsplit(cov,'\\+') %>% unlist
  new.dat=data.frame(ID=dat[complete.cases(dat[,c(coln,ls.cov)]),1],pheno=new.dat,stringsAsFactors=F)
  
  original.ID=data.frame(ID=dat$ID,stringsAsFactors=F)
  new.dat=merge(original.ID,new.dat,by='ID',all.x=T)
  
  output=data.frame(pheno=new.dat$pheno,stringsAsFactors=F)
  colnames(output)=coln
  return(output)
}

residualise_chr <- function(target.dat,cov.dat,linkage.dat){
  cov.dat = merge(cov.dat,linkage.dat,by='id') %>%
    dplyr::select(-id) %>%
    .[,c(ncol(.),1:(ncol(.)-1))]
  colnames(cov.dat)[1]='ID'
  meth.dat = t(target.dat) %>%
    data.frame()%>%
    dplyr::mutate(ID=rownames(.)) %>%
    merge(.,cov.dat,by='ID')

  input.for_resid = data.frame(pheno=colnames(meth.dat)[grep('^cg',colnames(meth.dat))],
                               covs=paste0(colnames(cov.dat)[!grepl('Sample_Sentrix_ID|ID',colnames(cov.dat))],collapse='+'),
                               stringsAsFactors=F)
                               
  cat('Start processing\n')
  tmp.dat.aftercorr <- pbapply(input.for_resid,1,residualise_dat,dat=meth.dat) %>% # A list created by residualise_dat
    bind_cols %>%  # Convert into a dataframe
    mutate(.,ID=meth.dat$ID) %>% # Add ID
    select('ID', everything()) # Move ID to the first column
  
}

chr = args[1]

tmp.dat = readRDS(paste0('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr/mvalues.chr',chr,'.rds'))
cat(paste0('Data for CHR ',chr, ' - Loaded\n'))
dat.forPCA=residualise_chr(target.dat = tmp.dat,cov.dat = covariates,linkage.dat = sentrix.id)

output.path = paste0('/exports/eddie/scratch/xshen33/LBC21_n_36_w1_DNAm_bychr/forPCA/mvalues.chr',chr,'.forPCA.rds')

saveRDS(dat.forPCA,file = output.path)

# Save ID list for methylation data
sentrix.id=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/sentrix.rds') %>%
  .[,c(4)] %>% data.frame(ID=.)

saveRDS(sentrix.id,'/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/methID.rds')

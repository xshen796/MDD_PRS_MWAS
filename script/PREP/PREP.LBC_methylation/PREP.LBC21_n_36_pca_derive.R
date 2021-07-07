# XS
# 21/11/2020
library(plyr)            # Data management
library(dplyr)           # Data management
library(tibble)          # Data management
library(pbapply)         # Add progress bar for apply functions
library(FactoMineR)      # PCA

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')

# Load data ---------------------------------------------------------------

ls.file=
  list.files('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr/forPCA',
  pattern='^mvalues',full.names=T)

mval.forPCA = pblapply(as.list(ls.file),FUN=function(x) readRDS(x)) %>%
  join_all(by='ID', type='left') %>%
  column_to_rownames(.,var = 'ID')

targetdata = mval.forPCA

# perform PCA -------------------------------------------------------------

fit.pca = PCA(targetdata, scale.unit = TRUE, ncp = 20, graph = F)

methPCA = fit.pca$ind$coord

# Reformat to fit ewas pipeline
colnames(methPCA)=paste0('PC',1:20)
methPCA = methPCA %>% data.frame %>%
      add_rownames(var='Sample_Sentrix_ID') %>% data.frame
      

saveRDS(methPCA, 
        file='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/mval_PC20.rds',version=2)
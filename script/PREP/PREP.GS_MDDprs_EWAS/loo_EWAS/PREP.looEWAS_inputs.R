library(dplyr)
library(data.table)
library(tidyverse)
library(pbapply)

setwd('/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS')


# Reformat data for EWAS ------------------------------------------------------------------------------------------

# File names
ls.prs.f = list.files(path = './',pattern = '.all_score') 
ls.scorename = gsub('loo_','',ls.prs.f) %>%
      gsub('.all_score','',.)

# Read and process files
ls.prs = as.list(ls.prs.f) %>%
      pblapply(.,FUN=fread,stringsAsFactors=F) %>%
      pblapply(.,FUN=dplyr::select,score=Pt_1)
names(ls.prs) = ls.scorename
all.prs = ls.prs %>%
      bind_cols 
colnames(all.prs)=ls.scorename
all.prs$id = fread(ls.prs.f[1])[,'IID']
all.prs = all.prs %>%
      select(id,everything()) %>%
      as.data.frame
all.prs[,grep('PRS',colnames(all.prs))]=scale(all.prs[,grep('PRS',colnames(all.prs))])

# Rename PRS based on which SNP was leading
target.ls = fread('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_genetic_meth/SNP_variantMWAS.txt',header=T)
target.ls = target.ls %>%
      filter(!is.na(R2)) %>% 
      mutate(looSNP.seq=paste0('PRS_',1:nrow(.))) %>%
      data.frame
rownames(target.ls)=target.ls$looSNP.seq

ls.col.tochange = colnames(all.prs)[2:ncol(all.prs)]
col.toreplace=target.ls[ls.col.tochange,'SNP']

colnames(all.prs)[2:ncol(all.prs)]=col.toreplace

# Write data for EWAS
write.table(all.prs,
            file='/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_forEWAS/MDDprs_for_looEWAS.txt',col.names=T,row.names=F,sep=',',quote=F)
# Write a list of phenotypes (PRS) to run EWAS
write.table(col.toreplace,
            file='/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_forEWAS/MDDprs_name.txt',col.names=F,row.names=F,sep='\n',quote=F)

library(dplyr)
library(pbapply)
library(readr)


# ID lists -------------------------------------------------------------

cov.wave1 = readRDS('data/GS_phenotype/covs_RosieData_wave1.rds')
cov.wave3 = readRDS('data/GS_phenotype/covs_RosieData_wave3.rds')

ls.f=list.files('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/',pattern = '.profile',full.names = T)

prs = as.list(ls.f) %>%
  lapply(.,fread,stringsAsFactors=F,header=T) %>%
  lapply(.,select,FID,IID,starts_with('SCORE')) %>%
  Reduce(function(x,y) left_join(x,y,by=c('IID','FID')),.) %>%
  mutate_at(vars(matches('SCORE')), ~scale(.)) %>% 
  as.data.frame 

colnames(prs) = colnames(prs) %>% gsub('SCORE','pT',.)

ID.wave1.genetic = prs %>% select(FID,IID) %>% .[.$IID %in% cov.wave1$id,]
ID.wave3.genetic = prs %>% select(FID,IID) %>% .[.$IID %in% cov.wave3$id,]


# Write files -------------------------------------------------------------

write_tsv(ID.wave1.genetic,col_names = F,file='data/GS_genetic_meth/GeneticID_wave1')
write_tsv(ID.wave3.genetic,col_names = F,file='data/GS_genetic_meth/GeneticID_wave3')



library(dplyr)              # Data management 
library(data.table)         # Data management
library(purrr)              # Data management
library(readr)              # Data management

setwd('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats')


# Load data ---------------------------------------------------------------
IDP=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/Processed_IM_dat/data/IDP_instance2_40kmerged.rds')
baseline=readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2018-10-phenotypes-ukb24262/BaselineCharacteristics.rds') 
bolt_covars <- read_tsv('/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/whitebritish_centre_array_flashpcs_457k.tsv.gz')


# gFA and gMD phenotypes --------------------------------------------------

gFA_gwas <- IDP %>%
  select(f.eid, g.FA.Total) %>%
  transmute(FID=f.eid, IID=f.eid, gFA=g.FA.Total)
write_tsv(gFA_gwas, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/gFA_forGWAS.tsv')


gMD_gwas <- IDP %>%
  select(f.eid, g.MD.Total) %>%
  transmute(FID=f.eid, IID=f.eid, gMD=g.MD.Total)
write_tsv(gMD_gwas, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/gMD_forGWAS.tsv')

# gTR and ATR phenotypes --------------------------------------------------
gFA_TR_gwas <- IDP %>%
  select(f.eid, g.FA.TR) %>%
  transmute(FID=f.eid, IID=f.eid, gFA_TR=g.FA.TR)
write_tsv(gFA_TR_gwas, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/gFA_TR_forGWAS.tsv')

gMD_TR_gwas <- IDP %>%
  select(f.eid, g.MD.TR) %>%
  transmute(FID=f.eid, IID=f.eid, gMD_TR=g.MD.TR)
write_tsv(gMD_TR_gwas, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/gMD_TR_forGWAS.tsv')


FA_ATR_gwas <- IDP %>%
  mutate(FA.ATR = (lh.FA.wm.anterior_thalamic_radiation+rh.FA.wm.anterior_thalamic_radiation)/2) %>%
  select(f.eid, FA.ATR) %>%
  transmute(FID=f.eid, IID=f.eid, FA_ATR=FA.ATR)
write_tsv(FA_ATR_gwas, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/FA_ATR_forGWAS.tsv')

MD_ATR_gwas <- IDP %>%
  mutate(MD.ATR = (lh.MD.wm.anterior_thalamic_radiation+rh.MD.wm.anterior_thalamic_radiation)/2) %>%
  select(f.eid, MD.ATR) %>%
  transmute(FID=f.eid, IID=f.eid, MD_ATR=MD.ATR)
write_tsv(MD_ATR_gwas, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/MD_ATR_forGWAS.tsv')

# Covariates --------------------------------------------------------------

baseline.sex = baseline %>% select(f.eid,sex=f.31.0.0)

covs <- IDP %>%
  select(f.eid,age=info.MRI_age.UKB,site=info.MRI_site,
         pos.x=info.Scanner_lateral_X_brain_position, 
         pos.y=info.Scanner_transverse_Y_brain_position,
         pos.z=info.Scanner_longitudinal_Z_brain_position,
         pos.table=info.Scanner_table_position) %>%
  merge(.,baseline.sex,by='f.eid')

covs$site = as.factor(covs$site)
covs$sex = as.factor(covs$sex)

covs.final=bolt_covars %>%
  left_join(covs, by=c('IID'='f.eid'))

write_tsv(covs.final, '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/bolt_covariates.tsv')

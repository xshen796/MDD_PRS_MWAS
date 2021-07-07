setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')
library(dplyr)
library(modEvA)
library(pbapply)


# **** GS sample **** -----------------------------------------------------

methPredic.MDDprsCpG=read.delim('data/methPredictor/predictor_osca_reml_MDDstatus_w2_wg.profile',
                                stringsAsFactors = F) %>%
  select(IID,SCORE) %>%
  dplyr::rename(ID=IID,all_CpG_predictor=SCORE)

methPredic.noassocCpG=read.delim('data/methPredictor/predictor_osca_reml_MDDstatus_w2_MDDprs_noassoc.profile',
                                 stringsAsFactors = F) %>%
  select(IID,SCORE) %>%
  dplyr::rename(ID=IID,MDDprs_noassoc_CpG_predictor=SCORE)

methPredic.assocCpG=read.delim('data/methPredictor/predictor_osca_reml_MDDstatus_w2_MDDprs_assoc.profile',
                               stringsAsFactors = F) %>%
  select(IID,SCORE) %>%
  dplyr::rename(ID=IID,MDDprs_assoc_CpG_predictor=SCORE)

mdd.gs=readRDS('data/methTraining/DNAm_training_pheno_MDD.rds')

mdd.gs=merge(mdd.gs,methPredic.MDDprsCpG,by='ID',all.x=T) %>%
  merge(.,methPredic.noassocCpG,by='ID',all.x=T) %>%
  merge(.,methPredic.assocCpG,by='ID',all.x=T)


# Load permutation MRS ----------------------------------------------------

ls.MRS.f = list.files('data/methTraining/permutation_reml/scores',pattern = '.profile$',full.names = T)

permu.MRS <- as.list(ls.MRS.f) %>%
  pblapply(.,FUN=read.delim,stringsAsFactors = F) %>%
  pblapply(.,FUN=select,SCORE) %>%
  bind_cols %>%

permu.MRS$ID = read.delim(ls.MRS.f[1],stringsAsFactors = F)$IID
colnames(permu.MRS)=gsub('\\.','',colnames(permu.MRS))

permu.MRS.dat = merge(permu.MRS,mdd.gs,by='ID',all.x=T)


# Association tests for MRS -----------------------------------------------

source('FUNs/reg_phewasStyle_r2.R')

# Define global vars 
targetdata=permu.MRS.dat
targetdata$ever_smoke = as.numeric(targetdata$ever_smoke)

# dependent variables
ls.dep.all=c('MDD_status')
# factors
ls.factor=colnames(targetdata)[grep('^SCORE',colnames(targetdata))]
# combine the two
ls.dep.factor.combo=expand.grid(ls.dep.all,ls.factor,stringsAsFactors = F)
# covs
ls.models=data.frame(dependent=ls.dep.factor.combo$Var1,
                     factor=ls.dep.factor.combo$Var2,
                     covs='',stringsAsFactors = F)
ls.models$covs=paste0(c('age','sex',paste0('meth.PC',1:20),'smoking_pack_year','ever_smoke'),collapse='+')#'involuntary_unemployment','household_income','ADHD_ever','recent_socialdprv',

# specify models
ls.models$model.est='glm'
ls.models$p_batch=1


# Analysis ----------------------------------------------------------------

result.permutation=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=NA,correctByFactor = T)

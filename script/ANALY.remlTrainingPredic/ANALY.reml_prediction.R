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



#  MDD analysis ------------------------------------

# Individual score
fit.mddprs.cpg=glm(MDD_status~age+sex+ever_smoke+smoking_pack_year+
                     scale(MDDprs_assoc_CpG_predictor)+scale(MDDprs_noassoc_CpG_predictor),data=mdd.gs,family='binomial') 
summary(fit.mddprs.cpg)
mddprs.cpg.r2=Dsquared(model = fit.mddprs.cpg, adjust = F)

fit.mddprs.cpg=glm(MDD_status~age+sex+Batch+ever_smoke+smoking_pack_year+
                     scale(MDDprs_noassoc_CpG_predictor),data=mdd.gs,family='binomial') 
summary(fit.mddprs.cpg)
mddprs.cpg.r2=Dsquared(model = fit.mddprs.cpg, adjust = F)


fit.all.cpg=glm(MDD_status~age+sex+ever_smoke+smoking_pack_year+
                  scale(all_CpG_predictor),data=mdd.gs,family='binomial')
summary(fit.all.cpg)
all.cpg.r2=Dsquared(model = fit.all.cpg, adjust = F)

groot=filter(mdd.gs,!is.na(all_CpG_predictor))
fit.h0=glm(MDD_status~age+sex+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10,data=groot,family='binomial')
summary(fit.h0)
h0.r2=Dsquared(model = fit.h0, adjust = F)

# MDDprs assoc CpGs vs no assoc CpGs
fit.mddprs.cpg=glm(MDD_status~age+sex+ever_smoke+smoking_pack_year+
                     scale(MDDprs_assoc_CpG_predictor)+scale(MDDprs_noassoc_CpG_predictor),data=mdd.gs,family='binomial') 
summary(fit.mddprs.cpg)



mddprs.cpg.varR2=mddprs.cpg.r2-h0.r2
all.cpg.varR2=all.cpg.r2-h0.r2

res.table=rbind(summary(fit.all.cpg)$coefficients[nrow(summary(fit.all.cpg)$coefficients),],
                  summary(fit.mddprs.cpg)$coefficients[nrow(summary(fit.mddprs.cpg)$coefficients),])
res.table=data.frame(res.table,R2=100*c(mddprs.cpg.varR2,all.cpg.varR2))
res.table$Phenotype=c('MRS_MDDprsCpG','MRS_allCpG')

save(res.table,file='result/MDDprs_assoc_mPrediction/MRS_predictMDD.RData')


# lifestyle factors   ------------------------------------
# Define functions 
source('FUNs/reg_phewasStyle_r2.R')

# Define global vars 
targetdata=mdd.gs
targetdata$ever_smoke = as.numeric(targetdata$ever_smoke)

# dependent variables
ls.dep.all=c('MDD_status','alcohol_units_week','smoking_pack_year','ever_smoke')
# factors
ls.factor=c('MDDprs_assoc_CpG_predictor','MDDprs_noassoc_CpG_predictor','all_CpG_predictor')
# combine the two
ls.dep.factor.combo=expand.grid(ls.dep.all,ls.factor,stringsAsFactors = F)
# covs
ls.models=data.frame(dependent=ls.dep.factor.combo$Var1,
                     factor=ls.dep.factor.combo$Var2,
                     covs='',stringsAsFactors = F)
ls.models$covs=paste0(c('age','sex',paste0('meth.PC',1:20)),collapse='+')#'involuntary_unemployment','household_income','ADHD_ever','recent_socialdprv',

# specify models
ls.models$model.est='glm'
ls.models$p_batch=1


# Analysis ----------------------------------------------------------------

result.lifestyle=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=NA,correctByFactor = T)

save(result.lifestyle,file='result/MDDprs_assoc_mPrediction/MRS_MDDprs_assoc_noassoc_lifestyle_MDD.RData')




# **** LBC **** -----------------------------------------------------------

methPredic.MDDprsCpG=readRDS('data/methPredictor/Predi_MDDresid_LBC_MDDprsRelatedCpG_covPCSmoking.rds')
methPredic.MDDprsCpG$MDDprs_assoc_CpG_predictor=rowSums(methPredic.MDDprsCpG[,grep('^chr',colnames(methPredic.MDDprsCpG))])

methPredic.MDDprs_noassocCpG=readRDS('data/methPredictor/Predi_MDDresid_LBC_MDDprsnoassocCpG_covPCSmoking.rds')
methPredic.MDDprs_noassocCpG$MDDprs_noassoc_CpG_predictor=rowSums(methPredic.MDDprs_noassocCpG[,grep('^chr',colnames(methPredic.MDDprs_noassocCpG))])

covariates=readRDS('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/ewas_covar.rds')

sentrix.ID=readRDS('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/sentrix.rds')

LBC1921.wave1.pheno=
  read.csv("/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/LBC_phenotype/LBC1921_XueyiShen_EWAS_DepressionDisorderPGR_AM_24AUG2020.csv",
           stringsAsFactors=F,header=T)
LBC1936.wave1.pheno=
  read.csv("/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/LBC_phenotype/LBC1936_XueyiShen_EWAS_DepressionDisorderPGR_AM_24AUG2020.csv",
           stringsAsFactors=F,header=T)

LBC1921.lifestyle.depre = LBC1921.wave1.pheno %>%
  rename(ID = studyno, alc=alcpw, bmi = bmi, smoke = smoker) %>%
  select(ID,alc,bmi,smoke,HADSD)

LBC1936.lifestyle.depre = LBC1936.wave1.pheno %>%
  rename(ID = lbc36no, alc=alcunitwk_w1, bmi = bmi_w1, smoke = smokcat_w1, HADSD = hadsd_w1) %>%
  select(ID,alc,bmi,smoke,HADSD)

LBCboth.lifestyle.depre = rbind(LBC1921.lifestyle.depre,LBC1936.lifestyle.depre)

targetdata = merge(methPredic.MDDprsCpG,methPredic.MDDprs_noassocCpG,by='ID') %>%
  select(ID,starts_with('MDDprs')) %>%                       # Get MRS
  merge(.,sentrix.ID,by.x='ID',by.y='Sample_Sentrix_ID') %>% # Change ID from meth to pheno ID
  select(id,starts_with('MDDprs')) %>% 
  rename(ID=id) %>%                  
  merge(.,covariates,by.x='ID',by.y='id') %>%                # Merge with covariates
  merge(.,LBCboth.lifestyle.depre,by='ID')                   # Merge with lifestyle and depression variables

# Analysis ----------------------------------------------------------------

# Define functions 
source('FUNs/reg_phewasStyle_r2.R')

# Define global vars 
targetdata=targetdata

# dependent variables
ls.dep.all=c('alc','bmi','smoke','HADSD')
# factors
ls.factor=c('MDDprs_assoc_CpG_predictor','MDDprs_noassoc_CpG_predictor')
# combine the two
ls.dep.factor.combo=expand.grid(ls.dep.all,ls.factor,stringsAsFactors = F)
# covs
ls.models=data.frame(dependent=ls.dep.factor.combo$Var1,
                     factor=ls.dep.factor.combo$Var2,
                     covs='',stringsAsFactors = F)
ls.models$covs=paste0(c('age','sex'),collapse='+')#'involuntary_unemployment','household_income','ADHD_ever','recent_socialdprv',

# specify models
ls.models$model.est='glm'
ls.models$p_batch=1


result.lifestyle=reg_phewasStyle(ls.models,dat_short=targetdata,dat_long=NA,correctByFactor = T)


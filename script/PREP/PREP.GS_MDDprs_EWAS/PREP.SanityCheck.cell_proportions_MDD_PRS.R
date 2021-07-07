library(dplyr)
library(pbapply)
library(data.table)
library(purrr)

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')


# Load data -------------------------------------------------------------------------------------------------------

# PRS GS
prs.ls=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/score_thresholds.txt',stringsAsFactors=F,header=F)
colnames(prs.ls)=c('fname','label.pT','pT')
prs.ls$label.pT=as.character(prs.ls$label.pT)
prs.ls$label.pT=format(prs.ls$pT,scientific=F,drop0trailing=T)

gs.prs = list.files('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/',
                    pattern = 'Meta_predict_FullGS.',full.names = T) %>%
      as.list %>%
      lapply(.,fread,header=T) %>%
      lapply(.,select,IID, starts_with('SCORE')) %>%
      Reduce(function(x,y) left_join(x,y,by='IID'),.)

      
# Cell proportions
cov.dat = list.files('data/GS_phenotype',pattern='covs_RosieData',full.names=T) %>%
      as.list %>%
      lapply(.,readRDS) %>%
      bind_rows

cell.dat = fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/w1_w3_celltypes.txt') %>%
      select(-Batch)

# MDD
link_file = read.csv('data/GS_genetic_meth/gwas_methid_link.csv')
MDD.dat = readRDS('data/methTraining/DNAm_training_pheno_MDD.rds') %>%
      select(ID,MDD_status) %>%
      right_join(.,link_file,by=c('ID'='meth_id')) %>%
      select(gwas_id,MDD_status)

targetdata = right_join(gs.prs,cov.dat,by=c('IID'='id')) %>%
      left_join(.,cell.dat,by=c('IID'='Sample_Name')) %>%
      left_join(.,MDD.dat,by=c('IID'='gwas_id'))

# Analysis --------------------------------------------------------------------------------------------------------

# LHS of formula for each outcomes
outcomes_formulae = c('CD8T', 'CD4T', 'NK', 'Bcell', 'Mono' , 'Gran') 
# predictors to use in the model (RHS formula)
prediction_formula <- colnames(targetdata)  %>%  .[grep('^SCORE|MDD',.)]  %>% 
      paste0('sex+age+ever_smoke+pack_years+Batch+wave+',.)  
# update the outcomes formulae to add in the RHS
model_formulae <- expand.grid(outcomes_formulae,prediction_formula,stringsAsFactors = F) %>%
      purrr::transpose(.) %>%
      lapply(.,paste0,collapse='~') %>%
      lapply(.,as.formula)
      

# run GLM on each formula using data.frame MyData
models <- pblapply(model_formulae, glm, data=targetdata) 

# Extract outcomes
summ.res = models %>%
      lapply(.,summary) %>%
      lapply(.,function(x) x$coefficients[nrow(x$coefficients),]) %>%
      bind_rows %>%
      data.frame(.,expand.grid(outcomes_formulae,prediction_formula,stringsAsFactors = F)) %>%
      group_by(Var2) %>%
      mutate(p.corrected = p.adjust(`Pr...t..`,method='fdr'))

save(summ.res,file='result/EWAS_MDDprs_Shen/Cell_proportion_MDDprs/Reg_res.RData')
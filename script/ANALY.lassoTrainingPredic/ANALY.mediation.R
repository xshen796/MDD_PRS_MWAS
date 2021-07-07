setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')
library(dplyr)
library(modEvA)
library(pbapply)

#methPredic.MDDprsCpG=readRDS('data/methPredictor/Predi_MDDresid_wave1_MDDprsRelatedCpG_covPCLifeStyle.rds')
methPredic.MDDprsCpG=readRDS('data/methPredictor/Predi_MDDresid_wave1_MDDprsRelatedCpG_covPCLifeStyle.rds')
methPredic.MDDprsCpG$MDDprs_assoc_CpG_predictor=rowSums(methPredic.MDDprsCpG[,grep('^chr',colnames(methPredic.MDDprsCpG))])

methPredic.MDDprs_noassocCpG=readRDS('data/methPredictor/Predi_MDDresid_wave1_MDDprs_noassoc_CpG_covPCLifeStyle.rds')
methPredic.MDDprs_noassocCpG$MDDprs_noassoc_CpG_predictor=rowSums(methPredic.MDDprs_noassocCpG[,grep('^chr',colnames(methPredic.MDDprs_noassocCpG))])

mdd.gs=readRDS('data/DNAm_training_pheno_MDD.rds')

mdd.gs=merge(mdd.gs,methPredic.MDDprsCpG[,c('ID','MDDprs_assoc_CpG_predictor')],by='ID',all.x=T)
mdd.gs=merge(mdd.gs,methPredic.MDDprs_noassocCpG[,c('ID','MDDprs_noassoc_CpG_predictor')],by='ID',all.x=T)

methPredic.MDDcpg=readRDS('data/methPredictor/Predi_MDDresid_wave1_all.rds')
colnames(methPredic.MDDcpg)[2]='all_CpG_predictor'
mdd.gs=merge(mdd.gs,methPredic.MDDcpg,by='ID',all.x=T)

gw.prs.wave1=read.csv('/gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/MDD_prs/related/MDDprs_e-8_phenoWave1_forEWAS.csv')
linkage=read.csv('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/gwas_methid_link.csv',header=T,stringsAsFactors=F)
gw.prs.wave1=merge(gw.prs.wave1,linkage,by.x='id',by.y='gwas_id',all.x=T)
mdd.gs=merge(mdd.gs,gw.prs.wave1[,c('meth_id','MDDprs_gwsig')],by.x='ID',by.y='meth_id',all.x=T)

# Mediation -----------------------------------------------------------------------------

vars=c('MDDprs_gwsig','MDDprs_assoc_CpG_predictor','MDD_status')
targetdata=mdd.gs[complete.cases(mdd.gs[,c('ID',vars)]),]

targetdata$sex=as.numeric(targetdata$sex)
targetdata$wave=as.numeric(targetdata$wave)
targetdata[,2:ncol(targetdata)]=scale(targetdata[,2:ncol(targetdata)])

set.seed(1234)
X <- targetdata$MDDprs_gwsig
M <- targetdata$MDDprs_assoc_CpG_predictor
Y <- targetdata$MDD_status
Data <- data.frame(X = X, Y = Y, M = M,targetdata)
model <- ' 
           # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # regression
             X ~~ age+sex+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10
             M ~~ age+sex
             Y ~~ age+sex
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
fit <- sem(model, data = Data)
summary(fit)
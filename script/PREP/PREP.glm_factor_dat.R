library(dplyr)

setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD')

# load DNAm predictors and change ID
predictor=readRDS('data/predictors_pack_year.rds')
colnames(predictor)[ncol(predictor)]='smoking_DNAm'
predictor=merge(predictor,readRDS('data/predictors_alc_units.rds'),by='ID',all=T)
colnames(predictor)[ncol(predictor)]='alc_DNAm'
predictor=merge(predictor,readRDS('data/predictors_bmi.rds'),by='ID',all=T)
colnames(predictor)[ncol(predictor)]='bmi_DNAm'

linkage_file=read.csv('data/stradl_genscot_methyl_id_link.csv')
predictor = merge(predictor,linkage_file,by.x='ID',by.y='meth_id',all.x=T)
predictor$ID = predictor$st_id
predictor=predictor[,!grepl('_id$',colnames(predictor))]


# load phenotype data and change ID
pheno = readRDS('data/DNAm_training_pheno.rds')
pheno = merge(pheno,linkage_file,by.x='ID',by.y='meth_id',all.y=T)
pheno$ID = pheno$st_id
pheno = pheno[,!grepl('_id$',colnames(pheno))]


# merge two datasets
factor.dat=merge(pheno,predictor,by='ID',all=T)

saveRDS(factor.dat,file='data/glm_factors.rds')

#!/bin/sh
#$ -N methPredictor
#$ -cwd
#$ -m beas
#$ -l h_vmem=16G
#$ -l h_rt=3:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load R

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD

# Top SNP assoc CpGs trained predictor: GS

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R \
--methdat /exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/chr/ \
--subID data/GS_phenotype/wave3_all_ID.rds \
--lassoCoef data/methTraining/weights_wave1_residMDD_EPIC_MDDprs_e-8assoc_CpG.txt \
--out data/methPredictor/Predi_MDDresid_GSwave3_MDDprsRelatedCpG_covPCSmoking.rds 


# Top SNP no assoc CpGs trained predictor: GS

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R \
--methdat /exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/chr/ \
--subID data/GS_phenotype/wave3_all_ID.rds \
--lassoCoef data/methTraining/weights_wave1_residMDD_EPIC_MDDprs_e-8_noassoc_CpG.txt \
--out data/methPredictor/Predi_MDDresid_GSwave3_MDDprs_noassoc_CpG_covPCSmoking.rds 


# WG CpGs trained predictor: GS

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R \
--methdat /exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/chr/ \
--subID data/GS_phenotype/wave3_all_ID.rds \
--lassoCoef data/methTraining/weights_wave1_residMDD_EPIC_wg_CpG.txt \
--out data/methPredictor/Predi_MDDresid_GSwave3_MDDprs_wg_CpG_covPCSmoking.rds 


# Top SNP assoc CpGs trained predictor: LBC

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R \
--methdat /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr \
--subID data/LBC_phenotype/methID.rds \
--lassoCoef data/methTraining/weights_wave1_residMDD_EPIC_MDDprs_e-8assoc_CpG.txt \
--out data/methPredictor/Predi_MDDresid_LBC_MDDprsRelatedCpG_covPCSmoking.rds


# Top SNP no assoc CpGs trained predictor: LBC

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R \
--methdat /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr \
--subID data/LBC_phenotype/methID.rds \
--lassoCoef data/methTraining/weights_wave1_residMDD_EPIC_MDDprs_e-8_noassoc_CpG.txt \
--out data/methPredictor/Predi_MDDresid_LBC_MDDprsnoassocCpG_covPCSmoking.rds


# WG CpGs trained predictor: GS

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R \
--methdat /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr \
--subID data/GS_phenotype/wave3_all_ID.rds \
--lassoCoef data/methTraining/weights_wave1_residMDD_EPIC_wg_CpG.txt \
--out data/methPredictor/Predi_MDDresid_LBC_MDDprs_wg_CpG_covPCSmoking.rds 


#!/bin/sh
#$ -N lassoMDD
#$ -cwd
#$ -m beas
#$ -l h_vmem=125G
#$ -l h_rt=12:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load R

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/

# Top SNP assoc CpGs training

Rscript script/FUNs_meth/bigLASSO_selectedCpG.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves/ data/GS_phenotype/wave1_all_ID.rds data/methTraining/cpgList_MDDprs_assoc.rds data/DNAm_training_pheno_MDD.rds resid.MDD_status data/methWeights/weights_wave3_residMDD_450k_MDDprs_e-8assoc_CpG.txt

# Top SNP no assoc CpGs training

Rscript script/FUNs_meth/bigLASSO_selectedCpG.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves/ data/GS_phenotype/wave1_all_ID.rds data/methTraining/cpgList_MDDprs_NOassoc.rds data/DNAm_training_pheno_MDD.rds resid.MDD_status data/methWeights/weights_wave3_residMDD_450k_MDDprs_e-8_noassoc_CpG.txt
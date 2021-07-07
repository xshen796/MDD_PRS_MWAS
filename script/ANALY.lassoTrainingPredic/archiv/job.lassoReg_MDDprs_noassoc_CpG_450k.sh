#!/bin/sh
#$ -N lassoMDD
#$ -cwd
#$ -m beas
#$ -l h_vmem=125G
#$ -l h_rt=4:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load R

cd /exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD

Rscript script/FUNs_meth/bigLASSO_selectedCpG.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/all_450mvalues.txt data/wave3_all_ID.rds data/cpgList_MDDprs_NOassoc.rds data/DNAm_training_pheno_MDD.rds resid.MDD_status data/methWeights/weights_wave3_residMDD_450k_MDDprs_e-8_noassoc_CpG.txt
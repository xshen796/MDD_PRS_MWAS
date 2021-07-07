#!/bin/sh
#$ -N methPredictor
#$ -cwd
#$ -m beas
#$ -l h_vmem=16G
#$ -l h_rt=3:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load R

cd /exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves data/wave1_all_ID.rds data/methWeights/weights_wave3_residMDD_450k_MDDprs_e-8_noassoc_CpG.txt data/methPredictor/Predi_MDDresid_wave1_MDDprs_noassoc_CpG_covPCLifeStyle.rds

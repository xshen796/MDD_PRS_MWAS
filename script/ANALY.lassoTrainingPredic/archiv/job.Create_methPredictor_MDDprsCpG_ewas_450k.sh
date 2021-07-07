#!/bin/sh
#$ -N methPredictor
#$ -cwd
#$ -m beas
#$ -l h_vmem=64G
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/

Rscript script/FUNs_meth/Create_DNAm_scores_byCHR.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves data/wave1_all_ID.rds data/weights_wave2_WholeGenome_WholeSample_450k/weights_mddprs_ewas.txt data/methPredictor/Predi_MDDprs_wave1.rds

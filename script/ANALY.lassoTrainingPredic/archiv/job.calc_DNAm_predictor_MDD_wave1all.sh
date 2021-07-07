#!/bin/sh
#$ -N methPredictor
#$ -cwd
#$ -m beas
#$ -l h_vmem=125G
#$ -l h_rt=6:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/

Rscript script/Create_DNAm_scores_Xshen.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/all_450mvalues.txt data/wave1_all_ID.rds data/weights_wave2_WholeGenome_WholeSample_450k/weights_wave2_residMDD_WholeG_WholeSample_450k.txt data/methPredictor/Predi_MDDresid_wave1_all_covagesex_noncovlifestyle.rds

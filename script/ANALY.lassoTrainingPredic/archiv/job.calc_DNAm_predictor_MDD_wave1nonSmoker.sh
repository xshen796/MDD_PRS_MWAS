#!/bin/sh
#$ -N methPredictor
#$ -cwd
#$ -m beas
#$ -l h_vmem=125G
#$ -l h_rt=4:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/eddie/scratch/xshen33/MR_meth_MDD

Rscript script/Create_DNAm_scores_Xshen.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/all_450mvalues.txt data/wave1_all_ID.rds data/weights_wave2_WholeGenome_nonSmoker_450k/weights_wave2_residMDD_WholeG_nonSmoker_450k.txt data/methPredictor/Predi_MDDresid_wave1_nonSmoker.rds

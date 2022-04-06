#!/bin/sh
#$ -N cor_mat_cpg
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -l h_rt=36:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile


Rscript /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/PREP/PREP.mQTL_MR/PREP.cor_mat_CpG.R

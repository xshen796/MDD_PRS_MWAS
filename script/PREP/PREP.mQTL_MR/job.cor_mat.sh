#!/bin/sh
#$ -N cor_mat_cpg
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -l h_rt=36:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load igmm/apps/R/3.5.1

cd /gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/script/PREP/

Rscript PREP.cor_mat_CpG.R $1
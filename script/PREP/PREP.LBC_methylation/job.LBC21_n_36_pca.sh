#!/bin/sh
#$ -cwd
#$ -M xueyi.shen@ed.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -l h_rt=36:00:00



. /etc/profile.d/modules.sh

#module unload R
module load igmm/apps/R/3.6.1

Rscript PREP.LBC21_n_36_pca.R $1

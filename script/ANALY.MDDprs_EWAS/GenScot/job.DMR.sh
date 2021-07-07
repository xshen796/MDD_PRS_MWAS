#!/bin/sh
#$ -N DMR
#$ -cwd
#$ -m beas
#$ -M xueyi.shen@ed.ac.uk
#$ -o logs
#$ -e logs
#$ -l h_vmem=64G
#$ -pe sharedmem 6
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

module load igmm/apps/R/3.6.1

R CMD BATCH ANALY.DMR.R

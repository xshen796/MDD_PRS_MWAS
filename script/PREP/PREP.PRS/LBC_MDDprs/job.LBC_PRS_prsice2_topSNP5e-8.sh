#!/bin/sh
#$ -N LBC_PRS
#$ -cwd
#$ -m beas
#$ -M xueyi.shen@ed.ac.uk
#$ -o logs
#$ -e logs
#$ -l h_vmem=32G
#$ -l h_rt=36:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile


module unload R
module load igmm/apps/R/3.5.1

bash LBC_PRS_prsice2_topSNP5e_08.sh -g $1 -p $2 -o $3

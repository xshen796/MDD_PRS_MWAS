#!/bin/sh
#$ -N GoDMC_dat_parse
#$ -cwd
#$ -M xueyi.shen@ed.ac.uk
#$ -m beas
#$ -l h_vmem=64G
#$ -l h_rt=48:00:00


. /etc/profile.d/modules.sh
source ~/.bash_profile

R CMD BATCH PREP.GoDMC_mQTL.R
#!/bin/sh
#$ -N lassoMDD
#$ -cwd
#$ -m beas
#$ -l h_vmem=125G
#$ -l h_rt=4:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD

Rscript script/CpG_whole_genome_bigLASSO.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/all_450mvalues.txt data/wave1_all_ID.rds data/DNAm_training_pheno_MDD.rds resid.MDD_status data/weights_wave1_residMDD_WholeG_WholeSample_450k.txt

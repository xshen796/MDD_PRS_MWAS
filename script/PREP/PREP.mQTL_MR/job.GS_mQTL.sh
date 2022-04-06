#!/bin/sh
#$ -N mQTL
#$ -m beas
#$ -l h_vmem=32G
#$ -pe sharedmem 3
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/eddie/scratch/xshen33/mQTL

#mkdir data

####################################################################################################################################
############################################         Prepare GS genetic data            ############################################
####################################################################################################################################

cp /exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA.* data/

plink --bfile data/GS20K_HRC_0.8_GCTA --keep data/GS_ID_withMeth.txt --make-bed --out data/GS10K_HRC_0.8_GCTA_withMeth

plink --bfile data/GS10K_HRC_0.8_GCTA_withMeth --update-ids data/updateID_withMeth.txt --make-bed --out data/GS10K_HRC_0.8_GCTA_MethID


####################################################################################################################################
############################################             methQTL analysis               ############################################
####################################################################################################################################

osca_Linux --eqtl \
--bfile data/GS10K_HRC_0.8_GCTA_MethID \
--befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w1 \
--extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/mQTL/cpg_ls.txt \
--thread-num 10 \
--task-num 1 \
--covar /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/methTraining/DNAm_training_dcov_osca \
--qcovar /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/methTraining/DNAm_training_qcov_osca \
--out /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/mQTL/mQTL_w1

osca_Linux --eqtl \
--bfile data/GS10K_HRC_0.8_GCTA_MethID \
--befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w2 \
--extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/mQTL/cpg_ls.txt \
--thread-num 10 \
--task-num 1 \
--covar /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/methTraining/DNAm_training_dcov_osca \
--qcovar /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/methTraining/DNAm_training_qcov_osca \
--out /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/mQTL/mQTL_w2


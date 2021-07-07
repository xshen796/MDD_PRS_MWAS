#!/bin/sh
#$ -N lassoMDD
#$ -cwd
#$ -m beas
#$ -pe sharedmem 4
#$ -l h_vmem=64G
#$ -l h_rt=24:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

#module load R

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/

# Top SNP assoc CpGs training  /exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/chr/

Rscript script/FUNs_meth/bigLASSO_selectedCpG.R \
--methFolder /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves \
--subID data/GS_phenotype/wave1_all_ID.rds \
--cpgList data/methTraining/cpgList_MDDprs_assoc.rds \
--phenoFile data/methTraining/DNAm_training_pheno_MDD.rds \
--phenoName resid.MDD_status \
--outputFile data/methTraining/weights_wave1_residMDD_EPIC_MDDprs_e-8assoc_CpG.txt



# Top SNP no assoc CpGs training

Rscript script/FUNs_meth/bigLASSO_selectedCpG.R \
--methFolder /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves \
--subID data/GS_phenotype/wave1_all_ID.rds \
--cpgList data/methTraining/cpgList_MDDprs_NOassoc.rds \
--phenoFile data/methTraining/DNAm_training_pheno_MDD.rds \
--phenoName resid.MDD_status \
--outputFile data/methTraining/weights_wave1_residMDD_EPIC_MDDprs_e-8_noassoc_CpG.txt


# All CpGs training

Rscript script/FUNs_meth/bigLASSO_selectedCpG.R \
--methFolder /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves \
--subID data/GS_phenotype/wave1_all_ID.rds \
--phenoFile data/methTraining/DNAm_training_pheno_MDD.rds \
--phenoName resid.MDD_status \
--outputFile data/methTraining/weights_wave1_residMDD_EPIC_wg_CpG.txt

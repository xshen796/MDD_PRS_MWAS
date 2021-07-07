#!/bin/sh
#$ -N remlMDD
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/

# REML
osca_Linux --reml \
--orm /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/orm_w1 \
--pheno /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/DNAm_training_pheno_MDD_osca \
--reml-pred-rand \
--covar /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/DNAm_training_dcov_osca \
--qcovar /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/DNAm_training_qcov_osca \
--out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1

####################################################################################################################################
########################################################   WG     ##################################################################
####################################################################################################################################

# Estimate probe effects
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w1 \
--blup-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1.indi.blp \
--out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1_wg_solution

# Create scores
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w2 \
 --score /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1_wg_solution.probe.blp \
 --score-has-header \
 --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methPredictor/predictor_osca_reml_MDDstatus_w2_wg


####################################################################################################################################
##################################################     MDDprs assoc     ############################################################
####################################################################################################################################

# Estimate probe effects
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w1 \
--extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/cpgList_MDDprs_assoc.txt \
--blup-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1.indi.blp \
--out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1_MDDprs_assoc_solution

# Create scores
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w2 \
  --extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/cpgList_MDDprs_assoc.txt \
 --score /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1_MDDprs_assoc_solution.probe.blp \
 --score-has-header \
 --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methPredictor/predictor_osca_reml_MDDstatus_w2_MDDprs_assoc



####################################################################################################################################
###############################################         MDDprs no assoc           ##################################################
####################################################################################################################################

# Estimate probe effects
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w1 \
--exclude-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/cpgList_MDDprs_assoc.txt \
--blup-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1.indi.blp \
--out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1_MDDprs_noassoc_solution

# Create scores
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w2 \
  --exclude-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/cpgList_MDDprs_assoc.txt \
 --score /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1_MDDprs_noassoc_solution.probe.blp \
 --score-has-header \
 --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methPredictor/predictor_osca_reml_MDDstatus_w2_MDDprs_noassoc



####################################################################################################################################
##############################################      MDDprs no assoc - permutation       ############################################
####################################################################################################################################

# Generate a randomised list of CpGs
mkdir data/methTraining/permutation_reml
mkdir data/methTraining/permutation_reml/cpg_list
mkdir data/methTraining/permutation_reml/scores
mkdir data/methTraining/permutation_reml/weight

echo "Permutation test start" > data/methTraining/permutation_reml/permu.log

for permu in {1..1000}
do

echo "Iteration: ${permu}" >> data/methTraining/permutation_reml/permu.log


# Generate CpG list
Rscript script/ANALY.remlTrainingPredic/assoc.remlReg_permutation.R --CpGList data/methTraining/cpgList_MDDprs_NOassoc.txt \
--outDir data/methTraining/permutation_reml/cpg_list --outFile list.$permu.txt

# Estimate probe effects
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w1 \
--extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/cpg_list/list.$permu.txt \
--blup-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1.indi.blp \
--out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/weight/blp_w1_solution_$permu

# Create scores
osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w2 \
 --extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/cpg_list/list.$permu.txt \
 --score /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/weight/blp_w1_solution_$permu.probe.blp \
 --score-has-header \
 --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/scores/predi_reml_w2_$permu

rm data/methTraining/permutation_reml/cpg_list/list.$permu.txt
rm /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/weight/blp_w1_solution_$permu*

done


#!/bin/sh
#$ -N remlMDD
#$ -cwd
#$ -m beas
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/

for permu in {1..10000}
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
 --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methPredictor/permutation_scores_CpG_MDDprs_assoc/predi_reml_w2_$permu

rm data/methTraining/permutation_reml/cpg_list/list.$permu.txt
rm /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/weight/blp_w1_solution_$permu*

done

# Clean up files
rm /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/scores/*.log
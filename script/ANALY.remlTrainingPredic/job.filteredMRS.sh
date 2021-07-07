#!/bin/sh
#$ -N remlMDD
#$ -cwd
#$ -m beas
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/

for permu in 0.00000005 0.0000001 0.0000005 0.000001 0.000005 0.00001 0.00005 0.0001 0.005 0.001 0.05 0.01 0.5 0.1 1
do

      echo "Iteration: ${permu}" >> data/methTraining/permutation_reml/permu.log
      
      
      # Generate CpG list
      Rscript script/ANALY.remlTrainingPredic/assoc.remlReg_filteredMRS.R --ewasSummstats result/EWAS_MDDprs_Shen/archiv/MDDprs_ewas_meta/MDDprs_pT5e_08_meta/mddprs_pT5e_08_meta_RosieData.toptable.txt1.tbl --pT $permu \
      --outDir data/methTraining/filteredMRS/CpG_list --outFile pT.$permu.txt
      
      # Estimate probe effects
      osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w1 \
      --extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/filteredMRS/CpG_list/pT.$permu.txt \
      --blup-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/blp_w1.indi.blp \
      --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/filteredMRS/weight/blp_w1_solution_$permu
      
      # Create scores
      osca_Linux --befile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves_Rosie_BOD/mvalue_w2 \
       --extract-probe /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/filteredMRS/CpG_list/pT.$permu.txt \
       --score /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/filteredMRS/weight/blp_w1_solution_$permu.probe.blp \
       --score-has-header \
       --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methPredictor/filteredMRS/predi_reml_w2_$permu
      
      rm data/methTraining/filteredMRS/CpG_list/list.$permu.txt
      rm /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/filteredMRS/weight/blp_w1_solution_$permu*

done

# Clean up files
rm /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/methTraining/permutation_reml/scores/*.log
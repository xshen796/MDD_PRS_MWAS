### Run MWAS on each wave separately. PCs included as covariates  =====================================================

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/MWAS_geneticPC

stradl_ewas --pdata /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/MDDprs_MWAS_input.txt --pheno pT_0.00000005 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave1.rds --no-prune --out  ewas_w1_geneticPC_pT_0.00000005

stradl_ewas_w3 --pdata /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/MDDprs_MWAS_input.txt --pheno pT_0.00000005 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave3.rds --no-prune --out  ewas_w3_geneticPC_pT_0.00000005


### Reformat sumstats and prepare Metal script for meta-analysis   ====================================================

R CMD BATCH /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/ANALY.MDDprs_EWAS/GenScot/GeneticPC_MWAS/Reformat_MetalScript.R

### Run Metal   =======================================================================================================

ls MetalScript_* > metal_script.list.txt

# Run meta analysis 
module load igmm/apps/metal

while read prs; do
  metal ${prs}
done < metal_script.list.txt

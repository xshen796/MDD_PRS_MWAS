# wiki page: https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=Psychiatry&title=STRADL%20EWAS%20pipeline&fbclid=IwAR13_bpLWyU-p6BXIb8XWpTR2VAxJIzXFg9cjGW8A57HlQBR4e1WxU3DMZI
# Pipeline location:  /exports/igmm/eddie/GenScotDepression/local/EWAS


####################################################      Run PRS      ###############################################################

# Run prsice
cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/PREP/PREP.PRS/GS_MDDprs
qsub job.GS_PRS_noMHC.sh

# Reformat PRS for EWAS
# check script: /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/PREP/PREP.GS_MDDprs_EWAS/PREP.MDD_prs_EWAS_input.R

#####################################################     Run EWAS      ##############################################################


cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/noMHC_and_gwsig

stradl_ewas --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC_forEWAS.txt --pheno Pt_5e.08 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave1.rds --no-prune --out ewas_w1_noMHC_5e_08

stradl_ewas_w3 --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC_forEWAS.txt --pheno Pt_5e.08 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave3.rds --no-prune --out ewas_w3_noMHC_5e_08


######################################################    Meta-analyse two waves     ##################################################################

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/noMHC_and_gwsig

#R CMD BATCH /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/ANALY.MDDprs_EWAS/loo_EWAS/metal_script.R

# Run meta analysis 
module load igmm/apps/metal
metal MetalScript_noMHC_5e_08.txt

# wiki page: https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=Psychiatry&title=STRADL%20EWAS%20pipeline&fbclid=IwAR13_bpLWyU-p6BXIb8XWpTR2VAxJIzXFg9cjGW8A57HlQBR4e1WxU3DMZI
# Pipeline location:  /exports/igmm/eddie/GenScotDepression/local/EWAS

stradl_ewas --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC_forEWAS.txt --pheno Pt_5e.08 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_phenotype/covs_RosieData_wave1.rds --no-prune --out ewas_w1_noMHC_5e_08

stradl_ewas_w3 --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC_forEWAS.txt --pheno Pt_5e.08 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_phenotype/covs_RosieData_wave3.rds --no-prune --out ewas_w3_noMHC_5e_08


######################################################    Meta-analyse two waves     ##################################################################

cd /exports/eddie/scratch/xshen33/GS_PRS_loo/
mkdir MetaWaves_EWAS
cd MetaWaves_EWAS

R CMD BATCH /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/ANALY.MDDprs_EWAS/loo_EWAS/metal_script.R

# Run meta analysis 
module load igmm/apps/metal
metal MetalScript_noMHC_5e_08.txt

cp ewas_BothWaves_noMHC_5e_08.metal.out1.tbl* /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_Shen/noMHC


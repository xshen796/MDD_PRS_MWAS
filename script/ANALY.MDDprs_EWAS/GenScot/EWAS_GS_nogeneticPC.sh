# wiki page: https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=Psychiatry&title=STRADL%20EWAS%20pipeline&fbclid=IwAR13_bpLWyU-p6BXIb8XWpTR2VAxJIzXFg9cjGW8A57HlQBR4e1WxU3DMZI
# Pipeline location:  /exports/igmm/eddie/GenScotDepression/local/EWAS


#########################################################     EWAS using pipeline    ##################################################################

# For all p thresholds - PRS calculated by Dave Howard (from the Nat Neuroscience paper) ---------------------------------------------------------------------------------
cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/MWAS_by_wave
while read prs; do

    stradl_ewas --pdata /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/MDDprs_MWAS_input.txt --pheno ${prs} --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave1.rds --no-prune --out  ewas_w1_${prs}

    stradl_ewas_w3 --pdata /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/MDDprs_MWAS_input.txt --pheno ${prs} --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave3.rds --no-prune --out ewas_w3_${prs}

done < /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/ls_MDDprs.txt

########################################################     EWAS using pipeline    ##################################################################

# pT=5e-8 - unrelated participants, PRS calculated by Dave Howard (from the Nat Neuroscience paper)
cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/MWAS_by_wave

stradl_ewas --pdata /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/MDDprs_MWAS_input_unrelated.txt --pheno pT_0.00000005 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave1.rds --no-prune --out  ewas_w1_unrelated_pT_0.00000005

cp ewas_w3_pT_0.00000005.toptable.txt ewas_w3_unrelated_pT_0.00000005.toptable.txt


# MWAS on PRS constructed from the leading variants (101 variants, one was missing because was the only sig SNP for that region and not present in GS data) --------------
# Calculate PRS pt=5e-8 so to obtain a SNP list after clumping
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_5e_08
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/data
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/

# variant PRS: check script/PREP/PREP.PRS/GS_MDDprs/pipeline.sh

# Reformat PRS output for MWAS: /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/PREP/PREP.GS_MDDprs_EWAS/loo_EWAS/PREP.looEWAS_inputs.R

# Run MWAS : ANALY.MDDprs_EWAS/GenScot/variantEWAS
stradl_ewas --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_5e_08_forEWAS.txt --pheno Pt_5e.08 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_RosieData_wave1.rds --no-prune --out ewas_w1_gwsig

stradl_ewas_w3 --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_5e_08_forEWAS.txt --pheno Pt_5e.08 --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_RosieData_wave3.rds --no-prune --out ewas_w3_gwsig

# MWAS on non-MHC PRS pt=5e-8   ------------------------------------------------------------------------------------------------------------------------------------------
# Check: noMHC_EWAS.sh



######################################################    Meta-analyse two waves     ##################################################################

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/MDDprs_ewas_meta

R CMD BATCH /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/ANALY.MDDprs_EWAS/GenScot/metal_script_dat_allpT.R

ls MetalScript_* > metal_script.list.txt

# Run meta analysis 
module load igmm/apps/metal

while read prs; do
  metal ${prs}
done < metal_script.list.txt


#mkdir /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_meta_LBC_ALSPAC

#cp REPmeta* /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_meta_LBC_ALSPAC/


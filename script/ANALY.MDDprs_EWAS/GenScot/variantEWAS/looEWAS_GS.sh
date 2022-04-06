# wiki page: https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=Psychiatry&title=STRADL%20EWAS%20pipeline&fbclid=IwAR13_bpLWyU-p6BXIb8XWpTR2VAxJIzXFg9cjGW8A57HlQBR4e1WxU3DMZI
# Pipeline location:  /exports/igmm/eddie/GenScotDepression/local/EWAS


#########################################################     EWAS using pipeline    ##################################################################
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/EWAS
cd /exports/eddie/scratch/xshen33/GS_PRS_loo/EWAS

while read snp; do

    stradl_ewas --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_forEWAS/MDDprs_for_looEWAS.txt --pheno ${snp} --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave1.rds --no-prune --out ewas_w1_${snp}

    stradl_ewas_w3 --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_forEWAS/MDDprs_for_looEWAS.txt --pheno ${snp} --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_phenotype/covs_geneticPC_RosieData_wave3.rds --no-prune --out ewas_w3_${snp}

done < ../PRS_forEWAS/MDDprs_name.txt



######################################################    Meta-analyse two waves     ##################################################################

cd /exports/eddie/scratch/xshen33/GS_PRS_loo/
rm -r MetaWaves_EWAS
mkdir MetaWaves_EWAS
cd MetaWaves_EWAS

R CMD BATCH /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/ANALY.MDDprs_EWAS/GenScot/variantEWAS/metal_script_dat.R

# Run meta analysis 
module load igmm/apps/metal

while read snp; do
  metal MetalScript_${snp}.txt
done < ../PRS_forEWAS/MDDprs_name.txt


#mkdir /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_meta_LBC_ALSPAC

#cp REPmeta* /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_meta_LBC_ALSPAC/


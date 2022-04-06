#!/bin/sh
#$ -N GS_PRS
#$ -cwd
#$ -m beas
#$ -l h_vmem=32G
#$ -l h_rt=24:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile



prsice_R_file='/exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R'
prsice_binary_file='/exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux'

# space/tab deliminated summary statistics in plain text
gwas_summstats='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.noMHC.forPRSice.gz' 

# prefix of the bed,fam,bim files
plink_files='/exports/eddie/scratch/xshen33/GS_PRS_loo/data/GS20K_HRC_0.8_GCTA'

# Phenotype for MDD diagnosis: space/tab deliminated file. Columns: FID, IID, phenotype (e.g. MDD_diagnosis)
phenotype_file='/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_loo_PRS/dummy.MDD.GS'
pheno_name='dummy_MDD'

# prefix for the output files
output_prefix='/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_allpT/MDDprs_noMHC'


# PRSice script:

Rscript $prsice_R_file \
    --prsice $prsice_binary_file \
    --base $gwas_summstats \
    --target $plink_files \
    --thread 3 \
    --stat BETA \
    --binary-target T \
    --allow-inter \
    --pheno $phenotype_file \
    --pheno-col $pheno_name \
    --clump-kb 500 \
    --bar-levels 1,0.5,0.1,0.05,0.01,0.005,0.001,0.0001,0.00001,0.00000005 \
    --extract /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_5e_08/MDDprs_5e_08.valid \
    --fastscore \
    --print-snp \
    --all-score \
    --out $output_prefix

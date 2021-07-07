# This script needs to run with R

# Settings ------------------------------------------------------------------------------

prsice_R_file='/exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R'
prsice_binary_file='/exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux'

# space/tab deliminated summary statistics in plain text
gwas_summstats='/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/daner_pgc_mdd_noUKBB_eur_hg19_v3.229.09_forPRS.txt' 

# prefix of the bed,fam,bim files
plink_files='/exports/eddie/scratch/xshen33/autosome.qc.maf01.hwe5e-6.geno02.mind02.snps.unrelatedcaucasian' 

# Phenotype for MDD diagnosis: space/tab deliminated file. Columns: FID, IID, phenotype (e.g. MDD_diagnosis)
phenotype_file='/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/ukb_mdd_phenotypes.txt'
pheno_name='MDD_diagnosis'

# prefix for the output files: directory can be added before the prefix to indicate a different directory to store generated PRS in a directory different from your current working directory
output_prefix='/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/mdd_prs/pgc3_meta_ukb_PRS' 



# PRSice script:

Rscript $prsice_R_file \
    --prsice $prsice_binary_file \
    --base $gwas_summstats \
    --target $plink_files \
    --thread 3 \
    --stat OR \
    --binary-target T \
    --pheno $phenotype_file \
    --pheno-col $pheno_name \
    --extract ls.102hits_MDDgwas.txt
    --dose-thres 0.9 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --bar-levels 1 \
    --fastscore \
    --all-score \
    --out /exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/data/mdd_prs/pgc3_meta_ukb_PRS
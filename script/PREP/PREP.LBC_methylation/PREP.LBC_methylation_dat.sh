# Methylation data preparation scripts copied from /exports/igmm/eddie/GenScotDepression/local/EWAS to /exports/igmm/eddie/GenScotDepression/shen/Tools/ewas_pipeline

# Divide Mvalue files by chromosome   ------------------------------------------------------------------

# LBC1921
mkdir /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_w1_DNAm_bychr

Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/ewas_pipeline/split_chr_mvalues.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_w1_DNAm_for_Shen_26Aug2020.rds /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_w1_DNAm_bychr 450k

# LBC1936
mkdir /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC36_w1_DNAm_bychr

Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/ewas_pipeline/split_chr_mvalues.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC36_w1_DNAm_for_Shen_26Aug2020.rds /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC36_w1_DNAm_bychr 450k

# Meta data - create m-values
mkdir /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr

Rscript /gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.LBC_methylation/PREP.combine_LBC21_n_36.R

# Split by chr
Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/ewas_pipeline/split_chr_mvalues.R /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/mval_LBC21_LBC36_w1_DNAm.rds /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr 450k

# new data for PCA: residualised against technical covariates (age+sex+plate+CD8T+CD4T+NK+Bcell+Mono+Gran+cohort)
for CHR in {1..22}
do
 qsub -N pca_${CHR} job.LBC21_n_36_pca.sh $CHR
done

cp -r /exports/eddie/scratch/xshen33/LBC21_n_36_w1_DNAm_bychr/forPCA /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_methylation/LBC21_n_36_w1_DNAm_bychr/
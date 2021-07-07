#$ -N BOLT_gMD
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 12
#$ -o gMD_bolt_o
#$ -e gMD_bolt_e
#$ -cwd

/exports/igmm/eddie/GenScotDepression/local/bin/bolt \
--phenoFile=/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/gMD_forGWAS.tsv \
--phenoCol=gMD \
--covarFile=/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/pheno/bolt_covariates.tsv \
--covarCol=sex \
--covarCol=site \
--qCovarCol=age \
--qCovarCol=pos.x \
--qCovarCol=pos.y \
--qCovarCol=pos.z \
--qCovarCol=pos.table \
--statsFile=/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/gMD.bolt.bed.stats.gz \
--statsFileBgenSnps=/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats/local_GWAS/gMD.bolt.bed.stats.gz \
--bed=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv2/bfile/autosome/autosome.bed \
--bim=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv2/bfile/autosome/autosome.bim \
--fam=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/ukb_chr1_v2.fixCol6.fam \
--remove=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/not_imputed.979.id \
--remove=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/non_whitebritish.30k.id \
--exclude=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/autosome_qc_maf001_geno1_remove.snplist \
--covarCol=centre \
--covarCol=genotyping_array \
--covarMaxLevels=30 \
--qCovarCol=PC{1:20} \
--LDscoresFile=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/gwas/BOLT/genetic_map_hg19_withX.txt.gz \
--lmmForceNonInf \
--numThreads=12 \
--bgenFile=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF=1e-3 \
--bgenMinINFO=0.3 \
--sampleFile=/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3/ukb4844_imp_chr1_v3_s487395.sample \
--verboseStats


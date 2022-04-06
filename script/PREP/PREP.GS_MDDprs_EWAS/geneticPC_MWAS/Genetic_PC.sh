###

cd /exports/eddie/scratch/xshen33/gs_genetic/

plink --bfile /exports/igmm/eddie/GenScotDepression/data/genscot/genetics/genotypes/GS20K_PLINK_files/QCd_data/QCdGS20K --keep /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_genetic_meth/GeneticID_wave1 --make-bed --out GS_genotype_w1

plink --bfile /exports/igmm/eddie/GenScotDepression/data/genscot/genetics/genotypes/GS20K_PLINK_files/QCd_data/QCdGS20K --keep /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_genetic_meth/GeneticID_wave3 --make-bed --out GS_genotype_w3


gcta64 --bfile GS_genotype_w1 --make-grm --autosome --out grm_w1 
gcta64 --bfile GS_genotype_w3 --make-grm --autosome --out grm_w3

gcta64  --grm grm_w1 --pca 20  --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_genetic_meth/Genetic_PC_w1
gcta64  --grm grm_w3 --pca 20  --out /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/GS_genetic_meth/Genetic_PC_w3
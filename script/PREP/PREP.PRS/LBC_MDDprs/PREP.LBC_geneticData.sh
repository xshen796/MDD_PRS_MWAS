#!/bin/sh
#$ -N LBC_PRS
#$ -cwd
#$ -m beas
#$ -M xueyi.shen@ed.ac.uk
#$ -o logs
#$ -e logs
#$ -l h_vmem=32G
#$ -l h_rt=36:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile


# Load necessary modules
module load igmm/apps/vcftools


# Copy data to scratch space

cd $myscratch
cp -r /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/ ./


#################################################################################################################
########################################### QC and convert to plink format ######################################
#################################################################################################################

# LBC1921  ----------------------------------------------------------------------

## Make plink files for each chromosome

cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/

cp ../LBC1921_HRCv1.1_INFO/* ./

for file in $(ls *.vcf.gz)
do
    input_filename=$file
    input_basename=$(basename ${input_filename} .vcf.gz)
    info_basename=$(basename ${input_basename} .dose)

    # Extract the INFO metric from VCF file for ALL variants (This has been done)
    # vcftools --gzvcf ${input_filename} --get-INFO INFO --out ${input_basename}
    # Unzip INFO metrics
    gunzip -c ${info_basename}.info.gz > ${info_basename}.INFO

    # Write a list of variants with INFO metric < 0.3
    awk '$7 < 0.8' ${info_basename}.INFO | awk '$5 < 0.005' ${info_basename}.INFO | cut -f 1 > ${info_basename}_ExcludePositions.txt

    # Call VCFtools to exclude those variants by position
    vcftools --gzvcf ${input_filename} --exclude ${info_basename}_ExcludePositions.txt --recode --out ${input_basename}_clean

    # Filter out multiallelic variants, set to missing calls with <90% posterior probability, then make plink binary file
    plink  --vcf ${input_basename}_clean.recode.vcf --double-id --biallelic-only strict list --vcf-min-gp 0.9 --make-bed --out ${input_basename}

    rm ${info_basename}.INFO ${info_basename}_ExcludePositions.txt ${input_basename}_clean.recode.vcf
done


# LBC1936   --------------------------------------------------------------------------

cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF
cp ../LBC1936_HRCv1.1_INFO/* ./

for file in $(ls *.vcf.gz)
do
    input_filename=$file
    input_basename=$(basename ${input_filename} .vcf.gz)
    info_basename=$(basename ${input_basename} .dose)

    # Extract the INFO metric from VCF file for ALL variants (This has been done)
    # vcftools --gzvcf ${input_filename} --get-INFO INFO --out ${input_basename}
    # Unzip INFO metrics
    gunzip -c ${info_basename}.info.gz > ${info_basename}.INFO

    # Write a list of variants with INFO metric < 0.3
    awk '$7 < 0.8' ${info_basename}.INFO | awk '$5 < 0.005' ${info_basename}.INFO | cut -f 1 > ${info_basename}_ExcludePositions.txt

    # Call VCFtools to exclude those variants by position
    vcftools --gzvcf ${input_filename} --exclude ${info_basename}_ExcludePositions.txt --recode --out ${input_basename}_clean

    # Filter out multiallelic variants, set to missing calls with <90% posterior probability, then make plink binary file
    plink  --vcf ${input_basename}_clean.recode.vcf --double-id --biallelic-only strict list --vcf-min-gp 0.9 --make-bed --out ${input_basename}

    rm ${info_basename}.INFO ${info_basename}_ExcludePositions.txt ${input_basename}_clean.recode.vcf
done


#################################################################################################################
###########################################      Merge 1-22 chrosomes      ######################################
#################################################################################################################

#  LBC1921
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/

for i in {1..22}
do
echo chr$i.dose >> mergelist.txt
done

plink --merge-list mergelist.txt --make-bed --out LBC1921_plink


# LBC1936
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF

for i in {1..22}
do
echo chr$i.dose >> mergelist.txt
done

plink --merge-list mergelist.txt --make-bed --out LBC1936_plink



#################################################################################################################
###########################################          Rename SNPs           ######################################
#################################################################################################################

# LBC1921  --------------------------------------------------------------------------------------------------------------------
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/

# rename SNPs
plink --bfile LBC1921_plink --update-name chr_bp_toRSname.txt --make-bed --out LBC1921_rsID_oldID

## remove duplicate SNPs
#plink --bfile LBC1921_plink --list-duplicate-vars suppress-first --out ls.LBC1921_plink_dupvar


# LBC1936  --------------------------------------------------------------------------------------------------------------------
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/

# rename SNPs
plink --bfile LBC1936_plink --update-name chr_bp_toRSname.txt --make-bed --out LBC1936_rsID

#################################################################################################################
####################################   Update participant ID for LBC 1921    ####################################
#################################################################################################################
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/
plink --bfile LBC1921_rsID_oldID --update-ids /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/subjID_genetic_pheno.txt --make-bed --out LBC1921_rsID


cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/Genotyped/
plink --bfile LBC21_clean_231009_liftover_ACGT_oldID --update-ids /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/subjID_genetic_pheno.txt --make-bed --out LBC21_clean_231009_liftover_ACGT


#################################################################################################################
####################################          Calculate GRM & PC           ######################################
#################################################################################################################

# Merge two cohorts
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC_both/
# genotyped
plink --merge-list mergelist.txt --make-bed --out LBC1921_LBC1936_genotyped_plink
# imputed
plink --merge-list mergelist_imputed.txt --make-bed --out LBC1921_LBC1936_imputed_plink

plink --bfile /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/LBC1921_rsID --exclude LBC1921_LBC1936_imputed_plink-merge.missnp --make-bed --out LBC1921_rsID_flipped
plink --bfile /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/LBC1936_rsID --exclude LBC1921_LBC1936_imputed_plink-merge.missnp --make-bed --out LBC1936_rsID_flipped

echo 'LBC1921_rsID_flipped' > mergelist_imputed_afterflip.txt
echo 'LBC1936_rsID_flipped' >> mergelist_imputed_afterflip.txt

plink --merge-list mergelist_imputed_afterflip.txt --make-bed --out LBC1921_LBC1936_imputed_plink


# Make GRM on genotyped data
gcta64  --bfile LBC1921_LBC1936_genotyped_plink  --autosome  --make-grm  --out LBC21_n_36

# PCA
gcta64  --grm LBC21_n_36 --pca 20  --out LBC21_n_36_PCs



#################################################################################################################
#################################          Remove duplicated var           ######################################
#################################################################################################################

plink2 --bfile LBC1921_LBC1936_genotyped_plink --rm-dup force-first --make-bed --out LBC1921_LBC1936_genotyped_nodupvar
plink2 --bfile LBC1921_LBC1936_imputed_plink --rm-dup force-first --make-bed --out LBC1921_LBC1936_imputed_nodupvar




#################################################################################################################
####################################        Back up processed data         ######################################
#################################################################################################################
cp /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/Genotyped/LBC21_clean_231009_liftover_ACGT* /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC1921/Genotyped/

cp /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/Genotyped/LBC36_clean_231009_liftover_ACGT* /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC1936/Genotyped/

cp /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/LBC1921_rsID* /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF

cp /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/LBC1936_rsID* /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF

mkdir /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC_both
cp /exports/eddie/scratch/xshen33/LBC_genetic/LBC_both/*genotyped* /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC_both/
cp /exports/eddie/scratch/xshen33/LBC_genetic/LBC_both/*imputed* /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC_both/


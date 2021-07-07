while getopts "g:p:o:s:" opt
do
   case "$opt" in
      g ) parameterG="$OPTARG" ;;
      p ) parameterP="$OPTARG" ;;
      o ) parameterO="$OPTARG" ;;
      s ) parameterS="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done




# This script needs to run with R

# Settings ------------------------------------------------------------------------------

prsice_R_file='/exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R'
prsice_binary_file='/exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux'

# space/tab deliminated summary statistics in plain text
gwas_summstats='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice' 

# prefix of the bed,fam,bim files
plink_files=$parameterG 

# Phenotype for MDD diagnosis: space/tab deliminated file. Columns: FID, IID, phenotype (e.g. MDD_diagnosis)
phenotype_file=$parameterP
pheno_name='dummy_MDD'

# prefix for the output files
output_prefix=$parameterO 

# SNPs to include
snp_include=$parameterS


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
    --extract $snp_include \
    --no-clump \
    --bar-levels 1 \
    --fastscore \
    --nonfounders \
    --all-score \
    --out $output_prefix
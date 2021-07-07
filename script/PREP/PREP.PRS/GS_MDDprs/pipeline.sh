
# Prepare imputed genetic data
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/data/
cp GS20K_HRC_0.8_GCTA.* /exports/eddie/scratch/xshen33/GS_PRS_loo/data/

cd /exports/eddie/scratch/xshen33/GS_PRS_loo/data/
plink2 --bfile GS20K_HRC_0.8_GCTA --rm-dup force-first --make-bed --out GS20K_HRC_0.8_GCTA_nodupvar

# Run PRS for pt=5e-8
qsub /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.PRS/GS_MDDprs/job.GS_PRS_5e_08_prsice2.sh


# Create a list of lead variants
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/data/ls_SNP/
R CMD BATCH /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.PRS/GS_MDDprs/PREP.102vars_GS.R

# Run loo PRS
cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.PRS/GS_MDDprs
qsub job.GS_PRS_prsice2_loo.sh

# After the script finishes, go to ANALY folder and run variant EWAS

# Entire MDD summstats, genotyped data   --------------------------------------------------------------------------------------
#bash LBC_PRS_prsice2.sh -g /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/Genotyped/LBC21_clean_231009_liftover_ACGT \
#-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBC1921 \
#-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBC1921.genotyped

#bash LBC_PRS_prsice2.sh -g /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/Genotyped/LBC36_clean_231009_liftover_ACGT \
#-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBC1936 \
#-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBC1936.genotyped


# Entire MDD summstats, imputed data   --------------------------------------------------------------------------------------
#bash LBC_PRS_prsice2_imputed1921.sh -g /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/LBC1921_rsID \
#-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBC1921 \
#-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBC1921.imputed

#bash LBC_PRS_prsice2_imputed1936.sh -g /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/LBC1936_rsID \
#-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBC1936 \
#-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBC1936.imputed

bash LBC_PRS_prsice2_imputed.sh -g /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC_both/LBC1921_LBC1936_imputed_nodupvar \
-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBCboth \
-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBCboth.imputed


# Top 102 hits, imputed data   -----------------------------------------------------------------------------------------------
#bash LBC_PRS_prsice2_topSNP5e_08.sh -g /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/LBC1921_rsID \
#-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBC1921 \
#-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBC1921.imputed.topSNP

#bash LBC_PRS_prsice2_topSNP5e_08.sh -g /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/LBC1936_rsID \
#-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBC1936 \
#-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.LBC1936.imputed.topSNP

bash LBC_PRS_prsice2_topSNP5e_08.sh -g /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC_both/LBC1921_LBC1936_imputed_nodupvar \
-p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/dummy.MDD.LBCboth \
-o /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/LBC_phenotype/MDDprs/MDDprs.imputed.topSNP

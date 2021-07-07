# wiki page: https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=Psychiatry&title=STRADL%20EWAS%20pipeline&fbclid=IwAR13_bpLWyU-p6BXIb8XWpTR2VAxJIzXFg9cjGW8A57HlQBR4e1WxU3DMZI
# Pipeline location:  /exports/igmm/eddie/GenScotDepression/local/EWAS

cd /exports/eddie/scratch/xshen33/GS_PRS_loo/EWAS

while read snp; do

    stradl_ewas --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_forEWAS/MDDprs_for_looEWAS.txt --pheno ${snp} --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_phenotype/covs_RosieData_wave1.rds --no-prune --out ewas_w1_${snp}

    stradl_ewas_w3 --pdata /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_forEWAS/MDDprs_for_looEWAS.txt --pheno ${snp} --cov /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_phenotype/covs_RosieData_wave3.rds --no-prune --out ewas_w3_${snp}


done < ../PRS_forEWAS/MDDprs_name.txt

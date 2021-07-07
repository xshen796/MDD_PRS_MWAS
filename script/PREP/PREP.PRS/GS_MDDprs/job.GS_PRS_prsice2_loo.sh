#!/bin/sh
#$ -N GS_PRS
#$ -cwd
#$ -m beas
#$ -M xueyi.shen@ed.ac.uk
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
. /etc/profile.d/modules.sh
source ~/.bash_profile

#cp /exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA.* /exports/eddie/scratch/xshen33/GS_PRS_loo/data/

#plink2 --bfile /exports/eddie/scratch/xshen33/GS_PRS_loo/data/GS20K_HRC_0.8_GCTA --rm-dup force-first --make-bed --out /exports/eddie/scratch/xshen33/GS_PRS_loo/data/GS20K_HRC_0.8_GCTA_nodupvar

cd /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.PRS/GS_MDDprs
ls /exports/eddie/scratch/xshen33/GS_PRS_loo/data/ls_SNP/ > ls_SNP_file
rm -r mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS/
mkdir /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS/

while read i; do
        
        result_tag=`echo $i | sed "s/.txt//" | sed "s/ls.SNP_//"`

	bash GS_PRS_prsice2_loo.sh -g /exports/eddie/scratch/xshen33/GS_PRS_loo/data/GS20K_HRC_0.8_GCTA_nodupvar \
	  -p /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_loo_PRS/dummy.MDD.GS \
	  -o  /exports/eddie/scratch/xshen33/GS_PRS_loo/PRS/loo_PRS_${result_tag} \
	  -s /exports/eddie/scratch/xshen33/GS_PRS_loo/data/ls_SNP/${i}

done < ls_SNP_file


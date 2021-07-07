module load igmm/apps/metal

cd /exports/eddie/scratch/xshen33/MDDprs_EWAS

for i in { '5e.08' '1e.06' '0.0001' '0.001' '0.01' '0.05' '0.1' '0.5' '1' }
do

  metal MetalScript_pT_EWAS_MDDprs_pT_${i}.txt

done


mkdir /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_meta_LBC_ALSPAC

cp REPmeta* /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_meta_LBC_ALSPAC/

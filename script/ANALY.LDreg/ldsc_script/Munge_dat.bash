##### set up the environment for LDSC : https://github.com/bulik/ldsc
# cd /exports/igmm/eddie/GenScotDepression/shen/LDSC/
# conda env create --file environment.yml
# conda activate ldsc
# sort -rn -k13,13 ampN10.bgenie.QC | head 

cd ../../../data/methQTL_GS/

if [ ! -d 'LDready.stats' ]; then
	 mkdir LDready.stats
fi

for fname in cg*.tsv
do

### prep the stats
## >0.9maximum N, >0.005MAF, info>0.1, phwe>1e-6  --->  for all imaging traits
Rscript /exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/script/ANALY.LDreg/ldsc_script/summstats_formating.R $fname ${fname}_forLDSC


### munge data
munge_sumstats.py \
--sumstats $fname'_forLDSC' \
--out $fname'_forLDSCmunged' \
--merge-alleles /exports/igmm/eddie/GenScotDepression/shen/Tools/LDSC/w_hm3.snplist

echo $fname done
 
done

mv *_forLDSC* LDready.stats/

if [ ! -d 'LDready.stats/logs' ]; then
	 mkdir LDready.stats/logs
fi

mv LDready.stats/*.log LDready.stats/logs
rm tmp.dat.*

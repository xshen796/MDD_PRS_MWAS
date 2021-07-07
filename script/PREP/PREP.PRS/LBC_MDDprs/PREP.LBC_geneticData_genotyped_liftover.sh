cd /exports/eddie/scratch/xshen33/LBC_genetic

# LBC 1921   -------------------------------------------------------------------------------------
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/Genotyped

plink --bfile LBC21_clean_231009 --recode --out LBC21_clean_231009

/exports/igmm/eddie/GenScotDepression/shen/Tools/liftover/liftOverPlink/liftOverPlink.py -m LBC21_clean_231009.map -p LBC21_clean_231009.ped -o LBC21_clean_231009_liftover -e /exports/igmm/eddie/GenScotDepression/shen/Tools/liftover/liftOver -c /exports/igmm/eddie/GenScotDepression/shen/Tools/liftover/hg18ToHg19.over.chain.gz

plink --file LBC21_clean_231009_liftover --make-bed --out LBC21_clean_231009_liftover_plink
plink --bfile LBC21_clean_231009_liftover_plink --alleleACGT --make-bed --out LBC21_clean_231009_liftover_oldID

# LBC 1936 ----------------------------------------------------------------------------------------
cd /exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/Genotyped

plink --bfile LBC36_clean_231009 --recode --out LBC36_clean_231009

/exports/igmm/eddie/GenScotDepression/shen/Tools/liftover/liftOverPlink/liftOverPlink.py -m LBC36_clean_231009.map -p LBC36_clean_231009.ped -o LBC36_clean_231009_liftover -e /exports/igmm/eddie/GenScotDepression/shen/Tools/liftover/liftOver -c /exports/igmm/eddie/GenScotDepression/shen/Tools/liftover/hg18ToHg19.over.chain.gz

plink --file LBC36_clean_231009_liftover --make-bed --out LBC36_clean_231009_liftover_plink
plink --bfile LBC36_clean_231009_liftover_plink --alleleACGT --make-bed --out LBC36_clean_231009_liftover_ACGT


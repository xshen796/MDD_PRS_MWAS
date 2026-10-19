
cd ../../../data/methQTL_GS/

if [ ! -d './LD.results' ]; then
	 mkdir LD.results
fi

fname_traitA='23andmePGCUKBnoMRI_3cohorts.sumstats.gz'
for fname_traitB in LDready.stats/*forLDSCmunged.sumstats.gz
do
	ldsc.py \
	--rg $fname_traitB,$fname_traitA \
	--ref-ld-chr /exports/igmm/eddie/GenScotDepression/shen/Tools/LDSC/eur_w_ld_chr/ \
	--w-ld-chr /exports/igmm/eddie/GenScotDepression/shen/Tools/LDSC/eur_w_ld_chr/ \
	--out $fname_traitB'_MDD_rgresult'

done

mv LDready.stats/*_rgresult* ./LD.results

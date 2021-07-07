plink --bfile /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/1000genomes/1000g_CEU_plink \
          --r2 \
          --ld-snp $1 \
          --ld-window-kb 1000 \
          --ld-window 99999 \
          --ld-window-r2 0 \
          --out $meth/data/ld_snp/$1
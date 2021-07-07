
# Simple command
osca_Linux --beqtl-summary mQTL_w1 --query 1 --out mQTL_summstats_w1
osca_Linux --beqtl-summary mQTL_w2 --query 1 --out mQTL_summstats_w2

####################################      Query a subset of probes      #####################################

cp /exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/mQTL/mQTL* /exports/eddie/scratch/xshen33/mQTL/data/
cd mQTL_byCpG

while read targetprobe; do

osca_Linux --beqtl-summary ../data/mQTL_w1 --query 1 --probe ${targetprobe} --out ${targetprobe}_w1
osca_Linux --beqtl-summary ../data/mQTL_w2 --query 1 --probe ${targetprobe} --out ${targetprobe}_w2

done < ../data/cpg_ls.txt

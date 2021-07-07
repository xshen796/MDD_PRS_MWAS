setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')

alspac.snp=read.delim('data/ALSPAC_genetic/shen-snps-alspac.txt',header=T,stringsAsFactors=F)
gs.snp=read.delim('/exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA.bim',header=F,stringsAsFactors=F)

alspac.snp.GSpresent=alspac.snp$RSID[alspac.snp$RSID %in% gs.snp$V2]
alspac.snp.GSpresent=data.frame(SNPrs=alspac.snp.GSpresent,stringsAsFactors=F)

write.table(alspac.snp.GSpresent,file='data/ALSPAC_genetic/snps-alspac-GSpresent.txt',col.names=T,row.names=F,quote=F)
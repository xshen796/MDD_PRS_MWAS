# Create a list of CpGs for the whole genome with CHR info

setwd("/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD")
library(dplyr)

ls.mvalue.file=list.files(path='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves',pattern='*mvalues450.chr',full.names=T)

for (i in 1:length(ls.mvalue.file)){
  meth.dat.loc=ls.mvalue.file[i]
  chr.dat=read.table(meth.dat.loc,row.names=1,header=T)
  chr.no=gsub('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/mvalues450.chr','',meth.dat.loc)
  chr.no=gsub('.txt','',chr.no)
  chr.no=as.numeric(chr.no)
  tmp.ls=data.frame(chr=chr.no,CpG=rownames(chr.dat),stringsAsFactors=F)
  if (i==1){CpG.ls=tmp.ls}else{CpG.ls=rbind(CpG.ls,tmp.ls)}
  cat(paste0(i,':  chr=',chr.no,'\n'))
}

CpG.ls=CpG.ls[order(CpG.ls$chr),]

saveRDS(CpG.ls,file='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/CpG_list.rds')
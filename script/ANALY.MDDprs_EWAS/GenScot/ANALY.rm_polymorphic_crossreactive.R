library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)   # Data management
library(here)
library(dplyr)
library(data.table)

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- anno %>% as.data.frame %>%
  dplyr::select(ID=Name, CHR=chr, MAPINFO=pos,ucsc_gene=UCSC_RefGene_Name,
         Islands_Name,Relation_to_Island)

ls.CrossReactive.CpG = fread(here::here('data/CpG_exclude/cross_hybridising_EPIC_GS.csv'),stringsAsFactors=F,header=F)
ls.polymorphic.CpG = fread(here::here('data/CpG_exclude/Polymorphic_CpG_EPIC_GS.csv'),stringsAsFactors=F,header=T) %>%
    filter(EUR_AF>0.05,EUR_AF<0.95)

ls.CpG_exclude = c(ls.CrossReactive.CpG$V1,ls.polymorphic.CpG$IlmnID) %>%
      unique

ls.pt=c('0.00000005','0.000001','0.0001','0.001','0.01','0.05','0.1','0.5','1')
for (pt in ls.pt){
    # Load data
    f.path=paste0('result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/ewas_meta_pT_',pt,'.metal.out1.tbl')
    ewasResults=read.delim(here::here(f.path),sep='\t',header=T,stringsAsFactors=F) %>%
      .[!.$MarkerName %in% ls.CpG_exclude,] 

    # Store data
    write.table(ewasResults,file=here::here(f.path),quote=F,sep='\t',col.names=T,row.names=F)
}

# noMHC
f.path=paste0('result/EWAS_MDDprs_Shen/noMHC_and_gwsig/ewas_BothWave_noMHC_5e_081.tbl')
ewasResults=read.delim(here::here(f.path),sep='\t',header=T,stringsAsFactors=F) %>%
      .[!.$MarkerName %in% ls.CpG_exclude,] 

# Store data
write.table(ewasResults,file=here::here(f.path),quote=F,sep='\t',col.names=T,row.names=F)


# top hits PRS
f.path=paste0('result/EWAS_MDDprs_Shen/noMHC_and_gwsig/ewas_BothWave_gwsig1.tbl')
ewasResults=read.delim(here::here(f.path),sep='\t',header=T,stringsAsFactors=F) %>%
  .[!.$MarkerName %in% ls.CpG_exclude,] 

# Store data
write.table(ewasResults,file=here::here(f.path),quote=F,sep='\t',col.names=T,row.names=F)


# unrelated, pt = 5e-8
f.path=paste0('result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/ewas_meta_unrelated_pT_0.00000005.metal.out1.tbl')
ewasResults=read.delim(here::here(f.path),sep='\t',header=T,stringsAsFactors=F) %>%
.[!.$MarkerName %in% ls.CpG_exclude,]

# Store data
write.table(ewasResults,file=here::here(f.path),quote=F,sep='\t',col.names=T,row.names=F)


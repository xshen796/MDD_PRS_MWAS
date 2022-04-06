library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)   # Data management
library(here)
library(dplyr)
library(data.table)
library(pbapply)

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- anno %>% as.data.frame %>%
  dplyr::select(ID=Name, CHR=chr, MAPINFO=pos,ucsc_gene=UCSC_RefGene_Name,
         Islands_Name,Relation_to_Island)

ls.CrossReactive.CpG = fread(here::here('data/CpG_exclude/cross_hybridising_EPIC_GS.csv'),stringsAsFactors=F,header=F)
ls.polymorphic.CpG = fread(here::here('data/CpG_exclude/Polymorphic_CpG_EPIC_GS.csv'),stringsAsFactors=F,header=T) %>%
    filter(EUR_AF>0.05,EUR_AF<0.95)

ls.CpG_exclude = c(ls.CrossReactive.CpG$V1,ls.polymorphic.CpG$IlmnID) %>%
      unique

# Adult replication
ls.f = list.files('result/EWAS_MDDprs_meta_LBC_ALSPAC/',full.name=T,pattern='REP') %>%
   .[!grepl('info',.)]
rm_cpg <- function(f.name){
       ewasResults=read.delim(f.name,sep='\t',header=T,stringsAsFactors=F) %>%
      .[!.$MarkerName %in% ls.CpG_exclude,] 
      write.table(ewasResults,file=f.name,quote=F,sep='\t',col.names=T,row.names=F)
}

ls.f %>% as.list %>%
   pblapply(rm_cpg)

# Adolescent replication
ls.f = list.files('result/EWAS_MDDprs_ALSPAC/',full.name=T,pattern='MDD3_k') %>%
   .[grepl('Rdata',.)]
rm_cpg <- function(f.name){
   load(f.name)
   ewasResults=ewas_res %>% 
      .[!.$probeID %in% ls.CpG_exclude,] 
   write.table(ewasResults,file=gsub('Rdata','csv',f.name),quote=F,sep='\t',col.names=T,row.names=F)
}

ls.f %>% as.list %>%
   pblapply(rm_cpg)
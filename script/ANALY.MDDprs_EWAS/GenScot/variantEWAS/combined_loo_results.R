library(dplyr)
library(data.table)
library(readr)
library(pbapply)
library(tidyverse)
library(superheat)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)   # Data management
setwd('/exports/eddie/scratch/xshen33/GS_PRS_loo/MetaWaves_EWAS')

# Reorganise MWAS summstats ---------------------------------------------------------------------------------------

ls.CrossReactive.CpG = fread(here::here('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/CpG_exclude/cross_hybridising_EPIC_GS.csv'),stringsAsFactors=F,header=F)
ls.polymorphic.CpG = fread(here::here('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/CpG_exclude/Polymorphic_CpG_EPIC_GS.csv'),stringsAsFactors=F,header=T) %>%
    filter(EUR_AF>0.05,EUR_AF<0.95)

ls.CpG_exclude = c(ls.CrossReactive.CpG$V1,ls.polymorphic.CpG$IlmnID) %>%
      unique

load_single_metaRes <- function(tmp.fname){
  variant.name = tmp.fname %>%
    gsub('ewas_BothWaves_','',.) %>%
    gsub('.metal.out1.tbl','',.)
  
  single.summstats = fread(tmp.fname,header = T,stringsAsFactors = F) %>%
    select(CpG=MarkerName, beta=Effect, se=StdErr,p=`P-value`)
  colnames(single.summstats)[2:ncol(single.summstats)]=
    paste0(variant.name,c('.beta','.se','.p'))
  
  return(single.summstats)
}

ls.f=list.files(path = './',pattern = 'metal.out1.tbl$')

all.res = as.list(ls.f) %>%
  pblapply(.,FUN=load_single_metaRes) 

all.res = all.res %>%
  Reduce(function(x,y) left_join(x,y,by='CpG'),.) %>%
  .[!.$CpG %in% ls.CpG_exclude,]


# Write to HPC 

write_tsv(all.res,
          '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_loo/variant_MWAS_Res.tsv',
          quote_escape='none')



# Draw heatmap ----------------------------------------------------------------------------------------------------

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS')
all.res=fread('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_loo/variant_MWAS_Res.tsv')

# position info for CpG sites and SNPs
ls.snp = colnames(all.res) %>% .[grep('.p$',.)] %>% gsub('.p','',.) %>% .[!. %in% 'CpG']
snp.ref = fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',stringsAsFactors=F) %>%
      .[.$V3 %in% ls.snp,]

#load('result/EWAS_MDDprs_fromKathryn/cpg_snp_mapping_pNbeta.RData')
ewasResult=read.delim('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/ewas_meta_pT_0.00000005.metal.out1.tbl')
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- anno %>% as.data.frame %>%
      dplyr::select(ID=Name, CHR=chr, MAPINFO=pos,ucsc_gene=UCSC_RefGene_Name,
                    Islands_Name,Relation_to_Island)
res.sig = filter(ewasResult,p.adjust(P.value,method='bonferroni')<0.05) %>%
      merge(.,ref.tomerge,by.x='MarkerName',by.y='ID',all.x=T) %>%
      mutate(is.MHC = if_else(MAPINFO>25000000&MAPINFO<35000000&CHR=='chr6','Yes','No')) %>%
      filter(.,is.MHC=='No'|P.value==1.78e-117)
p.mat=all.res %>%
  select(CpG, ends_with('.p')) %>%
  .[.$CpG %in% res.sig$MarkerName,] %>%
  column_to_rownames(var='CpG') %>%
  as.data.frame %>%
  mutate_all(., function(x) as.numeric(x))



# Reorder CpG by CHR
ls.CpG = res.sig %>%
  select(ID=MarkerName,CHR,MAPINFO) %>%
  mutate(.,CHR=gsub('chr','',CHR)) %>%
  mutate(.,CHR=as.numeric(CHR)) %>%
  .[order(.$CHR,.$MAPINFO,decreasing = F),] %>%
  mutate(ID=as.character(ID)) %>%
  .[.$ID %in% rownames(p.mat),] 

# Reorder SNP by CHR
colnames(snp.ref)=c('CHR','BP','RSID')
ls.SNP.102 = snp.ref %>%
  mutate(.,CHR=gsub('chr','',CHR)) %>%
  mutate(.,CHR=as.numeric(CHR)) %>%
  .[order(.$CHR,.$BP,decreasing = F),] %>%
  mutate(RSID=as.character(RSID)) %>%
  .[.$RSID %in% gsub('.p','',colnames(p.mat)),]

# Reorder p.mat and reduce size
colnames(p.mat)=colnames(p.mat) %>% gsub('.p','',.)
p.mat[p.mat==0]=1e-320  # Recover extremely low p-values to system threshold
p.mat[p.mat>0.05/102]=NA
p.mat.new = -log10(p.mat) %>%
  .[ls.CpG$ID,] %>%
  .[,ls.SNP.102$RSID]


ls.indep.cpg = ls.CpG[ls.CpG$ID %in% ls.cpg,]
p.mat.new=p.mat.new[ls.indep.cpg$ID,]

# Heatmap
superheat(p.mat.new,
          # row title
          row.title = "CpG",
          row.title.size = 6,
          # col title
          column.title = "SNP",
          column.title.size = 6,
          # scale the matrix columns
          scale = F,
          heat.col.scheme = "red",
          heat.na.col = "white",
          # cluster by gears
          membership.rows = paste('CHR',ls.CpG$CHR),
          membership.cols = paste('CHR',ls.SNP.102$CHR),
          # Text
          bottom.label.text.angle = 90,
          left.label.text.size = 4,
          bottom.label.text.size = 4)


k=rowSums(p.mat.new < -log10(0.05/525),na.rm=T)
process.corr = p.mat.new %>% mutate(Ninsig=k) %>%
  tibble::rownames_to_column(.,'CpG') %>%
  merge(.,ewasResult,by.x='CpG',by.y='MarkerName',all.x=T)

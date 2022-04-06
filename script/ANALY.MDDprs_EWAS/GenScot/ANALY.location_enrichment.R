setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/')

library(dplyr)
library(tibble)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)   # Data management
library(readr)

# Load data ---------------------------------------------------------------
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- anno %>% as.data.frame %>%
  dplyr::select(ID=Name, CHR=chr, MAPINFO=pos,ucsc_gene=UCSC_RefGene_Name,
                Islands_Name,Relation_to_Island)
pt='0.00000005'
f.path=paste0('result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/ewas_meta_pT_',pt,'.metal.out1.tbl')
ewasResults=read.delim(here::here(f.path),sep='\t',header=T,stringsAsFactors=F) %>%
  mutate(p.adj = p.adjust(P.value,method='bonferroni')) 
colnames(ewasResults)[1]='CpG'
# Add annotations
fig.dat.GS.pT.0.00000005=merge(ewasResults,ref.tomerge,by.x='CpG',by.y='ID',all.x=T)


# Summarise count data ----------------------------------------------------
location.counts = fig.dat.GS.pT.0.00000005 %>% 
  mutate(is.sig=p.adj<0.05) %>% 
  dplyr::count(is.sig,Relation_to_Island)

location.matrix = data.frame(location.counts %>% 
                               filter(is.sig==T) %>% 
                               .[,c('Relation_to_Island','n')],
                             location.counts %>% 
                               filter(is.sig==F) %>% 
                               .[,c('n')]) %>% 
  tibble::column_to_rownames(var='Relation_to_Island')
colnames(location.matrix)=c('sig','non_sig')
location.matrix.percentage = location.matrix %>%
  mutate(sig=sig/sum(sig),non_sig=non_sig/sum(non_sig))

test.shore = data.frame(sig=c(sum(location.matrix$sig[c(3,6)]),
                              sum(location.matrix$sig[c(1,2,4,5)])),
                        non_sig=c(sum(location.matrix$non_sig[c(3,6)]),
                                  sum(location.matrix$non_sig[c(1,2,4,5)])))
rownames(test.shore)=c('Shore','non_shore')


# chi-squre  tests --------------------------------------------------------

c_bytype <- function(tmp.relation,input.dat){
  new.dat=input.dat[tmp.relation,]
  other.rows = input.dat %>% .[!rownames(input.dat) %in% tmp.relation,] %>% 
    colSums
  new.dat=rbind(new.dat,other.rows)
  res=chisq.test(new.dat) 
  output=c(res$statistic,res$parameter,res$p.value) %>% data.frame
  return(output)
}

all.res = rownames(location.matrix) %>% as.list %>% 
  lapply(.,FUN=c_bytype,input.dat=location.matrix) %>% 
  bind_cols %>% 
  t %>% data.frame
colnames(all.res)[3]='p.value' 
rownames(all.res)=NULL

table.inpaper = data.frame(location.matrix.percentage %>% 
                             rownames_to_column(var='Relation_to_Island'),
                           all.res)
write.table(table.inpaper,file=here::here('result/location_mat.txt'),quote = F,row.names = T,col.names = T,sep = '\t')



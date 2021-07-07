library(dplyr)
library(data.table)
library(purrr)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(here)
library(pbapply)
source(here::here('MR_meth_MDD/script/FUNs_meth/clump_bychr.R')) 


setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')

# Load result matrix ------------------------------------------------------

all.res = fread('result/EWAS_loo/variant_MWAS_Res.tsv')
p.mat = all.res %>% column_to_rownames(var="CpG") 

# position info for CpG sites and SNPs
ls.snp = colnames(all.res) %>% .[grep('.p$',.)] %>% gsub('.p','',.) %>% .[!. %in% 'CpG']
snp.ref = fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',stringsAsFactors=F) %>%
  .[.$V3 %in% ls.snp,]

ewasResult=read.delim('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/ewas_meta_pT_0.00000005.metal.out1.tbl')
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ref.tomerge <- anno %>% as.data.frame %>%
  dplyr::select(ID=Name, CHR=chr, MAPINFO=pos,ucsc_gene=UCSC_RefGene_Name,
                Islands_Name,Relation_to_Island)
res.sig = filter(ewasResult,p.adjust(P.value,method='bonferroni')<0.05) %>%
  merge(.,ref.tomerge,by.x='MarkerName',by.y='ID',all.x=T) %>%
  mutate(is.MHC = if_else(MAPINFO>25000000&MAPINFO<35000000&CHR=='chr6','Yes','No')) %>%
  filter(.,is.MHC=='No'|P.value==1.78e-117)
p.mat=all.res %>%
  select(CpG, ends_with('.p')) %>%
  column_to_rownames(var='CpG') %>%
  as.data.frame %>%
  mutate_all(., function(x) as.numeric(x))
ls.CpG.filter = rowSums(p.mat[,2:ncol(p.mat)]<0.05/769525)
ls.cpg.tokeep = rownames(p.mat)[ls.CpG.filter>0]

# Reorder CpG by CHR
ewasResult = ewasResult %>%
  left_join(.,ref.tomerge,by=c('MarkerName'='ID')) %>% 
  select(ID=MarkerName,CHR,MAPINFO,P.value) %>%
  mutate(.,CHR=gsub('chr','',CHR)) %>%
  mutate(.,CHR=as.numeric(CHR)) 
ls.CpG = ewasResult %>%
  .[order(.$CHR,.$MAPINFO,decreasing = F),] %>%
  mutate(ID=as.character(ID)) %>%
  .[.$ID %in% ls.cpg.tokeep,] %>% 
  mutate(BP=MAPINFO)


# Reorder SNP by CHR
colnames(snp.ref)=c('CHR','BP','RSID')
ls.SNP.102 = snp.ref %>%
  mutate(.,CHR=gsub('chr','',CHR)) %>%
  mutate(.,CHR=as.numeric(CHR)) %>%
  .[order(.$CHR,.$BP,decreasing = F),] %>%
  mutate(RSID=as.character(RSID)) %>%
  .[.$RSID %in% gsub('.p','',colnames(p.mat)),] %>%
  .[.$CHR %in% unique(ls.CpG$CHR),]

# Lowest p-value for each SNP within 1Mbp
find_pval <- function(tmp.input,tmp.ewas_res){
  tmp.loc.chr = tmp.input$CHR
  tmp.loc.bp = tmp.input$BP
  block.res = tmp.ewas_res %>% filter(CHR==tmp.loc.chr) %>%
    filter((MAPINFO<(tmp.loc.bp+1000000))&(MAPINFO>(tmp.loc.bp-1000000)))
  tmp.p = min(block.res$P.value,na.rm=T)
  return(tmp.p)
}

input.p = ls.SNP.102 %>% transpose() %>%
  pblapply(.,find_pval,tmp.ewas_res=ewasResult) %>%
  unlist 
input.p[input.p>0.05/nrow(ewasResult)]=NA

cpg.input.p = ls.CpG %>% transpose() %>%
  pblapply(.,find_pval,tmp.ewas_res=ewasResult) %>%
  unlist 
cpg.input.p[cpg.input.p>0.05/nrow(ewasResult)]=NA

# Reorder p.mat and reduce size
colnames(p.mat)=colnames(p.mat) %>% gsub('.p','',.)
p.mat[p.mat==0]=1e-320  # Recover extremely low p-values to system threshold
p.mat[p.mat>0.05/96]=NA
p.mat.new = -log10(p.mat) %>%
  .[ls.CpG$ID,] %>%
  .[,ls.SNP.102$RSID]


# Calc number of trans CpGs associated with each SNP ----------------------
# result matrix: row = SNP, col = CpG
calc_trans <- function(tmp.snp.name,ref.snp,ref.cpg,tmp.mat){
  # Select results for the target SNP
  tmp.res = tmp.mat[tmp.snp.name,] %>% .[!is.na(.)] %>% 
    data.frame(log.p = .,CpG=names(.)) %>%
    # Add annotations
    left_join(.,ref.cpg,by=c('CpG'='ID')) %>% 
    select( -BP)
  # Local chromosome
  target.chr = ref.snp$CHR[ref.snp$RSID==tmp.snp.name]
  # Local BP
  target.bp = ref.snp$BP[ref.snp$RSID==tmp.snp.name] 
  
  # Summarise distal effects
  tmp.res = tmp.res %>% 
    mutate(trans.chr=ifelse(CHR!=target.chr,1,0)) %>% 
    mutate(trans.bp=ifelse(between(MAPINFO,target.bp-1000000,target.bp+1000000),0,1)) %>% 
    mutate(local=ifelse(CHR!=target.chr,0,
                        ifelse(between(MAPINFO,target.bp-1000000,target.bp+1000000),1,0)))
  trans.summstats = filter(tmp.res,trans.chr==1) 
  chrn_n.trans = unique(trans.summstats$CHR) %>% length
  #### N of trans effects in a distal chromosome
  output = data.frame(n.trans_chr_loci = as.vector(table(tmp.res$trans.chr)['1']),
                      n.trans_chr_chrn = chrn_n.trans,
                      n.trans_bp = as.vector(table(tmp.res$trans.bp)['1']),
                      n.cis=as.vector(table(tmp.res$local)['1']))
  return(output) 
}

summ.trans_n = ls.SNP.102$RSID %>% as.list %>% 
  pblapply(.,calc_trans,ref.snp=ls.SNP.102,ref.cpg=ls.CpG,tmp.mat=t(p.mat.new)) %>% 
  bind_rows %>% 
  data.frame(.,ls.SNP.102) %>% 
  .[order(.$n.trans_chr_loci,decreasing = T),] 

write.table(summ.trans_n,file='result/Trans_summary.txt',col.names = T,row.names = F,quote = F,sep = '\t')


# Save another version for FUMA -------------------------------------------

mddgwas.summstats=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice')

mddgwas.summstats = mddgwas.summstats %>% 
  .[.$SNP %in% summ.trans_n$RSID[summ.trans_n$n.trans_chr_chrn>0],]

write.table(mddgwas.summstats,file='data/TransEffect_MDD_gwas.txt',col.names = T,row.names = F,quote = F,sep = '\t')

# Aims:
#   1. Clump the HMC region to reduce the number CpGs tested
#   2. Retain the rest of significant CpGs to enter analysis
# Before the above two steps, the entire list of CpGs should be filtered with
# availability of mQTL data.

# Basic settings ----------------------------------------------------------

setwd('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GoDMC/mQTL_meta_nopgc/bychr')

library(dplyr)
library(data.table)
library(pbapply)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stringr)


# Create tables for no. of instruments per CpG ----------------------------

ls.mqtl=list.files(path = '.',pattern='^mQTL')

count_instrument <- function(fname,pT){
  dat = readRDS(fname) %>% filter(pval<pT)
  chr.name = gsub('mQTL.chr','',fname) %>% gsub('.rds','',.) %>% as.numeric
  output = dplyr::count(dat,cpg) %>% mutate(chr=chr.name)
}

summary.pT5e_08.instrument = pblapply(as.list(ls.mqtl),FUN = count_instrument,pT=5e-8) %>% bind_rows
summary.pT1e_06.instrument = pblapply(as.list(ls.mqtl),FUN = count_instrument,pT=1e-6) %>% bind_rows

summary.instrument = merge(summary.pT1e_06.instrument,summary.pT5e_08.instrument[,c('cpg','n')],by='cpg',all.x=T)
summary.instrument = summary.instrument[,c(1,3,2,4)]
colnames(summary.instrument)=c('cpg','chr','n.5e_08','n.1e_06')

setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')
#saveRDS(summary.instrument,file='data/mQTL/GoDMC.instrument.rds')


# Find CpGs for MR --------------------------------------------------------
# Basic settings
setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

# Define a function for clumping
clump_bychr <- function(dat,chr.no,tophits.bpwindow=250000,ref.cpgCorr){
  sig.rows=dat[dat$chr==chr.no,]
  if (nrow(sig.rows)>0){
    sig.rows=sig.rows[order(sig.rows$BP),]
    
    sig.rows$sigblock.no=99999
    
    for (i in 1:nrow(sig.rows)){
      if (i<=3){
        sig.rows$sigblock.no[i]=1
      }else{
        pair.r = abs(ref.cpgCorr[sig.rows$cpg[(i-3):i],sig.rows$cpg[(i-3):i]]) %>% 
          .[2:nrow(.),ncol(.)] # Correlation matrix of CpG(i) with three previous CpGs
        if ((sig.rows$BP[i]>sig.rows$BP[i-1]+tophits.bpwindow)|(sum(pair.r<0.3)==length(pair.r))){
          sig.rows$sigblock.no[i]=sig.rows$sigblock.no[i-1]+1
        }else{
          sig.rows$sigblock.no[i]=sig.rows$sigblock.no[i-1]
        }
      }
    }
    sig.rows$sigblock.no[sig.rows$sigblock.no==99999]=NA
    
    for (i in unique(sig.rows$sigblock.no)){
      tmp.block=filter(sig.rows,sigblock.no==i)
      tmp.top=tmp.block[tmp.block$p==min(tmp.block$p),]
      if (i==unique(sig.rows$sigblock.no)[1]){
        top.hits=tmp.top
      }else{
        top.hits=rbind(top.hits,tmp.top)
      }
    }
    
  }
  return(top.hits)
}


# Load summary info for GoDMC instruments
summary.instrument=readRDS('data/mQTL/GoDMC.instrument.rds')

# Load EWAS summary stats
ewas.summstat.tophitPRS=fread('result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.metal.out1.tbl',
                              header=T,stringsAsFactors=F)
ewas.summstat.tophitPRS = ewas.summstat.tophitPRS %>%
  mutate(p.adj=p.adjust(`P-value`,method = 'bonferroni'))
ls.cpg.ewas.sig = filter(ewas.summstat.tophitPRS,p.adj<0.05)

# In MHC region?
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- data.frame(ID=anno$Name, CHR=anno$chr, MAPINFO=anno$pos)
ls.cpg.ewas.sig = ls.cpg.ewas.sig %>%
  merge(.,ref.tomerge,by.x='MarkerName',by.y='ID',all.x=T) %>%
  mutate(is.MHC = if_else(MAPINFO>25000000&MAPINFO<35000000&CHR=='chr6','Yes','No')) %>%
  mutate(is.MHC = factor(is.MHC,levels=c('No','Yes')))

# Find suitable CpGs for MR (enough genetic instruments)
ls.cpg.sig_instru=ls.cpg.ewas.sig[ls.cpg.ewas.sig$MarkerName %in% summary.instrument$cpg[summary.instrument$n.5e_08>10],]

# Prepare for clumping
## Annotate CpGs and rename dataframe for clumping
anno.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno.tomerge = anno.450k[,c('chr','pos','Name')]
target.rows = merge(ls.cpg.ewas.sig,anno.tomerge,by.x='MarkerName',by.y='Name') %>%
  data.frame %>%
  dplyr::rename(cpg = MarkerName, p = .$`P-value` , BP = pos) %>%
  dplyr::rename(p = P.value) %>%
  filter(.,is.MHC=='Yes')
  
## Load CpG correlation data for clumping
ref.cpg_cormat = readRDS('data/Meth_Cor_Mat/cormat_MDDprs_tophits_CpG.rds')

# Clump CpGs significant in EWAS and available in the GoDMC data
cpg.forMR.MHC = as.list(unique(target.rows$chr)) %>%
  pblapply(.,FUN=clump_bychr,
           dat=target.rows,ref.cpgCorr = ref.cpg_cormat,tophits.bpwindow=1000000) %>%
  bind_rows %>%
  .[order(.$p.adj,decreasing = F),] 

# Check how many instruments they have
cpg.forMR.MHC = cpg.forMR.MHC %>% 
  merge(.,summary.instrument[,c('cpg','n.5e_08')],by='cpg',all.x=T) %>%
  select(-chr,-BP,-sigblock.no) %>%
  dplyr::rename(`P-value`=p)
# All non-MHC CpGs
cpg.forMR.all = merge(ls.cpg.sig_instru,summary.instrument[,c('cpg','n.5e_08')],by.x='MarkerName',by.y='cpg',all.x=T) %>%
  dplyr::rename(cpg=MarkerName) %>%
  filter(.,is.MHC=='No') %>%
  rbind(cpg.forMR.MHC) %>%
  dplyr::rename(chr=CHR)

# Load mQTL data for these CpGs with h19 annotation -----------------------
# Define function for extracting data
extract_mqtl <- function(target.dir,chr.no,target.ls,ref.snpName){
  ls.cpg = target.ls$cpg[target.ls$chr==chr.no]
  fname = paste0(target.dir,'/mQTL.',chr.no,'.rds')
  tmp.mqtl = readRDS(fname)
  
  target.mqtl = tmp.mqtl[tmp.mqtl$cpg %in% ls.cpg,] %>%
    data.frame(.,str_split(.$snp, ":", simplify = TRUE)) %>%
    mutate(gp.oldName=paste0(X1,':',X2)) %>%
    merge(ref.snpName,.,by='gp.oldName',all.y=T) %>%
    select(-matches('^X|V1|V2')) %>%
    dplyr::rename(SNP=V3)
  return(target.mqtl)
}

save_mqtl_forMR <- function(masterdat,dir.tosave,target.cpg){
  output = masterdat[masterdat$cpg==target.cpg,] %>%
     dplyr::filter(!is.na(SNP)) %>%
    .[!duplicated(.$SNP),]
  saveRDS(output,file = paste0(dir.tosave,'/mQTL.',target.cpg,'.rds'))
}

rs.chr_bp=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',
                header=F,sep=' ',stringsAsFactors=F)
rs.chr_bp$gp.oldName = paste0(rs.chr_bp$V1,':',rs.chr_bp$V2)
mqtl.dir='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GoDMC/mQTL_meta_nopgc/bychr'

mqtl.dat = as.list(unique(cpg.forMR.all$chr)) %>%
  pblapply(.,FUN=extract_mqtl,target.dir=mqtl.dir,target.ls=cpg.forMR.all,ref.snpName=rs.chr_bp) %>%
  bind_rows

as.list(unique(cpg.forMR.all$cpg)) %>%
  pblapply(.,FUN=save_mqtl_forMR,masterdat=mqtl.dat,dir.tosave='data/mQTL/mQTL_forMR/')


# Find a list of SNPs in linkage with maximum number of CpGs >23 ----------

nodup.mqtl.dat=mqtl.dat %>% mutate(check=paste0(SNP,':',cpg))
nodup.mqtl.dat=nodup.mqtl.dat[!is.na(nodup.mqtl.dat$SNP),]
nodup.mqtl.dat=nodup.mqtl.dat[!duplicated(nodup.mqtl.dat$check),]
overlap.snp = nodup.mqtl.dat %>% dplyr::count(.,SNP)
table(overlap.snp$n)

multivarMR.SNP = filter(overlap.snp,n>=23)
saveRDS(multivarMR.SNP,file='data/mQTL/mQTL_forMR/SNP_for_multivarMR.rds')


# Find a list of SNPs and CpGs for multi-variable MR ----------------------

# Process exposure data ---------------------------------------------------

# CpGs to include
res = as.list(paste0('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/GoDMC_MR_DNAm_to_MDD/Summary_DNAm_to_MDD.csv')) %>% 
  lapply(.,FUN=read.csv,header=T,sep='\t',stringsAsFactors = F) %>%
  .[[1]]

ls.cpg.multivarMR = res %>%
  group_by(exposure) %>%
  count(pval<0.05) %>%
  filter(`pval < 0.05`==T) %>%
  filter(n>=2)

write.table(ls.cpg.multivarMR$exposure,
            file='script/ANALY.MR/DNAm_to_MDD/ANALY_multivar/ls.cpg.multivar_MR.txt',
            sep = '\n',quote=F,col.names=F,row.names=F)

# SNPs to include
# Run exposure preparation for single-variant MR first
exposure.path = '/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/MR_InterFiles/exposure_stats/'
exposure.dat = list.files(path = exposure.path, pattern = '.exposure_dat',full.names=T) %>%
  .[grep(paste0(ls.cpg.multivarMR$exposure,collapse = '|'),.)] %>%
  as.list %>%
  lapply(.,FUN = fread,header=T) %>%
  bind_rows

ls.exposure.snp = unique(exposure.dat$SNP)

write.table(ls.exposure.snp,
            file='script/ANALY.MR/DNAm_to_MDD/ANALY_multivar/ls.snp.multivar_MR.txt',
            sep = '\n',quote=F,col.names=F,row.names=F)
library(dplyr)              # Data management 
library(data.table)         # Data management
library(purrr)              # Data management
library(readr)              # Data management
library(pbapply)            # Process management

setwd('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats')

# Load ref data
ref=fread('variants.txt',header = T,stringsAsFactors = F) %>%
  .[,c('rsid','af','info')]

# Reformat gFA and gMD summstats ------------------------------------------

ls.measure = c('gFA','gMD','gFA_TR','gMD_TR','FA_ATR','MD_ATR')

read_summstats <- function(tmp.measure){
  
  summstats = read.delim(paste0('local_GWAS/',tmp.measure,'.bolt.bed.stats.gz'),
                         stringsAsFactors = F)
  summstats.reformat = summstats %>%
    select(SNP,A1=ALLELE1,A2=ALLELE0,BETA,SE,P=P_BOLT_LMM) %>%
    mutate(Direction=NA) %>%
    merge(.,ref,by.x='SNP',by.y='rsid')
  
}

summstats.loaded = as.list(ls.measure) %>%
  pblapply(.,read_summstats)

names(summstats.loaded)=ls.measure

# Save data ---------------------------------------------------------------

pbmapply(write.table, summstats.loaded, 
         file=paste0('Processed/',names(summstats.loaded),'_BOLT_40k.InfoAdded.tbl'),
         col.names=T,row.names=F,sep='\t',quote=F)

# Rename data for MR analysis ---------------------------------------------

ls.summstats = list.files('Processed/',pattern = '.InfoAdded.tbl$',full.names = T)
ls.new.summstats = ls.summstats %>%
  gsub('aparc-Desikan_area_TotalSurface.metal','TotalSurface',.) %>%
  gsub('_BOLT_40k','',.) %>%
  gsub('IDP_T1_SIENAX_grey_normalised_volume','TotalVol',.) %>%
  gsub('aparc-Desikan_thickness_GlobalMeanThickness.metal','MeanThk',.)

for (i in 1:length(ls.summstats)){
  old.n = ls.summstats[i]
  new.n = ls.new.summstats[i]
  system(paste('mv',old.n,new.n))
}

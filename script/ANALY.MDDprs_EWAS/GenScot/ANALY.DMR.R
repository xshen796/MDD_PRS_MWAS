library('dplyr')
library('data.table')
library('pbapply')
library('dmrff')
library('tidyverse')
library('readr')
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)   # Data management

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')

# Load MWAS summstats -----------------------------------------------------

# Annotation file
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- anno %>% as.data.frame %>%
  dplyr::select(ID=Name, chr, pos,ucsc_gene=UCSC_RefGene_Name,
                Islands_Name,Relation_to_Island)

# Load ewas summstats and add annotations
ewas.res = fread('result/EWAS_MDDprs_Shen/archiv/MDDprs_ewas_meta/MDDprs_pT5e_08_meta/mddprs_pT5e_08_meta_RosieData.toptable.txt1.tbl') %>% 
  left_join(.,ref.tomerge,by=c('MarkerName'='ID')) %>% as.data.frame
rownames(ewas.res) = ewas.res$MarkerName

# Load methylation matrix

# wave 1

DNAm.w1 = list.files('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/chr/',
                  pattern='mvalues.chr',full.names = T) %>%
      as.list %>%
      pblapply(.,readRDS) %>%
      pblapply(.,function(x) x[,1:4000]) %>% 
      bind_rows %>%
      rownames_to_column('CpG') 
      

# DNAm.w3 = list.files('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/chr/',
#                    pattern='mvalues.chr',full.names = T) %>%
#        as.list %>%
#        pblapply(.,readRDS) %>%
#        pblapply(.,as.data.frame) %>% 
#        bind_rows%>%
#        rownames_to_column('CpG')

# DNAm = merge(DNAm.w1,DNAm.w3,by='CpG',all=T) %>%
DNAm = DNAm.w1 %>% 
       column_to_rownames(var="CpG") 


rm(DNAm.w1,DNAm.w3)

common <- intersect(rownames(DNAm),rownames(ewas.res))

DNAm = DNAm[common,] %>% as.matrix
ewas.res = ewas.res[common,] 
  

dmrs <- 
        dmrff(estimate=ewas.res$Effect,
              se=ewas.res$StdErr,
              p.value=ewas.res$`P-value`,
              methylation=DNAm,
              chr=ewas.res$chr,
              pos=ewas.res$pos,
              maxgap=500,
              verbose=T)


# Save results ----------------------------------------------------------------------------------------------------

write_tsv(dmrs,'result/EWAS_MDDprs_Shen/DMR_res.tsv')

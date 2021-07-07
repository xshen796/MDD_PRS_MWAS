library(ggplot2)
library(dplyr)
library(readr)
library(data.table)

source('/exports/igmm/eddie/GenScotDepression/shen/Tools/LocusZooms/functions/locus_zoom.R')

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')


# Load data ---------------------------------------------------------------

# New summstats that overlaps with GS variants
GS.SNP = fread('/exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA.bim')
GS.summstats=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice',
                        stringsAsFactors = F,header=T,sep='\t') %>% 
  .[.$SNP %in% GS.SNP$V2,] %>% 
  .[!duplicated(.$SNP),]

# gene list
UCSC_GRCh37_Genes_UniqueList.txt <- 
  read.delim("/exports/igmm/eddie/GenScotDepression/shen/Tools/LocusZooms/UCSC_GRCh37_Genes_UniqueList.txt", stringsAsFactors = FALSE, header = TRUE)


# Create LD  --------------------------------------------------------------
# Select SNP with global effect
trans_summary = fread('result/Trans_summary.txt',header=T) %>% 
  filter(n.trans_chr_chrn>21/2)
write.table(trans_summary$RSID,'script/Figs/FIG.zoomlocus/ls.snp',sep='\n',col.names = F,row.names = F,quote=F)

# Calculate LD
system('bash /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/Figs/FIG.zoomlocus/ld_snp.sh')
# ld files will be stored in data/ld_snp/rsIDxxx.ld


# Make zoom locus plots ---------------------------------------------------

# load necessary files into R
Example.assoc.linear <- GS.summstats

for (i in 1:nrow(trans_summary)){
  Example.ld <- read.table(paste0("data/ld_snp/",trans_summary$RSID[i],'.ld'), stringsAsFactors = FALSE, header = TRUE)
  
  
  # create a LocusZoom-like plot
  locus.zoom(data = Example.assoc.linear,                                    # a data.frame (or a list of data.frames) with the columns CHR, BP, SNP, and P
             snp = trans_summary$RSID[i],
             region = c(trans_summary$CHR[i], min(Example.ld$BP_B), max(Example.ld$BP_B)),                             # the chromosome region to be included in the plot
             offset_bp = 0,                                                  # how many basepairs around the SNP / gene / region of interest to plot
             ld.file = Example.ld,                                           # a file with LD values relevant to the SNP specified above
             genes.data = UCSC_GRCh37_Genes_UniqueList.txt,                  # a file of all the genes in the region / genome
             plot.title = trans_summary$RSID[i],        # the plot title
             file.name = paste0('Figs/zoomlocus/',trans_summary$RSID[i],'.jpg'),                        # the name of the file to save the plot to
             secondary.label = F)                                         # TRUE/FALSE whether to add rsIDs of secondary SNPs to plot
}

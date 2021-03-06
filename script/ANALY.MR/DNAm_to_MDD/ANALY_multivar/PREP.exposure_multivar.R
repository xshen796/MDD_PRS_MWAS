library(dplyr)
library(data.table)
library(TwoSampleMR)
library(optparse)
library(pbapply)
library(readr)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1

parse <- OptionParser()

option_list <- list(
   make_option('--SNPlist', type='character', help="SNP-CpG list to obtain", action='store'),
   make_option('--CpGlist', type='character', help="SNP-CpG list to obtain", action='store'),
   make_option('--mqtl', type='character', help="mQTL files in rds format", action='store'),
   make_option('--interexp', type='character', 
               help='Folder for intermediate files generated by processing exposure data', action='store')
)

args = commandArgs(trailingOnly=TRUE)

opt <- parse_args(OptionParser(option_list=option_list), args=args)
f.snp=opt$SNPlist
f.cpg=opt$CpGlist
exposure.path=opt$mqtl
output.path=opt$interexp

# f.snp='ls.snp.multivar_MR.txt'
# f.cpg='ls.cpg.multivar_MR.txt'
# exposure.path='/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/mQTL/mQTL_forMR_GS/'
# output.path='/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/data/MR_InterFiles/exposure_stats/'

# Basic settings ----------------------------------------------------------

dir.create(file.path(output.path), showWarnings = FALSE)

# Define functions --------------------------------------------------------
# Extract summstats
extract_summstats <- function(tmp.input){
   summ.stats.fname = as.character(tmp.input['fname'])
   traitA = as.character(tmp.input['exposure.name'])
   gwas_summarystat = readRDS(summ.stats.fname) %>% 
      dplyr::select(cpg,SNP,beta_a1,se,pval,samplesize,allele1,allele2,freq_a1) %>% 
      .[.$SNP %in% ls.snp,]
   colnames(gwas_summarystat) = c('CpG','SNP','Effect','se','p','N','Allele1','Allele2','Freq1')
   SigniSNP=gwas_summarystat
   
   SigniSNP$NewAllele1=SigniSNP$Allele1
   SigniSNP$NewAllele2=SigniSNP$Allele2
   SigniSNP$NewFreq1=SigniSNP$Freq1
   SigniSNP$NewEffect=SigniSNP$Effect
   
   lst=c(which(SigniSNP$Effect<0))
   SigniSNP$NewEffect[lst]=-1*SigniSNP$Effect[lst]
   SigniSNP$NewFreq1[lst]=1-1*SigniSNP$Freq1[lst]
   SigniSNP$NewAllele1[lst]=SigniSNP$Allele2[lst]
   SigniSNP$NewAllele2[lst]=SigniSNP$Allele1[lst]
   SigniSNP=data.frame(SigniSNP)
   SigniSNP=SigniSNP[,c('SNP','NewAllele1','NewAllele2','NewFreq1','NewEffect','se','p','N')]
   
   colnames(SigniSNP)=c("SNP","effect_allele","other_allele","eaf", "beta","se","pval",'samplesize')
   SigniSNP=data.frame(Phenotype=traitA,units='unit',SigniSNP)
   
   return(SigniSNP)
}

# Reformat data
format_exposure <- function(tmp.input,tmp.SNP.toinclude){
   exposure.fname = as.character(tmp.input['fname'])
   traitA = as.character(tmp.input['exposure.name'])
   
   exposure_dat <- read_exposure_data(exposure.fname,sep='\t') %>%
      .[.$SNP %in% tmp.SNP.toinclude,]
   
   #exposure_dat <- filter(exposure_dat,eaf.exposure>0.1,eaf.exposure<0.90,abs(beta.exposure)<0.038)
   exposure_dat <- exposure_dat[,c('SNP','effect_allele.exposure','other_allele.exposure',
                                   'eaf.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure',
                                   'units.exposure','exposure','mr_keep.exposure','pval_origin.exposure',
                                   'units.exposure_dat','id.exposure','data_source.exposure')]
   return(exposure_dat)
}

# Main process ------------------------------------------------------------
# Get inputs
ls.snp = fread(f.snp,header=F) %>% .$V1
ls.cpg = fread(f.cpg,header=F) %>% .$V1
ls.exposure.filename = list.files(exposure.path,pattern='^mQTL',full.names = T) %>%
   .[grep(paste0(ls.cpg,collapse = '|'),.)]
ls.exposure=ls.exposure.filename %>% gsub(exposure.path,'',.) %>% 
   gsub('.rds','',.) %>% gsub('mQTL.','',.) %>% gsub('/','',.)
ls.exposure.input = data.frame(fname=ls.exposure.filename,exposure.name=ls.exposure,
                               stringsAsFactors = F)
# Extract Summstats
exposure.dat.stp1=ls.exposure.input %>%
   pbapply(.,MARGIN = 1,FUN=extract_summstats)
names(exposure.dat.stp1) = ls.exposure.input$exposure.name

ls.exposure.input$exposure.name %>% as.list %>% 
   pblapply(.,FUN=function(x) write_tsv(exposure.dat.stp1[[x]],
                                      file = paste0(output.path,'/',x,'.tops'),
                                      col_names = T))

# Clump again using the first exposure
tmp.exposure_dat <- read_exposure_data(paste0(output.path,'/',ls.exposure.input$exposure.name[1],".tops"),
                                   sep='\t')
tmp.exposure_dat <- clump_data(tmp.exposure_dat,clump_kb=250,clump_r2=0.001)
ls.SNP.tokeep <- tmp.exposure_dat$SNP

# Reformat and re-select SNPs for other exposure data
ls.exposure.input = data.frame(fname=paste0(output.path,'/',ls.exposure,'.tops'),
                               exposure.name=ls.exposure,
                               stringsAsFactors = F)
exposure.dat.stp2=ls.exposure.input %>%
   pbapply(.,MARGIN = 1,FUN=format_exposure,
           tmp.SNP.toinclude=ls.SNP.tokeep)

# A summarised exposure data
summ.exposure.dat = exposure.dat.stp2 %>%
   bind_rows

write.table(summ.exposure.dat,file=paste0(output.path,'/multivar.exposure_dat'),
            col.names=T,row.names = F,sep="\t",quote=F)

rm(list=ls())
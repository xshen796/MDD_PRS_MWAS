# Create a randomised list of CpGs for reml regression and prediction

# Basic settings ----------------------------------------------------------

library(data.table)                                      # Data management   
library(dplyr)                                           # Data management 
library(optparse)                                        # Process control
library(tibble)                                          # Data management  

# Parse arguments
args <- commandArgs(trailingOnly = FALSE)
parse <- OptionParser()
option_list <- list(
  make_option('--ewasSummstats', type='character', help="CpGList", action='store'),
  make_option('--pT', type='numeric', help="p threshold", action='store'),
  make_option('--outDir', type='character', help="Output directory", action='store'),
  make_option('--outFile', type='character', help='Output file', action='store')
)
args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

CpG.loc = opt$ewasSummstats
pT = opt$pT
OUT_DIR = opt$outDir
OUT_FILE = opt$outFile


# Create a randomised list ------------------------------------------------
# Load EWAS summstats
ls.cpg=read.delim(CpG.loc,sep='\t',header=T,stringsAsFactors=F) %>%
  filter(.,P.value<pT) %>% # filter by a given p threshold
  select(CpG=MarkerName)

# Write CpG list ----------------------------------------------------------
write.table(ls.cpg,file=paste0(OUT_DIR,'/',OUT_FILE),
            append = F,quote = F,row.names = F,col.names = F)




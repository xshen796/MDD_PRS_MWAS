# Create a randomised list of CpGs for reml regression and prediction

# Basic settings ----------------------------------------------------------

library(data.table)
library(dplyr)
library(optparse)
library(tibble)

# Parse arguments
args <- commandArgs(trailingOnly = FALSE)
parse <- OptionParser()
option_list <- list(
  make_option('--CpGList', type='character', help="CpGList", action='store'),
  make_option('--outDir', type='character', help="Output directory", action='store'),
  make_option('--outFile', type='character', help='Output file', action='store')
)
args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

CpG.loc = opt$CpGList
OUT_DIR = opt$outDir
OUT_FILE = opt$outFile


# Create a randomised list ------------------------------------------------
# Load data 
ls.cpg.all = read.table(CpG.loc,header=F,stringsAsFactors = F)

# Create a randomised list 
target.cpg = sample(ls.cpg.all$V1,size=725,replace = F)

# Write CpG list ----------------------------------------------------------
write.table(target.cpg,file=paste0(OUT_DIR,'/',OUT_FILE),
            append = F,quote = F,row.names = F,col.names = F,sep = '\n')




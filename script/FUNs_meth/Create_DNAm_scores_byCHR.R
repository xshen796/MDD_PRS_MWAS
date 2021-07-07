## Create DNAm scores in an independent sample
## By X Shen
## 02/12/2020

## Info on file formats
## data is the pre_norm_beta_vals_clean450.Rdata
## pheno_file is a file with a column called ID which matches the methylation ID (colnames) in prenormalized beta file
## coef is a 2 column file with cpg and coefficient (coef.name) and (coef.value)
## output location is the name of the output file

# Basic settings ----------------------------------------------------------

library(data.table)
library(dplyr)
library(optparse)
library(tibble)

# Parse arguments
args <- commandArgs(trailingOnly = FALSE)
parse <- OptionParser()
option_list <- list(
   make_option('--methdat', type='character', help="Location for methylation data", action='store'),
   make_option('--subID', type='character', help="Subjects to include in the analysis", action='store'),
   make_option('--lassoCoef', type='character', help='Lasso coefficients', action='store'),
   make_option('--out', type='character', help='Output file', action='store')
)
args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

meth.dat.loc = opt$methdat
sub.ID = opt$subID
lasso.coef.datloc = opt$lassoCoef
output = opt$out


# Load data ---------------------------------------------------------------

## individuals to create scores in
  pheno_file <- readRDS(sub.ID)
## LASSO co-efficients
  coef.training <- read.table(lasso.coef.datloc, sep=" ", header=F, stringsAsFactors = F)
  colnames(coef.training)=c('coef.name','coef.value')
  
ls.meth.file.loc=list.files(path=meth.dat.loc,pattern='^mvalues.chr',full.names=T)


# Define functions --------------------------------------------------------
# Process managing functions:
# Log process
log_file <- gsub('.rds','.log',output)
logging <- function(str) { 
  cat(paste0(paste0(str, collapse=''), '\n'), file=log_file, append=TRUE) 
  cat(paste0(paste0(str, collapse=''), '\n')) 
}

# Calculate MRS by chromosome (main function)

calc_MRS <- function(F_mvalue,Obj_pheno,Obj_weight){
  
  chr.name=F_mvalue %>% strsplit('/') %>% unlist %>% tail(1) %>% 
            strsplit('\\.') %>% unlist %>% .[2]
  logging(c('Processing:\t', chr.name))
  
  # Load m-values
  data <- readRDS(F_mvalue)
  
  # Select participants with phenotype
  a = colnames(data) %in% Obj_pheno$ID
  meth = data[,a] %>% as.data.frame                           
  rm(data) # large files, remove after use  
  
  # Select CpGs with weights
  a = rownames(meth) %in% Obj_weight$coef.name
  meth.cpg = meth[a,]
  
  # Match weights and m-value data
  probes <- intersect(Obj_weight$coef.name, rownames(meth.cpg))
  rownames(Obj_weight) = Obj_weight$coef.name
  b = meth.cpg[probes,]
  p = Obj_weight[probes,]
  
  # Add weights
  for (i in probes) {
    b[i,]= b[i,]*p[i,"coef.value"]
  }
  
  # Calculate sum scores
  predicted_dep=colSums(b,na.rm=T)
  # if not trained on residuals then add the intercept: predicted_dep=colSums(b) + coef[1,2]
  pred_dep.chr = data.frame(predicted_dep)
  colnames(pred_dep.chr)=chr.name
  rownames(pred_dep.chr)=colnames(b)
  
  return(pred_dep.chr)
}



# Run analysis ------------------------------------------------------------

logging('Create MRS')
logging(c("Started: ", date()))
logging(c('M-value directory: ', meth.dat.loc))
logging(c('Testing sample: ', sub.ID))
logging(c('Lasso regression weights:', lasso.coef.datloc))
logging(c('Output file:', output))
logging(' ')


pred_dep.allCHR = as.list(ls.meth.file.loc) %>%
  lapply(.,FUN=calc_MRS,Obj_pheno=pheno_file,Obj_weight=coef.training) 

ref.ID=rownames(pred_dep.allCHR[[1]])
pred_dep.allCHR = pred_dep.allCHR %>%
  bind_cols %>%
  mutate(total.MRS=rowSums(.)) %>%
  mutate(ID=ref.ID)

saveRDS(pred_dep.allCHR, file=output)

logging(c('Output saved:', output))
logging(c('Ended:  ', date()))

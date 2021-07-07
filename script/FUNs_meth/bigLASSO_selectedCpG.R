# Lasso regression on high-dimensional data
# by XShen
# 02/12/2020


# Basic settings ----------------------------------------------------------

library(optparse)        # Parse arguments
library(biglasso)        # Lasso regression
library(dplyr)           # Data manipulation
library(bigmemory)       # For processing big matrices      
library(bigstatsr)       # For processing big matrices   

args = commandArgs(trailingOnly=TRUE)

parse <- OptionParser()

option_list <- list(
   make_option('--methFolder', type='character', help="Folder for methylation", action='store'),
   make_option('--subID', type='character', help='Data dictionary file', action='store'),
   make_option('--cpgList', type='character', help='Reference file for categorical variables', action='store'),
   make_option('--phenoFile', type='character', help='Output folder', action='store'),
   make_option('--phenoName', type='character', help='Whether to append or rewrite', action='store',default='yes'),
   make_option('--outputFile', type='character', help='Whether to append or rewrite', action='store',default='yes'))


args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

D_METH = opt$methFolder
F_subID = opt$subID
F_cpg = opt$cpgList
F_pheno = opt$phenoFile
phenoName = opt$phenoName
F_output = opt$outputFile
F_snp_ch_probe = '/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/SNP_CH_probes'

# Define process managing functions ---------------------------------------

log_file <- gsub('.txt','.log',F_output)
logging <- function(str) { 
   cat(paste0(paste0(str, collapse=''), '\n'), file=log_file, append=TRUE) 
   cat(paste0(paste0(str, collapse=''), '\n')) 
}



# Log arguments -----------------------------------------------------------

logging('Lasso regression')
logging(c("Started: ", date()))
logging(c('Methylation data folder: ', D_METH))
logging(c('Subject to include: ', F_subID))
logging(c('CpGs to include: ', F_cpg))
logging(c('Addition SNP CpGs to exclude: ', F_snp_ch_probe))
logging(c('Phenotype file: ', F_pheno))
logging(c('Phenotype: ', phenoName))
logging(c('Output file: ', F_output))
logging(' ')

# Load data ---------------------------------------------------------------

# M-values
ls.meth.f = list.files(path = D_METH,pattern = 'mvalues.chr',full.names = T)
meth <- as.list(ls.meth.f) %>%
   lapply(.,FUN=readRDS) %>%
   bind_rows

# CpGs to include in the analysis
if(!is.null(opt[['cpgList']])) {
   cpg.ID=readRDS(F_cpg)
   meth=meth[rownames(meth) %in% cpg.ID$cpg, ]
}

# SNP CpGs to exclude from analysis
probes_excl = read.table(file = F_snp_ch_probe,stringsAsFactors=F)
colnames(probes_excl)='cpg'
meth=meth[!rownames(meth) %in% probes_excl$cpg, ]

logging('Methylation data loaded')
logging(c('NCpG = ',nrow(meth)))

# depending on the structure of data.
# The purpose here is to set one participant per row.
meth = t(meth) 

# read linkage file
ID_file=readRDS(F_subID)
meth = meth[rownames(meth) %in% ID_file$ID,]

# load phenotypes
pheno.dat=readRDS(F_pheno)
rownames(pheno.dat)=pheno.dat$ID


# Prepare methylation and phenotype data ----------------------------------
# Two objects need to have matched IDs

# Find matching IDs
meth = meth[rownames(meth) %in% pheno.dat$ID,]
a = intersect(row.names(meth), row.names(pheno.dat))
logging(c('N=',length(a)))

# Intersect datasets
meth_intersected <- meth[a,]
pheno_intersected <- pheno.dat[a,]
# Check ordering is correct
if (!all.equal(row.names(meth_intersected), row.names(pheno_intersected))){
   stop('Methylation and phenotype data need to have matching IDs',call.=F)
}

tmp.y=pheno_intersected[,phenoName]
y.input = tmp.y

rm(meth) # Remove large methylation file to clean up memory

logging('Phenotype loaded')


# Analysis ----------------------------------------------------------------

# Note: run on ecdflibrary (biglasso). Preferrably h_vmem>150G
# Rcpp needs to be used for very large matrices.

## => Set up file structure so there is a folder 'data' for 
## storing temporary matrix data

if(!dir.exists('data')){
   system('mkdir data')
}


if(file.exists('data/test.bk')){
   system('rm data/test.*')
}

## => Transform dataframe to big matrix
big_meth <- bigstatsr::as_FBM(meth_intersected, backingfile = "data/test")
(desc <- sub("\\.bk$", ".desc", big_meth$backingfile))
dput(big_meth$bm.desc(), desc)

X.bm <- attach.big.matrix(desc)
#X.bm = as.big.matrix(meth_intersected,backingfile='')
logging('Big.matrix transformed')

## => Lasso regression
cvfit <- cv.biglasso(X.bm, y.input, seed = 1234, nfolds = 5,ncores= 3) 
logging('Lasso finished')

## => Summarise results
# Note: coef(cvfit) is not a dataframe object
out <- data.frame(Coefs = coef(cvfit,type='coefficients')[which(coef(cvfit,type='coefficients') != 0)])
out$Vars = c('Intercept',colnames(meth_intersected))[which(coef(cvfit) != 0)]
weights <- out[-1,] %>% ## no need to use the intercept if trained on residuals
   select(Vars,Coefs)
out.bakup <- coef(cvfit,type='coefficients')

## => Write results    
write.table(weights, file=F_output, quote=F, row.names=F, col.names=F)
save(out.bakup,cvfit, file='data/tmp.res.RData')


logging('Lasso results saved')
logging(c("Finished: ", date()))
logging('\n')
gc()      


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

meth.dat.loc = args[1]
sub.ID = args[2]
cpg.ID = args[3]
pheno.datloc = args[4]
pheno.name = args[5]
output.name = args[6]

### load methylation data (for each chromosome)

probes_excl = read.table('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/SNP_CH_probes',stringsAsFactors=F)
colnames(probes_excl)='cpg'

ls.meth.f = list.files(path = meth.dat.loc,pattern = 'mvalues.chr',full.names = T)

meth <- as.list(ls.meth.f) %>%
   lapply(.,FUN=readRDS) %>%
   bind_rows
   
cpg.ID=readRDS(cpg.ID)
meth=meth[rownames(meth) %in% cpg.ID$cpg, ]
meth=meth[!rownames(meth) %in% probes_excl$cpg, ]

meth = t(meth) # depending on the structure of data. The purpose here is to set one participant in a row

print('Methylation data loaded')
print(paste0('NCpG = ',ncol(meth)))


### 100G node on ecdflibrary(biglasso)

# read linkage file
ID_file=readRDS(sub.ID)
meth = meth[rownames(meth) %in% ID_file$ID,]


# load phenotypes
pheno.dat=readRDS(pheno.datloc)
rownames(pheno.dat)=pheno.dat$ID

# read in phenotype
row.names(pheno.dat) <- pheno.dat$ID 

meth = meth[rownames(meth) %in% pheno.dat$ID,]
a = intersect(row.names(meth), row.names(pheno.dat))
print(paste0('N=',length(a)))
meth_intersected <- meth[a,]
pheno_intersected <- pheno.dat[a,]
# Check ordering is correct
if (!all.equal(row.names(meth_intersected), row.names(pheno_intersected))){
      stop('Methylation and phenotype data need to have matching IDs',call.=F)
}

rm(meth)

print('Phenotype loaded')

if(file.exists('data/test.bk')){
   system('rm data/test.*')
}

big_meth <- bigstatsr::as_FBM(meth_intersected, backingfile = "data/test")
(desc <- sub("\\.bk$", ".desc", big_meth$backingfile))
dput(big_meth$bm.desc(), desc)

X.bm <- attach.big.matrix(desc)
           
print('Big.matrix transformed')
      tmp.y=pheno_intersected[,pheno.name]
      y.input = tmp.y
      cvfit <- cv.biglasso(X.bm, y.input, seed = 1234, nfolds = 10, ncores = 2) 
print('Lasso finished')
      out <- as.data.frame(cbind(rownames(coef(cvfit))[which(coef(cvfit) != 0)], coef(cvfit)[which(coef(cvfit) != 0)]))
      
      #out1 <- out[-1,]## no need to use the intercept if trained on residuals
      write.table(out, file=output.name, quote=F, row.names=F, col.names=F)
print('Lasso results saved')
      gc()      


args = commandArgs(trailingOnly=TRUE)
#if (length(args)<4) {
#      stop("At least four arguments must be supplied: methylation data position, file for sub ID, phenotype file and column name for phenotype", 
#           call.=FALSE)
#} 

library('biglasso')

meth.dat.loc = args[1]
sub.ID = args[2]
pheno.datloc = args[3]
pheno.name = args[4]
output.name = args[5]

### load methylation data (for each chromosome)

meth <- read.table(meth.dat.loc,row.names=1,header=T)
colnames(meth)=gsub('^X','',colnames(meth))

meth = t(meth) # depending on the structure of data. The purpose here is to set one participant in a row

print('Methylation data loaded')


### 100G node on ecdflibrary(biglasso)

# read linkage file
ID_file=readRDS(sub.ID)
meth = meth[rownames(meth) %in% ID_file$ID,]

# load phenotypes
pheno.dat=readRDS(pheno.datloc)
rownames(pheno.dat)=pheno.dat$ID

# read in alcohol phenotype
row.names(pheno.dat) <- pheno.dat$ID 
a = intersect(row.names(meth), row.names(pheno.dat))
print(paste0('N=',length(a)))
X <- meth[a,]
y <- pheno.dat[a,]
# Check ordering is correct
if (!all.equal(row.names(X), row.names(y))){
      stop('Methylation and phenotype data need to have matching IDs',call.=F)
}

rm(meth)

print('Phenotype loaded')


      tmp.y=y[,pheno.name]
      X.bm <- as.big.matrix(X[!is.na(tmp.y),])
print('Big.matrix transformed')
      y.input = tmp.y[!is.na(tmp.y)]
      cvfit <- cv.biglasso(X.bm, y.input, seed = 1234, nfolds = 10, ncores = 2) 
print('Lasso finished')
      out <- as.data.frame(cbind(rownames(coef(cvfit))[which(coef(cvfit) != 0)], coef(cvfit)[which(coef(cvfit) != 0)]))
      
      out1 <- out[-1,]## no need to use the intercept if trained on residuals
      write.table(out1, file=output.name, quote=F, row.names=F, col.names=F)
print('Lasso results saved')
      gc()      


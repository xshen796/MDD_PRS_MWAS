## Create DNAm scores in independent sample
## Uses pre-normalized beta values

## Info on file formats
## data is the pre_norm_beta_vals_clean450.Rdata
## pheno_file is a file with a column called ID which matches the methylation ID (colnames) in prenormalized beta file
## coef is a 2 column file with cpg and coefficient (coef.name) and (coef.value)
## output location is the name of the output file

args=commandArgs(TRUE)

  if (length(args)==0) {
        stop("At least four arguments must be supplied including (betavals), (phenofile), (lasso-coefficients).n and output location", call.=FALSE)
  }

meth.dat.loc = args[1]
sub.ID = args[2]
lasso.coef.datloc = args[3]
output = args[4]

## Load data
## methylation data
  data <- read.table(meth.dat.loc,row.names=1,header=T)
  colnames(data)=gsub('^X','',colnames(data))
  data=t(data)
## individuals to create scores in
  pheno_file <- readRDS(sub.ID)
## LASSO co-efficients
  coef <- read.table(lasso.coef.datloc, sep=" ", header=F)
  colnames(coef)=c('coef.name','coef.value')


a = which(rownames(data) %in% pheno_file$ID)
meth = data[a,]                           
rm(data) # large files, remove after use                                  

meth1 = as.data.frame(meth)
meth1$id = as.character(rownames(meth1))


a = which(names(meth1) %in% coef$coef.name)
meth2 = meth1[,a]
meth3 = t(meth2)
probes <- intersect(coef$coef.name, rownames(meth3))
rownames(coef) = coef$coef.name

b = meth3[probes,]
p = coef[probes,]

for (i in probes) {
      b[i,]= b[i,]*p[i,"coef.value"]
}

predicted_dep=colSums(b)
# if not trained on residuals then add the intercept: predicted_dep=colSums(b) + coef[1,2]
pred_dep = as.data.frame(predicted_dep)
pred_dep$ID = rownames(pred_dep)

dep = merge(pheno_file, pred_dep, by="ID")

saveRDS(dep, file=output)
cat(paste0('Output saved as: ',output))

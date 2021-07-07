## Create DNAm scores in an independent sample

library(data.table)
library(dplyr)
library(optparse)

## Info on file formats
## data is the pre_norm_beta_vals_clean450.Rdata
## pheno_file is a file with a column called ID which matches the methylation ID (colnames) in prenormalized beta file
## coef is a 2 column file with cpg and coefficient (coef.name) and (coef.value)
## output location is the name of the output file

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

## Load data
## individuals to create scores in
  pheno_file <- readRDS(sub.ID)
## LASSO co-efficients
  coef <- read.table(lasso.coef.datloc, sep=" ", header=F)
  colnames(coef)=c('coef.name','coef.value')
  
ls.meth.file.loc=list.files(path=meth.dat.loc,pattern='^mvalues',full.names=T)

for (mvalue.file in ls.meth.file.loc){
      
      chr.name=gsub(meth.dat.loc,'',mvalue.file) %>%
                  strsplit(split='\\.') %>% unlist %>% .[2]
      ## methylation data
      data <- readRDS(mvalue.file)
      data=t(data)
  
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
      
      predicted_dep=colSums(b,na.rm=T)
      # if not trained on residuals then add the intercept: predicted_dep=colSums(b) + coef[1,2]
      pred_dep.tmp = as.data.frame(predicted_dep)
      colnames(pred_dep.tmp)=chr.name
      
      if (mvalue.file==ls.meth.file.loc[1]){pred_dep.allCHR=pred_dep.tmp}else{pred_dep.allCHR=data.frame(pred_dep.allCHR,pred_dep.tmp)}

}

pred_dep.allCHR$ID = rownames(pred_dep.allCHR)

dep = merge(pheno_file, pred_dep.allCHR, by="ID")

saveRDS(dep, file=output)
cat(paste0('Output saved as: ',output))

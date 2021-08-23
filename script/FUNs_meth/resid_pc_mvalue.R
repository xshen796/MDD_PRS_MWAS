# XShen
# 23/08/2021

library(plyr)            # Data management
library(dplyr)           # Data management
library(tibble)          # Data management
library(FactoMineR)      # PCA
library(optparse)

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
args = commandArgs(trailingOnly=TRUE)

parse <- OptionParser()

option_list <- list(
  make_option('--meth', type='character', help="Methylation data", action='store'),
  make_option('--covar', type='character', help='Covariates (rds format)', action='store'),
  make_option('--residualise', type='character', help='y/n', action='store'),
  make_option('--makePC', type='character', help='y/n', action='store'),
  make_option('--out', type='character', help='Output directory', action='store'))


opt <- parse_args(OptionParser(option_list=option_list), args=args)

D_METH = opt$meth
F_covar = opt$covar
D_out = opt$out


if (grepl('/',D_METH)){
  METH_dir = D_METH %>% strsplit(.,'/') %>% unlist %>% .[length(.)-1]
}else{
  METH_dir = '.'
}

METH_pattern = D_METH %>% strsplit(.,'/') %>% unlist %>% tail(1)

resid.y = opt$residualise
makepc.y = opt$makePC

# Define functions  -------------------------------------------------------
residualise_dat <- function(col.input,dat) {
  coln=as.character(col.input[1])
  cov=as.character(col.input[2])
  colnames(dat)[1]='ID'
  
  eq=paste0(coln,'~',cov)
  fit=glm(as.formula(eq),dat=dat)
  new.dat=resid(fit)
  ls.cov=strsplit(cov,'\\+') %>% unlist
  new.dat=data.frame(ID=dat[complete.cases(dat[,c(coln,ls.cov)]),1],pheno=new.dat,stringsAsFactors=F)
  
  original.ID=data.frame(ID=dat$ID,stringsAsFactors=F)
  new.dat=merge(original.ID,new.dat,by='ID',all.x=T)
  
  output=data.frame(pheno=new.dat$pheno,stringsAsFactors=F)
  colnames(output)=coln
  return(output)
}

residualise_chr <- function(target.dat,cov.dat,linkage.dat){

    meth.dat = t(target.dat) %>%
    data.frame()%>%
    dplyr::mutate(ID=rownames(.)) %>%
    merge(.,cov.dat,by='ID')
  
  input.for_resid = data.frame(pheno=colnames(meth.dat)[grep('^cg',colnames(meth.dat))],
                               covs=paste0(colnames(cov.dat)[!grepl('ID',colnames(cov.dat))],collapse='+'),
                               stringsAsFactors=F)
  
  cat('Start processing\n')
  tmp.dat.aftercorr <- apply(input.for_resid,1,residualise_dat,dat=meth.dat) %>% # A list created by residualise_dat
    bind_cols %>%  # Convert into a dataframe
    mutate(.,ID=meth.dat$ID) %>% # Add ID
    select('ID', everything()) # Move ID to the first column
  
}


# Load data ---------------------------------------------------------------

covariates=readRDS(F_covar) 
colnames(covariates)[1]='ID'
covariates$ID = as.character(covariates$ID)


# Residualise -------------------------------------------------------------

meth.dat = list.files(path = METH_dir, pattern = METH_pattern,full.names = T,include.dirs = T) %>% 
  lapply(readRDS) 

if (resid.y=='y'){
  resid.meth.dat = meth.dat %>% 
    lapply(.,FUN=residualise_chr,cov.dat=covariates)
}else if (resid.y=='n'){
  resid.meth.dat = meth.dat %>% 
    lapply(.,as.data.frame) %>% 
    lapply(.,FUN=tibble::rownames_to_column,var='ID')
}


# Create methylation PC ---------------------------------------------------

if (makepc.y=='y'){
  mval.forPCA = resid.meth.dat %>%
    join_all(by='ID', type='left') %>%
    tibble::column_to_rownames(.,var = 'ID')
  
  # Perform PCA
  fit.pca = PCA(mval.forPCA, scale.unit = TRUE, ncp = 20, graph = F)
  
  methPCA = fit.pca$ind$coord
  
  # Reformat to fit ewas pipeline
  colnames(methPCA)=paste0('PC',1:20)
  methPCA = methPCA %>% data.frame %>%
    tibble::rownames_to_column(var='ID') %>% data.frame
}


# Save data ---------------------------------------------------------------

if (resid.y=='y'){
  names(resid.meth.dat) = list.files(path = METH_dir, pattern = METH_pattern,full.names = T,include.dirs = T) %>%
    strsplit(.,'/') %>% 
    lapply(.,tail,n=1) %>% 
    unlist %>% 
    gsub('.rds','_resid.rds',fixed = T) 
  
  names(output.dat) %>% as.list %>% 
    lapply(.,FUN=function(x) saveRDS(output.dat[[x]],paste0(D_out,'/',x)))
}

if (makepc.y=='y'){
  saveRDS(methPCA,paste0(D_out,'/','methPCs.rds'))
}


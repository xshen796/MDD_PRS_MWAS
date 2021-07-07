library(dplyr)
library(data.table)

setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD')

# Load data -------------------------------------------------------------------------------------------------------

# Load PRS
prs.ls=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/score_thresholds.txt',
                  stringsAsFactors=F,header=F)
colnames(prs.ls)=c('fname','label.pT','pT')
prs.ls$label.pT=as.character(prs.ls$label.pT)
prs.ls$label.pT=format(prs.ls$pT,scientific=F,drop0trailing=T)

prs.all = as.list(1:10) %>%
      lapply(.,function(x) read.table(paste0('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/Meta_predict_FullGS.S',
                                             x,'.profile'),stringsAsFactors=F,header=T)) %>%
      lapply(.,function(x) select(x,ID=IID,starts_with('SCORE_'))) %>%
      Reduce(function(x,y) left_join(x,y,by='ID'),.)

# Load MDD and other covariates
load("/exports/igmm/eddie/GenScotDepression/data/genscot/phenotypes/masterDB.Rdata")
load('data/GS_phenotype/STRADL.Rdata')

totaldata = merge(totaldata,smoking,by='id',all.x=T)
totaldata = merge(totaldata,x[,c('id','CIDI_MDD')],by='id',all.x=T)

dat.training$MDD_status=totaldata$dep_status
dat.training$MDD_status[(totaldata$CIDI_MDD!=0)&(dat.training$MDD_status==0)]=NA
dat.training$MDD_status[is.na(totaldata$CIDI_MDD)&(dat.training$MDD_status==0)]=NA
dat.training$MDD_status[dat.training$MDD_status==3]=NA

PC=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/HM3mds.mds',sep='',header=T)
PC=PC[,c(2,4:13)]
PC=PC[!duplicated(PC$IID),]


# Regression ------------------------------------------------------------------------------------------------------


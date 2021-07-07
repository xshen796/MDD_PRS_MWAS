# By Shen
# Use igmm/apps/R/3.5.1

setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')
library(dplyr)

# Load data  ----------------------------------------------------------------------------------------------

link_file = read.csv('data/gwas_methid_link.csv',stringsAsFactors=F)
# sample mvalue data: groot=readRDS('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/chr/mvalues.chr22.rds')
load("/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/masterDB.Rdata")

demographic=totaldata[,c('id','sex','age','famid')]
demographic=merge(demographic,link_file,by.x='id',by.y='gwas_id',all.y=T)

# Age and sex ---------------------------------------------------------------------------------------------
# overall sample
targetdata=demographic
sex_age=data.frame(sample='All',age.mean=mean(targetdata$age,na.rm=T),age.sd=sd(targetdata$age,na.rm=T),male=sum(targetdata$sex=='M',na.rm=T)/nrow(targetdata),N=nrow(targetdata))
# wave1
targetdata=filter(demographic,wave=='Wave1')
sex_age=rbind(sex_age,
          data.frame(sample='wave1',age.mean=mean(targetdata$age,na.rm=T),age.sd=sd(targetdata$age,na.rm=T),male=sum(targetdata$sex=='M',na.rm=T)/nrow(targetdata),N=nrow(targetdata)))
# wave2
targetdata=filter(demographic,wave=='Wave3')
sex_age=rbind(sex_age,
          data.frame(sample='wave2',age.mean=mean(targetdata$age,na.rm=T),age.sd=sd(targetdata$age,na.rm=T),male=sum(targetdata$sex=='M',na.rm=T)/nrow(targetdata),N=nrow(targetdata)))

sex_age

# N changes after data merging -----------------------------------------------------------------------------
# Merging w PRS data
prs=read.table('data/MDD_prs/related/Meta_predict_FullGS.S1.profile',stringsAsFactors=F,header=T)
prs=prs[,c('IID','SCORE_0.00000005')]

targetdata=demographic
targetdata=merge(targetdata,prs,by.x='id',by.y='IID',all.x=T)
targetdata=filter(targetdata,!is.na(SCORE_0.00000005))
N=data.frame(merge.with='PRS',Nwave1=sum(targetdata$wave=='Wave1'),Nwave2=sum(targetdata$wave=='Wave3'),Nmeta=nrow(targetdata))
N

# Merging w technical covariates
covs=read.table('data/EWAS_MDDprs_fromKathryn/wave1/mdd-pgrs-wave1-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',')
covs=covs[,c('Sample_Name','sex','age','ever_smoke','pack_years','Sentrix_Position','Gran','Mono','Bcell','NK','CD4T','CD8T','Batch')]
covs$Batch=as.factor(covs$Batch)
covs.wave1=covs
covs=read.table('data/EWAS_MDDprs_fromKathryn/wave3/mdd-pgrs-wave3-snp_ch/mdd_SNPCH_prs_5e_08/dat.csv',header=T,sep=',')
covs=covs[,c('Sample_Name','sex','age','ever_smoke','pack_years','Sentrix_Position','Gran','Mono','Bcell','NK','CD4T','CD8T','Batch')]
covs$Batch=as.factor(covs$Batch)
covs.wave3=covs
covs=rbind(covs.wave1,covs.wave3)

targetdata=merge(targetdata,covs,by.x='id',by.y='Sample_Name',all.x=T)
targetdata=targetdata[complete.cases(targetdata[,c('Sentrix_Position','Gran','Mono','Bcell','NK','CD4T','CD8T','Batch')]),]

N=rbind(N,data.frame(merge.with='tech.covs',Nwave1=sum(targetdata$wave=='Wave1'),Nwave2=sum(targetdata$wave=='Wave3'),Nmeta=nrow(targetdata)))
N

# Merging w lifestyle factors, age and sex
targetdata=targetdata[complete.cases(targetdata[,c('sex.x','age.x','ever_smoke','pack_years')]),]
N=rbind(N,data.frame(merge.with='behav.covs',Nwave1=sum(targetdata$wave=='Wave1'),Nwave2=sum(targetdata$wave=='Wave3'),Nmeta=nrow(targetdata)))
N

# N for lifestyle factors -----------------------------------------------------------------------------
load("/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/masterDB.Rdata")
load('data/STRADL.Rdata')
totaldata = merge(totaldata,x[,c('id','CIDI_MDD')],by='id',all.x=T)

totaldata = totaldata[,c('id','dep_status','CIDI_MDD','alcohol_units_week','bmi')]
targetdata = merge(targetdata,totaldata,by='id',all.x=T)

# MDD
targetdata$MDD_status=targetdata$dep_status
targetdata$MDD_status[(targetdata$CIDI_MDD!=0)&(targetdata$MDD_status==0)]=NA
targetdata$MDD_status[is.na(targetdata$CIDI_MDD)&(targetdata$MDD_status==0)]=NA
targetdata$MDD_status[targetdata$MDD_status==3]=NA

table(targetdata$MDD_status[targetdata$wave=='Wave1'])
table(targetdata$MDD_status[targetdata$wave=='Wave3'])
table(targetdata$MDD_status)

# Alc
sum(!is.na(targetdata$))
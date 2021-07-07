setwd("/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD")
library(dplyr)
load("/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/masterDB.Rdata")
smoking = read.csv('data/GS_phenotype/smoking_pack_years.csv')
load('data/GS_phenotype/STRADL.Rdata')
totaldata = merge(totaldata,smoking,by='id',all.x=T)
totaldata = merge(totaldata,x[,c('id','CIDI_MDD')],by='id',all.x=T)

# Covariates --------------------------------------------------------------------------

# BMI
dat.training=totaldata[,c('id','bmi')]
dat.training$bmi[dat.training$bmi<15] = NA
dat.training$bmi[dat.training$bmi>60] = NA
#dat.training$log_bmi = log(dat.training$bmi)

# alcohol consumption
dat.training$alcohol_units_week = totaldata$units
dat.training$alcohol_units_week[dat.training$alcohol_units_week>300] = NA 
dat.training$alcohol_units_week[totaldata$drink_status==2|totaldata$drink_status==3]=NA
#dat.training$log_alcohol_units_week = log(dat.training$alcohol_units_week+1)

# pack years
dat.training$smoking_pack_year = totaldata$pack_years
#dat.training$smoking_pack_year[totaldata$ever_smoke.y==2]=NA
#dat.training$smoking_pack_year[totaldata$ever_smoke.y==3]=NA
dat.training$smoking_pack_year[(totaldata$ever_smoke.y==4)]=0
#dat.training$smoking_pack_year[totaldata$age_started<=10]=NA
#dat.training$log_pack_year = log(dat.training$smoking_pack_year+1)

# ever smoke
dat.training$ever_smoke=totaldata$ever_smoke.y
dat.training$ever_smoke[dat.training$ever_smoke==4]=0
dat.training$ever_smoke[dat.training$ever_smoke==1]=999
dat.training$ever_smoke[dat.training$ever_smoke==3]=1
dat.training$ever_smoke[dat.training$ever_smoke==999]=3
dat.training$ever_smoke=as.factor(dat.training$ever_smoke)

dat.training=dat.training[,c('id','bmi','alcohol_units_week','smoking_pack_year','ever_smoke')]
dat.training$sex=totaldata$sex
dat.training$age=totaldata$age.x

# add genetic PCs
PC=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/HM3mds.mds',sep='',header=T)
PC=PC[,c(2,4:13)]
PC=PC[!duplicated(PC$IID),]

dat.training=merge(dat.training,PC,by.x='id',by.y='IID',all.x=T)


# MDD -------------------------------------------------------------------------------
# Case = SCID.MDD==1
# Control = SCID.MDD==0 & CIDI.MDD==0
dat.training$MDD_status=totaldata$dep_status
dat.training$MDD_status[(totaldata$CIDI_MDD!=0)&(dat.training$MDD_status==0)]=NA
dat.training$MDD_status[is.na(totaldata$CIDI_MDD)&(dat.training$MDD_status==0)]=NA
dat.training$MDD_status[dat.training$MDD_status==3]=NA

# Add batch ---------------------------------------------------------------
batch.cell=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/w1_w3_celltypes.txt',
                      header=T,stringsAsFactors = F) %>%
   select(-wave)

dat.training=merge(dat.training,batch.cell,by.x='id',by.y='Sample_Name')


# Change into methylation IDs -------------------------------------------------------
link_file = read.csv('data/GS_genetic_meth/gwas_methid_link.csv')

dat.DNAm.pheno = merge(dat.training,link_file,by.x='id',by.y='gwas_id')
dat.DNAm.pheno = dat.DNAm.pheno[,c('meth_id','Batch','wave','age','sex',paste0('C',1:10),
                                  'smoking_pack_year','ever_smoke','MDD_status','alcohol_units_week')]
dat.DNAm.pheno$meth_id = as.character(dat.DNAm.pheno$meth_id) 
colnames(dat.DNAm.pheno)[1]='ID'


# Add methylation PCs  -------------------------------------------------------------

# add methylation PCs
PC.meth.wave1 = readRDS('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/Factominer_PCs_5087_resid_PC100.rds') 
PC.meth.wave3 = readRDS('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/wave3_pcs.rds')

PC.meth.wave1[,paste0('PC',1:100)]=scale(PC.meth.wave1[,paste0('PC',1:100)])
PC.meth.wave3[,paste0('PC',1:100)]=scale(PC.meth.wave3[,paste0('PC',1:100)])

PC.meth = rbind(PC.meth.wave1,PC.meth.wave3) %>% select(Sample_Sentrix_ID,paste0('PC',1:20))
colnames(PC.meth)[grep('PC',colnames(PC.meth))]=paste0('meth.',colnames(PC.meth)[grep('PC',colnames(PC.meth))])

dat.DNAm.pheno=merge(dat.DNAm.pheno,PC.meth,by.x='ID',by.y='Sample_Sentrix_ID',all.x=T)


# Residualise -----------------------------------------------------------------------

residulise_dat <- function(col.input1,col.input2,dat) {
      coln=as.character(col.input1)
      cov=as.character(col.input2) 
      
      eq=paste0(coln,'~',paste0(cov,collapse='+'))
      fit=glm(as.formula(eq),dat=dat)
      new.dat=resid(fit)
      new.dat=data.frame(ID=dat[complete.cases(dat[,c(coln,cov)]),1],pheno=new.dat)
      colnames(new.dat)[2]=paste0('resid.',coln)
      return(new.dat)
}


covs = c('age','wave','sex',paste0('C',1:10),'Batch')
# If use Carmen's preprocessed data then no need to control for age and sex 
# 'bmi','alcohol_units_week','smoking_pack_year','ever_smoke',paste0('meth.PC',1:20),paste0('C',1:10),paste0('meth.PC',1:20),

# residualise
#dat.DNAm.pheno.wave1=merge(dat.DNAm.pheno,residulise_dat(c('MDD_status',covs),dat=filter(dat.DNAm.pheno,wave=='Wave1')),by='ID',all=T)
#dat.DNAm.pheno.wave3=merge(dat.DNAm.pheno,residulise_dat(c('MDD_status',covs),dat=filter(dat.DNAm.pheno,wave=='Wave3')),by='ID',all=T)
#dat.DNAm.pheno.all=rbind(dat.DNAm.pheno.wave1,dat.DNAm.pheno.wave3) %>%
#	filter(.,!is.na(resid.MDD_status))
dat.DNAm.pheno.all=merge(dat.DNAm.pheno,residulise_dat('MDD_status',covs,dat=dat.DNAm.pheno),by='ID',all=T) %>%
   filter(.,!is.na(resid.MDD_status))

saveRDS(dat.DNAm.pheno.all,file='data/methTraining/DNAm_training_pheno_MDD.rds')

#dat.DNAm.pheno.nonSmoker=merge(dat.DNAm.pheno,residulise_dat(c('MDD_status',covs),dat=dat.DNAm.pheno),by='ID',all=T)
#dat.DNAm.pheno.nonSmoker=filter(dat.DNAm.pheno.nonSmoker,!is.na(resid.MDD_status))
#saveRDS(dat.DNAm.pheno.nonSmoker,file='data/methTraining/DNAm_training_pheno_MDD_nonSmoker.rds')


# create subj IDs ---------------------------------------------------------
wave3.ID.all=data.frame(ID=dat.DNAm.pheno$ID[dat.DNAm.pheno$wave=='Wave3'],stringsAsFactors = F)
wave3.nonsmoker=dat.DNAm.pheno$ID[(dat.DNAm.pheno$wave=='Wave3')&(dat.DNAm.pheno$ever_smoke==0)]
wave3.ID.nonsmoker=data.frame(ID=wave3.nonsmoker[!is.na(wave3.nonsmoker)],stringsAsFactors = F)

wave1.ID.all=data.frame(ID=dat.DNAm.pheno$ID[dat.DNAm.pheno$wave=='Wave1'],stringsAsFactors = F)

saveRDS(wave3.ID.all,file='data/wave3_all_ID.rds')
saveRDS(wave3.ID.nonsmoker,file='data/wave3_nonsmoker_ID.rds')
saveRDS(wave1.ID.all,file='data/wave1_all_ID.rds')



# Create osca format phenotype and covariate files ------------------------

osca.pheno = dat.DNAm.pheno %>%
   mutate(FID=ID) %>%
   select(ID,FID,MDD_status)
osca.cov.discrete = dat.DNAm.pheno %>%
   mutate(FID=ID) %>%
   select(ID,FID,ever_smoke,sex)
osca.cov.quant = dat.DNAm.pheno %>%
   mutate(FID=ID) %>%
   select(ID,FID,age,smoking_pack_year,starts_with('meth.PC'))

write.table(osca.pheno,file='data/methTraining/DNAm_training_pheno_MDD_osca',
            sep = '\t',row.names = F,col.names = F,quote = F)
write.table(osca.cov.discrete,file='data/methTraining/DNAm_training_dcov_osca',
            sep = '\t',row.names = F,col.names = F,quote = F)
write.table(osca.cov.quant,file='data/methTraining/DNAm_training_qcov_osca',
            sep = '\t',row.names = F,col.names = F,quote = F)


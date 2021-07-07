setwd("/sdata/images/projects/UKBIOBANK/users/Shen/ii.PGRS/IDP_4thrls_20k/pheWAS/result/x.supplementary_materials/LD.results")

##### load sig hits
ls.files=list.files(path = '/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/methQTL_GS/LD.results',pattern = 'exposure_dat$')
ls.input=data.frame(obj.name=gsub('\\.exposure_dat','',ls.files),file.name=ls.files)
N_sig_hits <- function(ls.input,path.name){
      ls.input[1]->obj.name
      ls.input[2]->file.name
      f.name=paste0(path.name,file.name)
      im_gwas=read.table(f.name,header=F,sep='\t',)
      n.gwas.hits=nrow(im_gwas)
      
      #output=data.frame(phenotype=obj.name,n_hits=n.gwas.hits)
      output=c(obj.name,n.gwas.hits)
      return(output)
}

summary.GWAS.Nhits=t(apply(X = ls.input,MARGIN = 1, N_sig_hits, path.name='IM_GWAS_results/exposure_stats/'))
summary.GWAS.Nhits=data.frame(phenotype=summary.GWAS.Nhits[,1],Nhits=as.numeric(summary.GWAS.Nhits[,2]))

ls.files=list.files(path = 'IM_GWAS_results/outcome_stats/',pattern = 'MDDoutcome_match$')
ls.input=data.frame(obj.name=gsub('\\.MDDoutcome_match','',ls.files),file.name=ls.files)
summary.hits.forMR=t(apply(X = ls.input[c(1:6,8:22,24:27),],MARGIN = 1, N_sig_hits, path.name='IM_GWAS_results/outcome_stats/'))
summary.hits.forMR=data.frame(phenotype=summary.hits.forMR[,1],N.SNPs.MR=as.numeric(summary.hits.forMR[,2]))
summary.GWAS.Nhits=merge(summary.GWAS.Nhits,summary.hits.forMR,by='phenotype',all.x=T)
summary.GWAS.Nhits$N.SNPs.MR[is.na(summary.GWAS.Nhits$N.SNPs.MR)]=0

##### load LD heritability analyses
ls.files=list.files(path='/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/methQTL_GS/LD.results/',pattern = 'MDD_rgresult',full.names=T)

for (fname in ls.files){
  obj.name=gsub('.tsv_forLDSCmunged.sumstats.gz_MDD_rgresult.log','',fname)
  tmp.obj=read.table(fname,header=F,sep='\t',stringsAsFactors=F)
  
  # LDregression result
  anchor.row=grep('Summary of Genetic Correlation Results',tmp.obj$V1)
  tmp.rgresult=tmp.obj$V1[(anchor.row+1):(anchor.row+2)]
  write.table(tmp.rgresult,file=paste0(obj.name,'.input.txt'),col.names = F, row.names = F,
              quote=F)
  
  # heritability summary of IM traits
#  n.anchor=grep('Heritability of phenotype 1',tmp.obj$V1)
#  n.snps=n.anchor-1
#  ns.h2.traitA=(n.anchor+2):(n.anchor+5)
#  if (fname == ls.files[1]){write.table(tmp.obj$V1[n.snps],file='SNP.N.used.txt',
#                                        col.names = F, row.names = F, quote = F)}else{
#    write.table(tmp.obj$V1[n.snps],file='SNP.N.used.txt',append = T,
#                col.names = F,row.names = F,quote=F)
#  }
  
#  tmp.h2=tmp.obj$V1[ns.h2.traitA]
#  tmp.h2=gsub(': ','\t',tmp.h2)
#  write.table(tmp.h2,file=paste0(obj.name,'_h2.txt'),col.names=F,row.names=F,quote=F)
}

rm(list=ls(pattern = '^tmp.|^n.|^ls.|.name|expre'))

ls.files=list.files(path='/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/methQTL_GS/LD.results/',pattern = 'input.txt$',full.names=T)

for (fname in ls.files[c(1,3:4)]){
  obj.name=gsub('.input.txt','',fname)
  tmp.obj=read.table(fname,header=T)
  
  if (fname == ls.files[1]){
    LDrg.results=tmp.obj
  }else{
    LDrg.results=rbind(LDrg.results,tmp.obj)
  }
}

LDrg.results$p1=gsub('./LDready.stats/','',LDrg.results$p1)
LDrg.results$p1=gsub('.bgenie.QC_forLDSCmunged.sumstats.gz','',LDrg.results$p1)
LDrg.results$p2='MDD_3cohorts'
LDrg.results=LDrg.results[,c(1:6,11:12)]


ls.files=list.files(pattern = '_h2.txt$')

for (fname in ls.files){
  obj.name=gsub('_h2.txt','',fname)
  expre=paste0('tmp.obj=read.table(\'',fname,'\',header=F,sep=\'\t\',stringsAsFactors = F)')
  eval(parse(text=expre))
  
  tmp.obj=c(tmp.obj$V2,obj.name)
  if(fname==ls.files[1]){h2.table=tmp.obj}else{h2.table=rbind(h2.table,tmp.obj)}
  
}
rownames(h2.table)=NULL
h2.table=data.frame(h2.table)
colnames(h2.table)=c('Total Observed scale h2','Lambda GC','Mean Chi^2','Intercept','p1')


LDrg.results=merge(LDrg.results,h2.table,by='p1')
rm(list=ls(pattern = '^tmp.|^n.|^ls.|.name|expre'))



###### merge the two together
h2.table=merge(h2.table,summary.GWAS.Nhits,by.x='p1',by.y='phenotype')
colnames(h2.table)[grep('p1',colnames(h2.table))]='Neuroimaging phenotype'
colnames(h2.table)[grep('Nhits',colnames(h2.table))]='Number of genome-wide significant hits'
colnames(h2.table)[grep('N.SNPs.MR',colnames(h2.table))]='Number of SNPs available for MR'


h2.table$`Neuroimaging phenotype`=gsub('^gFA_','gFA-',h2.table$`Neuroimaging phenotype`)
h2.table$`Neuroimaging phenotype`=gsub('^FA_','FA in ',h2.table$`Neuroimaging phenotype`)
h2.table$`Neuroimaging phenotype`=gsub('^gMD_','gMD-',h2.table$`Neuroimaging phenotype`)
h2.table$`Neuroimaging phenotype`=gsub('^MD_','MD in ',h2.table$`Neuroimaging phenotype`)
h2.table$`Neuroimaging phenotype`=gsub('^amp','Amplitude in ',h2.table$`Neuroimaging phenotype`)

setwd("/sdata/images/projects/UKBIOBANK/users/Shen/ii.PGRS/IDP_4thrls_20k/pheWAS/")
#write.table(h2.table,file='table/IM_GWAS_summary.txt',col.names = T, row.names = F, quote = F,sep='\t')


# add LDreg (genetic correlation) --------------------------------------------------

# LDreg
labels_category=c('^gFA','^gICVF','^gMD','^rsamp','^tFA','^tICVF','^tMD')
targetdata=LDrg.results
for (i in labels_category){
      tmp.p.cor=p.adjust(targetdata$p[grep(i,targetdata$p1)],method='fdr')
      if (i==labels_category[1]){p.cor=tmp.p.cor}else{p.cor=c(p.cor,tmp.p.cor)}
}
LDrg.results$p.corrected=p.cor
LDreg.table=data.frame(h2.table,
                       rg_se=paste0(LDrg.results$rg,' (',LDrg.results$se,')'),
                       rg_p_corr=LDrg.results$p.corrected)


save(LDreg.table,file='result/x.supplementary_materials/LDreg_summary.RData')

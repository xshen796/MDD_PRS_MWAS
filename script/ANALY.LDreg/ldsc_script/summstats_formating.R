args=commandArgs(TRUE)

  if (length(args)==0) {
        stop("At least four arguments must be supplied including (betavals), (phenofile), (lasso-coefficients).n and output location", call.=FALSE)
  }

library(dplyr)

dat.loc = args[1]
output.loc = args[2]

tmp.midfile.loc = gsub('.tsv','.txt',dat.loc)

dat=read.delim(dat.loc,fill=T,skip=1,header=F)


col.ls=c('SNP','snpid','chr','pos','A1','A2','N','af','BETA','SE','p','strand','info_type','info','r2')
colnames(dat)=col.ls
dat=filter(dat,snpid!='snpid',SNP!='.')
dat=dat[,!grepl('snpid',colnames(dat))]

write.table(dat,file=tmp.midfile.loc,sep='\t',quote=F,col.names=T,row.names=F)
rm(dat)
dat=read.table(tmp.midfile.loc,header=T,stringsAsFactors=F)

dat=filter(dat,info>0.1,af>0.005,af<0.995)

write.table(dat,file=output.loc,sep='\t',quote=F,col.names=T,row.names=F)

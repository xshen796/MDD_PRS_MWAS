# Run under igmm/apps/R/3.6.1
# by XShen

# Basic settings
setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')

library(ggrepel)
library(dplyr)
library(qqman)
options(bitmapType='cairo')
source('FUNs/manhattan_plot_XShen.R')

ref=read.delim('data/EWAS_MDDprs_Shen/MDDprs_ewas_wave1/MDDprs_pT1_wave1/mddprs_pT1_wave1_RosieData.toptable.txt',sep='\t',header=T,stringsAsFactors=F)
ref.tomerge=ref[,c('ID','CHR','MAPINFO')]


ls.pt=c('5e_08','0.000001','0.0001','0.01','0.1','1')


for (pt in ls.pt){
    f.path=paste0('data/EWAS_MDDprs_Shen/MDDprs_ewas_meta/MDDprs_pT',pt,'_meta/mddprs_pT',pt,'_meta_RosieData.toptable.txt1.tbl')
    ewasResults=read.delim(f.path,sep='\t',header=T,stringsAsFactors=F)
    colnames(ewasResults)[1]='CpG'
    ewasResults=merge(ewasResults,ref.tomerge,by.x='CpG',by.y='ID',all.x=T)

    qman.dat=ewasResults
    qman.dat$CHR=gsub('chr','',qman.dat$CHR)
    qman.dat$CHR=as.numeric(qman.dat$CHR)
    qman.dat=filter(qman.dat,!is.na(CpG),!is.na(CHR),!is.na(MAPINFO),!is.na(P.value))

    tmp.fig=m_plot(dat=qman.dat,chr="CHR", bp="MAPINFO", snp="CpG", p="P.value",plot_title=paste0('EWAS: MDDprs pT<',pt),
                 add_category_name=F,fig_size=c(22,13),tophits_annot=T,y_lim=c(0,120),
                 outputpath=paste0('Figs/ewas_MDDprs_pt',pt,'.png'))
                 
                 
    expre=paste0('fig.pT.',pt,'=tmp.fig')
    eval(parse(text=expre))

}




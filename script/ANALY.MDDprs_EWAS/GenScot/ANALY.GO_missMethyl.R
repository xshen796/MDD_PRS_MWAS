setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')
library('dplyr')
library('tibble')
library('missMethyl')
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

ann <- getAnnotation('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
ref.tomerge=ann %>% data.frame %>%
    rownames_to_column(.,'ID') %>%
    dplyr::rename(.,CHR=chr,MAPINFO=pos) %>%
    select(ID,CHR,MAPINFO)

ls.pt=c('5e_08','0.000001','0.0001','0.01','0.1','1')


for (pt in ls.pt){
    f.path=paste0('result/EWAS_MDDprs_Shen/archiv/MDDprs_ewas_meta/MDDprs_pT',
                  pt,'_meta/mddprs_pT',pt,'_meta_RosieData.toptable.txt1.tbl')
    ewasResults=read.delim(f.path,sep='\t',header=T,stringsAsFactors=F)
    colnames(ewasResults)[1]='CpG'
    ewasResults=merge(ewasResults,ref.tomerge,by.x='CpG',by.y='ID',all.x=T)

    qman.dat=ewasResults %>%
        mutate(is.MHC = if_else(MAPINFO>25000000&MAPINFO<35000000&CHR==6,'Yes','No'))
    qman.dat$CHR=gsub('chr','',qman.dat$CHR)
    qman.dat$CHR=as.numeric(qman.dat$CHR)
    qman.dat = filter(qman.dat,is.MHC=='No') %>%
        mutate(P.adj=p.adjust(P.value,method='bonferroni'))
    
    
    sigCpGs=qman.dat$CpG[qman.dat$P.adj<0.05&qman.dat$is.MHC=='No']
    gst <- gometh(sig.cpg=sigCpGs, all.cpg=qman.dat$CpG, collection="GO",array.type = "EPIC", anno = anno)
    sig.gst=filter(gst,FDR<0.05)
    
    sig.gst=sig.gst[order(sig.gst$P.DE),]
    
    gst.region=goregion()
    
    write.table(sig.gst,file=paste0('result/GeneOntology/gstTable.MDDprs.pt',pt,'.txt'),col.names=T,row.names=F,quote=F,sep='\t')

}

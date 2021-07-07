
setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

library(dplyr)
library(stringr)


# Generation Scotland -----------------------------------------------------

discovery.summstats = list.files('result/EWAS_MDDprs_Shen/archiv/MDDprs_ewas_meta/',
                                 recursive = T,pattern = '_meta_RosieData.toptable.txt1.tbl',full.names = T) %>% 
    .[!grepl('.info$',.)] %>% 
    as.list %>% 
    pblapply(.,read.delim,sep='\t',header=T,stringsAsFactors=F) %>% 
    pblapply(.,rename,CpG=MarkerName) 


names(discovery.summstats) = list.files('result/EWAS_MDDprs_Shen/archiv/MDDprs_ewas_meta/',
                                        recursive = T,pattern = '_meta_RosieData.toptable.txt1.tbl',full.names = T) %>% 
    .[!grepl('.info$',.)] %>% 
    strsplit(.,'/') %>% 
    lapply(.,function(x) tail(x,1)) %>% 
    unlist %>% as.vector %>% 
    gsub('mddprs_pT','',.) %>% gsub('_meta_RosieData.toptable.txt1.tbl','',.)

# calculate lambda
discovery.lambda = discovery.summstats %>% 
    pblapply(.,mutate,chisq=qchisq(P.value,1,lower.tail=F)) %>% 
    pblapply(.,function(x) median(x$chisq)/qchisq(0.5,1)) %>% 
    bind_rows

for (pt in ls.pt){
    f.path=paste0('data/EWAS_MDDprs_Shen/MDDprs_ewas_meta/MDDprs_pT',pt,'_meta/mddprs_pT',pt,'_meta_RosieData.toptable.txt1.tbl')
    ewasResults=read.delim(f.path,sep='\t',header=T,stringsAsFactors=F)
    colnames(ewasResults)[1]='CpG'
    chisq <- qchisq(ewasResults$P.value, 1, lower.tail = F)
    tmp.lambda <- median(chisq) / qchisq(0.5,1)

    if (pt==ls.pt[1]){lambda.all=tmp.lambda}else{lambda.all=c(lambda.all,tmp.lambda)}

}


png(file=paste0(root, "-QQ-fine.png"))
qqman::qq(dat$P, main=title, pch='.', col="blue4", cex=1.0, las=1, mai=c(1.02,1,0.82,0.52))
dev.off()

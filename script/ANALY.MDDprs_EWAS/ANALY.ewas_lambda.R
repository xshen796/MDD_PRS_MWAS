
setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

library(dplyr)
library(stringr)
library(qqman)
library(pbapply)

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
    bind_rows %>% 
    t %>% as.data.frame %>% rownames_to_column('pT') %>% 
    mutate(pT=gsub('e_','e-',pT)) %>% mutate(pT=as.numeric(pT)) %>% rename(GS.lambda=V1)


# Replication (adults) ----------------------------------------------------

replication.summstats = list.files('result/EWAS_MDDprs_meta_LBC_ALSPAC/',
                                 recursive = T,pattern = 'REPmeta',full.names = T) %>% 
    .[!grepl('.info$',.)] %>% 
    as.list %>% 
    pblapply(.,read.delim,sep='\t',header=T,stringsAsFactors=F) %>% 
    pblapply(.,rename,CpG=MarkerName) 


names(replication.summstats) = list.files('result/EWAS_MDDprs_meta_LBC_ALSPAC/',
                                        recursive = T,pattern = 'REPmeta',full.names = T) %>% 
    .[!grepl('.info$',.)] %>% 
    strsplit(.,'/') %>% 
    lapply(.,function(x) tail(x,1)) %>% 
    unlist %>% as.vector %>% 
    gsub('REPmeta_EWAS_MDDprs_pT_','',.) %>% gsub('.metal.out1.tbl','',.)

# calculate lambda
replication.lambda = replication.summstats %>% 
    pblapply(.,mutate,chisq=qchisq(P.value,1,lower.tail=F)) %>% 
    pblapply(.,function(x) median(x$chisq)/qchisq(0.5,1)) %>% 
    bind_rows %>% 
    t %>% as.data.frame %>% rownames_to_column('pT') %>% 
    mutate(pT=gsub('_','',pT)) %>% mutate(pT=gsub('e.','e-',pT)) %>%
    mutate(pT=as.numeric(pT)) %>% rename(LBC_ALSPACadult.lambda=V1)


# Replication (adolescents) -----------------------------------------------

load_youthdat <- function(tmp.fname){
    load(tmp.fname)
    return(ewas_res)
}

youth.summstats = list.files('result/EWAS_MDDprs_ALSPAC/',
                                   recursive = T,pattern = 'MDD3_k_',full.names = T) %>% 
    .[!grepl('.info$|MDD3_k_102',.)] %>% 
    as.list %>% 
    pblapply(.,load_youthdat) %>% 
    pblapply(.,rename,CpG=probeID) 


names(youth.summstats) = list.files('result/EWAS_MDDprs_ALSPAC/',
                                          recursive = T,pattern = 'MDD3_k_',full.names = T) %>% 
    .[!grepl('.info$|MDD3_k_102',.)] %>% 
    strsplit(.,'/') %>% 
    lapply(.,function(x) tail(x,1)) %>% 
    unlist %>% as.vector %>% 
    gsub('MDD3_k_','',.) %>% gsub('_std_sex,age,child_smoke,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10_houseman_2021-02-12.Rdata','',.)

# calculate lambda
youth.lambda = youth.summstats %>% 
    pblapply(.,mutate,chisq=qchisq(P_VAL,1,lower.tail=F)) %>% 
    pblapply(.,function(x) median(x$chisq)/qchisq(0.5,1)) %>% 
    bind_rows %>% 
    t %>% as.data.frame %>% rownames_to_column('pT') %>% 
    mutate(pT=gsub('e','e-',pT)) %>% mutate(pT=as.numeric(pT)) %>% rename(ALSPACyouth.lambda=V1)


# Summarise lambda results ------------------------------------------------

all.lambda = left_join(discovery.lambda,replication.lambda,by='pT') %>% 
    left_join(.,youth.lambda,by='pT') %>% 
    .[order(.$pT,decreasing = F),] %>% 
    filter(pT!=0.2)

write.table(all.lambda,file='result/lambda_MWAS_discovery_repAdult_repYouth.txt',quote = F,col.names = T,row.names = F,sep='\t')



#qqman::qq(dat$P.value, main=title, pch='.', col="blue4", cex=1.0, las=1, mai=c(1.02,1,0.82,0.52))


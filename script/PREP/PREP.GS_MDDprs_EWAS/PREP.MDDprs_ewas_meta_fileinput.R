library(dplyr)
library(dplyr)
library(tidyverse)
library(pbapply)
library(data.table)
library(readr)

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_Shen/MDDprs_ewas_meta')

# Reformat files for metal ------------------------------------------------

wave1.fname = list.files('../MDDprs_ewas_wave1/',pattern = 'RosieData.toptable.txt',
                         recursive = T,full.names = T)
wave3.fname = list.files('../MDDprs_ewas_wave3/',pattern = 'RosieData.toptable.txt',
                         recursive = T,full.names = T)

wave1.newname = wave1.fname %>% 
    str_split(.,'/',simplify = T) %>% 
    .[,ncol(.)] 

wave3.newname = wave3.fname %>% 
    str_split(.,'/',simplify = T) %>% 
    .[,ncol(.)]

ls.input = data.frame(targetf=c(wave1.fname,wave3.fname),
                      outputf=c(wave1.newname,wave3.newname))

reformat_summstats <- function(ls.input){
    old.file = as.character(ls.input[1])
    output.file = as.character(ls.input[2])
    
    summstas = fread(old.file,header=T) %>% 
        select(ID,beta,se,P.Value,N)
    write_tsv(summstas,output.file)
}

ls.input %>% 
    pbapply(.,FUN=reformat_summstats,MARGIN = 1)


# Generate metal scripts --------------------------------------------------
setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

func=read.table('script/FUNs_meth/MDDprs_ewas_metal_func.txt',sep='\n',blank.lines.skip=F)

ls.pt=c('0.00000005','0.000001','0.0001','0.01','0.05','0.1','0.2','0.5','1')
ls.w1.file=paste0('mddprs_pT',ls.pt,'_wave1_RosieData.toptable.txt')
ls.w3.file=paste0('mddprs_pT',ls.pt,'_wave3_RosieData.toptable.txt')

ls.meta.path=paste0('result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/')
ls.meta.output=paste0('mddprs_pT',ls.pt,'_meta_RosieData')

ls.metal.script = paste0(ls.meta.path,'/MetalScript_',ls.pt,'.txt')

ls.forfunc=data.frame(ls.w1.file,ls.w3.file,ls.meta.output,ls.metal.script,stringsAsFactors=F)

generate_metal_input <- function(function.input,ls.file.input){
    anchor.func=function.input
    input.1=as.character(ls.file.input[1])
    input.2=as.character(ls.file.input[2])
    input.3=as.character(ls.file.input[3])
    output.path=as.character(ls.file.input[4])
    anchor.func$V1=gsub('INPUT1',input.1,anchor.func$V1)
    anchor.func$V1=gsub('INPUT2',input.2,anchor.func$V1)
    anchor.func$V1=gsub('INPUT3',input.3,anchor.func$V1)
    
    write.table(anchor.func,file=output.path,sep='\n',col.names=F,row.names=F,quote=F)
}

apply(X = ls.forfunc,MARGIN = 1,FUN = generate_metal_input,function.input=func)


#  Make a script to run metal script

ls.wd = list.files(ls.meta.path,pattern='MetalScript')

write.table(ls.wd,file='result/EWAS_MDDprs_Shen/metal_script.list.txt',sep='\n',col.names=F,row.names=F,quote=F)

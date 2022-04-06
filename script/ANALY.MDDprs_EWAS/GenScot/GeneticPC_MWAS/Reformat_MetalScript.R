library(dplyr)
library(data.table)
library(readr)
library(pbapply)
setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/result/EWAS_MDDprs_Shen/MWAS_geneticPC/')

# Reformat summstats for METAL --------------------------------------------

reformat_summstats <- function(tmp.file,is.pc=F){
  if (is.pc==F){
    summstats=fread(tmp.file,header=T) %>% 
      select(ID,beta,P.Value,se)
    write_tsv(summstats,file = paste0('./',gsub('toptable.txt','forMetal.gz',tmp.file)))
  }else if(is.pc==T){
    sumstats=fread(tmp.file,header=T)
    sumstats.pc = paste0('C',1:8) %>% as.list %>% 
      lapply(.,FUN = function(x){
        output = dplyr::filter(sumstats,param==x) %>% 
        select(ID,beta,P.Value=p,se) 
        return(output)
          })
    names(sumstats.pc)=paste0(gsub('linear.txt','',tmp.file),'geneticPC',1:8,'.forMetal.gz')
    names(sumstats.pc) %>%  
      sapply(.,function(x) write_tsv(sumstats.pc[[x]], 
                                       file = x) )
  }
}

# MDD PRS MWAS
ls.f = list.files('./',pattern = 'toptable')

ls.f %>% as.list %>% 
  pblapply(.,reformat_summstats)

# Genetic PCs' MWAS
ls.f.pc = list.files('./', pattern = 'linear.txt')

pc.dat = ls.f.pc %>% as.list %>% 
  pblapply(.,reformat_summstats,is.pc=T)
  

# Load and create metal scripts -------------------------------------------

metal.ref = 
  read.table('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS/script/ANALY.MDDprs_EWAS/GenScot/variantEWAS/summstats_metal_template',
             sep = '\n',stringsAsFactors = F,header = F,blank.lines.skip = F)

mk_metal <- function(ref.file,tmp.input,filename.1,filename.2,outputname){
  
  filename.1=tmp.input[1]
  filename.2=tmp.input[2]
  outputname=tmp.input[3]
  
  new.file = ref.file$V1 %>%
    gsub('FILE1',filename.1,.) %>%
    gsub('FILE2',filename.2,.) %>%
    gsub('OUTNAME',outputname,.) %>%
    as.data.frame
  return(new.file)
}

# Prepare files to process

ls.w1.f = list.files(path = './',pattern = '_w1',full.names = T) %>%
  .[grep('forMetal.gz',.)]
ls.w3.f = gsub('w1','w3',ls.w1.f)
ls.output.f = gsub('ewas_w1','ewas_meta',ls.w1.f) %>% gsub('.forMetal.gz','.metal.out',.)

input.ls = data.frame(file1=ls.w1.f,file2=ls.w3.f,outputfile=ls.output.f,stringsAsFactors = F)

metal.script = input.ls %>%
  apply(.,FUN=mk_metal,MARGIN = 1,
        ref.file=metal.ref)

names(metal.script)= ls.w1.f %>%
  gsub('ewas_w1_','',.) %>%
  gsub('.forMetal.txt','',.) %>%
  gsub('.//','',.)

lapply(names(metal.script),function(x) write.table(metal.script[[x]],
                                                   file=paste0('MetalScript_',x,'.txt'),
                                                   col.names=F,row.names=F,sep=' ',quote=F,append=F))

# Prepare a list of metal scripts to process --------------------------


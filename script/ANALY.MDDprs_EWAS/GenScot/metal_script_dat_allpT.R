library(dplyr)
library(data.table)
library(readr)
library(pbapply)
setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/')

# Reformat summstats for METAL --------------------------------------------

reformat_summstats <- function(tmp.file){
  summstats=fread(paste0('../MWAS_by_wave/',tmp.file),header=T) %>% 
    select(ID,beta,P.Value,se)
  write_tsv(summstats,path = paste0('./',gsub('toptable.txt','forMetal.txt',tmp.file)))
}

ls.f = list.files('../MWAS_by_wave',pattern = 'toptable')

ls.f %>% as.list %>% 
  pblapply(.,reformat_summstats)

# Load and create metal scripts -------------------------------------------

metal.ref = 
  read.table('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/ANALY.MDDprs_EWAS/GenScot/variantEWAS/summstats_metal_template',
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
  .[grep('forMetal.txt',.)]
ls.w3.f = gsub('w1','w3',ls.w1.f)
ls.output.f = gsub('ewas_w1','ewas_meta',ls.w1.f) %>% gsub('.forMetal.txt','.metal.out',.)

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


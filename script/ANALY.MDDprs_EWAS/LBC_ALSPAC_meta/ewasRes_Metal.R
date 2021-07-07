library(data.table)
library(dplyr)
library(pbapply)

setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

# Reformat ALSPAC results ---------------------------------------------------------------------------------------------

reformat_alspac_res <- function(tmp.pt){
      f.path=paste0('MR_meth_MDD/result/EWAS_MDDprs_ALSPAC/MDD3_m_',tmp.pt,
                    '_std_age,ante_ever_smoke,mpc1,mpc2,mpc3,mpc4,mpc5,mpc6,mpc7,mpc8,mpc9,mpc10_houseman_2021-02-12.Rdata')
      load(here::here(f.path))
      ewasResults=ewas_res 

      # Reformat data for meta-analysis
      reformatted.dat = ewasResults %>% 
            dplyr::select(ID=probeID,beta=BETA,se=SE,p=P_VAL,N=N_for_probe) %>%
            mutate(ID=as.character(ID)) 
      return(reformatted.dat)
}

ls.pt=c('5e08','1e06','1e04','1e03','1e02','5e02','1e01','5e01','1')

res.new = as.list(ls.pt) %>%
      pblapply(.,FUN=reformat_alspac_res)

names(res.new)=c('5e.08','1e.06','0.0001','0.001','0.01','0.05','0.1','0.5','1')

pblapply(as.list(names(res.new)), function(x) write.table( data.frame(res.new[[x]]), 
                                                           file = paste0('/exports/eddie/scratch/xshen33/MDDprs_EWAS/ALSPAC_EWAS_MDDprs_pT_',x,'_formeta.txt'), 
                                                           sep='\t', col.names=T, row.names=F,quote=F))


# Reformat LBC results --------------------------------------------------------------------------------------------

reformat_lbc_res <- function(tmp.pt){
      f.path=paste0('MR_meth_MDD/result/EWAS_MDDprs_LBC/MDDprs_Pt_',tmp.pt,
                    '_ewas.toptable.txt')
      ewasResults=read.delim(here::here(f.path),sep='\t',header=T,stringsAsFactors=F) 
      
      # Reformat data for meta-analysis
      reformatted.dat = ewasResults %>% 
            dplyr::select(ID=ID,beta=beta,se=se,p=P.Value,N=N)  
      return(reformatted.dat)
}

ls.pt=c('5e.08','1e.06','0.0001','0.001','0.01','0.05','0.1','0.5','1')

res.new = as.list(ls.pt) %>%
      pblapply(.,FUN=reformat_lbc_res)

names(res.new)=c('5e.08','1e.06','0.0001','0.001','0.01','0.05','0.1','0.5','1')

pblapply(as.list(names(res.new)), function(x) write.table( data.frame(res.new[[x]]), 
                                                           file = paste0('/exports/eddie/scratch/xshen33/MDDprs_EWAS/LBC_EWAS_MDDprs_pT_',x,'_formeta.txt'), 
                                                           sep='\t', col.names=T, row.names=F,quote=F))


# Make Metal scripts ----------------------------------------------------------------------------------------------
setwd('/exports/eddie/scratch/xshen33/MDDprs_EWAS')

metal.ref = 
      read.table('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/ANALY.MDDprs_EWAS/LBC_ALSPAC_meta/summstats_metal_template',
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

ls.f = list.files(path = './',full.names = F)
ls.alspac.f = ls.f[grep('ALSPAC',ls.f)]
ls.lbc.f = gsub('ALSPAC','LBC',ls.alspac.f)
ls.output.f = gsub('ALSPAC','REPmeta',ls.alspac.f) %>% gsub('_formeta.txt','.metal.out',.)

sum(!ls.alspac.f %in% ls.f[grep('ALSPAC',ls.f)]) # Check if name patterns changed 
# between rh and lh (need to be 0)

input.ls = data.frame(file1=ls.lbc.f,file2=ls.alspac.f,outputfile=ls.output.f,stringsAsFactors = F)

metal.script = input.ls %>%
      apply(.,FUN=mk_metal,MARGIN = 1,
            ref.file=metal.ref)

names(metal.script)= ls.alspac.f %>%
      gsub('ALSPAC_','',.) %>%
      gsub('_formeta.txt','',.)

lapply(names(metal.script),function(x) write.table(metal.script[[x]],
                                                   file=paste0('MetalScript_pT_',x,'.txt'),
                                                   col.names=F,row.names=F,sep=' ',quote=F,append=F))


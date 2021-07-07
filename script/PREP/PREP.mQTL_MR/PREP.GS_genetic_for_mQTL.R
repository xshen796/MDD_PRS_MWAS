library(dplyr)
library(data.table)
library(pbapply)
library(Hmisc)
library(readr)

setwd('/exports/eddie/scratch/xshen33/mQTL')

# Load data ---------------------------------------------------------------
ID.wave1=readRDS('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/sentrix.rds')
ID.wave3=readRDS('/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/wave3_sentrix.rds')

ID.all = rbind(ID.wave1,ID.wave3)

ID.plink = read.table('/exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA.fam')
ID.plink = ID.plink %>%
  select(id = V2, fid=V1)

# ID in genetic data to keep ----------------------------------------------

ID.genetic.tokeep = ID.all %>% merge(.,ID.plink,by='id',all.x=T) %>% 
  select(fid=fid,iid=id)


write.table(ID.genetic.tokeep,file='data/GS_ID_withMeth.txt',
            sep='\t',quote=F,row.names = F,col.names = F)

# Update ID ---------------------------------------------------------------
ID.toupdate = ID.all %>% merge(.,ID.plink,by='id',all.x=T) %>% 
  select(old.fid=fid,old.iid=id,new.fid=Sample_Sentrix_ID,new.iid=Sample_Sentrix_ID)

write.table(ID.toupdate,file='data/updateID_withMeth.txt',
            sep='\t',quote=F,row.names = F,col.names = F)


# Probe IDs for mQTL analysis ---------------------------------------------
ls.cpg = c('cg07519229','cg08116408','cg17925084',
           'cg03270340','cg14345882','cg08344181','cg17862947')

write.table(ls.cpg,file='data/cpg_ls_sigBrain.txt',
            sep='\n',quote=F,row.names = F,col.names = F)

ls.cpg = list.files('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/mQTL/mQTL_forMR',pattern='^mQTL') %>%
        gsub('mQTL.','',.) %>%
	gsub('.rds','',.)
write.table(ls.cpg,file='data/cpg_ls.txt',sep='\n',quote=F,row.names = F,col.names = F)


# Run mQTL analysis -------------------------------------------------------

system('qsub /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.mQTL_MR/job.GS_mQTL.sh')
# Extract results
system('bash /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.mQTL_MR/extract_mQTL.sh')


# Parse mQTL summstats ----------------------------------------------------

source.path = '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/mQTL/'
target.path = '/exports/eddie/scratch/xshen33/mQTL/mQTL_result/'

extract_qtl_bycpg <- function(path1,path2,fname_summ,tmp.suffix,manual_N){
  groot = fread(paste0(path1,'/',fname_summ),header=T,stringsAsFactors=F) %>%
    mutate(N=manual_N)
  groot.subset =
    as.list(unique(groot$Probe)) %>%
    pblapply(.,FUN=function(x) filter(groot,Probe==x))
  names(groot.subset)=unique(groot$Probe)
  as.list(names(groot.subset)) %>%
    pblapply(., function(x) write_tsv(groot.subset[[x]], 
                                      file = paste0(path2,x,'_',tmp.suffix,'.txt'), sep='\t',
                                      quote=F,row.names = F))
}

extract_qtl_bycpg(source.path,target.path,'mQTL_summstats_w1','w1',manual_N = 4890)
extract_qtl_bycpg(source.path,target.path,'mQTL_summstats_w2','w2',manual_N = 4323)


# Set up metal inputs for meta analysis -----------------------------------
setwd('/exports/eddie/scratch/xshen33/mQTL/mQTL_result/')

# Metal scripts
metal.template= 
  read.table('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/mQTL/metal_input_template',
             sep='\n',blank.lines.skip = F)

ls.txt = list.files(path='.') 
ls.input = data.frame(input1 = ls.txt[grep('_w1',ls.txt)],
                      input2 = ls.txt[grep('_w2',ls.txt)],
                      stringsAsFactors = F) %>%
  mutate(input3 = gsub('_w1.txt','_bothWaves',input1),
         script.file = gsub('_w1.txt','.MetalInput',input1))

make_metal_input <- function(ls.input,tmp.tmplate){
  summstats.1=as.character(ls.input[1])
  summstats.2=as.character(ls.input[2])
  outfile=as.character(ls.input[3])
  scriptfile=as.character(ls.input[4])
  
  new.tmplate = gsub('INPUT1',summstats.1,tmp.tmplate$V1) %>% 
    gsub('INPUT2',summstats.2,.) %>% 
    gsub('INPUT3',outfile,.)
  
  write.table(new.tmplate,file = scriptfile,quote = F,
              row.names = F,col.names = F,sep = '\n')
}

apply(ls.input,MARGIN = 1,FUN = make_metal_input,tmp.tmplate=metal.template)

# Reference to metal scripts
write.table(ls.input$script.file,file = 'ls.metal.script.txt',
            sep = '\n',row.names = F,col.names = F,quote=F)


# Run metal ---------------------------------------------------------------
system('qsub job.run_metal.sh')


# Reformat metal outputs for MR -------------------------------------------
setwd('/exports/eddie/scratch/xshen33/mQTL/mQTL_meta/')
ls.meta = list.files('.',patter='.tbl$')

reformat_metal <- function(tmp.fname,target.path){
  cpg = gsub('_bothWaves1.tbl','',tmp.fname)
  
  summstats = read.delim(tmp.fname,stringsAsFactors = F)
  colnames(summstats) = c('SNP','allele1','allele2','freq_a1','freq_se','beta_a1','se','pval',
                          'direction','hetisq','hetchisq','hetdf','hetpval')
  summstats = summstats %>%
    mutate(allele1=capitalize(allele1),allele2=capitalize(allele2),cpg=cpg,samplesize=9213)
  saveRDS(summstats,file=paste0(target.path,'/mQTL.',cpg,'.rds'))
}

as.list(ls.meta) %>%
  pblapply(.,reformat_metal,
           target.path='/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/mQTL/mQTL_forMR_GS/')

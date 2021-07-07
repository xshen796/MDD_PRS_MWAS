library(dplyr)              # Data management 
library(data.table)         # Data management
library(purrr)              # Data management

setwd('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/UKB_IM_summstats')


# Load data ---------------------------------------------------------------

# Load file names
ls.f = list.files('.') %>%
  .[!. %in% c('Processed','README','variants.txt')]

if (!file.exists('Processed')){
    dir.create('Processed')
}

# Load data to workspace
summstats = 
  as.list(ls.f) %>%
  lapply(.,FUN=fread,header = T,stringsAsFactors = F)

# Load ref data
ref=fread('variants.txt',header = T,stringsAsFactors = F) %>%
  .[,c('rsid','af','info')]


# Process -----------------------------------------------------------------

# Add info
summstats = summstats %>%
  lapply(.,FUN=merge,y=ref,by='rsid') %>%
  lapply(.,FUN=mutate, pval=10^(-`pval(-log10)`))
names(summstats)=ls.f

# Write new summstats
sapply(names(summstats), function(x) write.table(summstats[[x]],file=gsub('.txt','_InfoAdded',x),col.names=T,row.names=F,sep=' ',quote=F))

system('mv *_InfoAdded Processed')

# Make metal script
# Metal script template
metal.ref = 
  read.table('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.IM_summstats/PREP.IM_summstats_metal_template',
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

ls.f = list.files('Processed/',pattern = '_InfoAdded$',full.names = T)
ls.lh.f = ls.f[grep('_lh_',ls.f)]
ls.rh.f = gsub('_lh_','_rh_',ls.lh.f)
ls.output.f = gsub('_lh_','_',ls.lh.f) %>% gsub('_InfoAdded','.metal.out',.)

sum(!ls.rh.f %in% ls.f[grep('_rh_',ls.f)]) # Check if name patterns changed 
                                           # between rh and lh (need to be 0)

input.ls = data.frame(file1=ls.lh.f,file2=ls.rh.f,outputfile=ls.output.f,stringsAsFactors = F)

metal.script = input.ls %>%
  apply(.,FUN=mk_metal,MARGIN = 1,
         ref.file=metal.ref)

lapply(metal.script,write.table,file='/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.IM_summstats/IM_metal.input',col.names=F,row.names=F,sep=' ',quote=F,append=T)


# Run meta analysis -------------------------------------------------------

system('module load igmm/apps/metal')
system('metal /exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/script/PREP/PREP.IM_summstats/IM_metal.input')



# Add allele info ---------------------------------------------------------

ls.f = list.files('Processed/',pattern = 'metal.out1.tbl$',full.names = T)

summstats = 
  as.list(ls.f) %>%
  lapply(.,FUN=fread,header = T,stringsAsFactors = F) %>%
  lapply(.,FUN=merge,y=ref,by.x='MarkerName',by.y='rsid')

names(summstats)=ls.f
sapply(names(summstats), function(x) write.table(summstats[[x]],file=gsub('metal.out1.tbl','metal.InfoAdded.tbl',x),col.names=T,row.names=F,sep=' ',quote=F))



# Unilateral measures: rename columns to match metal outputs --------------
ls.f = list.files('Processed/',pattern = 'InfoAdded$',full.names = T) %>%
  .[!grepl('_lh_|_rh_',.)]

summstats = 
  as.list(ls.f) %>%
  lapply(.,FUN=fread,header = T,stringsAsFactors = F) %>%
  lapply(.,FUN=dplyr::rename,SNP=rsid,A1=a1,A2=a2,BETA=beta,SE=se,P=pval) %>%
  lapply(.,FUN=dplyr::mutate,Direction=NA) %>%
  lapply(.,FUN=dplyr::select,SNP,A1,A2,BETA,SE,P,Direction,af,info)

names(summstats)=ls.f
sapply(names(summstats), function(x) write.table(summstats[[x]],file=gsub('InfoAdded','InfoAdded.tbl',x),col.names=T,row.names=F,sep=' ',quote=F))

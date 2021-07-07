setwd('/exports/eddie/scratch/xshen33/GoDMC')
library('data.table')
library('dplyr')
library('pbapply')

# Download and unzip files 
# GoDMC mQTL data before removing any overlapping sites
# http://mqtldb.godmc.org.uk/downloads
#system('wget http://fileserve.mrcieu.ac.uk/mqtl/README')
#system('wget http://fileserve.mrcieu.ac.uk/mqtl/assoc_meta_all.csv.gz')
#system('gunzip assoc_meta_all.csv.gz')

mQTL.GoDMC = fread('assoc_meta_all_no_pgc.csv',showProgress=T,stringsAsFactors=F)
#cpg.450k = readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/450k_bothWaves/CpG_list.rds')
cpg.EPIC = readRDS('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/GS_methylation/EPIC_bothWaves/EPIC_CpGs.rds')

print('Master data loaded')


extract_file <- function(target.chr,masterdata,cpg.ref){
      tmp.cpg.ls = cpg.ref$CpG[cpg.ref$chr==target.chr]
      tmp.dat = masterdata[masterdata$cpg %in% tmp.cpg.ls,]
      tmp.savename = paste0('mQTL.chr',target.chr,'.rds')
      saveRDS(tmp.dat,file=tmp.savename)
}

print('File parsing starts')
chr.input = as.list(1:22)

lapply(chr.input,FUN=extract_file,masterdata=mQTL.GoDMC,cpg.ref=cpg.EPIC)

#setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')

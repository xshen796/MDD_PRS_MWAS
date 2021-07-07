library(data.table)
library(dplyr)
library(pbapply)
library(LDlinkR)
setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/')

# Load data -------------------------------------------------------------------------------------------------------

paper.102=fread('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_genetic_meth/ls.102SNP.GS.txt',header=T)
clumped.GS.summstats=fread('/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS_5e_08/MDDprs_5e_08.snp') %>%
      filter(P<5e-8)
GS.SNP = fread('/exports/eddie/scratch/xshen33/GS_PRS_loo/data/GS20K_HRC_0.8_GCTA.bim')
clumped.GS.summstats=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice') %>% 
  .[.$SNP %in% GS.SNP$V2,] %>% 
  filter(P<5e-8)


# Whether in MHC region --------------------------------------------------------------------------------
clumped.GS.summstats = clumped.GS.summstats %>%
      mutate(is.MHC = if_else(BP>25000000&BP<35000000&CHR==6,'Yes','No')) 


# Remove ambiguous variants -----------------------------------------------

old.snp.ambig = fread('/exports/eddie/scratch/xshen33/GS_PRS_loo/data/ambig_snp.txt',header=F) %>% .$V1

ls.f.allsnp = list.files('/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS/',pattern = '.log') %>% 
  gsub('loo_PRS','ls.SNP',.) %>% 
  gsub('.log','.txt',.)
ls.f.validsnp = list.files('/exports/eddie/scratch/xshen33/GS_PRS_loo/PRS/',pattern = '.all_score') %>% 
  gsub('loo_PRS','ls.SNP',.) %>% 
  gsub('.all_score','.txt',.)

ls.f.rerun = ls.f.allsnp %>% .[!. %in% ls.f.validsnp]

ls.snp.ambig = ls.f.rerun %>% paste0('/exports/eddie/scratch/xshen33/GS_PRS_loo/data/ls_SNP/',.) %>% as.list %>% 
  pblapply(.,fread,header=F) %>% 
  bind_rows %>% 
  .$V1 %>% 
  c(.,old.snp.ambig)

write.table(ls.snp.ambig,file='/exports/eddie/scratch/xshen33/GS_PRS_loo/data/ambig_snp.txt',sep='\n',row.names = F,col.names = F,quote=F)

# Find proxies in LD with the tops SNPs in paper --------------------------
find_proxy <- function(tmp.rs,tmp.summstats,omit.ls=''){
  tmp.proxy = LDproxy(tmp.rs,pop='CEU',r2d='r2',token = '81423ba617e2', file = FALSE) %>% 
    .[!as.character(.$RS_Number) %in% omit.ls,] %>% 
    filter(R2>0.1) #%>% select(RS_Number,Coord,Distance,R2)
  tmp.summstats = tmp.summstats %>%
        merge(.,tmp.proxy,by.x='SNP',by.y='RS_Number') %>%
        .[order(-log10(.$P),.$R2,decreasing=T),] %>%
        .[1,]
  return(tmp.summstats)
}

k = as.list(paper.102$MarkerName) %>%
pblapply(.,find_proxy,tmp.summstats=clumped.GS.summstats,omit.ls=ls.snp.ambig) %>%
bind_rows %>% 
  mutate(paperSNP=paper.102$MarkerName)


write.table(k$SNP[!is.na(k$R2)],file='data/GS_genetic_meth/SNP_gwsig.txt',sep='\n',row.names=F,col.names=F,quote=F)
write.table(k,file='data/GS_genetic_meth/SNP_variantMWAS.txt',sep='\t',row.names=F,col.names=T,quote=F)

# Load SNP list for PRS 5e-8  -------------------------------------------------------
# Run PRS for 5e-8 first so to select SNP list 
target.ls = k %>% filter(!is.na(R2))

create_newls <- function(tmp.target.ls,tmp.i,output.path){
   snp.tokeep = tmp.target.ls$SNP[tmp.i]
   new.ls = tmp.target.ls$SNP %>%
            .[. %in% snp.tokeep]
   new.fname = paste0('ls.SNP_',tmp.i,'.txt') %>%
               paste0(output.path,'/',.)
   write.table(new.ls,file=new.fname,sep='\n',col.names=F,row.names=F,quote=F)
}

as.list(1:nrow(target.ls)) %>%
   pblapply(.,FUN=create_newls,tmp.target.ls=target.ls,output.path='/exports/eddie/scratch/xshen33/GS_PRS_loo/data/ls_SNP/')

# Phenotype files                   --------------------------------------------------
fam.dat.GS=fread('/exports/igmm/eddie/GenScotDepression/data/genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA.fam',
                      header=F,stringsAsFactors=F)

pheno.GS <- fam.dat.GS %>%
  .[,1:2] %>%
  mutate(.,dummy_MDD=round(runif(nrow(.))))
colnames(pheno.GS)[1:2]=c('FID','IID')

write.table(pheno.GS,
            file='data/GS_loo_PRS/dummy.MDD.GS',
            row.names=F,col.names=T,sep='\t',quote=F)


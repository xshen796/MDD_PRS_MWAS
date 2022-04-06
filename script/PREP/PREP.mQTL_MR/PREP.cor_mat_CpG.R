

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MDD_PRS_MWAS')
library(dplyr)
library(data.table)
library(propagate)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# Target CpG list
hitsEwas=fread('result/EWAS_MDDprs_Shen/MDDprs_ewas_meta/no_geneticPC/ewas_meta_pT_0.00000005.metal.out1.tbl')
A <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
genes <- data.frame(ID=A$Name, geneSymbol=A$UCSC_RefGene_Name, CHR=A$chr, MAPINFO=A$pos, FEATURE=A$UCSC_RefGene_Group, CpGISLAND=A$Relation_to_Island)
hitsEwas = merge(hitsEwas,genes,by.x='MarkerName',by.y='ID',all.x=T)
hitsEwas$CHR = gsub('chr','',hitsEwas$CHR) %>% as.numeric
hits.MDDprsEwas = hitsEwas %>%
  mutate(.,p.corrected=p.adjust(`P.value`,method='bonferroni')) %>%
  filter(.,p.corrected<0.05)


# set DNAm paths
path.wave1='/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/Mvalues/chr/'  # data format: mvalues.chr1.rds
path.wave3='/exports/igmm/eddie/GenScotDepression/data/genscot/methylation/STRADL/wave3-final/chr/' # data format: as above

for (chr.no in 6){

    dat.1=readRDS(paste0(path.wave1,'mvalues.chr',chr.no,'.rds'))
    dat.2=readRDS(paste0(path.wave3,'mvalues.chr',chr.no,'.rds'))
    dat.1=data.frame(dat.1)
    dat.2=data.frame(dat.2)
    
    dat.1$CpG.id=rownames(dat.1)
    dat.2$CpG.id=rownames(dat.2)
    
    dat.both=merge(dat.1,dat.2,by='CpG.id',all=T) %>%
         .[.$CpG.id %in% hitsEwas$MarkerName,]
             
    
    #if (chr.no==1){
        dat.both.targetcpg=dat.both
    #}else{
    #    dat.both.targetcpg=rbind(dat.both.targetcpg,dat.both)
    #}
    cat(paste0(chr.no),'\n')
    rm(dat.1,dat.2)
            
}

    dat.for_cormat=dat.both.targetcpg %>% 
      dplyr::select(-CpG.id) %>% t(.)
    # RSC
    cormat=cor(dat.for_cormat,method='pearson',use='pairwise.complete.obs')
    rownames(cormat)<-colnames(cormat)<-dat.both.targetcpg$CpG.id

    output.file=paste0('data/Meth_Cor_Mat/cormat_MDDprs_tophits_CpG_R2.rds')
    saveRDS(cormat,file=output.file)

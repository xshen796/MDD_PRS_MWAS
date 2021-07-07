setwd("/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD")
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(data.table)


# Create a set of significant CpGs from meta summary stats  ----------------------------------------------
mddprs.ewas=fread('result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.metal.out1.tbl',stringsAsFactors=F,fill=T)
mddprs.ewas=mddprs.ewas[,c('MarkerName','Effect','StdErr','P-value','Direction')]
colnames(mddprs.ewas)=c('CpG','beta.meta_prs','se.meta_prs','P.Value.meta_prs','Direction.meta_prs')
mddprs.ewas$p.corrected = p.adjust(mddprs.ewas$P.Value.meta_prs,method='bonferroni')
mddprs.ewas.sig=filter(mddprs.ewas,p.corrected<0.05)
MDDprs_assoc_cpg=data.frame(cpg=as.character(mddprs.ewas.sig$CpG),stringsAsFactors=F)

saveRDS(MDDprs_assoc_cpg,file='data/methTraining/cpgList_MDDprs_assoc.rds')

# Create a list of CpGs free of MDDprs assoc   -----------------------------------------------------------

# Read ewas results
ewasResults=read.delim('result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.metal.out1.tbl',sep='\t',header=T,stringsAsFactors=F)
colnames(ewasResults)[1]='CpG'
# Annotate location
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ref.tomerge <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
# Tidy up format
ewasResults=merge(ewasResults,ref.tomerge,by.x='CpG',by.y='ID',all.x=T) 
ewasResults$CHR=gsub('chr','',ewasResults$CHR)
ewasResults$CHR=as.numeric(ewasResults$CHR)

input.ewas.sig <- ewasResults %>%  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(MAPINFO)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(ewasResults, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, MAPINFO) %>%
  mutate( BPcum=MAPINFO+tot)
  
BP_window=250000

# Select significant regions
input.ewas.sig$p.corrected = p.adjust(input.ewas.sig$P.value,method='bonferroni')
ewas.sig=filter(input.ewas.sig,p.corrected<0.05)


generate_BP_toavoid <- function(fulldat,sigdat,map_window){
      upper.range=max(fulldat$BPcum)
      lower.range=min(fulldat$BPcum)
      
      prc.BP.avoid=sigdat$BPcum
      range.BP.avoid <- 
          sapply(prc.BP.avoid,function(x) c((x-map_window):(x+map_window))) %>%
          as.vector %>%
          unique %>%
          (function(x) x[(x<=upper.range)&(x>=lower.range)])
      
      return(range.BP.avoid)   
}

for (i in unique(ewas.sig$CHR)){
    dat.1=filter(input.ewas.sig,CHR==i)
    dat.2=filter(ewas.sig,CHR==i)
    
    if(i==unique(ewas.sig$CHR)[1]){
        total.BP.avoid=generate_BP_toavoid(dat.1,dat.2,BP_window)
    }else{
        total.BP.avoid=c(total.BP.avoid,generate_BP_toavoid(dat.1,dat.2,BP_window))
    }
}

# MDDprs-free CpG
MDDprs_free.cpg=data.frame(cpg=input.ewas.sig$CpG[!input.ewas.sig$BPcum %in% total.BP.avoid],stringsAsFactors=F)
saveRDS(MDDprs_free.cpg,file='data/methTraining/cpgList_MDDprs_NOassoc.rds')


# Transform summary stats of MDDprs EWAS to weights format  ------------------------------------------------------------
mddprs.ewas=fread('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/EWAS_MDDprs_fromKathryn/wave3/mdd-pgrs-wave3-snp_ch/mdd_SNPCH_prs_5e_08/toptables/mdd_SNPCH_prs_5e_08.PR=SNP_CH.PC=20.MV=mvals-noST-TT.txt',stringsAsFactors=T,fill=T)
mddprs.ewas=mddprs.ewas[,c('ID','beta_SE')]
colnames(mddprs.ewas)=c('coef.name','coef.value')

write.table(mddprs.ewas,file='/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/data/weights_wave2_WholeGenome_WholeSample_450k/weights_mddprs_ewas.txt',col.names=F,row.names=F,quote=F,sep=' ')
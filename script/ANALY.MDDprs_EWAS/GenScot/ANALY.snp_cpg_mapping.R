library(dplyr)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(TwoSampleMR)

setwd('/gpfs/igmmfs01/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')

# Find hits ------------------------------------------------------------------------------------------------

get_hits <- function(dat,chr="CHR", bp="BP", snp="SNP", p="P",tophits.bpwindow=3000000){
      
    data_for_fig=dat
    colnames(data_for_fig)[colnames(data_for_fig)==chr]='CHR'
    colnames(data_for_fig)[colnames(data_for_fig)==bp]='BP'
    colnames(data_for_fig)[colnames(data_for_fig)==p]='p'
    colnames(data_for_fig)[colnames(data_for_fig)==snp]='SNP'

    # Prepare data   -----------
    don <- data_for_fig %>% 
      
      # Compute chromosome size
      group_by(CHR) %>% 
      summarise(chr_len=max(BP)) %>% 
      
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      
      # Add this info to the initial dataset
      left_join(data_for_fig, ., by=c("CHR"="CHR")) %>%
      
      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate(BPcum=BP+tot)
    
     
    # Prepare X axis   ------------
    axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    
    # find top hits to annotate  ----------
    sig.rows = filter(don,p<5e-8)
    if (nrow(sig.rows)>0){
        sig.rows=sig.rows[order(sig.rows$BPcum),]
    
        sig.rows$sigblock.no=99999
        
        for (i in 1:nrow(sig.rows)){
            if (i==1){
                sig.rows$sigblock.no[i]=1
            }else{
                if (sig.rows$BPcum[i]>sig.rows$BPcum[i-1]+tophits.bpwindow){
                      sig.rows$sigblock.no[i]=sig.rows$sigblock.no[i-1]+1
                    }else{
                      sig.rows$sigblock.no[i]=sig.rows$sigblock.no[i-1]
                    }
            }
        }
        
        sig.rows$sigblock.no[sig.rows$sigblock.no==99999]=NA
        
        for (i in unique(sig.rows$sigblock.no)){
            tmp.block=filter(sig.rows,sigblock.no==i)
            tmp.top=tmp.block[tmp.block$p==min(tmp.block$p),]
            if (i==unique(sig.rows$sigblock.no)[1]){
                top.hits=tmp.top
            }else{
                top.hits=rbind(top.hits,tmp.top)
            }
        
        }
        
        cpg.ls.return=top.hits$SNP
    }else{cpg.ls.return=NA}
    
    return(cpg.ls.return)
}



ref=read.delim('result/EWAS_MDDprs_Shen/MDDprs_ewas_wave1/MDDprs_pT1_wave1/mddprs_pT1_wave1_RosieData.toptable.txt',sep='\t',header=T,stringsAsFactors=F)
ref.tomerge=ref[,c('ID','CHR','MAPINFO')]

ls.SNP=list.files(pattern='^mdd_SNPCH_rs',path='result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd')

cpg.tokeep=NA
for (snp in ls.SNP){
    f.path=paste0('result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd/',snp,'/',snp,'.metal.out1.tbl')
    ewasResults=read.delim(f.path,sep='\t',header=T,stringsAsFactors=F)
    colnames(ewasResults)[1]='CpG'
    ewasResults=merge(ewasResults,ref.tomerge,by.x='CpG',by.y='ID',all.x=T)

    qman.dat=ewasResults
    qman.dat$CHR=gsub('chr','',qman.dat$CHR)
    qman.dat$CHR=as.numeric(qman.dat$CHR)
    qman.dat=filter(qman.dat,!is.na(CpG),!is.na(CHR),!is.na(MAPINFO),!is.na(P.value))

    tmp.cpg=get_hits(dat=qman.dat,chr="CHR", bp="MAPINFO", snp="CpG", p="P.value")
    
    cpg.tokeep=c(cpg.tokeep,tmp.cpg)

}

save(cpg.tokeep,file='result/EWAS_MDDprs_fromKathryn/cpg_snp_mapping.RData')


# Association matrix for the hits ------------------------------------------------------------------------------------------------

load('result/EWAS_MDDprs_fromKathryn/cpg_snp_mapping.RData')
#cpg.tokeep=cpg.tokeep[!is.na(cpg.tokeep)]
cpg.tokeep=list.files('data/mQTL/mQTL_forMR/') %>%
  gsub('mQTL.','',.)%>%
  gsub('.rds','',.)

#  cpg.tokeep=fig.dat.pT.5e_08$CpG[p.adjust(fig.dat.pT.5e_08$P.value,method='bonferroni')<0.05]

ls.SNP=list.files(pattern='^mdd_SNPCH_rs',path='result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd') %>%
  strsplit(.,split = '_') %>% unlist %>% .[grep('^rs',.)]

rs.chr_bp=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',
  header=F,sep=' ',stringsAsFactors=F)

ls.SNP.102 = rs.chr_bp[rs.chr_bp$V3 %in% ls.SNP,] %>% .[!duplicated(.$V3),]

ls.SNP=list.files(pattern='^mdd_SNPCH_rs',path='result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd') 

for (snp in ls.SNP){
    f.path=paste0('result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd/',snp,'/',snp,'.metal.out1.tbl')
    ewasResults=read.delim(f.path,sep='\t',header=T,stringsAsFactors=F) %>%
                mutate(P.adj=p.adjust(P.value,method='bonferroni'))
    colnames(ewasResults)[1]='CpG'
    
    tmp.p.corr = ewasResults[ewasResults$CpG %in% cpg.tokeep,c('CpG','P.adj')]
    tmp.p = ewasResults[ewasResults$CpG %in% cpg.tokeep,c('CpG','P.value')]
    tmp.beta = ewasResults[ewasResults$CpG %in% cpg.tokeep,c('CpG','Effect')]
    
    if(snp==ls.SNP[1]){
        p.mat.corr=tmp.p.corr
        p.mat=tmp.p
        beta.mat=tmp.beta
    }else{
        p.mat.corr=merge(p.mat.corr,tmp.p.corr,by='CpG',all=T)
        p.mat=merge(p.mat,tmp.p,by='CpG',all=T)
        beta.mat=merge(beta.mat,tmp.beta,by='CpG',all=T)
    }
    cat(paste0(match(snp,ls.SNP),'\n'))
}

#colnames(p.mat.corr)[2:ncol(p.mat.corr)]=gsub('mdd_SNPCH_','',ls.SNP)
colnames(p.mat)[2:ncol(p.mat)]=gsub('mdd_SNPCH_','',ls.SNP)
colnames(beta.mat)[2:ncol(beta.mat)]=gsub('mdd_SNPCH_','',ls.SNP)

# Add refs
ls.102hits=read.table('data/GS_genetic_meth/ls.102hits_MDDgwas.txt',
                      header = F,stringsAsFactors=F)
colnames(ls.102hits)=c('CHR','MarkerName','BP')

A <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
genes <- data.frame(ID=A$Name, CHR=A$chr, MAPINFO=A$pos)
ls.CpG = genes[genes$ID %in% beta.mat$CpG,]


save(ls.CpG,ls.SNP.102,p.mat,beta.mat,file='result/EWAS_MDDprs_fromKathryn/cpg_snp_mapping_pNbeta.RData')


# Top hits PRS: cpg location and snp location ----------------------------------------------------------------------------
# Load data
MDDgwas=fread('/exports/igmm/eddie/GenScotDepression/data/ukb/summary_stats/PGC/MDD_Howard2019/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta')
ls.102hits=read.table('data/GS_genetic_meth/ls.102hits_MDDgwas.txt',
                      header = F,stringsAsFactors=F)
colnames(ls.102hits)=c('CHR','MarkerName','BP')
hitsEwas=fread('result/EWAS_MDDprs_fromKathryn//meta-pgrs-mdd/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.metal.out1.tbl')
A <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
genes <- data.frame(ID=A$Name, geneSymbol=A$UCSC_RefGene_Name, CHR=A$chr, MAPINFO=A$pos, FEATURE=A$UCSC_RefGene_Group, CpGISLAND=A$Relation_to_Island)
hitsEwas = merge(hitsEwas,genes,by.x='MarkerName',by.y='ID',all.x=T)

# pos for gwas and MDDprs ewas
hits.MDDgwas = merge(MDDgwas,ls.102hits,by='MarkerName',all.y=T)

hitsEwas$CHR = gsub('chr','',hitsEwas$CHR) %>% as.numeric

# Log distance to the nearest gwas hits
distance_nearest_hit <- function(target.cpg,gwas.hits){
  tmp.chrcpg=target.cpg$CHR
  tmp.BPcpg=target.cpg$MAPINFO
  
  tmp.gwas.hits = filter(gwas.hits,CHR==tmp.chrcpg)
  nearest.distance = min(abs(tmp.gwas.hits$BP-tmp.BPcpg),na.rm = T)
  return(nearest.distance)
}

cpg.list = hitsEwas %>% transpose
cpg.distance =
  lapply(cpg.list,FUN=distance_nearest_hit,gwas.hits=hits.MDDgwas) %>%
  unlist
hitsEwas = mutate(hitsEwas,
                  distance_to_MDDgwas_hits=cpg.distance)
hitsEwas$distance_to_MDDgwas_hits[hitsEwas$distance_to_MDDgwas_hits==Inf]=NA

hits.MDDprsEwas = hitsEwas %>%
  mutate(.,p.corrected=p.adjust(.$`P-value`,method='bonferroni')) %>%
  filter(.,p.corrected<0.05)

# Plot distance vs p-value/beta
fig.p_distance = ggplot(hitsEwas, aes(x=distance_to_MDDgwas_hits, 
                                             y=-log10(`P-value`))) +
  geom_point(alpha=0.5)+
  geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "red")+
  xlab('Minimal distance to MDD GWAS loci')+
  ylab('-log10(p)')

ggsave(fig.p_distance,file='Figs/p_distance.png',height=8,width=8,dpi=200)

# Interval of 95% distribution
quantile(hits.MDDprsEwas$distance_to_MDDgwas_hits, q = c(0.025, 0.5, 0.975))/1000

quantile(hitsEwas$distance_to_MDDgwas_hits, q = c(0.025, 0.5, 0.975))/1000



# Number of sig CpGs for each SNP -----------------------------------------

ls.SNP=list.files(pattern='^mdd_SNPCH_rs',path='result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd') %>%
  .[grep('^mdd_SNPCH_rs',.)]

rs.chr_bp=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',
                header=F,sep=' ',stringsAsFactors=F)
A <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
cpg.ref <- data.frame(ID=A$Name, CHR=A$chr, MAPINFO=A$pos,stringsAsFactors = F)

tmp.snp.name = ls.SNP %>% strsplit(.,split = '_') %>% 
  unlist %>% .[grep('^rs',.)]
ls.SNP.102 =rs.chr_bp[rs.chr_bp$V3 %in% tmp.snp.name,] %>% .[!duplicated(.$V3),]
colnames(ls.SNP.102)=c('CHR','BP','SNP')


for (snp in ls.SNP){
  f.path=paste0('result/EWAS_MDDprs_fromKathryn/meta-pgrs-mdd/',snp,'/',snp,'.metal.out1.tbl')
  ewasResults=read.delim(f.path,sep='\t',header=T,stringsAsFactors=F) %>%
    mutate(P.adj=p.adjust(P.value,method='bonferroni'))
  colnames(ewasResults)[1]='CpG'
  
  tmp.snp.name = snp %>% strsplit(.,split = '_') %>% 
    unlist %>% .[grep('^rs',.)]
  info.snp = ls.SNP.102[grep(tmp.snp.name,ls.SNP.102$SNP),]
  
  tmp.n.sig = sum(ewasResults$P.adj<0.05)
  tmp.sig.result = filter(ewasResults,P.adj<0.05,!is.na(P.adj)) %>%
    merge(.,cpg.ref,by.x='CpG',by.y='ID')
  snp.cis.range = seq(from=info.snp$BP-3000000,to=info.snp$BP+3000000,by=1)
  tmp.n.cis_sig = filter(tmp.sig.result,CHR==info.snp$CHR)
    tmp.n.cis_sig = sum(tmp.n.cis_sig$MAPINFO %in% snp.cis.range)
  tmp.n.trans_sig = tmp.n.sig - tmp.n.cis_sig
  tmp.n.trans_chr = length(unique(tmp.sig.result$CHR))-1
  tmp.trans_chrs = paste0(unique(tmp.sig.result$CHR[tmp.sig.result$CHR!=info.snp$CHR]),collapse = ';')

  tmp.summary = data.frame(info.snp,Nsig=tmp.n.sig,
                           NsigTrans=tmp.n.trans_sig,NchrTrans=tmp.n.trans_chr,
                           chrTrans=tmp.trans_chrs)
  if(snp==ls.SNP[1]){
    Nsig.summary=tmp.summary
  }else{
    Nsig.summary=rbind(Nsig.summary,tmp.summary)
  }
  cat(paste0(match(snp,ls.SNP),'\n'))
}

#rm(list=ls(pattern = '^tmp.'))
save(Nsig.summary,file='result/EWAS_MDDprs_fromKathryn/NsigCpG_perSNP.RData')


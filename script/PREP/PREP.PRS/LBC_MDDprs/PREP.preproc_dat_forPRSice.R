library(data.table)
library(dplyr)

setwd('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD')

# This summary stats contains 9765 variants that passed genome-wide significant threshold
# Caution! - some other versions may contain only two cohorts
# This summary stats was processed by Dave Howard
# QC: MAF>0.005, INFO>0.1; Ncase>0.8

# load summary stats ------------------------------------------------------
summstats.Howard2019.3cohorts=
  fread('/exports/igmm/eddie/GenScotDepression/data/ukb/summary_stats/PGC/MDD_Howard2019/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta',
  header=T,sep='\t',stringsAsFactors=F)

rs.chr_bp=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',
  header=F,sep=' ',stringsAsFactors=F)

# hg18=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg18/snp130_chr_bp_rs.txt',
#            header=F,sep=' ',stringsAsFactors=F)



# Reformat for PRSice 2.0 -------------------------------------------------
colnames(summstats.Howard2019.3cohorts)=c('SNP','A1','A2','Freq1','FreqSE','MinFreq','MaxFreq','BETA','SE','P','Direction')

summstats.Howard2019.3cohorts <- 
  summstats.Howard2019.3cohorts %>%
  mutate(.,A1=toupper(A1)) %>%
  mutate(.,A2=toupper(A2))

# Add chr and bp info
colnames(rs.chr_bp)=c('CHR','BP','SNP')
rs.chr_bp$CHR=gsub('chr','',rs.chr_bp$CHR)
summstats.Howard2019.3cohorts=merge(summstats.Howard2019.3cohorts,rs.chr_bp,by='SNP',all.x=T)


# Save summary stats ------------------------------------------------------
write.table(summstats.Howard2019.3cohorts,
            file='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/summstats/23andme_PGCNoGS_UKB_Aug8_FinalGCldsc_3cohort1.meta.forPRSice',
            row.names=F,col.names=T,sep='\t',quote=F)

# Create dummy phenotypes for LBC -----------------------------------------
fam.dat.LBC1921=fread('/exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/LBC1921_rsID.fam',
                      header=F,stringsAsFactors=F)
fam.dat.LBC1936=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/LBC_genetic/LBC1936/Genotyped/LBC36_clean_231009.fam',
                      header=F,stringsAsFactors=F)

pheno.LBC1921 <- fam.dat.LBC1921 %>%
  .[,1:2] %>%
  mutate(.,dummy_MDD=round(runif(nrow(.))))
colnames(pheno.LBC1921)[1:2]=c('FID','IID')

pheno.LBC1936 <- fam.dat.LBC1936 %>%
  .[,1:2] %>%
  mutate(.,dummy_MDD=round(runif(nrow(.))))
colnames(pheno.LBC1936)[1:2]=c('FID','IID')

pheno.both = rbind(pheno.LBC1921,pheno.LBC1936)


# Save phenotype files for LBC --------------------------------------------
write.table(pheno.LBC1921,
            file='data/LBC_phenotype/dummy.MDD.LBC1921',
            row.names=F,col.names=T,sep='\t',quote=F)
write.table(pheno.LBC1936,
            file='data/LBC_phenotype/dummy.MDD.LBC1936',
            row.names=F,col.names=T,sep='\t',quote=F)
write.table(pheno.both,
            file='data/LBC_phenotype/dummy.MDD.LBCboth',
            row.names=F,col.names=T,sep='\t',quote=F)



# Prepare LBC genetic data ------------------------------------------------

## Imputed data

# Create a mapping file from chr:bp to rsname
LBC1921.bim=fread('/exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/LBC1921_plink.bim',header=F,stringsAsFactors=F)

LBC1921.bim$gp.oldName = paste0(LBC1921.bim$V1,':',LBC1921.bim$V4)

LBC1936.bim=fread('/exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/LBC1936_plink.bim',header=F,stringsAsFactors=F)

LBC1936.bim$gp.oldName = paste0(LBC1936.bim$V1,':',LBC1936.bim$V4)


rs.chr_bp=fread('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ucsc_annotation/hg19/snp151Common_chr_bp_rs.txt',
  header=F,sep=' ',stringsAsFactors=F)
rs.chr_bp$gp.oldName = paste0(gsub('chr','',rs.chr_bp$V1),':',rs.chr_bp$V2)

LBC1921.merged = merge(LBC1921.bim,rs.chr_bp,by='gp.oldName') %>%
	.[!duplicated(.$gp.oldName),] %>%
	.[,c(3,10)]

LBC1936.merged = merge(LBC1936.bim,rs.chr_bp,by='gp.oldName') %>%
	.[!duplicated(.$gp.oldName),] %>%
	.[,c(3,10)]


write.table(LBC1921.merged,
            file='/exports/eddie/scratch/xshen33/LBC_genetic/LBC1921/HRCv1.1/LBC1921_HRCv1.1_VCF/chr_bp_toRSname.txt',
            row.names=F,col.names=F,sep='\t',quote=F)

write.table(LBC1936.merged,
            file='/exports/eddie/scratch/xshen33/LBC_genetic/LBC1936/HRCv1.1/LBC1936_HRCv1.1_VCF/chr_bp_toRSname.txt',
            row.names=F,col.names=F,sep='\t',quote=F)

# For both genotyped and imputed data: create a mapping file from old subj ID to IDs in phenotype file -- LBC1921
LBC1921.blood_pheno=fread('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/LBC_phenotype/LBC1921_XueyiShen_EWAS_DepressionDisorderPGR_AM_12NOV2020.csv',
                          header=T,stringsAsFactors=F)

linkage.dat=LBC1921.blood_pheno[,2:1] %>% filter(bloodnum!='',studyno!)
linkage.dat=linkage.dat[,c(1,1,2,2)]

write.table(linkage.dat,file='data/LBC_phenotype/subjID_genetic_pheno.txt',
            row.names=F,col.names=F,sep='\t',quote=F)
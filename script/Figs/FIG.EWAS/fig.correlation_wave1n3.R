# Run under igmm/apps/R/3.6.1
# by XShen

# Basic settings
setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/STRADL/personal_folder/DNAm_projects/MR_meth_MDD/')

library(ggrepel)
library(dplyr)
library(ggpubr)
options(bitmapType='cairo')

# Plot function  ------------------------------------------------------------------------------------------

generate_corr_plot <- function(ewasRes.w1,ewasRes.w3,target.col,input.lab){

    w1.tmp=ewasRes.w1[,c('ID',target.col)]
    w3.tmp=ewasRes.w3[,c('ID',target.col)]
    
    total.dat=merge(w1.tmp,w3.tmp,by='ID',all.x=T)
    colnames(total.dat)=c('ID','w1.stats','w3.stats')
    
    fig.dat=data.frame(w1.stats=total.dat$w1.stats,w3.stats=total.dat$w3.stats)
    cor.dat=cor(fig.dat,use="pairwise.complete.obs", method="pearson")
    r.beta=round(cor.dat,digits = 3)
    lab.tmp=input.lab
    
    fig=ggplot(fig.dat, aes(x=fig.dat[,1], 
                                      y=fig.dat[,2])) +
    geom_point(alpha=0.7,size=1)+
    annotate(geom="text",label=paste0('r = ',r.beta), x=-0.02, y=0.01)+
    xlab(paste0(lab.tmp,': Wave 1'))+
    ylab(paste0(lab.tmp,': Wave 2'))+
    geom_smooth(method=lm,size=0.3)
    
    return(fig)
}

# Make plots -----------------------------------------------------------------------------------------------

ls.pt=c('1','0.1','0.01','0.0001','0.000001')

for (pt in ls.pt){
      file.w1=paste0('data/EWAS_MDDprs_Shen/MDDprs_ewas_wave1/MDDprs_pT',pt,'_wave1/mddprs_pT',pt,'_wave1_RosieData.toptable.txt')
      file.w3=paste0('data/EWAS_MDDprs_Shen/MDDprs_ewas_wave3/MDDprs_pT',pt,'_wave3/mddprs_pT',pt,'_wave3_RosieData.toptable.txt')
      ewasResults.w1=read.delim(file.w1,sep='\t',header=T,stringsAsFactors=F)
      ewasResults.w1$p.log=-log10(ewasResults.w1$P.Value)
      ewasResults.w3=read.delim(file.w3,sep='\t',header=T,stringsAsFactors=F)
      ewasResults.w3$p.log=-log10(ewasResults.w3$P.Value)

      fig.tmp.beta=generate_corr_plot(ewasResults.w1,ewasResults.w3,'t','T-value')
      fig.tmp.p=generate_corr_plot(ewasResults.w1,ewasResults.w3,'p.log','-log10(p-value)')

      eval(parse(text=paste0('fig.beta.pT.',pt,'=fig.tmp.beta')))
      eval(parse(text=paste0('fig.p.pT.',pt,'=fig.tmp.p')))
}

# ewas results from Kathryn

ewasResults.w1=read.delim('data/EWAS_MDDprs_fromKathryn/meta_wave1_wave3/mdd-pgrs-snp_ch/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.w1.txt',sep=' ',header=T,stringsAsFactors=F)
ewasResults.w1$p.log=-log10(ewasResults.w1$p)
ewasResults.w3=read.delim('data/EWAS_MDDprs_fromKathryn/meta_wave1_wave3/mdd-pgrs-snp_ch/mdd_SNPCH_prs_5e_08/mdd_SNPCH_prs_5e_08.w3.txt',sep=' ',header=T,stringsAsFactors=F)
ewasResults.w3$p.log=-log10(ewasResults.w3$p)

colnames(ewasResults.w1)[1]='ID'
colnames(ewasResults.w3)[1]='ID'

fig.beta.pT.5e_08=generate_corr_plot(ewasResults.w1,ewasResults.w3,'beta','Beta')
fig.p.pT.5e_08=generate_corr_plot(ewasResults.w1,ewasResults.w3,'p.log','-log10(p-value)')


# Arrange plots -------------------------------------------------------------------------------------------
labs.plot=paste0('pT<',c('1','0.1','0.01','0.0001','0.000001','5e-08'))

fig.beta=ggarrange(fig.beta.pT.1,fig.beta.pT.0.1,fig.beta.pT.0.01,fig.beta.pT.0.0001,fig.beta.pT.0.000001,fig.beta.pT.5e_08,legend = NULL,widths = 3,heights = 2,labels=labs.plot)
fig.p=ggarrange(fig.p.pT.1,fig.p.pT.0.1,fig.p.pT.0.01,fig.p.pT.0.0001,fig.p.pT.0.000001,fig.p.pT.5e_08,legend = NULL,widths = 3,heights = 2,labels=labs.plot)

png(file = "Figs/Correlation.tvalue.png", bg = "transparent", width = 2180, height = 1580, units = "px",res=100)
fig.beta
dev.off()

png(file = "Figs/Correlation.p.png", bg = "transparent", width = 2180, height = 1580, units = "px",res=100)
fig.p
dev.off()
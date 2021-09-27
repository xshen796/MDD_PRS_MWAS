library(data.table)
library(dplyr)
library(knitr)
library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)
library(optparse)
library(ggpubr)
library(harrypotter)
options(bitmapType='cairo') # Specific to Eddie - to enable writing figures

# Load arguments ----------------------------------------------------------
# Environment: R 3.6.1
parse <- OptionParser()

option_list <- list(
   make_option('--MRexposure', type='character', help="Intermediate files of exposure data for MR", action='store'),
   make_option('--MRoutcome', type='character', help='Intermediate files of outcome data for MR', action='store'),
   make_option('--outcomeName', type='character', help='Name of outcome phenotype', action='store'),
   make_option('--outputfig', type='character', help='Folder for QC figures', action='store'),
   make_option('--outputtable', type='character', help='Summary file for MR analysis', action='store'),
   make_option('--saveHarmonisedData', type='character', help='Save harmonised data for analysis', action='store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

exposure.path=opt$MRexposure
outcome.path=opt$MRoutcome
outcome.name=opt$outcomeName
ofig.path=opt$outputfig
otable.path=opt$outputtable
saveDat=opt$saveHarmonisedData

# Basic settings ----------------------------------------------------------
dir.create(file.path(ofig.path), showWarnings = FALSE)
#dir.create(file.path(otable.path), showWarnings = FALSE)


# MR analysis -------------------------------------------------------------
ls.exposure=list.files(path = exposure.path, pattern = '.exposure_dat') %>% 
   gsub('.exposure_dat','',.) 
ls.exposure.filename = list.files(path = exposure.path, pattern = '.exposure_dat',full.names=T)
ls.outcome=list.files(path = outcome.path, pattern = '.outcome_dat') %>% 
   gsub('.outcome_dat','',.) 
ls.outcome.filename = list.files(path = outcome.path, pattern = '.outcome_dat',full.names=T)

# Check if exposure and outcome match
if (length(ls.exposure) != length(ls.outcome)){
   cat('Numbers of exposure and outcome variables do NOT match\n')
   cat('Analysing those that match\n\n')
   
   ls.exposure=ls.exposure[ls.exposure %in% ls.outcome]
   ls.exposure.filename=ls.exposure.filename[grep(paste(ls.exposure,collapse = '|'),
                                                  ls.exposure.filename)]
}

ls.exposure = ls.exposure[order(ls.exposure)]
ls.exposure.filename = ls.exposure.filename[order(ls.exposure.filename)]
ls.outcome = ls.outcome[order(ls.outcome)]
ls.outcome.filename = ls.outcome.filename[order(ls.outcome.filename)]


for (f in 1:length(ls.exposure)){
      # file names
      exposure.fname=ls.exposure.filename[f]
      traitA=ls.exposure[f]
      traitB=outcome.name
      
      ## exposure dat
      exposure_dat=read.table(exposure.fname,header=T,sep='\t')
      ## outcome dat
      f.outcome = grep(traitA,ls.outcome.filename) %>% ls.outcome.filename[.]
      outcome_dat <- read_outcome_data(file=f.outcome,sep='\t')
      
      # Harmonise the exposure and outcome data
      dat <- harmonise_data(exposure_dat, outcome_dat)
      
      if (as.logical(saveDat)){
            write.table(dat,file=paste0(ofig.path,'/',traitA,'_',traitB,"_HarmonisedDat.csv"),col.names=T,quote=F,row.names = F,sep='\t')
      }
      
      
      cat(paste0('\nHarmonised data: Nsnp = ',nrow(dat),'\n\n'))
      if (nrow(dat)>3){
         # analyse
         MR=mr(dat = dat,method_list =c("mr_ivw","mr_egger_regression","mr_weighted_median"))
         het = mr_heterogeneity(dat)
         plt= mr_pleiotropy_test(dat)
         res_single = mr_singlesnp(dat)
         res_loo = mr_leaveoneout(dat)
         res.direction = directionality_test(dat)
         
         tmp.mr_res=data.frame(MR,mr_pleiotropy_test(dat)[5:7],
                               mr_heterogeneity(dat)[2,6:8],
                               res.direction[,5:8])
         
         # summarise stats
         if (f==1){
            MRsummary=tmp.mr_res
         }else{
            MRsummary=rbind(MRsummary,tmp.mr_res)
         }
         ## visualise
         mr_report(dat,output_path=ofig.path)
         
         fig.scatter = mr_scatter_plot(MR, dat) %>% .[[1]]+scale_colour_hp_d(option = "ronweasley2", name = "MR method")
         fig.forest = mr_forest_plot(res_single) %>% .[[1]]+scale_colour_hp_d(option = "ronweasley2", name = "MR method")
         fig.loo = mr_leaveoneout_plot(res_loo) %>% .[[1]]
         fig.funnel = mr_funnel_plot(res_single) %>% .[[1]]
        
         fig.total = ggarrange(fig.scatter,fig.funnel,fig.forest,fig.loo,ncol = 2,nrow=2,
                               labels = c('Scatter plot','Funnel plot','Forest plot','Leave-one-out plot'),
                               widths = c(1,1),heights = c(1,1.5),align = 'v')
         fig.scatter=fig.scatter
	 save(fig.scatter,file=paste0(ofig.path,'/',traitA,'_',traitB,"_scatterplot.RData"))
         ggsave(fig.total, file=paste0(ofig.path,'/',traitA,'_',traitB,"_plots.png"), width=10, height=15,dpi=300)
	 
         system(paste0('echo ',traitA,' against ',traitB,'   analysis  done >> XS.log'))
      }

}

write.table(MRsummary,file=otable.path,col.names = T,row.names=F,sep="\t",quot=F)
rm(list=ls())

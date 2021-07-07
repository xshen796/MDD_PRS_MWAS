library(dplyr)
library(pbapply)
library(performance)


# **** GS sample **** -----------------------------------------------------

mdd.gs=readRDS('data/methTraining/DNAm_training_pheno_MDD.rds') %>%
      mutate(MDD_status=as.factor(MDD_status)) %>%
      mutate(sex=as.factor(sex))

# p+T method: pT < 0.1
prs.ls=read.table('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/score_thresholds.txt',stringsAsFactors=F,header=F)
colnames(prs.ls)=c('fname','label.pT','pT')
prs.ls$label.pT=as.character(prs.ls$label.pT)
prs.ls$label.pT=format(prs.ls$pT,scientific=F,drop0trailing=T)

ls.f=list.files('/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/gs_genetic/MDD_prs_GS/related/',pattern = '.profile',full.names = T)

prs = as.list(ls.f) %>%
      lapply(.,fread,stringsAsFactors=F,header=T) %>%
      lapply(.,select,IID,starts_with('SCORE')) %>%
      Reduce(function(x,y) left_join(x,y,by='IID'),.) %>%
      mutate_at(vars(matches('SCORE')), ~scale(.))


linkage.f = read.csv('/exports/igmm/eddie/GenScotDepression/shen/ActiveProject/Genetic/MR_meth_MDD/data/GS_genetic_meth/gwas_methid_link.csv',
                     stringsAsFactors=F) %>%
      select(ends_with('id'))

target.dat = left_join(mdd.gs,linkage.f,by=c('ID'='meth_id')) %>%
      left_join(.,prs,by=c('gwas_id'='IID'))

model.input = paste0('C',1:10,collapse = '+') %>% 
      paste0('MDD_status~age+sex+',.,'+',colnames(target.dat)[grep('SCORE',colnames(target.dat))]) %>%
      .[order(.)]

fit = model.input %>%
      as.list %>%
      lapply(.,as.formula) %>%
      lapply(.,glm,data=target.dat,family='binomial') 
fit.0 = glm(MDD_status~age+sex+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10,data=target.dat,family='binomial')


summ.r2 = fit %>%
      lapply(.,function(x) r2(x)[[1]]-r2(fit.0)[[1]]) %>%
      bind_rows
summ.res = fit %>%
      lapply(.,summary) %>%
      lapply(.,function(x) tail(x$coefficients,1)%>%as.data.frame) %>%
      bind_rows %>%
      mutate(R2=summ.r2$`Tjur's R2`)


# Save table ------------------------------------------------------------------------------------------------------

write.table(summ.res,file='result/MDDprs_predicion.txt',quote=F,row.names=F,col.names=T,sep='\t')



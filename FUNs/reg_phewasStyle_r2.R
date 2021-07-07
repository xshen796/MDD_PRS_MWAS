reg_phewasStyle <- function (ls.models,dat_short,dat_long,correctByFactor=F){

library(modEvA)
      
      run_model <- function(ls.mod,dat_short,dat_long){
            # define vars
            dep = as.character(ls.mod[1])
            factor = as.character(ls.mod[2])
            covs = as.character(ls.mod[3])
            mod.type = as.character(ls.mod[4])
            
            # run model
            if (mod.type=='lme'){
                  # model.expression
                  fh_r=1:(nrow(dat_long)/2)
                  sh_r=(nrow(dat_long)/2+1):nrow(dat_long)
                  
                  dat_long[fh_r,dep]=scale(dat_long[fh_r,dep])
                  dat_long[sh_r,dep]=scale(dat_long[sh_r,dep])
                  if(is.numeric(dat_long[,factor])){
                        dat_long[fh_r,factor]=scale(dat_long[fh_r,factor])
                        dat_long[sh_r,factor]=scale(dat_long[sh_r,factor])
                  }
                  
                  mod=paste0(dep,'~',covs,'+',factor)
                  
                  fit=lme(as.formula(as.character(mod)),data=dat_long,na.action=na.exclude,random=~1|f.eid,control=lmeControl(opt = "optim"))
                  table = summary(fit)$tTable
                  tarv = nrow(table)
                  stats = table[tarv,c(1,2,4,5)]
                  mod_result = data.frame(mod_name=paste0(dep,'~',factor),t(stats))
                  colnames(mod_result)[2:5]=c('beta','std','t.value','p.value')
                  
            }else{
                  all.var=c(dep,factor,strsplit(covs,split='+',fixed=T)[[1]])
                  dat_short=dat_short[complete.cases(dat_short[,all.var]),]
                  N=nrow(dat_short)
                  dep.dat=dat_short[,dep]            
                  if (length(table(dep.dat))==2){
                        mod=paste0(dep,'~',covs,'+scale(',factor,')')
                        mod.h0=paste0(dep,'~',covs)
                        fit=glm(as.formula(as.character(mod)),data=dat_short,na.action=na.exclude,family = 'binomial')
                        fit.h0=glm(as.formula(as.character(mod.h0)),data=dat_short,na.action=na.exclude,family = 'binomial')
                  }else{
                        mod=paste0('scale(',dep,')~',covs,'+scale(',factor,')')
                        mod.h0=paste0(dep,'~',covs)
                        fit=glm(as.formula(as.character(mod)),data=dat_short,na.action=na.exclude)
                        fit.h0=glm(as.formula(as.character(mod.h0)),data=dat_short,na.action=na.exclude)
                  }            
                  
                  table = summary(fit)$coefficients
                  tarv = nrow(table)
                  stats = table[tarv,c(1:4)]
                  
                  r2=Dsquared(model = fit, adjust = F)-Dsquared(model = fit.h0, adjust = F)

                  mod_result = data.frame(t(stats),r2,N)
                  #colnames(mod_result)[2:7]=c('beta','std','t.value','p.value','R2','N')
            }
            
            return(mod_result)
      }

      tmp.result = pbapply(X = ls.models,MARGIN = 1,FUN = run_model,dat_short=targetdata,dat_long=dat_long)
      result.table =  matrix(unlist(tmp.result),ncol=6,byrow = T)
      result.table = data.frame(result.table)
      colnames(result.table) = c('beta','std','t.value','p.value','R2','N')
      #result.table$mod_name = lapply(tmp.result, function(l) as.character(l[[1]]))
      result.table = data.frame(ls.models[,1:2],result.table,stringsAsFactors = F)
      
      ls.factor=unique(ls.models$p_batch)
      result.table$p.corrected=99999
      if (correctByFactor==T){
            for (f in ls.factor){
                  loc=grep(f,ls.models$p_batch)
                  result.table$p.corrected[loc]=p.adjust(result.table$p.value[loc],method='fdr')
            }
      }else{
            result.table$p.corrected=p.adjust(result.table$p.value,method='fdr')
            
      }
      
      return(result.table)
}
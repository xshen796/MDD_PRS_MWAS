UKB_lme_reg = function(targetdata_long,targetdata_short,
                       no_fvar_long,no_fvar_short,no_dvar_long,no_dvar_short,vinteraction,targetvar){

# input data in the format of: 1 - f.eid, 
#                              2 to 1+no_fvar_xx - fvar,
#                              2+no_fvar_xx to 1+no_fvar_xx+no_dvar_xx - dvar,
#                              2+no_fvar_xx+no_dvar_xx to end - covar

      
#################################################################################
#########################       both lme and lm        ##########################
#################################################################################
  
  
   ############################      lme model       #############################
  if (sum(!is.na(targetdata_long))>0){
    # set up parameters for all the models
    fvars=colnames(targetdata_long)[2:(no_fvar_long+1)]
    dvars=colnames(targetdata_long)[(no_fvar_long+2):(no_fvar_long+no_dvar_long+1)]
    if ((no_fvar_long+no_dvar_long+2)>ncol(targetdata_long)){covars=''}else{
          covars=colnames(targetdata_long)[(no_fvar_long+no_dvar_long+2):ncol(targetdata_long)]}   
    
    mod_build = cbind(rep(fvars,each=length(dvars)),dvars)
    covar_build = rep(covars, each=nrow(mod_build))
    covar_build = matrix(covar_build,ncol=length(covars),nrow=length(covar_build)/length(covars))
    mod_build = cbind(mod_build,covar_build)
    colnames(mod_build)[1:2]=c('fvar','dvar')
    colnames(mod_build)[3:ncol(mod_build)]='covar'
    
    cat('Done: setting up lme models \n')
    
    # scale numeric dependent variables and factors
    fh_r=1:(nrow(targetdata_long)/2)
    sh_r=(nrow(targetdata_long)/2+1):nrow(targetdata_long)
    
    targetdata_long[fh_r,dvars]=scale(targetdata_long[fh_r,dvars])
    targetdata_long[sh_r,dvars]=scale(targetdata_long[sh_r,dvars])
    for (i in 1:length(fvars)){
          col_n=fvars[i]
          if(is.numeric(targetdata_long[,col_n])){
                targetdata_long[fh_r,col_n]=scale(targetdata_long[fh_r,col_n])
                targetdata_long[sh_r,col_n]=scale(targetdata_long[sh_r,col_n])
          }
    }
    
    cat('Done: scaling numeric variables \n')
        
    # build up regression expressions and run
    n_mod = nrow(mod_build)
    vn_mod = ncol(mod_build)
    i_mod = 0
    
    cat(paste0(sprintf('%d lme models are being tested', n_mod),'\n\n'))
    
    for (i in 1:n_mod){
      if (!vinteraction){
          f = mod_build[i,1]
          d = mod_build[i,2]
          c = paste0(mod_build[i,3:vn_mod],collapse='+')
      }else{
          f = paste0(mod_build[i,1],'*',mod_build[i,3])
          d = mod_build[i,2]
          c = paste0(mod_build[i,4:vn_mod],collapse='+')
      }
      
      if (i==1){
            cat(paste0('Your first lme model is:\n',d,'~',c,'+',f,',random=~1|f.eid,data=targetdata_long \n\n'))
      }
      
      exp_mod = paste0('fit=lme(',d,'~',c,'+',f,',random=~1|f.eid,data=targetdata_long,na.action=na.exclude,control=lmeControl(opt = "optim"))')
      eval(parse(text=exp_mod))
      table = summary(fit)$tTable
      tarv = nrow(table)-targetvar+1
      beta = table[tarv,1]
      std = table[tarv,2]
      t.value = table[tarv,4]
      p.value = table[tarv,5]
      
      coef.tarv = rownames(table)[tarv]
      if (i==1){
            cat(paste0('The coefficients of ',coef.tarv,' is going to be reported\n'))
      }
      mod_result = data.frame(mod_name=paste0(f,'-',d),beta=beta,std=std,t.value=t.value,p.value=p.value)
          if (i==1){results.long=mod_result}else{
                results.long=rbind(results.long,mod_result)
          }
      
      # report progression
      i_mod=i_mod+1
      progres_mod = data.frame(perct=paste(seq(10,100,by=10),'%'),
                               perct.count=seq(0.1,1,0.1),
                               n.count=round(seq(0.1,1,0.1)*n_mod))
      if (sum(i_mod==progres_mod$n.count)==1){cat(paste0(sprintf('%s', progres_mod$perct[i_mod==progres_mod$n.count]),' - '))}   
      
    }
    cat('Done: lme models \n\n\n')
    
  }
    
    ############################       lm model       #############################
  
  if (sum(!is.na(targetdata_short))>0){
    # set up parameters for all the models
    fvars=colnames(targetdata_short)[2:(no_fvar_short+1)]
    dvars=colnames(targetdata_short)[(no_fvar_short+2):(no_fvar_short+no_dvar_short+1)]
    if ((no_fvar_short+no_dvar_short+2)>ncol(targetdata_short)){covars=''}else{
    covars=colnames(targetdata_short)[(no_fvar_short+no_dvar_short+2):ncol(targetdata_short)]}
    
    mod_build = cbind(rep(fvars,each=length(dvars)),dvars)
    covar_build = rep(covars, each=nrow(mod_build))
    covar_build = matrix(covar_build,ncol=length(covars),nrow=length(covar_build)/length(covars))
    mod_build = cbind(mod_build,covar_build)
    colnames(mod_build)[1:2]=c('fvar','dvar')
    colnames(mod_build)[3:ncol(mod_build)]='covar'
    
    cat('Done: setting up glm models \n')
        
    # scale numeric dependent variables and factors
    #targetdata_short[,dvars]=scale(targetdata_short[,dvars])
    for (i in 1:length(dvars)){
          col_n=dvars[i]
          if(is.numeric(targetdata_short[,col_n])){
                targetdata_short[,col_n]=scale(targetdata_short[,col_n])
          }
    }
    for (i in 1:length(fvars)){
          col_n=fvars[i]
          if(is.numeric(targetdata_short[,col_n])){
                targetdata_short[,col_n]=scale(targetdata_short[,col_n])
          }
    }
    
    cat('Done: scaling numeric variables \n')
        
    # build up regression expressions and run
    n_mod = nrow(mod_build)
    vn_mod = ncol(mod_build)
    cat(paste0(sprintf('%d glm models are being tested', n_mod),'\n\n'))
    i_mod=0
    for (i in 1:n_mod){
      if (!vinteraction){
            f = mod_build[i,1]
            d = mod_build[i,2]
            c = paste0(mod_build[i,3:vn_mod],collapse='+')
      }else{
            f = paste0(mod_build[i,1],'*',mod_build[i,3])
            d = mod_build[i,2]
            c = paste0(mod_build[i,4:vn_mod],collapse='+')
      }
 
      if (i==1){
            cat(paste0('Your first lm model is:  \n',d,'~',c,'+',f,',data=targetdata_short \n\n'))
      }
      
      if (is.numeric(targetdata_short[,d])){
            exp_mod = paste0('fit=glm(',d,'~',c,'+',f,',data=targetdata_short,na.action=na.exclude)')
      }else if(is.factor(targetdata_short[,d])){
            exp_mod = paste0('fit=glm(',d,'~',c,'+',f,',data=targetdata_short,na.action=na.exclude,family=binomial(link=\'logit\'))')
      }
      
      
      eval(parse(text=exp_mod))
      table = summary(fit)$coefficients
      tarv = nrow(table)-targetvar+1
      if (is.numeric(targetdata_short[,d])){
            beta = table[tarv,1]
      }else if(is.factor(targetdata_short[,d])){
            beta = exp(table[tarv,1])
      }
      #beta = table[tarv,1]
      std = table[tarv,2]
      t.value = table[tarv,3]
      p.value = table[tarv,4]
      
      coef.tarv = rownames(table)[tarv]
      if (i==1){
            cat(paste0('The coefficients of ',coef.tarv,' is going to be reported\n'))
      }
      
      mod_result = data.frame(mod_name=paste0(f,'-',d),beta=beta,std=std,t.value=t.value,p.value=p.value)
            if (i==1){results.short=mod_result}else{
              results.short=rbind(results.short,mod_result)
            }
      
      # report progression
      i_mod=i_mod+1
      progres_mod = data.frame(perct=paste(seq(10,100,by=10),'%'),
                               perct.count=seq(0.1,1,0.1),
                               n.count=round(seq(0.1,1,0.1)*n_mod))
      if (sum(i_mod==progres_mod$n.count)==1){cat(paste0(sprintf('%s', progres_mod$perct[i_mod==progres_mod$n.count]),' - '))}   
      
    }   
    cat('Done: lm models \n\n\n')
        
  }
    
    
    #########################     Clean up results      ###########################
   
    # based on the factors, split results.all into different objects, store into a list
    if (sum(!is.na(targetdata_long))>0&sum(!is.na(targetdata_short))>0){
      results.all=rbind(results.long,results.short)
    }else if(sum(!is.na(targetdata_long))==0&sum(!is.na(targetdata_short))>0){
      results.all=results.short
    }else if(sum(!is.na(targetdata_long))>0&sum(!is.na(targetdata_short))==0){
      results.all=results.long
    }
  
    n_dvar=length(dvars)

    for (i in 1:length(fvars)){
      if(i==1){loc_results_f=grep(fvars[i],results.all[,1])}
      else{loc_results_f=rbind(loc_results_f,grep(fvars[i],results.all[,1]))}
    }
    
    if (length(fvars)==1){
          list_loc=1:length(loc_results_f)
          for (i in list_loc){
                var_name = fvars[i]
                eval(parse(text=paste0('result.',var_name,'=results.all')))
                eval(parse(text=paste0('result.',var_name,'$p.corrected=p.adjust(result.',
                                       var_name,'[,5], method=\'fdr\')')))
          }
    }else{list_loc=1:nrow(loc_results_f)
          for (i in list_loc){
            loc = loc_results_f[i,]
            var_name = fvars[i]
            eval(parse(text=paste0('result.',var_name,'=results.all[loc,]')))
            eval(parse(text=paste0('result.',var_name,'$p.corrected=p.adjust(result.',
                                   var_name,'[,5], method=\'fdr\')')))
          }
    }
    result_names = paste0('result.',fvars)
          # remove row names
          temp.expre = paste0('rownames(',result_names,')=NULL')
          eval(parse(text=temp.expre))
    temp.expre = paste0(result_names,'=',result_names,collapse=', ')
    eval(parse(text=paste0('output=list(',temp.expre,')')))


      cat('Results are ready')

    
return(output)
}
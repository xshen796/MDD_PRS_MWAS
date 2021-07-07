add_p_batch <- function (ls.models,ls.factor,ls.dep){
      
      ls.models$p_batch=99999
      target.model=ls.models
      cate.no = 1
      for (fac in ls.factor){
            for (dep in ls.dep){
                  loc = grepl(dep,target.model$dependent)&grepl(fac,target.model$factor)
                  target.model$p_batch[loc]=cate.no
                  cate.no=cate.no+1
            }
      }
      return(target.model)

}
clump_bychr <- function(dat,chr.no,tophits.bpwindow=1000000){
  sig.rows=dat %>% filter(.,CHR==chr.no) %>% mutate(BP=MAPINFO)
  if (nrow(sig.rows)>0){
    sig.rows=sig.rows[order(sig.rows$BP),]
    
    sig.rows$sigblock.no=99999
    
    for (i in 1:nrow(sig.rows)){
      if (i==1){
        sig.rows$sigblock.no[i]=1
      }else{
        if (sig.rows$BP[i]>sig.rows$BP[i-1]+tophits.bpwindow){
          sig.rows$sigblock.no[i]=sig.rows$sigblock.no[i-1]+1
        }else{
          sig.rows$sigblock.no[i]=sig.rows$sigblock.no[i-1]
        }
      }
    }
    sig.rows$sigblock.no[sig.rows$sigblock.no==99999]=NA
    
    for (i in unique(sig.rows$sigblock.no)){
      tmp.block=filter(sig.rows,sigblock.no==i)
      tmp.top=tmp.block[tmp.block$P.value==min(tmp.block$P.value),]
      if (i==unique(sig.rows$sigblock.no)[1]){
        top.hits=tmp.top
      }else{
        top.hits=rbind(top.hits,tmp.top)
      }
    }
    
  }
  return(top.hits)
}


clump_bychr_withcorr <- function(dat,chr.no,tophits.bpwindow=250000,ref.cpgCorr){
  sig.rows=dat[dat$chr==chr.no,]
  if (nrow(sig.rows)>0){
    sig.rows=sig.rows[order(sig.rows$BP),]
    
    sig.rows$sigblock.no=99999
    
    for (i in 1:nrow(sig.rows)){
      if (i<=3){
        sig.rows$sigblock.no[i]=1
      }else{
        pair.r = abs(ref.cpgCorr[sig.rows$cpg[(i-3):i],sig.rows$cpg[(i-3):i]]) %>% 
          .[2:nrow(.),ncol(.)] # Correlation matrix of CpG(i) with three previous CpGs
        if ((sig.rows$BP[i]>sig.rows$BP[i-1]+tophits.bpwindow)|(sum(pair.r<0.3)==length(pair.r))){
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
    
  }
  return(top.hits)
}

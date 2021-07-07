m_plot<-function(dat,chr="CHR", bp="BP", snp="SNP", p="P",shape_sig=T,
                 category.input,labels_annot=F,outputpath=NA,plot_title=NA,tophits_annot=T,man_annot=NA,tophits.bpwindow=3000000,
                 add_category_name=F,fig_size=c(13,8),y_lim=NULL,sig_line='adj'){

# Reformat: change col names  ----------------------------------------------------------
data_for_fig=dat
colnames(data_for_fig)[colnames(data_for_fig)==chr]='CHR'
colnames(data_for_fig)[colnames(data_for_fig)==bp]='BP'
colnames(data_for_fig)[colnames(data_for_fig)==p]='p'
colnames(data_for_fig)[colnames(data_for_fig)==snp]='SNP'

# Prepare data   -----------------------------------------------------------------------
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
  mutate( BPcum=BP+tot) %>%
  
  # Add adjusted p values
  mutate( p.bonf = p.adjust(p,method='bonferroni')) 


# Find p thresholds -------------------------------------------------------

if (sum(don$p.bonf<0.05)>0){
  pT.bonf = max(don$p[don$p.bonf<0.05], na.rm = T)  
}else{
  pT.bonf = 0.05/nrow(don)
}

p.annot = if_else(sig_line=='gw',5e-8,
                  if_else(sig_line=='adj',pT.bonf,5e-8))

# Prepare X axis   ---------------------------------------------------------------------
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# find top hits to annotate  -----------------------------------------------------------
sig.rows = filter(don,p.bonf<0.05,!is.na(BPcum))
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
    
    don$is_annotate=''
    if (tophits_annot==T){
      # Add highlight and annotation information
      #mutate( is_highlight=ifelse(SNP %in% SNPOfInterest, "yes", "no")) %>%
       don$is_annotate[don$SNP %in% top.hits$SNP]="yes"
       don$is_annotate[!don$SNP %in% top.hits$SNP]="no"
    }else{don$is_annotate='no'}

}else{don$is_annotate='no'}

# Annotate manually marked CpGs  -------------------------------------------------------
if(sum(!is.na(man_annot))>0){
   don = don %>%
	mutate( is_highlight=ifelse(SNP %in% man_annot, "yes", "no"))
}


# Make the plot  -----------------------------------------------------------------------
fig<-ggplot(don, aes(x=BPcum, y=-log10(p))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.3) +
    scale_color_manual(values = rep(c("dodgerblue2", "gray66"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0),limits=y_lim) +     # remove space between plot area and x axis
    
    # axis titles
    ylab('-log10(p)')+
    xlab('Chromosome')+

    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size=5)
    )
    
if (!is.na(plot_title)){
    fig=fig+
    ggtitle(plot_title)
}else{}


fig=fig+
geom_hline(yintercept=-log10(p.annot), linetype="dashed", 
            color = "red", size=0.2)

if ((labels_annot==T)|(tophits_annot==T)){
    fig=fig+
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=1)
}else if(sum(!is.na(man_annot))>0){
    fig = fig +
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_highlight=="yes"), aes(label=SNP), size=1)
}
    

    if (!is.na(outputpath)){
           ggsave(plot = fig,filename = outputpath,
           width = fig_size[1], height = fig_size[2], 
           units = 'cm', dpi = 100)
    }else{}

return(fig)

}

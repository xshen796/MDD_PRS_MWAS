library(ggraph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)


circular_graph_mr <- function(target.res,output.file=NA,target.col='p.adj'){
  target.result = target.res
  # create a data frame giving the hierarchical structure of your individuals
  set.seed(1234)
  d1 <- data.frame(from="origin", to=c('Brain structures',
                                       'CpG'))
  d2 <- data.frame(from='CpG', to=unique(target.result$exposure))
  d3 <- data.frame(from='Brain structures', to=unique(target.result$outcome))
  hierarchy <- rbind(d1, d2, d3)
  
  
  connect <- target.result %>%
    select(from=exposure,to=outcome,value=log.p.val)
  
  connect$value[target.result[,target.col]>0.05]=NA
  connect = filter(connect,!is.na(value))
  
  # create a vertices data.frame. One line per object of our hierarchy
  vertices  <-  data.frame(
    name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) , 
    value = 1,
    stringsAsFactors = F
  ) 
  rownames(vertices)=vertices$name
  new.size  <-  rbind(plyr::count(connect$from),plyr::count(connect$to))
  colnames(new.size) <- c('name','value')
  new.size$name=as.character(new.size$name)
  vertices[new.size$name,'value']=vertices[new.size$name,'value']+new.size$value
  
  # Let's add a column with the group of each name. It will be useful later to color points
  vertices$group  <-  hierarchy$from[match(vertices$name, hierarchy$to)]
  
  #Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
  #calculate the ANGLE of the labels
  vertices$id <- NA
  myleaves <- which(is.na( match(vertices$name, hierarchy$from) ))
  nleaves <- length(myleaves)
  vertices$id[ myleaves ] <- seq(1:nleaves)
  vertices$angle <- 160 -1.3*nrow(target.result) - 360 * vertices$id / nleaves
  
  # calculate the alignment of labels: right or left
  # If I am on the left part of the plot, my labels have currently an angle < -60
  vertices$hjust <- ifelse( (vertices$angle < -90)&(vertices$angle > -270), 1.2, -0.2)
  
  # flip angle BY to make them readable
  vertices$angle <- ifelse((vertices$angle < -90)&(vertices$angle > -270), vertices$angle+180, vertices$angle)
  
  # Create a graph object
  mygraph <- graph_from_data_frame(hierarchy, vertices=vertices)
  
  # The connection object must refer to the ids of the leaves:
  from  <-  match( connect$from, vertices$name)
  to  <-  match( connect$to, vertices$name)
  
  p=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    geom_conn_bundle(data = get_con(from = from, to = to), 
                     width=0.6, alpha=0.4, colour="lightslategrey") +
    geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, 
                       label=name, angle=angle, hjust=hjust,colour=group), 
                   size=3.5, alpha=1) +
    theme_void() 
    #theme(legend.position = "none")
  p.label = p + 
    geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, size=value), alpha=0.5) +
    scale_colour_manual(values= c('lightcoral','dodgerblue3'),name = 'Exposure/Outcome') +
    scale_size_continuous(range = c(1,9) ,name = 'No. of sig effects') +
    expand_limits(x = c(-2, 2), y = c(-2, 2))
  
  if (!is.na(output.file)){
    ggsave(p.label,filename = output.file,width = 23,height = 18,
           units = 'cm',dpi = 300,device = 'png')
  }
  
  return(p.label)
} 
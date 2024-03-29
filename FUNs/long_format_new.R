long_format <- function(targetdata,cols_nonimg,cols_img){

    #=====================================================================================
    # prep non-imaging phenotypes
    targetdata_nonimg = targetdata[,c('f.eid',cols_nonimg)] ## non-imaging phenotypes
    targetdata_img=targetdata[,dat_colnames]
    length_long = nrow(targetdata_nonimg)
    data_long <- rbind(targetdata_nonimg,targetdata_nonimg) 
    data_long <- data.frame(data_long,hemi = c(rep(1,length_long),rep(2,length_long)))
    
    # add up the tracts
    lh_data = targetdata_img[,grep(paste0('lh.'),colnames(targetdata_img),fixed = T)]
    rh_data = targetdata_img[,grep(paste0('rh.'),colnames(targetdata_img),fixed = T)]
       
    colnames_longformat = gsub(pattern = 'lh\\.',replacement = '',colnames(lh_data))
    colnames_longformat = unlist(colnames_longformat)
    colnames(lh_data) = colnames(rh_data) = colnames_longformat
    
    long_data_clean = rbind(lh_data,rh_data)
    data_long = cbind(data_long[,1],long_data_clean,data_long[,2:ncol(data_long)])
    colnames(data_long)[1]='f.eid'


  return(data_long)
}
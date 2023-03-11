library(dplyr)
library(purrr)

pathogen_pca <- function(data_pathogen,
                         n.pc,
                         filter_threshold = 0.05,
                         id.list = NULL)
{
  data_pathogen_mat <- data_pathogen %>% dplyr::select(-ID)
  
  PC <- prcomp(t(data_pathogen_mat),scale = F)
  ## initiate out
  out <- data.frame(PID = data_pathogen$ID)
  
  ## return ppca
  ppca_out <- as.data.frame(PC$rotation) %>% dplyr::select(paste0('PC',1:n.pc))
  names(ppca_out) <- paste0("pPC", 1:n.pc)
  
  out <- cbind(out, ppca_out)
  return(out)
  
  
}

args <- commandArgs(trailingOnly = TRUE)

tbls <- args[sapply(args,function(x) grepl(x=x,pattern = 'AA_Table'))]
data_pathogen <- lapply(tbls,function(x) data.table::fread(x) %>% dplyr::rename(ID = PID)) %>% purrr::reduce(full_join, by = "ID")

n.pc <- args[!sapply(args,function(x) grepl(x=x,pattern = 'AA_Table'))][[1]]
out_path <- args[!sapply(args,function(x) grepl(x=x,pattern = 'AA_Table'))][[2]]

if(n.pc == 0){
  pPC <- data.frame(PID = data_pathogen$ID)
  
}else{
  pPC <- pathogen_pca(data_pathogen=data_pathogen,n.pc=as.numeric(n.pc))
}

data.table::fwrite(pPC,out_path,sep = ' ',col.names = T,row.names = F)
library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
library(data.table,quietly=TRUE)
setwd('/home/')
args <- commandArgs(trailingOnly = TRUE)
covar_file_path <- args[[1]]
pPCA_file_path <- args[[2]]
PCA_file_path <- args[[3]]
out_path <- args[[4]]

#Read pPCA file
ppca_file = data.table::fread(pPCA_file_path)
#Read provided covar file
covar_file = data.table::fread(covar_file_path)

#Read PCA file, rename headers
pca_file <- data.table::fread(PCA_file_path)
if(nrow(pca_file) == 0){
  pca_file <- data.frame(IID = covar_file$IID)
}else{
  pca_file = pca_file[,-1]
  colnames(pca_file) <- c('IID',paste0('PC',seq(1,ncol(pca_file)-1,1)))
}

#Join all covar files
merged_covar_file = dplyr::left_join(covar_file,pca_file) %>% dplyr::left_join(ppca_file) 
data.table::fwrite(merged_covar_file,
                   file = out_path,col.names = T,row.names = F,sep = ' ',na = 'NA',quote = F)
data.table::fwrite(merged_covar_file %>% dplyr::select(-PID),
                   file = gsub(out_path,pattern = '.txt',replacement = '_plink.txt'),col.names = T,row.names = F,sep = ' ',na = 'NA',quote = F)

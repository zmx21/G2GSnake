library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
library(data.table,quietly=TRUE)

args <- commandArgs(trailingOnly = TRUE)
covar_file_path <- args[[1]]

#Read pPCA file
ppca_file = data.table::fread("../results/tmp/pPCA.txt")
#Read PCA file, rename headers
pca_file = data.table::fread("../results/tmp/G2G_QC.eigenvec")[,-1]
colnames(pca_file) <- c('#IID',paste0('PC',seq(1,ncol(pca_file)-1,1)))
#Read provided covar file
covar_file = data.table::fread(covar_file_path)

#Join all covar files
merged_covar_file = dplyr::left_join(covar_file,pca_file) %>% dplyr::left_join(ppca_file) 
data.table::fwrite(merged_covar_file %>% dplyr::rename(IID = `#IID`),file = "../results/tmp/merged_covar.txt",col.names = T,row.names = F,sep = ' ',na = 'NA',quote = F)

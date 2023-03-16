library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
args <- commandArgs(trailingOnly = TRUE)
setwd('/home/')
covar_file_path <- args[[1]]
aa_tbl_path <- args[[2]]
aa_info_path <- args[[3]]

pathogen_MAC <- as.numeric(args[[4]])
pathogen_missing <- as.numeric(args[[5]])
out_path_tbl <- args[[6]]
out_path_variants <- args[[7]]
out_path_info <- args[[8]]

#Read provided covar file
covar_file = data.table::fread(covar_file_path)

#Read AA Tbl
aa_tbl = data.table::fread(aa_tbl_path) 

#Check if variant pass QC
NA_freq_by_snps <- apply(is.na(aa_tbl[,-'PID']), 2, sum, na.rm = TRUE) / nrow(aa_tbl)
counts_by_snps_0 <- apply(aa_tbl[,-'PID'] == 0, 2, sum, na.rm = TRUE)
counts_by_snps_1 <- apply(aa_tbl[,-'PID'] != 0, 2, sum, na.rm = TRUE)
ind_to_keep <- which(NA_freq_by_snps < pathogen_missing & counts_by_snps_0 >pathogen_MAC & counts_by_snps_1 > pathogen_MAC)
aa_tbl_filt <- cbind(aa_tbl[,'PID'],aa_tbl[,-'PID'][,..ind_to_keep])

#Join AA table
aa_tbl_jned = dplyr::left_join(covar_file,aa_tbl_filt,by=c('PID'='PID')) %>% dplyr::select(-'PID')
data.table::fwrite(aa_tbl_jned,col.names = T,row.names = F,quote = F,na = 'NA',sep = '\t',file = out_path_tbl)

aa_tbl_jned_PLINK = dplyr::left_join(covar_file %>% dplyr::select(IID,PID),aa_tbl_filt,by=c('PID'='PID')) %>% dplyr::select(-'PID') %>% dplyr::relocate(IID)
data.table::fwrite(aa_tbl_jned_PLINK,col.names = T,row.names = F,quote = F,na = 'NA',sep = '\t',file = gsub(out_path_tbl,pattern = '.txt',replacement = '_plink.txt'))

write(colnames(aa_tbl_filt[,-'PID']),out_path_variants)

#Write out info table
info_tbl <- data.table::fread(aa_info_path) %>% dplyr::filter(ID %in% colnames(aa_tbl_filt))
data.table::fwrite(info_tbl,col.names = T,row.names = F,quote = F,na = 'NA',sep = '\t',file = out_path_info)

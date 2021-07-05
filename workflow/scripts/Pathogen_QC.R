library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
args <- commandArgs(trailingOnly = TRUE)
covar_file_path <- args[[1]]
aa_tbl_path <- args[[2]]
pathogen_MAC <- as.numeric(args[[3]])
pathogen_missing <- as.numeric(args[[4]])
cur_variant <- args[[5]]

#Read provided covar file
covar_file = data.table::fread(covar_file_path)

#Join AA table
aa_tbl = data.table::fread(aa_tbl_path,select = c('PID',cur_variant)) 
aa_tbl_jned = dplyr::left_join(covar_file,aa_tbl,by=c('PID'='PID')) %>% dplyr::select(-'PID')

#Check if variant pass QC
NA_freq_by_snps <- apply(is.na(aa_tbl_jned[,..cur_variant]), 2, sum, na.rm = TRUE) / nrow(aa_tbl_jned)
counts_by_snps_0 <- apply(aa_tbl_jned[,..cur_variant] == 0, 2, sum, na.rm = TRUE)
counts_by_snps_1 <- apply(aa_tbl_jned[,..cur_variant] != 0, 2, sum, na.rm = TRUE)

#If variant did not pass QC, write out empy file
if(NA_freq_by_snps > pathogen_missing | counts_by_snps_0 < pathogen_MAC | counts_by_snps_1 < pathogen_MAC){
  system(glue::glue("touch ../results/tmp/pathogen_variants/{cur_variant}"))
}else{
  data.table::fwrite(aa_tbl_jned,
                     glue::glue("../results/tmp/pathogen_variants/{cur_variant}"),sep = ' ',row.names = F,na = 'NA',quote = F)
}
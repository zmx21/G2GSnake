library(glue,quietly=TRUE)
library(filematrix,quietly=TRUE)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[sapply(args,function(x) grepl(x=x,pattern = 'log.out'))]
results_dir <- sapply(results_dir,function(x) gsub(x=x,pattern = 'log.out',replacement = ''))

host_variants_path <- args[!sapply(args,function(x) grepl(x=x,pattern = 'log.out'))][[1]]
out_path <- args[!sapply(args,function(x) grepl(x=x,pattern = 'log.out'))][[2]]


host_variants <- data.table::fread(host_variants_path,header = F)
g2g_file_paths <- lapply(results_dir,function(x) dir(x)[sapply(dir(x),function(x) grepl(pattern = '.allchr.txt',x = x))])
pathogen_variants <- lapply(g2g_file_paths,function(x) sapply(x,function(y) strsplit(x=y,'.allchr.txt')[[1]][1],USE.NAMES = F))

#Initialize file matrix
fm = fm.create(filenamebase = "../results/G2G_Results", nrow = nrow(host_variants),
               ncol = length(unlist(pathogen_variants)), type = "double")
colnames(fm) = unlist(pathogen_variants)
rownames(fm) = host_variants$V1

for(i in 1:length(g2g_file_paths)){
  cur_gene_results <- paste0(results_dir[i],g2g_file_paths[[i]])
  cur_gene_pathogen_variants <- pathogen_variants[[i]]
  for(j in 1:length(cur_gene_results)){
    df <- data.table::fread(cur_gene_results[j])
    if(nrow(df) == 0){
      fm[,which(colnames(fm) == cur_gene_pathogen_variants[j])] <- rep(NA,nrow(fm))
    }else{
      fm[,which(colnames(fm) == cur_gene_pathogen_variants[j])] <- rep(NA,nrow(fm))
      fm[match(df$SNPID,rownames(fm)),which(colnames(fm) == cur_gene_pathogen_variants[j])] <- df$p.value
    }
  }
}
close(fm)

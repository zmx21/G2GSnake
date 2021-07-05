library(glue,quietly=TRUE)
library(filematrix,quietly=TRUE)

host_variants <- data.table::fread('../results/tmp/G2G_QC.variants',header = F)
g2g_file_path <- dir('../results/')[sapply(dir('../results/'),function(x) grepl(pattern = '.allchr.txt',x = x))]
pathogen_variants <- sapply(g2g_file_path,function(x) strsplit(x=x,'.allchr.txt')[[1]][1],USE.NAMES = F)

#Initialize file matrix
fm = fm.create(filenamebase = "../results/G2G_Results", nrow = nrow(host_variants), ncol = length(g2g_file_path), type = "double")
colnames(fm) = pathogen_variants
rownames(fm) = host_variants$V1

for(i in 1:length(g2g_file_path)){
  df <- data.table::fread(paste0('../results/',g2g_file_path[i]))
  if(nrow(df) == 0){
    fm[,i] <- rep(NA,nrow(fm))
  }else{
    fm[match(df$SNPID,rownames(fm)),i] <- df$p.value
  }
}
close(fm)

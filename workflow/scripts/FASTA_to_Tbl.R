library(seqinr)
library(dplyr)
library(parallel)

FastaToTbl <- function(input_fasta_path){
  #Read in AA fasta for protein
  input_fasta <- seqinr::read.fasta(input_fasta_path)
  gene_name <- strsplit(strsplit(input_fasta_path,split = 'gene.')[[1]][2],split = '.fasta')[[1]][1]
  seq_length <- unique(sapply(input_fasta,length))
  if(length(seq_length) > 1){
    stop('error: sequences with different length')
  }
  
  #Get unique aa residues at each position, exclude x
  aa_residues <- lapply(1:seq_length,function(i) unique(sapply(1:length(input_fasta),function(y) input_fasta[[y]][i])))
  aa_residues <- lapply(1:seq_length,function(i) setdiff(toupper(aa_residues[[i]]),c('X','-')))
  variant_names <- sapply(1:seq_length,function(i) paste0('gene_',gene_name,'_pos_',i,'_',ifelse(aa_residues[[i]]=='*','STOP',aa_residues[[i]])))
  
  #Construct Table
  samples <- sapply(names(input_fasta),function(x) strsplit(x=x,split = ';ambiguous')[[1]][1])
  uniq_samples <- unique(samples)
  aa_tbl <- matrix(0,nrow = length(uniq_samples),ncol = length(unlist(variant_names)))
  rownames(aa_tbl) <- uniq_samples
  colnames(aa_tbl) <- unlist(variant_names)
  
  for(i in 1:seq_length){
    cur_pos_residues <- lapply(uniq_samples,function(x) unique(sapply(input_fasta[samples == x],function(y) y[[i]])))
    cur_pos_result <- do.call(rbind,lapply(cur_pos_residues,function(cur_sample_residues){
      if(all(cur_sample_residues == 'x')){
        tmp <- rep(NA,length(variant_names[[i]]))
        names(tmp) <- variant_names[[i]]
      }else if(all(cur_sample_residues == '-')){
        tmp <- rep(0,length(variant_names[[i]]))
        names(tmp) <- variant_names[[i]]
      }else{
        cur_sample_residues <- toupper(setdiff(cur_sample_residues,'x'))
        tmp <- rep(0,length(variant_names[[i]]))
        names(tmp) <- variant_names[[i]]
        tmp[match(cur_sample_residues,aa_residues[[i]])] <- 1
      }
      return(tmp)
    }))
    aa_tbl[,colnames(cur_pos_result)] <- cur_pos_result
  }
  return(as.data.frame(aa_tbl))
}
args <- commandArgs(trailingOnly = T)

#Write out variants
tbl_out = FastaToTbl(paste0(args[[1]]))
tbl_out$PID <- rownames(tbl_out)
tbl_out <- tbl_out %>% dplyr::relocate(PID) %>% dplyr::filter(PID %in% data.table::fread(args[[2]])$PID)

data.table::fwrite(x = tbl_out,
                   col.names = T,row.names = F,quote = F,sep = '\t',file =  args[[3]],na = 'NA')

#Write out variant info
variants <- colnames(tbl_out %>% dplyr::select(-PID))
genes <- sapply(variants,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
pos <- sapply(variants,function(x) strsplit(x=strsplit(x=x,split = 'pos_')[[1]][2],split = '_')[[1]][1])

info_tbl <- data.frame(ID = variants,Gene = genes,Pos = pos)
data.table::fwrite(x = info_tbl,
                   col.names = T,row.names = F,quote = F,sep = '\t',file = args[[4]],na = 'NA')

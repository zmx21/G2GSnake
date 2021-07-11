library(seqinr)
library(R.utils)
library(parallel)

GetAmbigiousCode <- function(nucl){
  result = switch(  
    tolower(nucl),  
    "r"= c('a','g'),  
    "y"= c('c','t'),  
    "s"= c('g','c'),  
    "w"= c('a','t'),  
    "k"= c('g','t'),  
    "m"= c('a','c'),  
    "b"= c('c','g','t'),  
    "d"= c('a','g','t'),  
    "h"= c('a','c','t'),  
    "v"= c('a','c','g'),  
    "n"= 'n',
    "t"= "t", 
    "a"= "a", 
    "c"= "c", 
    "g"= "g",
    "n"= "n"
  )  
  return(result)
}

ExpandFasta <- function(input_fasta_path,n_cores){
  input_fasta <- seqinr::read.fasta(input_fasta_path)
  final_output <- mclapply(1:length(input_fasta),function(i){
    output_fasta <- list()
    counter <- 1
    cur_seq <- input_fasta[[i]]
    if(any(!tolower(cur_seq) %in% c('a','c','t','g','n'))){
      ambigious_pos <- sort(which(!tolower(cur_seq) %in% c('a','c','t','g','n')),decreasing = F)
      
      #Get consecutive positions in potentially the same codons (difference less than two)
      ambigious_pos_diff <- which(diff(ambigious_pos) <= 2)
      if(length(ambigious_pos_diff) > 0){
        ambigious_pos_diff <- as.data.frame(seqToIntervals(ambigious_pos_diff))
        ambigious_pos_diff$to <- ambigious_pos_diff$to + 1

        ambigious_pos_list <- lapply(1:nrow(ambigious_pos_diff),function(k){
          from = ambigious_pos[ambigious_pos_diff$from[k]]
          to = ambigious_pos[ambigious_pos_diff$to[k]]
          if(to - from + 1 <= 3){
            return(list(seq(from,to,1)))
          }else{
            return(lapply(1:(to - from + 1 - 2),function(q) seq(from,to,1)[q:(q+2)]))
          }
        })
        ambigious_pos_list <- c(unlist(ambigious_pos_list,recursive = F),as.list(setdiff(ambigious_pos,unlist(ambigious_pos_list))))
        
      }else{
        ambigious_pos_list <- as.list(ambigious_pos)
      }
      
      for(k in 1:length(ambigious_pos_list)){
        cur_pos <- ambigious_pos_list[[k]]
        if(length(cur_pos) == 1){
          ambigious_nucl <- GetAmbigiousCode(cur_seq[[cur_pos]])
          tmp_seq <- lapply(ambigious_nucl,function(x) {
            tmp <- cur_seq
            tmp[cur_pos] <- x
            return(tmp)
          })
        }else{
          ambigious_nucl <- lapply(cur_pos,function(x) GetAmbigiousCode(cur_seq[[x]]))
          nucl_comb <- expand.grid(ambigious_nucl)
          tmp_seq <- lapply(1:nrow(nucl_comb),function(x) {
            tmp <- cur_seq
            tmp[cur_pos] <- as.vector(t(nucl_comb[x,]))
            return(tmp)
          })
        }
        names(tmp_seq) <- paste0(names(input_fasta)[i],';ambiguous=',seq(counter,counter + length(tmp_seq) - 1,1))
        output_fasta <- c(output_fasta,tmp_seq)
        counter <- counter + length(tmp_seq)
      }
    }else{
      tmp_seq <- list(cur_seq)
      names(tmp_seq) <- paste0(names(input_fasta)[i],';ambiguous=1')
      output_fasta <- c(output_fasta,tmp_seq)
    }
    return(output_fasta)
  },mc.cores = n_cores)
  return(unlist(final_output,recursive = F))
}
args <- commandArgs(trailingOnly = T)
fasta <- ExpandFasta(args[[1]],as.numeric(args[[3]]))
seqinr::write.fasta(fasta,names(fasta),args[[2]])


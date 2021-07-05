library(dplyr)

pathogen_pca <- function(path_tree,
                         path_pathogen,
                         n.pc,
                         filter_threshold = 0.05,
                         id.list = NULL)
{
  
  ## 1) TREE --------------------------------------
  ## ----------------------------------------------
  tree <- ape::read.tree(path_tree)
  #plot(tree)
  
  ## check if matrix is singular
  #solve(ape::vcv.phylo(tree))
  
  ## 2) Y matrix ----------------------------------
  ## ----------------------------------------------
  
  data_pathogen <- data.table::fread(path_pathogen) %>% dplyr::rename(ID = PID)
  
  if (!is.null(id.list))
  {
    data_pathogen <- data_pathogen %>%
      filter(ID %in% id.list)
    
  }
  
  tree_subset <-
    ape::keep.tip(tree, as.character(data_pathogen$ID))
  
  ## reorder rows, same tiplabels
  id_tree <- tibble(ID = tree_subset$tip.label)
  
  ## turn data into matrix and reorder rows according to tips
  Y <- left_join(id_tree,data_pathogen)
  
  stopifnot(identical((tree_subset$tip.label), as.character(Y$ID)))
  
  
  id_pathogen <- Y$ID
  
  ## remove ID
  Y_no_id <- Y %>% dplyr::select(-ID) %>% as.matrix()
  
  ## remove columns with only wild type
  Y_no_id <-
    Y_no_id[, which(apply(Y_no_id, 2, function(x)
      mean(x, na.rm = TRUE)) != 0)]
  
  ## Remove variant between 0.05 and 0.95 allele frequency
  Y_no_id.freq <-
    parallel::mclapply(1:ncol(Y_no_id), function(x) {
      mean(Y_no_id[, x], na.rm = T)
    }, mc.cores = 20) %>% unlist()
  Y_no_id <- Y_no_id[, which(Y_no_id.freq > filter_threshold & Y_no_id.freq < (1-filter_threshold)) ]
  
  ## only non varying ones
  Y_no_id.var <-
    parallel::mclapply(1:ncol(Y_no_id), function(x) {
      var(Y_no_id[, x], na.rm = T)
    }, mc.cores = 20) %>% unlist()
  Y_no_id <- Y_no_id[, which(Y_no_id.var != 0)]
  
  
  
  ##  pPCA --------------------------------------
  ## ----------------------------------------------
  
  ## initiate out
  out <- data.frame(PID = id_pathogen)
  
  #print(sessionInfo())
  
  ## from Nimisha
  vir_4d <- phylobase::phylo4d(tree_subset, Y_no_id)
  
  vir_pca <- adephylo::ppca(
    vir_4d,
    scale = TRUE,
    ## if scaled, very similar to pca
    scannf = FALSE,
    nfposi = 16,
    method = "Abouheif"
  )
  
  ## return ppca
  ppca_out <- vir_pca$li %>% dplyr::select(1:n.pc)
  names(ppca_out) <- paste0("pPC", 1:n.pc)
  
  out <- cbind(out, ppca_out)
  return(out)
  
}

args <- commandArgs(trailingOnly = TRUE)
pPC <- pathogen_pca(path_pathogen=args[[1]],
                    path_tree=args[[2]],n.pc=as.numeric(args[[3]]))
data.table::fwrite(pPC,args[[4]],sep = ' ',col.names = T,row.names = F)
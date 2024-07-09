library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)
library(genomation)
library(sva)

# Get sample type (E or T) for sample annotation for removing batch effects
get_sample_location <- function(sample.ids){
  sample_location <- sample_ids %>% 
    lapply(strsplit, '_') %>% 
    lapply(function(x) x[[1]][[2]]) %>% 
    lapply(substr, 1, 1) %>% 
    unlist
  return(sample_location)
}

# Replace chrIDs of a methylRawKit object from refseq to UCSC
replace_chr_ids <- function(obj){
  mkit_obj <- obj
  # read in alias file
  chrom_alias <- read.csv('data/reference/bosTau9.chromAlias.txt', sep='\t')
  # replace in each sub-object
  for(i in seq(1,length(obj))){
    mkit_obj[[i]]$chr <- chrom_alias$X..ucsc[match(mkit_obj[[i]]$chr, chrom_alias$refseq)]
  }
  return(mkit_obj)
}

# convert CSV to format usable by moGCN (not working as a function yet)
to_mogcn_format <- function(csvfile){
  mkit_df <- read.csv(csvfile, header=T)
  n_samples <- as.integer(str_match(colnames(mkit_df)[length(colnames(mkit_df))-7], '\\d+'))
  
  output_matrix <- matrix(NA, dim(mkit_df)[1], n_samples)
  for(i in 1:n_samples){
    output_matrix[,i] <- mkit_df[,3*i+4] / mkit_df[,3*i+3] * 1000   # if I end up removing the rownames when exporting, just need
  }                                                                 # to change the +4 and +3 here to +3 and +2
  
  output_df <- as.data.frame(output_matrix)
  colnames(output_df) <- c('IVP10','IVP11','IVP14','IVPF10','IVPF15','VE24','VE26','VE28','VE29','VE4','VE5','VE6','VE3')
  colnames(output_df) <- c('IVP10','IVP11','IVP14','IVPF9','VE24','VE26','VE28','VE29','VE4','VE5','VE3')
  
  output_transposed <- t(output_df)
  colnames(output_transposed) <- mkit_df$feature.name
  
  return(output_transposed)
}

# ComBat Batch Correction (using output df from to_mogcn_format function)
if(F){
  batched_df <- ComBat(output_df, batch=c(batch_ordered))
  
  
}


###### code for rearranging tables from RNAseq
if(F){
  ED_rnaseq <- read.csv('data/miscellaneous/D15_ED_AdjData.txt')

  ED_rnaseq_reordered <- ED_rnaseq[,c(1,2,3,4,5,8,10,11,12,13,15,16,17,14)]
  ED_rnaseq_reordered <- ED_rnaseq_reordered[rowSums(ED_rnaseq_reordered[,-1]) != 0, ]
  sapply(strsplit(colnames(ED_rnaseq_reordered),'_'), `[`, 1)
  rownames(ED_rnaseq_reordered) <- ED_rnaseq_reordered[,1]
  ED_rnaseq_reordered <- ED_rnaseq_reordered[,-1]
  
  ED_rnaseq_transposed <- t(ED_rnaseq_reordered)
  colnames(ED_rnaseq_transposed) <- sapply(strsplit(colnames(ED_rnaseq_transposed),':'), `[`, 2)
  
  write.csv(ED_rnaseq_transposed,'data/miscellaneous/rnaseq_ED_mogcn.csv')
  
  # same for TE
  TE_rnaseq <- read.csv('data/miscellaneous/D15_TE_AdjData.txt', sep='\t')
  
  TE_rnaseq_reordered <- TE_rnaseq[,c(1,2,3,4,7,9,10,11,12,14,15,13)]
  TE_rnaseq_reordered <- TE_rnaseq_reordered[rowSums(TE_rnaseq_reordered[,-1]) != 0, ]
  colnames(TE_rnaseq_reordered) <- sapply(strsplit(colnames(TE_rnaseq_reordered),'_'), `[`, 1)
  rownames(TE_rnaseq_reordered) <- TE_rnaseq_reordered[,1]
  TE_rnaseq_reordered <- TE_rnaseq_reordered[,-1]
  
  TE_rnaseq_transposed <- t(TE_rnaseq_reordered)
  colnames(TE_rnaseq_transposed) <- sapply(strsplit(colnames(TE_rnaseq_transposed),':'), `[`, 2)
  
  write.csv(TE_rnaseq_transposed,'data/miscellaneous/rnaseq_TE_mogcn.csv')
}

# comparing variances
if(F){
  IV_pm <- pm[,1:8]
  VE_pm <- pm[,9:16]
  
  IV_var <- apply(IV_pm, 1, var, na.rm=T)
  VE_var <- apply(VE_pm, 1, var, na.rm=T)
  
  t.test(IV_var, VE_var, paired=F)
}

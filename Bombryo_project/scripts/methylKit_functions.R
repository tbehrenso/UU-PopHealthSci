library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)
library(genomation)

test_function1 <- function(...){
  return(sum(...))
}

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



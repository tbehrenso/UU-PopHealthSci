setwd('C:/Users/tbehr/Desktop/UU/WD')

library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)
library(genomation)

file_list <- as.list(list.files('./Outputs_Bismark/'))

# extract sample IDs from file names, and make sure its a list of character vectors
sample_ids <- file_list %>% 
  lapply(str_match, '^[^\\W_]+_[^\\W_]+') %>% 
  lapply(as.character)

# Define which samples to exclude
ED_samples <- c("IVP10_ED", "IVP11_ED", "IVP14_ED", "IVP9_EC", "VE24_ED", "VE26_ED", "VE28_ED", "VE29_ED")
T1D_samples <- c("IVP10_T1D", "IVP11_T1D", "IVP14_T1D", "IVP9_T1C", "VE24_T1D", "VE26_T1D", "VE28_T1D", "VE29_T1D")
exclude_list <- c(ED_samples)

# Find indeces in list to remove
exclude_index <- which(sample_ids %in% exclude_list)
# Remove samples from file_list and sample_ids
if(length(exclude_list) > 0){
  file_list <- file_list[-exclude_index]
  sample_ids <- sample_ids[-exclude_index]
}

# Get sample type (E or T) for sample annotation for removing batch effects
sample_location <- sample_ids %>% 
  lapply(strsplit, '_') %>% 
  lapply(function(x) x[[1]][[2]]) %>% 
  lapply(substr, 1, 1) %>% 
  unlist

file_paths <- lapply(file_list, function(x) paste0('./Outputs_Bismark/', x))

# create binary vector for treatment argument (in this case, in vivo (VE) is 0, in vitro cryo (IVP) is 1)
treatment_binary <- as.numeric(str_detect(unlist(sample_ids), 'IVP'))

mkit_obj <- methRead(file_paths, sample.id=sample_ids, assembly='UCD1.3', treatment=treatment_binary, pipeline='bismarkCoverage')

### Convert accession_value to chrN to match with UCSC table for annotation
chrom_alias <- read.csv("ncbi_dataset/bosTau9.chromAlias.txt", sep='\t')
# replace in each file
for(i in seq(1,length(file_list))){
  mkit_obj[[i]]$chr <- chrom_alias$X..ucsc[match(mkit_obj[[i]]$chr, chrom_alias$refseq)]
}

# Coverage and Methylation Statistics Plots
if(F){
  # Plot all MethylationStats Histograms with same axis (need to manually specify)
  for(i in seq(1,length(file_list))){
    hist(getData(mkit_obj[[i]])$numCs / (getData(mkit_obj[[i]])$numTs + getData(mkit_obj[[i]])$numCs),
         col = '#6495ed', main = sample_ids[[i]], ylim=c(0,4000000), xlab='% methylation per base')
  }
  
  # Plot all CoverageSats Histograms
  for(i in seq(1,length(file_list))){
    hist(log10(getData(mkit_obj[[i]])$coverage),
          col = 'darkgreen', main = sample_ids[[i]], xlab='log10 of read coverage per base')
  }
  
  getMethylationStats(mkit_obj[[1]], both.strands=F, plot=T)
  getMethylationStats(mkit_obj[[1]], both.strands=F, plot=F)
  
  getCoverageStats(mkit_obj[[1]], both.strands=F, plot=T)
  getCoverageStats(mkit_obj[[1]], both.strands=F, plot=F)
}

# filter by coverage
## NOTE: never seems to be below 10. Has this filtering been done before?
mkit_obj_filt <- filterByCoverage(mkit_obj, lo.count=10, hi.perc=99)

# Normalization
## (if multiple samples)
mkit_obj_norm <- normalizeCoverage(mkit_obj_filt, method = "median")  # can also use 'mean' method. Should investigate differences

# Merge Data -> extracts sites common in all samples
mkit_merged <- unite(mkit_obj_norm)

## Further Filtering (remove sites with little variation)
# get percent methylation matrix
pm=percMethylation(mkit_merged)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation to determine a suitable cutoff
hist(sds, breaks = 100)
# keep only CpG with standard deviations larger than 2%
mkit_merged <- mkit_merged[sds > 2]


# Correlation
getCorrelation(mkit_merged,plot=F)
# Cluster
clusterSamples(mkit_merged, dist="correlation", method="ward.D", plot=TRUE)
# PCA
PCASamples(mkit_merged, screeplot = F)

# Batch Effects (if using samples combined between ED and T1D)
if(F){
  # Batch Effects
  sample_annotation <- data.frame(batch.id=sample_location)
  component_association <- assocComp(mkit_merged, sample_annotation)
  
  # remove components with association p-value <0.01
  associated_componenets <- which(component_association$association[1,] < 0.01)
  cat('Removing Principal Components:', associated_componenets)
  mkit_batched <- removeComp(mkit_merged, comp=associated_componenets)
}

#### Differential Methylation
# create data frame with sample location as factor for covariate
covariate_df <- data.frame(covariate = as.factor(sample_location))
diff_methylation <- calculateDiffMeth(mkit_merged, overdispersion = 'MN')

diff_methylation_25p <- getMethylDiff(diff_methylation,difference=25,qvalue=0.05)


# Gene Annotation

refseq_features <- readTranscriptFeatures("ncbi_dataset/UCD1.2.gz")

diff_methylation_annotated <- annotateWithGeneParts(as(diff_methylation_25p,"GRanges"),refseq_features)











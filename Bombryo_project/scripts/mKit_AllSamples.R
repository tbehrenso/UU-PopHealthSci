library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)
library(genomation)
library(biomaRt)
library(ggrepel)

source('scripts/methylKit_functions.R')

SAMPLE_LOC <- 'ED'         # should be either 'ED' or 'T1D'
SAMPLE_EXCLUDE <- NA        # NA or name of ONE sample

MINIMUM_COVERAGE <- 5
COVERAGE_HI_PERC <- 99.9
METH_DIFF_PERC <- 10
METH_DIFF_Q <- 0.1
MIN.PER.GROUP <- 7    # integer or NA
TILED <- T
OVERDISPERSION_CORRECTION <- F
MORPHOLOGY_BATCH_CORRECTION <- F
BATCH_PVAL <- 0.01
ONLY_MATCHING <- F

# -------------------------------------------------------
#   Define and Extract Filepaths and Samples, 
# -------------------------------------------------------
file_list <- as.list(list.files('data/raw/Outputs_Coverage/'))

# extract sample IDs from file names, and make sure its a list of character vectors
sample_ids <- file_list %>% 
  lapply(str_match, '^[^\\W_]+_[^\\W_]+') %>% 
  lapply(as.character)

# Define which samples to exclude
ED_samples <- str_match(unlist(sample_ids), '.*E[D|C]')[!is.na(str_match(unlist(sample_ids), '.*E[D|C]'))]
T1D_samples <- str_match(unlist(sample_ids), '.*T1[D|C]')[!is.na(str_match(unlist(sample_ids), '.*T1[D|C]'))]
if(SAMPLE_LOC=='ED'){
  exclude_list <- c(T1D_samples)
} else if(SAMPLE_LOC=='T1D'){
  exclude_list <- c(ED_samples)
}
if(!is.na(SAMPLE_EXCLUDE)){
  exclude_list <- c(exclude_list, SAMPLE_EXCLUDE)
}

# Find indeces in list to remove
exclude_index <- which(sample_ids %in% exclude_list)
# Remove samples from file_list and sample_ids
if(length(exclude_list) > 0){
  file_list <- file_list[-exclude_index]
  sample_ids <- sample_ids[-exclude_index]
}

file_paths <- lapply(file_list, function(x) paste0('data/raw/Outputs_Coverage/', x))

## Pre-defined sample list (that coincide with matching RNAseq samples)
if(ONLY_MATCHING){
  ED_samples <- c('IVP10_ED', 'IVP11_ED', 'IVP14_ED', 'IVPF10_ED', 'IVPF15_ED',
                  'VE24_ED', 'VE26_ED', 'VE28_ED', 'VE29_ED', 'VE4_ED', 'VE5_ED', 'VE6_ED', 'VE3_EC')
  T1D_samples <- c('IVP10_T1D', 'IVP11_T1D', 'IVP14_T1D', 'IVPF9_T1C',
                   'VE24_T1D', 'VE26_T1D', 'VE28_T1D', 'VE29_T1D', 'VE4_T1D', 'VE5_T1D', 'VE3_T1C')
  if(SAMPLE_LOC=='ED'){
    sample_ids <- ED_samples
    file_paths <- lapply(ED_samples, function(x) paste0('data/raw/Outputs_Coverage/', x, '_1_bismark_bt2_pe.bismark.cov.gz'))
  } else if (SAMPLE_LOC=='T1D'){
    sample_ids <- T1D_samples
    file_paths <- lapply(T1D_samples, function(x) paste0('data/raw/Outputs_Coverage/', x, '_1_bismark_bt2_pe.bismark.cov.gz'))
  }
}

# Get sample information (morphology + sex) from file
sample_info <- read.table('data/reference/SampleInfo.csv', sep=',', header = T)
batch_ordered <- sapply(sample_ids, function(x) {sample_info$Size[which(sample_info$H.Folders==x)]})
sex_ordered <- sapply(sample_ids, function(x) {sample_info$Sex[which(sample_info$H.Folders==x)]})

# create binary vector for treatment argument (in this case, in vivo (VE) is 0, in vitro cryo (IVP) is 1)
treatment_binary <- as.numeric(str_detect(unlist(sample_ids), 'IVP'))

# -------------------------------------------------------
#   Load and Process methylRawList object
# -------------------------------------------------------
# Code for creating new object instead of loading
if(F){
  mkit_obj <- methRead(file_paths, header=F, sample.id=as.list(sample_ids), mincov = MINIMUM_COVERAGE,
                       assembly='UCD1.3', treatment=treatment_binary, pipeline='bismarkCoverage')
  mkit_obj <- replace_chr_ids(mkit_obj)
}

# load a saved mkit object
if(TILED){
  if(ONLY_MATCHING){
    switch(SAMPLE_LOC,
           T1D = load('data/processed/AA_Data/mkit_tiled1000_T1D_matching_mincov10.RData'), 
           ED = load('data/processed/AA_Data/mkit_tiled1000_ED_matching_mincov10.RData'))  # have a mincov5 variant of this one
  } else {
    switch(SAMPLE_LOC,
           T1D = load('data/processed/AA_Data/mkit_tiled1000_T1D_All_mincov10.RData'),
           ED = load('data/processed/AA_Data/mkit_tiled1000_ED_All_mincov10.RData'))
  }
  mkit_obj <- mkit_tiled
} else {
  if(ONLY_MATCHING){
    switch(SAMPLE_LOC,
           T1D = load('data/processed/AA_Data/mkit_obj_T1D_matching_mincov10.RData'), 
           ED = load('data/processed/AA_Data/mkit_obj_ED_matching_mincov10.RData'))
  } else {
    switch(SAMPLE_LOC,
           T1D = load('data/processed/AA_Data/mkit_obj_T1D_raw_mincov10.RData'),
           ED = load('data/processed/AA_Data/mkit_obj_ED_raw_mincov10.RData'))
  }
}

mkit_obj_norm <- filterByCoverage(mkit_obj, lo.count=MINIMUM_COVERAGE, hi.perc=COVERAGE_HI_PERC)

mkit_obj_norm <- normalizeCoverage(mkit_obj_norm, method = 'median')

# -------------------------------------------------------
#   Create and examine merged methylBase object
# -------------------------------------------------------

if(is.na(MIN.PER.GROUP)){
  mkit_merged <- methylKit::unite(mkit_obj_norm)
}else{
  mkit_merged <- methylKit::unite(mkit_obj_norm, min.per.group=as.integer(MIN.PER.GROUP))
}

## Further Filtering (remove sites with little variation)
# get percent methylation matrix
pm=percMethylation(mkit_merged)
# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)
# Visualize the distribution of the per-CpG standard deviation to determine a suitable cutoff
hist(sds, breaks = 100)
# keep only CpG with standard deviations larger than 2%. Keep NA for when using min.per.group argument
mkit_merged <- mkit_merged[sds > 2 | is.na(sds)]

# Correlation
getCorrelation(mkit_merged,plot=F)
# Cluster
clusterSamples(mkit_merged, dist='correlation', method='ward.D', plot=TRUE)
# PCA
PCASamples(mkit_merged, screeplot = F, transpose = T, adj.lim = c(0.3,0.1))

# can only do this part if not using the min.per.group argument
if(is.na(MIN.PER.GROUP)){
  
  # PCA Custom
  pm=percMethylation(mkit_merged)
  pca_stats <- prcomp(pm)
  pca_df <- data.frame(PC1=pca_stats$rotation[,1], PC2=pca_stats$rotation[,2], type=as.factor(treatment_binary),
                       morph=batch_ordered, sex=sex_ordered)
  ggplot(pca_df, aes(x=PC1, y=PC2, col=type, label=sapply(strsplit(row.names(pca_df),'_'), `[`, 1), shape=morph)) + 
    geom_point(size=3) + geom_text_repel()
  ggplot(pca_df, aes(x=PC1, y=PC2, col=morph, label=sapply(strsplit(row.names(pca_df),'_'), `[`, 1), shape=type)) + 
    geom_point(size=3) + geom_text_repel() +
    scale_shape_manual(values=c(1, 16), labels=c('VE', 'IVP')) +
    scale_color_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
    #stat_ellipse()
  
  ### Batch Effects
  sample_annotation <- data.frame(batch.id=batch_ordered)
  component_association <- assocComp(mkit_merged, sample_annotation)


  # remove components with association p-value <0.01
  associated_componenets <- which(component_association$association[1,] < BATCH_PVAL)
  cat('Removing Principal Components:', associated_componenets)
  mkit_batched <- removeComp(mkit_merged, comp=associated_componenets)
} else {
  # else if min.per.group is active, just continue with mkit_merged
  mkit_batched <- mkit_merged
}

# Remove those mapped to chrUn before differential methylation calculation
mkit_batched <- mkit_batched[substring(mkit_batched$chr, 1, 5) != 'chrUn']

# -------------------------------------------------------
#   Differential Methylation and Annotation
# -------------------------------------------------------

# create data frame with sample location as factor for covariate
covariate_df_both <- data.frame(morph = as.factor(batch_ordered), sex = as.factor(sex_ordered))
covariate_df_sex <- data.frame(sex = as.factor(sex_ordered))


diff_methylation <- calculateDiffMeth(mkit_batched, overdispersion = ifelse(OVERDISPERSION_CORRECTION, 'MN', 'none'),
                                      covariates = covariate_df_both)

diff_methylation_25p <- getMethylDiff(diff_methylation,difference=METH_DIFF_PERC,qvalue=METH_DIFF_Q)

hist(diff_methylation_25p$meth.diff)

#### Gene Annotation
## Examine percentage of differentially methylated sites are in introns/exons/promoters/intergenic
if(!exists('refseq_features')){refseq_features <- readTranscriptFeatures('data/reference/UCD1.2_Genes.gz')}

diff_meth_annotated <- annotateWithGeneParts(as(diff_methylation_25p,'GRanges'),refseq_features)

# View the distance to the nearest Transcription Start Site
# the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(diff_meth_annotated)

# get location of each CpG 
gene_members <- getMembers(diff_meth_annotated)

plotTargetAnnotation(diff_meth_annotated, main = 'Differential Methylation Annotation', cex.legend = 0.7)

# Summarize methylation over specific regions (using original methylKit obj)
promoters <- regionCounts(mkit_obj, refseq_features$promoters)

## Examine percentage in CpG Islands
if(!exists('refseq_cpgislands')){refseq_cpgislands <- readFeatureFlank('data/reference/UCD1.2_CpG.gz')}

diff_meth_cpg <- annotateWithFeatureFlank(as(diff_methylation_25p,'GRanges'), refseq_cpgislands$features, refseq_cpgislands$flanks,
                                          feature.name='CpGi',flank.name='shores')

island_members <- getMembers(diff_meth_cpg)

plotTargetAnnotation(diff_meth_cpg, col=c('green','gray','white'), main='differential methylation annotation')

# get percentage of intron/exon/promoters that overlap with differentially methylated bases
getFeatsWithTargetsStats(diff_meth_annotated,percentage=F)
getFeatsWithTargetsStats(diff_meth_cpg,percentage=F)


# extract relevant info into one dataframe (only for sites on known chromosomes)
known_chromosome <- substring(diff_methylation_25p$chr, 1, 5) != 'chrUn'
dmc_transcript_assoc <- apply(gene_members, 1, function(x) {ifelse(sum(x)==0, 'intergenic', colnames(gene_members)[which(x==1)])})
dmc_island_assoc <- apply(island_members, 1, function(x) {ifelse(sum(x)==0, 'other', colnames(island_members)[which(x==1)])})
all_known_features <- cbind(methylKit::getData(diff_methylation_25p)[known_chromosome,], dist_tss, 
                            dmc_transcript_assoc[known_chromosome], dmc_island_assoc[known_chromosome])
colnames(all_known_features)[12:13] <- c('transcript.assoc', 'island.assoc')

# output file with DMCs, adjacent feature, and details
if(F){
  output_filepath <- paste0('output/AllAnalysis/dmclist_',
                            SAMPLE_LOC,
                            ifelse(is.na(SAMPLE_EXCLUDE),'', paste0('_m',strsplit(SAMPLE_EXCLUDE, '_')[[1]])),
                            ifelse(ONLY_MATCHING, '_matching', '_All'),
                            ifelse(TILED,'_tiled', ''),
                            ifelse(is.na(MIN.PER.GROUP),'', paste0('_mpg',MIN.PER.GROUP)),
                            ifelse(OVERDISPERSION_CORRECTION, '', '_noOC'),
                            ifelse(MORPHOLOGY_BATCH_CORRECTION, paste0('_morphbatch',BATCH_PVAL), ''),
                            '_mincov', MINIMUM_COVERAGE,
                            '_dif', METH_DIFF_PERC,
                            '_q', METH_DIFF_Q,
                            '.csv')
  if(!file.exists(output_filepath)){
    write.csv(all_known_features, output_filepath, row.names=F, quote=F)
  }else{
    print('Warning: Output file with that name already exists')
  }
}







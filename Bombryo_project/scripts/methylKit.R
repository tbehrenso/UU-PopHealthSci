library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)
library(genomation)
library(biomaRt)

source('scripts/methylKit_functions.R')

MINIMUM_COVERAGE <- 5
COVERAGE_HI_PERC <- 99.9
METH_DIFF_PERC <- 10
METH_DIFF_Q <- 0.1

file_list <- as.list(list.files('data/raw/Outputs_Coverage/'))

# extract sample IDs from file names, and make sure its a list of character vectors
sample_ids <- file_list %>% 
  lapply(str_match, '^[^\\W_]+_[^\\W_]+') %>% 
  lapply(as.character)

# Define which samples to exclude
ED_samples <- c('IVP10_ED', 'IVP11_ED', 'IVP14_ED', 'IVP9_EC', 'VE24_ED', 'VE26_ED', 'VE28_ED', 'VE29_ED')
T1D_samples <- c('IVP10_T1D', 'IVP11_T1D', 'IVP14_T1D', 'IVP9_T1C', 'VE24_T1D', 'VE26_T1D', 'VE28_T1D', 'VE29_T1D')
exclude_list <- c(T1D_samples, 'IVP9_EC')

# Find indeces in list to remove
exclude_index <- which(sample_ids %in% exclude_list)
# Remove samples from file_list and sample_ids
if(length(exclude_list) > 0){
  file_list <- file_list[-exclude_index]
  sample_ids <- sample_ids[-exclude_index]
}


file_paths <- lapply(file_list, function(x) paste0('data/raw/Outputs_Coverage/', x))

# create binary vector for treatment argument (in this case, in vivo (VE) is 0, in vitro cryo (IVP) is 1)
treatment_binary <- as.numeric(str_detect(unlist(sample_ids), 'IVP'))

# create methylKit object. (Note: by default, removes reads with coverage < 10)
mkit_obj <- methRead(file_paths, header=F, sample.id=sample_ids, mincov = MINIMUM_COVERAGE,
                     assembly='UCD1.3', treatment=treatment_binary, pipeline='bismarkCoverage')

# Create tiled objected to evaluate differential methylation over regions rather than bases (SLOW)
if(F){
  mkit_tiled <- tileMethylCounts(mkit_obj,win.size=1000,step.size=1000)
  
  # load tiled object if saved
  load('data/processed/mkit_tiled_ED.RData')
  # replace mkit_obj with tiled
  mkit_obj <- mkit_tiled
}

# create methylRawObj grouped based on regions of interest (promoter, intron, etc. / CpG Islands, shores)
if(F){
  ## for introns/exons/promoters
  refseq_features <- readTranscriptFeatures('data/reference/UCD1.2_Genes.gz')
  mkit_features <- regionCounts(mkit_obj, refseq_features$promoters)
  
  mkit_obj <- mkit_features
  # for CpG Islands
  refseq_cpgislands <- readFeatureFlank('data/reference/UCD1.2_CpG.gz')
  mkit_cpg <- regionCounts(mkit_obj, refseq_cpgislands$features)
  
  mkit_obj <- mkit_cpg
}


#Convert accession_value to chrN to match with UCSC table for annotation
mkit_obj <- replace_chr_ids(mkit_obj)

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
## NOTE: never seems to be below 10. Low count filtering redundent, because already done when creating mkit_obj
mkit_obj_filt <- filterByCoverage(mkit_obj, lo.count=MINIMUM_COVERAGE, hi.perc=COVERAGE_HI_PERC)

# Normalization
## (if multiple samples)
mkit_obj_norm <- normalizeCoverage(mkit_obj_filt, method = 'median')  # can also use 'mean' method. Should investigate differences

# Merge Data -> extracts sites common in all samples
mkit_merged <- unite(mkit_obj_norm)
#mkit_merged <- unite(mkit_obj_norm, min.per.group=as.integer(2))

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
clusterSamples(mkit_merged, dist='correlation', method='ward.D', plot=TRUE)
# PCA
PCASamples(mkit_merged, screeplot = F, transpose = T)

# Batch Effects
if(F){
  # Batch Effects
  sample_annotation <- data.frame(batch.id=c('a','a','B','a','a','a','a','a'))
  component_association <- assocComp(mkit_merged, sample_annotation)
  
  # remove components with association p-value <0.01
  associated_componenets <- which(component_association$association[1,] < 0.01)
  cat('Removing Principal Components:', associated_componenets)
  mkit_batched <- removeComp(mkit_merged, comp=associated_componenets)
}

#### Differential Methylation
# create data frame with sample location as factor for covariate

#covariate_df <- data.frame(covariate = as.factor(sample_location))
diff_methylation <- calculateDiffMeth(mkit_merged, overdispersion = 'MN')

diff_methylation_25p <- getMethylDiff(diff_methylation,difference=METH_DIFF_PERC,qvalue=METH_DIFF_Q)

hist(diff_methylation_25p$meth.diff)

#### Gene Annotation
## Examine percentage of differentially methylated sites are in introns/exons/promoters/intergenic
refseq_features <- readTranscriptFeatures('data/reference/UCD1.2_Genes.gz')

diff_meth_annotated <- annotateWithGeneParts(as(diff_methylation_25p,'GRanges'),refseq_features)

# View the distance to the nearest Transcription Start Site
# the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(diff_meth_annotated)

# get location of each CpG 
gene_members <- getMembers(diff_meth_annotated)

plotTargetAnnotation(diff_meth_annotated, main = 'Differential Methylation Annotation')

# Summarize methylation over specific regions (using original methylKit obj)
promoters <- regionCounts(mkit_obj, refseq_features$promoters)

## Examine percentage in CpG Islands
refseq_cpgislands <- readFeatureFlank('data/reference/UCD1.2_CpG.gz')

diff_meth_cpg <- annotateWithFeatureFlank(as(diff_methylation_25p,'GRanges'), refseq_cpgislands$features, refseq_cpgislands$flanks,
                                    feature.name='CpGi',flank.name='shores')

island_members <- getMembers(diff_meth_cpg)

plotTargetAnnotation(diff_meth_cpg, col=c('green','gray','white'), main='differential methylation annotation')

# get percentage of intron/exon/promoters that overlap with differentially methylated bases
getFeatsWithTargetsStats(diff_meth_cpg,percentage=F)

# extract relevant info into one dataframe (only for sites on known chromosomes)
known_chromosome <- substring(diff_methylation_25p$chr, 1, 5) != 'chrUn'
dmc_transcript_assoc <- apply(gene_members, 1, function(x) {ifelse(sum(x)==0, 'intergenic', colnames(gene_members)[which(x==1)])})
dmc_island_assoc <- apply(island_members, 1, function(x) {ifelse(sum(x)==0, 'other', colnames(island_members)[which(x==1)])})
all_known_features <- cbind(getData(diff_methylation_25p)[known_chromosome,], dist_tss, 
                      dmc_transcript_assoc[known_chromosome], dmc_island_assoc[known_chromosome])
colnames(all_known_features)[12:13] <- c('transcript.assoc', 'island.assoc')

### biomaRt to convert RefSeq accession to Gene Symbols
mart <- useDataset('btaurus_gene_ensembl', mart = useMart('ensembl'))
refseq_test <- list('XR_236776', 'NM_001083778', 'NM_001001156')
# remove version number from feature names
refseq_feature_names <- sapply(strsplit(all_known_features$feature.name, '\\.'), `[`, 1)

refseq_feature_names_NM <- refseq_feature_names[grepl('^NM', refseq_feature_names)]
refseq_feature_names_NR <- refseq_feature_names[grepl('^NR', refseq_feature_names)]
refseq_feature_names_NP <- refseq_feature_names[grepl('^NP', refseq_feature_names)]
refseq_feature_names_XM <- refseq_feature_names[grepl('^XM', refseq_feature_names)]
refseq_feature_names_XR <- refseq_feature_names[grepl('^XR', refseq_feature_names)]
refseq_feature_names_XP <- refseq_feature_names[grepl('^XP', refseq_feature_names)]

geneid_mrna <- getBM(filters=c('refseq_mrna'), attributes=c('refseq_mrna', 'ensembl_gene_id'),
                     values=refseq_feature_names_NM, mart=mart)
geneid_mrna_pred <- getBM(filters=c('refseq_mrna_predicted'), attributes=c('refseq_mrna_predicted', 'ensembl_gene_id'),
                          values=refseq_feature_names_NM, mart=mart)
geneid_ncrna <- getBM(filters=c('refseq_ncrna'), attributes=c('refseq_ncrna', 'ensembl_gene_id'),
                      values=refseq_feature_names_NM, mart=mart)
geneid_ncrna_pred <- getBM(filters=c('refseq_ncrna_predicted'),
                           attributes=c('refseq_ncrna_predicted', 'ensembl_gene_id'),
                           values=refseq_feature_names_NM, mart=mart)
geneid_peptide <- getBM(filters=c('refseq_peptide'), attributes=c('refseq_peptide', 'ensembl_gene_id'),
                        values=refseq_feature_names_NM, mart=mart)
geneid_peptide_pred <- getBM(filters=c('refseq_peptide_predicted'),
                             attributes=c('refseq_peptide_predicted', 'ensembl_gene_id'),
                             values=refseq_feature_names_NM, mart=mart)

all_gene_ids <- rbind(geneid_mrna, geneid_mrna_pred, geneid_ncrna, geneid_ncrna_pred, geneid_peptide, geneid_peptide_pred)

all_gene_ids_ordered <- sapply(refseq_feature_names, function(x) ifelse(x %in% all_gene_ids$refseq_mrna,
                                                all_gene_ids$ensembl_gene_id[which(all_gene_ids$refseq_mrna==x)], 'NA'))

all_known_features_geneid <- cbind(all_known_features, all_gene_ids_ordered)

# output file with DMCs, adjacent feature, and details
if(F){
  output_filepath <- paste0('output/dmclist_',
                            'ED_tiled1k_mIVP9',           # <- REMEMBER to change this value, not set automatically
                            '_mincov', MINIMUM_COVERAGE,
                            '_dif', METH_DIFF_PERC,
                            '_q', METH_DIFF_Q,
                            '.csv')
  write.csv(all_known_features_geneid, output_filepath, row.names=F, quote=F)
}




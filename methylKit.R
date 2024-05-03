setwd('C:/Users/tbehr/Desktop/UU/WD')

library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)

file_list <- as.list(list.files('./Outputs_Bismark/'))

# extract sample IDs from file names, and make sure its a list of character vectors
sample_ids <- file_list %>% 
  lapply(str_match, '^[^\\W_]+_[^\\W_]+') %>% 
  lapply(as.character)

file_paths <- lapply(file_list, function(x) paste0('./Outputs_Bismark/', x))

# create binary vector for treatment argument (in this case, in vivo (VE) is 0, in vitro cryo (IVP) is 1)
treatment_binary <- as.numeric(str_detect(unlist(sample_ids), 'IVP'))

mkit_obj <- methRead(file_paths, sample.id=sample_ids, assembly='UCD1.3', treatment=treatment_binary, pipeline='bismarkCoverage')

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

# filter by coverage
## NOTE: never seems to be below 10. Has this filtering been done before?
mkit_obj_filt <- filterByCoverage(mkit_obj, lo.count=10, hi.perc=99.9)

# Normalization
## (if multiple samples)
mkit_obj_norm <- normalizeCoverage(mkit_obj_filt, method = "median")  # can also use 'mean' method. Should investigate differences

# Merge Data
IVPF_merged <- unite(mkit_obj_norm)

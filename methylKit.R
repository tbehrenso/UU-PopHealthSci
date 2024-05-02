setwd('C:/Users/tbehr/Desktop/UU/WD')

library('methylKit')
library('ggplot2')
library(stringr)

file_list <- as.list(list.files('./Outputs_Bismark/'))

# extract sample IDs from file names, and make sure its a list of character vectors
sample_ids <- file_list %>% 
  lapply(str_match, '^[^\\W_]+_[^\\W_]+') %>% 
  lapply(as.character)

file_paths <- lapply(file_list, function(x) paste0('./Outputs_Bismark/', x))

# create binary vector for treatment argument (in this case, in vivo (VE) is 0, in vitro cryo (IVP) is 1)
treatment_binary <- as.numeric(str_detect(unlist(sample_ids), 'IVP'))

IVPF_obj <- methRead(file_paths, sample.id=sample_ids, assembly='UCD1.3', treatment=treatment_binary, pipeline='bismarkCoverage')

getMethylationStats(IVPF_obj[[1]], both.strands=F, plot=T)
getMethylationStats(IVPF_obj[[1]], both.strands=F, plot=F)

getCoverageStats(IVPF_obj[[1]], both.strands=F, plot=T)
getCoverageStats(IVPF_obj[[1]], both.strands=F, plot=F)

# filter by coverage
## NOTE: never seems to be below 10. Has this filtering been done before?
IVPF_obj_filt <- filterByCoverage(IVPF_obj, lo.count=10)

# Normalization
## (if multiple samples)
IVPF_obj_norm <- normalizeCoverage(IVPF_obj_filt, method = "median")  # can also use 'mean' method. Should investigate differences

# Merge Data
IVPF_merged <- unite(IVPF_obj_norm)

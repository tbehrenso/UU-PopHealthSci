setwd('C:/Users/tbehr/Desktop/UU/WD')

library('methylKit')
library('ggplot2')

file_list <- list(
  './IVPF15_MethylExt_Paired/IVPF15_ED_1_bismark_bt2_pe.bismark.cov.gz',
  './IVPF11_MethylExt_Paired/IVPF11_T1D_1_bismark_bt2_pe.bismark.cov.gz'
)

IVPF_obj <- methRead(file_list, sample.id=list('IVPF15','IVPF11'), assembly='UCD1.3', treatment=c(0,0), pipeline='bismarkCoverage')

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

setwd('C:/Users/tbehr/Desktop/UU/WD')

library('methylKit')
library('ggplot2')

IVPF15_obj <- methRead('./IVPF15_MethylExt_Ign10/bedGraph/IVPF15_bismark.cov.gz', sample.id='IVPF15', assembly='UCD1.3', pipeline='bismarkCoverage')

getMethylationStats(IVPF15_obj, both.strands=F, plot=T)
getMethylationStats(IVPF15_obj, both.strands=F, plot=F)


getCoverageStats(IVPF15_obj, both.strands=F, plot=T)
getCoverageStats(IVPF15_obj, both.strands=F, plot=F)

# filter by coverage
## NOTE: never seems to be below 10. Has this filtering been done before?
IVPF_obj_filt <- filterByCoverage(IVPF15_obj, lo.count=10)

# Normalization
## (if multiple samples)


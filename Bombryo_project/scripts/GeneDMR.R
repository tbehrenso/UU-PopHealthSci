library(GeneDMRs)
library(ffbase)

ED_controls <- c('C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE24_ED_1_bismark_bt2_pe.bismark.cov.gz',
              'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE26_ED_1_bismark_bt2_pe.bismark.cov.gz',
              'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE28_ED_1_bismark_bt2_pe.bismark.cov.gz',
              'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE29_ED_1_bismark_bt2_pe.bismark.cov.gz')
ED_cases <- c('C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP10_ED_1_bismark_bt2_pe.bismark.cov.gz',
           'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP11_ED_1_bismark_bt2_pe.bismark.cov.gz',
           'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP14_ED_1_bismark_bt2_pe.bismark.cov.gz',
           'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP9_EC_1_bismark_bt2_pe.bismark.cov.gz')
T1D_controls <- c('C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE24_T1D_1_bismark_bt2_pe.bismark.cov.gz',
                 'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE26_T1D_1_bismark_bt2_pe.bismark.cov.gz',
                 'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE28_T1D_1_bismark_bt2_pe.bismark.cov.gz',
                 'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/VE29_T1D_1_bismark_bt2_pe.bismark.cov.gz')
T1D_cases <- c('C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP10_T1D_1_bismark_bt2_pe.bismark.cov.gz',
              'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP11_T1D_1_bismark_bt2_pe.bismark.cov.gz',
              'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP14_T1D_1_bismark_bt2_pe.bismark.cov.gz',
              'C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs/IVP9_T1C_1_bismark_bt2_pe.bismark.cov.gz')

# Read in relevant files
inputmethfile_ED <- Methfile_read(control_paths = ED_controls, case_paths = ED_cases, WGBS=T)
inputmethfile_T1D <- Methfile_read(control_paths = T1D_controls, case_paths = T1D_cases, WGBS=T)

inputrefseqfile <- Bedfile_read(paths = 'GeneDMRs', bedfile = "refseq", suffix = '.txt')

inputgenebodyfile <- Bedfile_read(paths = 'GeneDMRs', bedfile = "refseq", feature = TRUE, featurewrite = TRUE)
inputcpgifeaturefile <- Bedfile_read(paths = 'GeneDMRs', bedfile = "cpgi", feature = TRUE, featurewrite = FALSE)

### (from here, easier to create 'inputmethfile' file from one of the 'inputmethfile_XX' files)

#replace Chr IDs
chrom_alias <- read.csv('C:/Users/tbehr/Desktop/UU/Bombryo_project/data/reference/bosTau9.chromAlias.txt', sep='\t')
inputmethfile$chr <- chrom_alias$X..ucsc[match(inputmethfile$chr, chrom_alias$refseq)]


# step-by-step
# QC
inputmethfile_QC <- Methfile_QC(inputmethfile)
# Methylation mean (had 60919 unmatched ID)
regiongeneall <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "all")

regiongeneall_Qvalue <- Logic_regression(regiongeneall)

regiongeneall_significant <- Significant_filter(regiongeneall_Qvalue, qvalue=0.01, methdiff=0)


# Circos Plot
cytofile <- Cytofile_read(paths='C:/Users/tbehr/Desktop/UU/Bombryo_project/GeneDMRs', cytofile='cytoUCD1.2', suffix='.txt')
Circos_plot(cytofile, inputmethfile_QC, inputrefseqfile, inputcpgifeaturefile)

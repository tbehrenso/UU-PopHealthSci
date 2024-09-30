library(methylKit)
library(ggplot2)
library(stringr)
library(gridExtra)
library(genomation)
library(sva)
library(VennDiagram)

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

# Quadrant Plot
create_quadrant_plot <- function(){
  dmclist <- read.csv('output/AllAnalysis/dmclist_ED_All_tiled_mpg7_noOC_bothcovar_mincov10_dif10_q0.1.csv')
  ED_conv <- read.csv('C:/Users/tbehr/Desktop/ED_matching_toEntrez.csv')
  # create entrezid column
  dmclist$feature.name <- sapply(strsplit(dmclist$feature.name, '\\.'), `[`, 1)
  dmclist$entrezid <- sapply(dmclist$feature.name, function(x) 
    ifelse(x %in% ED_conv$RefSeq, ED_conv$Entrez[ED_conv$RefSeq==x], NA))
  dmclist <- dmclist[!is.na(dmclist$entrezid),]
  
  entrez_features_unique <- unique(dmclist$entrezid)
  entrez_features_unique <- entrez_features_unique[!is.na(entrez_features_unique)]
  
  # create dataframe that contains relevant info for each entrezid
  dmclist_unique_df <- data.frame(entrezid=entrez_features_unique, meth.diff=0, direction=NA)
  for(i in 1:length(entrez_features_unique)){
    entrezid <- entrez_features_unique[i]
    dmclist_subset <- dmclist[dmclist$entrezid==entrezid,]
    
    if(dim(dmclist_subset)[1]==1){
      dmclist_unique_df$meth.diff[i] <- dmclist_subset$meth.diff
      dmclist_unique_df$direction[i] <- ifelse(dmclist_subset$meth.diff > 0, '+', '-')
    } else {
      dmclist_unique_df$meth.diff[i] <- mean(dmclist_subset$meth.diff)
      dmclist_unique_df$direction[i] <- ifelse(all(dmclist_subset$meth.diff>0), '+', 
                                               ifelse(all(dmclist_subset$meth.diff<0), '-', 'Mix'))
    }
  }
  
  # import DEG data
  deglist <- read.csv('data/processed/D15_ED_AdjData.txt', sep='\t')
  deglist$X <- sapply(strsplit(deglist$X, ':'), `[`, 2)
  
  venn.diagram(x=list(deglist$X, dmclist_entrez_df$entrezid), category.names=c('DEG','DMR'), filename='ED_DEG_DMR_Venn.png')
  
  # find common gene ids and extract
  shared_ids <- intersect(deglist$X, dmclist_unique_df$entrezid)
  combined_df <- data.frame(entrezid=shared_ids, foldchange=NA, methdiff=NA, direction=NA)
  for(i in 1:length(shared_ids)){
    id <- shared_ids[i]
    combined_df$foldchange[i] <- mean(as.numeric(deglist[deglist$X==id, 9:17])) / mean(as.numeric(deglist[deglist$X==id, 2:8]))
    combined_df$methdiff[i] <- dmclist_unique_df$meth.diff[dmclist_unique_df$entrezid==id]
    combined_df$direction[i] <- dmclist_unique_df$direction[dmclist_unique_df$entrezid==id]
  }
  combined_df$logfc <- log(combined_df$foldchange)
  
  # filter out Inf and NaN values
  combined_df <- combined_df[is.finite(combined_df$logfc),]
  
  combined_df$category <- ifelse(combined_df$logfc > 0 & combined_df$methdiff > 0, 'HyperM-UpReg',
                                 ifelse(combined_df$logfc > 0 & combined_df$methdiff < 0, 'HypoM-UpReg',
                                        ifelse(combined_df$logfc < 0 & combined_df$methdiff > 0, 'HyperM-DownReg',
                                               ifelse(combined_df$logfc < 0 & combined_df$methdiff < 0, 'HypoM-DownReg',
                                                      NA))))
  
  ggplot(combined_df, aes(x=methdiff, y=logfc, color=category)) +
    geom_point()
  
}


PLS_integration <- function(){
  library(mixOmics)
  library(data.table)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  
  ED_matching_methyl <- as.data.frame(
    data.table::fread('data/processed/AA_Data/ED_matching_TILED_mincov5_comBatched_mogFormat.csv',header = T))
  methyl_colnames_noaccession <- sapply(strsplit(colnames(ED_matching_methyl), '\\.'), `[`, 1)
  methyl_refseq_unique <- unique(methyl_colnames_noaccession)[-1]
  
  pls.result <- pls(X, Y)
  plotIndiv(pls.result)
  plotVar(pls.result)
}

# create venn diagram to find common features from moGCN output
moGCN_output_comparison <- function(){
  topn_A <- read.csv('C:/Users/tbehr/Desktop/UU/MoGCN/MoGCN/result/ED_Results_e1000/topn_omics_2.csv')
  topn_B <- read.csv('C:/Users/tbehr/Desktop/UU/MoGCN/MoGCN/result/ED_Results_e1000/topn_omics_3.csv')
  
  # remove _x or _y
  epoch1000_A <- sapply(strsplit(topn_A$epoch_1000, '_'), `[`, 1)
  epoch1000_B <- sapply(strsplit(topn_B$epoch_1000, '_'), `[`, 1)
  
  venn.diagram(list(A=epoch1000_A, B=epoch1000_B), filename='mogcn_compAB.png')
  venn.diagram(list(A=topn_B$epoch_10, B=topn_B$epoch_20), filename='mogcn_TE_20v1000.png')
}

belen_id_conversion <- function(){
  library(mygene)
  library(data.table)
  
  dataAll <- as.data.frame(data.table::fread("data/processed/AA_Data/ED_matching_TILED_merged_annotated_mincov5.csv",header = T))
  
  genes <- queryMany(dataAll$feature.name, scopes="refseq.rna", fields=c("symbol","name","entrezgene","type_of_gene"), species="9913")
  geneList <- data.frame(genes)
  dataAll <- cbind(dataAll,geneList); colnames(dataAll)
  
  dataEpigenome <- data.frame(row.names=gsub(" ","",paste(dataAll$chr,dataAll$start,sep="_")),
                              Refseq=dataAll$feature.name, Entrez=dataAll$entrezgene,
                              Symbol=dataAll$symbol, GeneName=dataAll$name,
                              Region=dataAll$transcript.assoc,Island=dataAll$island.assoc,
                              IVP10=dataAll$numCs1*100/dataAll$coverage1,
                              IVP11=dataAll$numCs2*100/dataAll$coverage2,
                              IVP14=dataAll$numCs3*100/dataAll$coverage3,
                              IVPF10=dataAll$numCs4*100/dataAll$coverage4,
                              IVPF15=dataAll$numCs5*100/dataAll$coverage5,
                              VE24=dataAll$numCs6*100/dataAll$coverage6,
                              VE26=dataAll$numCs7*100/dataAll$coverage7,
                              VE28=dataAll$numCs8*100/dataAll$coverage8,
                              VE29=dataAll$numCs9*100/dataAll$coverage9,
                              VE4=dataAll$numCs10*100/dataAll$coverage10,
                              VE5=dataAll$numCs11*100/dataAll$coverage11,
                              VE6=dataAll$numCs12*100/dataAll$coverage12,
                              VE3=dataAll$numCs13*100/dataAll$coverage13)
  
  
  
}


# temp code for subsetting data 
tempcode <- function(){

  diff_methylation_25p <- diff_methylation[diff_methylation$chr=='chr10']
  
  all_known_features[all_known_features$start>53000000 & all_known_features$end<54000000,]
  
  
}


get_gene_ids <- function(){
  library(mygene)
  library(data.table)
  
  genes <- queryMany(all_known_features$feature.name, scopes="refseq.rna", fields=c("symbol","name","entrezgene","type_of_gene"), species="9913")
  geneList <- data.frame(genes)
  
  all_known_features_geneid_belen <- cbind(all_known_features, geneList$symbol)
  
  all_known_features_geneid_belen[grepl('CDH3', all_known_features_geneid_belen$`geneList$symbol`),]
}


make_volcano_plot <- function(){
  
  
}


if(F){
  # just testing converting DEG list from Entrez to gene symbol
  library(biomaRt)
  library(mygene)
  
  deglist <- read.csv('data/processed/D15_ED_AdjData.txt', sep='\t')
  
  deg_entrez_ids <- sapply(strsplit(deglist[,1], ':'), `[`, 2)
  
  genes <- queryMany(deg_entrez_ids, scopes="entrezgene", fields=c("symbol","name","entrezgene","type_of_gene"), species="9913")
  geneList <- data.frame(genes)
  
  
  
  # mart <- useDataset('btaurus_gene_ensembl', mart = useMart('ensembl'))
  # 
  # deg_gene_symbols <- getBM(filters=c('entrezgene_id'), attributes=c('entrezgene_id', 'external_gene_name'), values=deg_entrez_ids, mart=mart)

}

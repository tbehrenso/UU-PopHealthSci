library(stringr)
library(ggplot2)
library(ggrepel)

### ED PCA 

ED_adj <- read.csv('data/D15_ED_AdjData.txt', sep='\t')

pca_stats <- prcomp(ED_adj[-1])

# getting morphology, sex, and treatment information
treatment_binary <- as.numeric(str_detect(colnames(ED_adj)[-1], 'IVP'))
sample_ids <- colnames(ED_adj)[-1]
sample_info <- read.table('data/reference/SampleInfo.csv', sep=',', header = T)
batch_ordered <- c('Short', 'Short', 'Short', 'Short', 'Round', 'Round', 'Round', 
                   'Long', 'Short', 'Short', 'Short', 'Short', 'Long', 'Long', 'Long', 'Long')
sex_ordered <- c('Male', 'Male', 'Male', 'Male', 'Male', 'Male', 'Male',
                 'Male','Female','Female','Female','Female','Male','Female','Male','Female')


pca_df <- data.frame(PC1=pca_stats$rotation[,1], PC2=pca_stats$rotation[,2], type=as.factor(treatment_binary),
                     morph=batch_ordered, sex=sex_ordered)
ggplot(pca_df, aes(x=PC1, y=PC2, col=type, label=sapply(strsplit(row.names(pca_df),'_'), `[`, 1), shape=morph)) + 
  geom_point(size=3) + geom_text_repel()
ggplot(pca_df, aes(x=PC1, y=PC2, col=morph, label=sapply(strsplit(row.names(pca_df),'_'), `[`, 1), shape=type)) + 
  geom_point(size=3) + geom_text_repel() +
  scale_shape_manual(values=c(1, 16), labels=c('VE', 'IVP')) +
  scale_color_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))

### TE PCA 

TE_adj <- read.csv('data/D15_TE_AdjData.txt', sep='\t')

pca_stats <- prcomp(TE_adj[-1])

# getting morphology, sex, and treatment information
treatment_binary <- as.numeric(str_detect(colnames(TE_adj)[-1], 'IVP'))
sample_ids <- colnames(TE_adj)[-1]
sample_info <- read.table('data/reference/SampleInfo.csv', sep=',', header = T)
batch_ordered <- c('Short', 'Short', 'Short', 'Round', 'Round', 'Round', 
                   'Long', 'Short', 'Short', 'Short', 'Short', 'Long', 'Long', 'Long')
sex_ordered <- c('Male', 'Male', 'Male', 'Male', 'Male', 'Female',
                 'Male','Female','Female','Female','Female','Male','Female','Male')


pca_df <- data.frame(PC1=pca_stats$rotation[,1], PC2=pca_stats$rotation[,2], type=as.factor(treatment_binary),
                     morph=batch_ordered, sex=sex_ordered)
ggplot(pca_df, aes(x=PC1, y=PC2, col=type, label=sapply(strsplit(row.names(pca_df),'_'), `[`, 1), shape=morph)) + 
  geom_point(size=3) + geom_text_repel()
ggplot(pca_df, aes(x=PC1, y=PC2, col=morph, label=sapply(strsplit(row.names(pca_df),'_'), `[`, 1), shape=type)) + 
  geom_point(size=3) + geom_text_repel() +
  scale_shape_manual(values=c(1, 16), labels=c('VE', 'IVP')) +
  scale_color_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))

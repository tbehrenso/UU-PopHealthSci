library(tidyverse)
library(ggplot2)
library(stringr)
library(reshape2) 

#FILENAME <- './IVPF15_Paired_IgnBoth10/IVPF15_ED_1_bismark_bt2_pe.M-bias.txt'
FILENAME <- 'data/raw/MBias/IVP10_T1D_NoIgn_1_bismark_bt2_pe.M-bias.txt'
COLNAMES <- c('position', 'count', 'methylated', '%methylation', 'coverage')

sample_id <- str_sub(str_match(FILENAME, '[/]\\w+\\d+_\\w+[/]'), 2, -2)

mbias_data_raw <- read.table(FILENAME, sep = '\n', as.is = T)

rawtext <- readLines(FILENAME)

context_matches <- na.omit(str_match(rawtext,".+context.+"))[,1]

#contexts <- lapply(context_matches, function(x) substr(x, 1, 3))
contexts <- context_matches

split_text <- split(rawtext, cumsum(!grepl('[^,\t]', rawtext)))

# remove empty lines
split_text <- lapply(split_text, function(x) x[nzchar(x)])

# remove empty list elements
split_data <- split_text[lengths(split_text)>0]

methylation_data <- split_data %>% 
  lapply(function(x) tail(x, -3) %>% 
           strsplit('\t') %>% 
           as.data.frame() %>% 
           t() %>% 
           type.convert(as.is = TRUE) %>% 
           `rownames<-`(NULL) %>% 
           `colnames<-`(COLNAMES) %>% 
           as.data.frame()
  )

names(methylation_data) <- context_matches

# add column to each that indicates context
methylation_data[[1]]$context <- contexts[[1]]
methylation_data[[2]]$context <- contexts[[2]]
methylation_data[[3]]$context <- contexts[[3]]
methylation_data[[4]]$context <- contexts[[4]]
methylation_data[[5]]$context <- contexts[[5]]
methylation_data[[6]]$context <- contexts[[6]]


methylation_combined_R1 <- rbind(methylation_data[[1]], methylation_data[[2]], methylation_data[[3]])
methylation_combined_R2 <- rbind(methylation_data[[4]], methylation_data[[5]], methylation_data[[6]])


# choose if plotting forward or reverse read
##### should make it not manual later #########################
methylation_combined <- methylation_combined_R1


# scaling the methylation counts and combining to the same column as percent methylation
# This is allows the legend to be generated automatically
methylation_percmeth_as_countscaled <- methylation_combined
methylation_percmeth_as_countscaled$`%methylation` <- methylation_percmeth_as_countscaled$count / max(methylation_combined$count) * 100
# add columns to ID values after being combined
methylation_combined$measurement <- 'percentmethyl'
methylation_percmeth_as_countscaled$measurement <- 'methylcount'

methylation_combined_stacked <- rbind(methylation_combined, methylation_percmeth_as_countscaled)
methylation_combined_stacked$measurement <- as.factor(methylation_combined_stacked$measurement)

ggplot(data = methylation_combined_stacked) +
  geom_line(aes(x = position, y = `%methylation`, group=interaction(context,measurement), colour=context, linetype=measurement), linewidth=1) +
  scale_y_continuous(limits = c(0,100), sec.axis = sec_axis(~ . * max(methylation_combined$count)/100, name = 'Methylated Sites Count')) +
  #geom_vline(xintercept = 10, linetype='dashed') +
  ylab('Percent Methylation') + xlab('Position along read (bp)') +
  labs(colour='Context', linetype='Value') +
  scale_linetype_manual(labels = c('Methylated\nSites Count','Percent\nMethylation'), values = c('dashed','solid')) +
  scale_color_manual(values=c('#377eb8', '#4daf4a', '#e41a1c')) +
  ggtitle(paste('M-Bias for', sample_id)) +
  theme_classic(base_size = 20) +
  geom_vline(xintercept = 10, linetype='dashed',linewidth=1)



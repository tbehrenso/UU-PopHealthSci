library(ggplot2)

data<-read.csv2("C:/Users/tbehr/Desktop/enrichment_DownBM_graphCSV.csv") #it might be read.csv
data$EnrichmentFDR<-as.numeric(data$EnrichmentFDR)
data$logFDR<--log10(data$EnrichmentFDR)
colnames(data)

tiff(filename="Plot-UpRegBioMark.tiff",units="in", width=8, heigh=4, res=300)
ggplot(data, aes(x=Pathway, y=logFDR, colour=logFDR)) +
  geom_segment( aes(x=reorder(Pathway,logFDR), xend=Pathway, y=0, yend=logFDR)) +
  geom_point(aes(size=nGenes)) +
  labs(col="-log10(FDR)",size="Genes#") +
  scale_colour_continuous(low = 'blue', high = 'red') + #color of the lollipop
  coord_flip() + #flip the graph
  xlab("") + ylab("") + #labels for the axis
  theme(
    axis.text.y = element_text(color="black",size = 12),
    axis.text.x = element_text(color="black",size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank())

dev.off()

data<-read.csv("C:/Users/tbehr/Desktop/PLOTTING_TE_ALL_CLUSTERED.csv", header=TRUE)
data$logFDR <- -log(data$FDR, 10)
data <- data[!duplicated(data$Term), ]

data$termsize <- rescale(data$logFDR, to = c(68,70))
data$Category <- factor(data$Category, levels=unique(data$Category))

cols <- c("royalblue",'green','yellow','purple')
cols <- c("royalblue",'green')
cols <- c('yellow','purple')
cols <- c('red')

pdf("BubblePlot.pdf",width=7, height=4)
require("ggrepel")
ggplot(data, aes(Term, logFDR, fill = Category, size = Count))+
  labs(title = NULL, x = NULL, y = '-logFDR')+
  geom_point(shape = 21, col = 'black', alpha = 0.8)+
  scale_size(range = c(1,12), guide='none') +
  facet_grid(.~Category, space = 'free_x', scales = 'free_x')+ 
  scale_fill_manual(values = cols, guide ='none') +
  geom_text_repel (aes(label = Term,size = termsize),color = 'black',box.padding = unit(0.4, "lines"), 
                   point.padding = unit(0.3, "lines"),arrow = arrow(length = unit(0.01, 'npc'))) +
  ylab("-log10 FDR") +
  theme(axis.text.y = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=10),
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'grey30'), 
        panel.border = element_rect(fill = 'transparent', colour = 'grey30'),
        panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank()) 
dev.off()






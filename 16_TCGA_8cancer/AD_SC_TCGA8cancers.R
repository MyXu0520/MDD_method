#shared_MDD cluster
TCGA8_sample_CMDD <- matrix(nrow = length(shared_MDD_cg))
for (i in 1:length(all_data)) {
  temp_sample <- all_data[[i]][shared_MDD_cg,]
  TCGA8_sample_CMDD <- cbind(TCGA8_sample_CMDD,temp_sample)
}
TCGA8_sample_CMDD <- TCGA8_sample_CMDD[,-1]
set.seed(123)
# gsub("\\.","_",colnames(merged_data))
merged_data_CMDD <- t(TCGA8_sample_CMDD)
tsne_result_CMDD <- Rtsne(merged_data_CMDD, dims = 2, perplexity = 30, verbose = FALSE,theta = 0.0,max_iter = 1000)
tsne_out_CMDD = as.data.frame(tsne_result_CMDD$Y)
colnames(tsne_out_CMDD) <- c("TSNE_1","TSNE_2")
library(plotly)
# library(Seurat)
library(cols4all)
library(tidydr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(knitr)
mycol <- c4a('10',8)
mytheme <- theme_void() + theme(plot.margin = margin(5.5,15,5.5,5.5))
p <- ggplot(tsne_out_CMDD, aes(x = TSNE_1, y = TSNE_2)) +
  # stat_ellipse(aes(fill = group), level = 0.95, linetype = 1,
  #              show.legend = F, geom = 'polygon', alpha = 0.1) +
  geom_point(aes(colour  = group), size =1.0, alpha = 1) +
  mytheme +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) +
  theme(panel.grid = element_blank())+
  # geom_text(data = group, aes(x = TSNE_1, y = TSNE_2, label = group), fontface = "bold", color = 'black', size = 4) +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_fill_manual(values = mycol) +
  # scale_fill_brewer(palette = "Spectral")+
  # scale_color_brewer(palette = "Spectral")
  scale_color_manual(values = mycol)
# scale_color_gradientn(values = seq(0,1,0.2),
# colours = c('cyan','blue','green','orange','red'))
p
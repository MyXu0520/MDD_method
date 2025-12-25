#shared_MDD
setwd("D:/MDD_result/9/ESCA_MDD_meth")
file <- list.files()
file <- file[-length(file)]
temp_file <- read.table(file[1])
chr_start_end <- paste0(temp_file$V1,"_",temp_file$V2,"_",temp_file$V3)
tsne_data <- matrix(ncol = length(chr_start_end))
for (i in 1:45) {
  temp <- read.table(file[i])
  temp_pattern <- paste0(temp$V1,"_",temp$V2,"_",temp$V3)
  rownames(temp) <- temp_pattern
  temp <- temp[chr_start_end,]
  tsne_data <- rbind(tsne_data,temp$V4)
}
tsne_data <- tsne_data[-1,]
rownames(tsne_data) <- file
colnames(tsne_data) <- chr_start_end
library(tsne)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(Rtsne)
tsne_out <- Rtsne(
  tsne_data,
  dims = 2,
  pca = FALSE,
  perplexity = 10,
  theta = 0.0,
  max_iter = 1000
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("TSNE_1","TSNE_2")
group <- c(rep("EAC",5),rep("ESCC",21),rep("ESCC_N",5),rep("EAC",7),rep("EAC_N",7))
group <- as.factor(group)
library(ggplot2)
library(pals)
library(dplyr)
define_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_text(color='black', size=18),
  axis.ticks.length = unit(0.4, "lines"),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size=18),
  legend.key = element_blank(),
  legend.key.size = unit(1, 'cm')
)
# define_color1 <- pals::glasbey()[1:4]
# define_color2 <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
#                    "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
#                    "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
#                    "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
#                    "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
#                    "#aeae5c","#1e90ff","#00bfff","#BBDD78","#ffff00",
#                    "#F3AE63","#93CBaE","#73558B")[1:4]
# color <- c("#FFCD44","#EE7C79", "#008F91", "#4489C8")
# df_label <- summarise(group_by(data, cluster),
#                       x_pos = mean(TSNE_1),
#                       y_pos = mean(TSNE_2))
p <- ggplot(tsne_result, aes(x = TSNE_1, y = TSNE_2)) +
  geom_point(aes(color = group), alpha = .8,size = 3) +
  scale_color_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))+theme_bw()
p+stat_ellipse(data=tsne_result,
                geom = "polygon",level=0.8,
                linetype = 1,size=0.8,
                aes(fill=group),
                alpha=0.1,
                show.legend = T)+
  scale_fill_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))
#ESCC_SPEC
setwd("D:/MDD_result/9/ESCA_ESCC_MDD_meth")
file <- list.files()
file <- file[-length(file)]
temp_file <- read.table(file[1])
chr_start_end <- paste0(temp_file$V1,"_",temp_file$V2,"_",temp_file$V3)
tsne_data <- matrix(ncol = length(chr_start_end))
for (i in 1:45) {
  temp <- read.table(file[i])
  temp_pattern <- paste0(temp$V1,"_",temp$V2,"_",temp$V3)
  rownames(temp) <- temp_pattern
  temp <- temp[chr_start_end,]
  tsne_data <- rbind(tsne_data,temp$V4)
}
tsne_data <- tsne_data[-1,]
rownames(tsne_data) <- file
colnames(tsne_data) <- chr_start_end
library(tsne)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(Rtsne)
tsne_out <- Rtsne(
  tsne_data,
  dims = 2,
  pca = FALSE,
  perplexity = 10,
  theta = 0.0,
  max_iter = 1000
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("TSNE_1","TSNE_2")
group <- c(rep("EAC",5),rep("ESCC",21),rep("ESCC_N",5),rep("EAC",7),rep("EAC_N",7))
group <- as.factor(group)
library(ggplot2)
library(pals)
library(dplyr)
define_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_text(color='black', size=18),
  axis.ticks.length = unit(0.4, "lines"),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size=18),
  legend.key = element_blank(),
  legend.key.size = unit(1, 'cm')
)
# define_color1 <- pals::glasbey()[1:4]
# define_color2 <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
#                    "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
#                    "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
#                    "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
#                    "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
#                    "#aeae5c","#1e90ff","#00bfff","#BBDD78","#ffff00",
#                    "#F3AE63","#93CBaE","#73558B")[1:4]
# color <- c("#FFCD44","#EE7C79", "#008F91", "#4489C8")
# df_label <- summarise(group_by(data, cluster),
#                       x_pos = mean(TSNE_1),
#                       y_pos = mean(TSNE_2))
p <- ggplot(tsne_result, aes(x = TSNE_1, y = TSNE_2)) +
  geom_point(aes(color = group), alpha = .8,size = 3) +
  scale_color_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))+theme_bw()
p+stat_ellipse(data=tsne_result,
               geom = "polygon",level=0.8,
               linetype = 1,size=0.8,
               aes(fill=group),
               alpha=0.1,
               show.legend = T)+
  scale_fill_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))
#CMI
setwd("D:/MDD_result/9/ESCA_CMI_meth")
file <- list.files()
file <- file[-length(file)]
temp_file <- read.table(file[1])
chr_start_end <- paste0(temp_file$V1,"_",temp_file$V2,"_",temp_file$V3)
tsne_data <- matrix(ncol = length(chr_start_end))
for (i in 1:45) {
  temp <- read.table(file[i])
  temp_pattern <- paste0(temp$V1,"_",temp$V2,"_",temp$V3)
  rownames(temp) <- temp_pattern
  temp <- temp[chr_start_end,]
  tsne_data <- rbind(tsne_data,temp$V4)
}
tsne_data <- tsne_data[-1,]
rownames(tsne_data) <- file
colnames(tsne_data) <- chr_start_end
library(tsne)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(Rtsne)
tsne_out <- Rtsne(
  tsne_data,
  dims = 2,
  pca = FALSE,
  perplexity = 10,
  theta = 0.0,
  max_iter = 1000
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("TSNE_1","TSNE_2")
group <- c(rep("EAC",5),rep("ESCC",21),rep("ESCC_N",5),rep("EAC",7),rep("EAC_N",7))
group <- as.factor(group)
library(ggplot2)
library(pals)
library(dplyr)
define_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_text(color='black', size=18),
  axis.ticks.length = unit(0.4, "lines"),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size=18),
  legend.key = element_blank(),
  legend.key.size = unit(1, 'cm')
)
# define_color1 <- pals::glasbey()[1:4]
# define_color2 <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
#                    "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
#                    "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
#                    "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
#                    "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
#                    "#aeae5c","#1e90ff","#00bfff","#BBDD78","#ffff00",
#                    "#F3AE63","#93CBaE","#73558B")[1:4]
# color <- c("#FFCD44","#EE7C79", "#008F91", "#4489C8")
# df_label <- summarise(group_by(data, cluster),
#                       x_pos = mean(TSNE_1),
#                       y_pos = mean(TSNE_2))
p <- ggplot(tsne_result, aes(x = TSNE_1, y = TSNE_2)) +
  geom_point(aes(color = group), alpha = .8,size = 3) +
  scale_color_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))+theme_bw()
p+stat_ellipse(data=tsne_result,
               geom = "polygon",level=0.8,
               linetype = 1,size=0.8,
               aes(fill=group),
               alpha=0.1,
               show.legend = T)+
  scale_fill_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))
# EAC_SPEC
#EAC_SPEC
setwd("D:/MDD_result/9/ESCA_EAC_MDD_meth")
file <- list.files()
file <- file[-length(file)]
temp_file <- read.table(file[1])
chr_start_end <- paste0(temp_file$V1,"_",temp_file$V2,"_",temp_file$V3)
tsne_data <- matrix(ncol = length(chr_start_end))
for (i in 1:45) {
  temp <- read.table(file[i])
  temp_pattern <- paste0(temp$V1,"_",temp$V2,"_",temp$V3)
  rownames(temp) <- temp_pattern
  temp <- temp[chr_start_end,]
  tsne_data <- rbind(tsne_data,temp$V4)
}
tsne_data <- tsne_data[-1,]
rownames(tsne_data) <- file
colnames(tsne_data) <- chr_start_end
library(tsne)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(Rtsne)
tsne_out <- Rtsne(
  tsne_data,
  dims = 2,
  pca = FALSE,
  perplexity = 10,
  theta = 0.0,
  max_iter = 1000
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("TSNE_1","TSNE_2")
group <- c(rep("EAC",5),rep("ESCC",21),rep("ESCC_N",5),rep("EAC",7),rep("EAC_N",7))
group <- as.factor(group)
library(ggplot2)
library(pals)
library(dplyr)
define_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_text(color='black', size=18),
  axis.ticks.length = unit(0.4, "lines"),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size=18),
  legend.key = element_blank(),
  legend.key.size = unit(1, 'cm')
)
# define_color1 <- pals::glasbey()[1:4]
# define_color2 <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
#                    "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
#                    "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
#                    "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
#                    "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
#                    "#aeae5c","#1e90ff","#00bfff","#BBDD78","#ffff00",
#                    "#F3AE63","#93CBaE","#73558B")[1:4]
# color <- c("#FFCD44","#EE7C79", "#008F91", "#4489C8")
# df_label <- summarise(group_by(data, cluster),
#                       x_pos = mean(TSNE_1),
#                       y_pos = mean(TSNE_2))
p <- ggplot(tsne_result, aes(x = TSNE_1, y = TSNE_2)) +
  geom_point(aes(color = group), alpha = .8,size = 3) +
  scale_color_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))+theme_bw()
p+stat_ellipse(data=tsne_result,
               geom = "polygon",level=0.8,
               linetype = 1,size=0.8,
               aes(fill=group),
               alpha=0.1,
               show.legend = T)+
  scale_fill_manual(values = c("#FFCD44","#EE7C79", "#008F91", "#4489C8"))
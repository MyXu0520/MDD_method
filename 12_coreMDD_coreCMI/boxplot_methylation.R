###Code for WGBS project######
options(scipen = 20)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(plyr)
library(DelayedMatrixStats)
library(tibble)
library(valr)
library(ggbeeswarm)
library(ggpubr)
library(VennDiagram)
library(dmrseq)
library(bsseq)
library(gridExtra)
library(readr)
library(vipor)
library(dplyr)
library(formattable)
library(stringr)
# library(reshape2)
library(BSgenome)
library(MutationalPatterns)
################################################meta_info################################################
EAC_Tumor_sampleList=c(paste0("EAC_", c(1:4,6)), paste0("GEJ_", 1:7))
ESCC_Tumor_sampleList=paste0("ESCC_", c(1:17,19:22))
EAC_Nonmalignant_sampleList=paste0("GEJ_Nonmalignant_", 1:7)
ESCC_Nonmalignant_sampleList=paste0("ESCC_Nonmalignant_", c(1:3, "53F", "54M"))
types=c(rep("ESCC_Tumor",21),rep("ESCC_Nonmalignant", 5), rep("EAC/GEJ_Tumor",12), rep("GEJ_Nonmalignantl",7))
fileList=c(ESCC_Tumor_sampleList, ESCC_Nonmalignant_sampleList, EAC_Tumor_sampleList, EAC_Nonmalignant_sampleList)
chrListTarget=paste("chr",1:22,sep="")
annotation_row=data.frame(Type=c(rep("ESCC_Tumor", 21), rep("ESCC_Nonmalignant", 5), rep("EAC/GEJ_Tumor", 12), rep("GEJ_Nonmalignant", 7)), stringsAsFactors = F)
rownames(annotation_row)=fileList
annotation_row$Type=factor(annotation_row$Type, levels=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor"))
annotation_row=annotation_row[order(annotation_row$Type), ,drop=F]
annotation_col=annotation_row
annotation_col$Sample=rownames(annotation_col)
ann_colors = list(Type = c("#4489C8", "#EE7C79", "#008F91", "#FFCD44"))
names(ann_colors[[1]])=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor")
my_theme=theme_classic()+theme(axis.text = element_text(color="black", size=10),
                               axis.title.y=element_text(color="black", size=12),
                               plot.title=element_text(hjust = 0.5, face="bold", size=14))
shared_MDD <- read.table("D:/MDD_result/10/ESCA_sharedMDDs.txt")
shared_CMI <- read.table("D:/MDD_result/10/intersect_CMI.txt")
methylation_level3 <- cbind(shared_MDD,shared_CMI)
colnames(methylation_level3) <- c("shared_MDD","CMI")
core_MDD <- read.table("D:/MDD_result/12/core_MDD.txt")
core_CMI <- read.table("D:/MDD_result/12/core_CMI.txt")
methylation_level2 <- cbind(core_MDD,core_CMI)
colnames(methylation_level2) <- c("core_MDD","core_CMI")
methylation_level <- cbind(methylation_level,methylation_level2)
library(tidyverse)
library(gaCMInder)
library(ggpubr)
library(rstatix)
library(ggprism)
library(tidyverse)
library(reshape2)
methylation_level$Group <- c(rep("EAC/GEJ_Tumor",5),rep("ESCC_Tumor",21),rep("ESCC_Nonmalignant",5),rep("EAC/GEJ_Tumor",7),rep("GEJ_Nonmalignant",7))
methylation_level$sample <- rownames(methylation_level)
Plotdata <- methylation_level
Plotdata=melt(Plotdata, id.vars = c("sample", "Group"))
Plotdata$value <- as.numeric(Plotdata$value)
Plotdata <- Plotdata[order(Plotdata$Group),]
Plotdata <- Plotdata[c(1:24,77:90,35:76,25:34),]
# data$group <- factor(data$group,levels = c("EAC/GEJ_Tumor", "GEJ_Nonmalignant","ESCC_Tumor","ESCC_Nonmalignant"))
Plotdata$Group <- factor(Plotdata$Group,levels = c("EAC/GEJ_Tumor", "GEJ_Nonmalignant","ESCC_Tumor","ESCC_Nonmalignant"))
mytheme <- theme_prism() +
  theme(strip.text = element_text(size = 18),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 18),
        axis.text.x = element_text(color = "black",size = 16),
        axis.title = element_text(color = "black",size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
        legend.position = "none")
p <- ggplot(Plotdata,aes(x = Group,y = value)) +
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.4),width = 0.1) +
  geom_boxplot(position = position_dodge(width = 0.4),
               outlier.shape = NA) +
  geom_signif(comparisons =  list(c("EAC/GEJ_Tumor", "GEJ_Nonmalignant"),c("ESCC_Tumor","ESCC_Nonmalignant")),y_position = c(0.9,0.95), map_signif_level = T)+
  geom_point(aes(fill = Group),
             pch = 21,size = 2,
             position = position_jitter(0.2))+
  facet_wrap(.~ variable,scales = "free_x",nrow = 1) +
  scale_x_discrete(guide = "prism_bracket") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25),
                     # minor_breaks = seq(0,1,0.25),
                     guide = "prism_offset_minor") +
  # scale_fill_discrete()+
  scale_fill_manual(values = c("#FFCD44","#EE7C79","#008F91","#4489C8") )+
  # scale_fill_discrete()+
  # scale_fill_manual(values = ann_colors)
  labs(x = NULL,y = "Methylation_level") +
  # theme(axis.text.x = element_blank())+
  mytheme
p <- p + theme(legend.position='right',axis.text.x = element_blank(),axis.text.y = element_text(size = 10),axis.title = element_text(size = 10),plot.title = element_text(size = 10))
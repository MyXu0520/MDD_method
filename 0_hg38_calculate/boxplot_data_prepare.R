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
###########################################global_methylation_plot
# data <- read.table("C:/Users/MyXu/Desktop/global.txt")
setwd("D:/MDD_result/0/")
# methylation_level <- matrix(ncol = length(list.files("C:/Users/MyXu/Desktop/result")))
methylation_level <- matrix(nrow = 45)
methylation_level <- as.data.frame(methylation_level)
order_name <- c(EAC_Tumor_sampleList,EAC_Nonmalignant_sampleList,ESCC_Tumor_sampleList,ESCC_Nonmalignant_sampleList)
b <- list.files("D:/MDD_result/0/",pattern = ".txt")
for (i in 1:length(b)) {
  if(b[i] == "global.txt"){
    temp_file <- read.table(b[i])
    temp_file <- as.data.frame(temp_file)
    colnames(temp_file) <- c("sample","beta")
    rownames(temp_file) <- temp_file$sample
    temp_file <- temp_file[order_name,]
    # methylation_level <- cbind(methylation_level,temp_file[,2])
    methylation_level[,i] <- temp_file[,2]
    colnames(methylation_level)[i] <- b[i]
  }
  else if(b[i] == "CGI.txt"){
    temp_file <- read.table(b[i])
    temp_file3 <- as.data.frame(temp_file)
    colnames(temp_file3)<- "beta"
    temp_file3$sample <- rownames(temp_file3)
    temp_file3 <- temp_file3[order_name,]
    # methylation_level <- cbind(methylation_level,temp_file3$beta)
    methylation_level[,i] <- temp_file3$beta
    colnames(methylation_level)[i] <- b[i]
  }
  else{
    temp_file <- read.table(b[i])
    temp_file2 <- data.frame()
    temp_file2 <- cbind(rownames(temp_file),temp_file[[1]])
    temp_file2 <- as.data.frame(temp_file2)
    colnames(temp_file2) <- c("sample","beta")
    rownames(temp_file2) <- temp_file2$sample
    temp_file2 <- temp_file2[order_name,]
    # methylation_level <- cbind(methylation_level,temp_file2[,2])
    methylation_level[,i] <- temp_file2[,2]
    colnames(methylation_level)[i] <- b[i]
  }
}
rownames(methylation_level) <- order_name
save(methylation_level,file = "methylation_level_differentRegion.RData")
write.table(methylation_level,"methylation_level_differentRegion.txt")
################################################Figure1B################################################
#bash Shell/Figure1B/getGlobalMeanBetaValues.sh global
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/commonMDDs_hg38.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/gencode.v31.basic.promoter.Takai_Jones.CGI.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/hg38_repeat.change.LINE.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/hg38_repeat.change.LTR.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/hg38_repeat.change.SINE.bed
###combine the result and get plotdata "Figure1B/Methylaton_sixTypes.txt"
# plotdata=read.table("Figure1B/Methylaton_sixTypes.txt", sep="\t", stringsAsFactors = F, header=T)
# plotdata <- methylation_level
# plotdata=plotdata[plotdata$Type%in%c("ESCC_Nonmalignant", "GEJ_Nonmalignant","ESCC_Tumor","EAC/GEJ_Tumor"),]
# plotdata$Sample=factor(plotdata$Sample, levels=fileList)
# plotdata$Type=factor(plotdata$Type, levels=c("ESCC_Nonmalignant", "GEJ_Nonmalignant","ESCC_Tumor","EAC/GEJ_Tumor"))
# plotdata=melt(plotdata, id.vars = c("Sample", "Type"))
# colnames(plotdata)=c("Sample", "Type", "DataType", "Methylation")
# plotdata$DataType=factor(plotdata$DataType, levels=unique(plotdata$DataType))
# plotDotPlot=function(myPlotData, saveFile){
#   p<-ggplot(myPlotData, aes(x=Type, y=Methylation, fill=Type, color=Type)) +
#     # geom_dotplot(binaxis='y', stackdir='center', dotsize=1.2) + facet_wrap(~ DataType, nrow=1)
#     geom_boxplot(binaxis='y', stackdir='center', dotsize=1.2) + facet_wrap(~ DataType, nrow=1)
#   # p=p+scale_color_manual(values=c("#EA3323", "#0000F5", "#EA3323", "#0000F5"))
#   p=p+scale_color_manual(values=c("#4489C8", "#EE7C79", "#008F91", "#FFCD44"))
#
#   p=p+scale_fill_manual(values=ann_colors$Type[1:4])
#   p=p+ylim(0,1)+xlab("")+ylab("Mean methylation")
#   # p=p+ stat_summary(fun=mean, geom="point", shape=18, size=3, color="black")
#   p=p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
#                      geom="pointrange", color="black")
#   p=p+theme_bw()+theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),
#                        axis.text = element_text(color="black", size=10),
#                        axis.title.y=element_text(color="black", size=12))
#   print(p)
#   pdf(saveFile, width=20, height=3.5)
#   print(p)
#   dev.off()
# }
#
# plotDotPlot(plotdata, "Figure1B/methylation_sixType2.pdf")
library(tidyverse)
library(gaCMInder)
library(ggpubr)
library(rstatix)
library(ggprism)
library(tidyverse)
library(reshape2)
methylation_level <- methylation_level[,-c(1,6)]
methylation_level$Group <- c(rep("EAC/GEJ_Tumor",12),rep("GEJ_Nonmalignant",7),rep("ESCC_Tumor",21),rep("ESCC_Nonmalignant",5))
# Plotdata$Group <- as.factor(Plotdata$Group)
methylation_level$sample <- rownames(methylation_level)
Plotdata <- methylation_level
Plotdata=melt(Plotdata, id.vars = c("sample", "Group"))
Plotdata$value <- as.numeric(Plotdata$value)
Plotdata$Group <- as.factor(Plotdata$Group)
Plotdata <- Plotdata[c(1:84,267:315,85:119,120:266),]
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
  scale_fill_manual(values = c("#FFCD44","#EE7C79", "#008F91","#4489C8")) +
  # scale_fill_discrete()+
  # scale_fill_manual(values = ann_colors)
  labs(x = NULL,y = "Methylation_level") +
  # theme(axis.text.x = element_blank())+
  mytheme
p <- p + theme(legend.position='right',axis.text.x = element_blank(),axis.text.y = element_text(size = 10),axis.title = element_text(size = 10),plot.title = element_text(size = 10))
# my_comparisons = list(c("EAC/GEJ_Tumor", "GEJ_Nonmalignant"), c("ESCC_Tumor" ,"ESCC_Nonmalignant"))
# p + stat_compare_means(comparisons = my_comparisons,
#                        label = "p.signif",
#                        method = "t.test")
# p <- p+ stat_compare_means(method = "anova")
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size = rel(1.6)))
# p+theme(axis.text = element_text(size = 20,angle = 50))
  # ann_colors = list(Type = c( "#FFCD44","#EE7"#FFCD44","#EE7C79", "#008F91","#4489C8"C79", "#008F91","#4489C8"))
# p+ theme(legend.title =element_text(size= 20 ),
         # legend.text =element_text(size= 14 ))
# p+scale_fill_discrete()
# df_p_val1 <- Plotdata %>% group_by(variable) %>%
#   wilcox_test(value  ~ Group) %>%
#   adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   add_significance(p.col = "p.adj") %>%
#   add_xy_position(step.increase = 0.18)
p3 <- p +
  geom_smooth(method = "lm",
              size = 1,se = T,
              color = "black",
              linetype = "dashed",
              aes(group = 1))+
  stat_cor(label.y = -2.2,size = 6,
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),group = 1),
           color = "black",method = "spearman",
           label.x.npc = "left") +
  stat_regline_equation(label.y = -2.75,size = 6,
                        aes(group = 1),color = "red")
save.image(file = "boxplot_differentRegion.RData")
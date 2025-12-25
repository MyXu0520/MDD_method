# heatmap
setwd("C:/Users/MyXu/Desktop/")
data <- read.csv("D:/MDD_result/2/merged_methylation_data.csv")
rownames(data) <- data$Region_ID
gdata <- data[,-1]
cdata <- colnames(data[,-1])
cdata <- as.data.frame(cdata)
cdata$group <- c(rep("EAC_T",5),rep("ESCC_T",21),rep("ESCC_N",5),rep("EAC_N",7),rep("EAC_N",7))
cdata$state <- c(rep(1,26),rep(0,19))
gdata <- t(gdata)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
# Heatmap(gdata[,grep("chr1",colnames(gdata))],show_column_names = F,clustering_method_rows = "complete",
#         cluster_rows = T,cluster_columns = F)
#         row_split = nrow(gdata),
#         name = "heat",
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         show_column_names = FALSE,
# Heatmap(gdata, clustering_method_columns = "complete",
#         cluster_columns = TRUE,show_column_names = F,column_names_rot = 90, show_row_names = F, cluster_rows = F,
#         row_split = 3, column_split = 3,name = "heat",
# Heatmap(gdata,col = colorRamp2(c(0,2,4), c("#6395C7", "white", "firebrick3")),
#         name = "Expression",
#         column_split = celltype_info,
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         show_column_names = FALSE,
#         show_row_names = TRUE)
# pheatmap(gdata,
#          )
cdata$group <- as.factor(cdata$group)
cdata$state <- as.factor(cdata$state)
annotation_row <- data.frame(cdata[,c(2,3)])
pdf("heatmap.pdf")
pheatmap(gdata[,grep("chr2",colnames(gdata))],
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         fontsize_row=11,
         show_colnames = F,
         annotation_row = annotation_row
         )
dev.off()
chr_methy <- read.csv("C:/Users/MyXu/Desktop/result/2/merged_methylation_data.csv")
rownames(chr_methy) <- chr_methy[,1]
Region <- strsplit(chr_methy$Region_ID,"_")
temp <- matrix(ncol = 3,nrow = nrow(chr_methy))
for (i in 1:nrow(chr_methy)) {
  if(length(Region[[i]])>3){temp[i,] <- NA}
  else{temp[i,] <- Region[[i]]}
}
temp <- na.omit(temp)
colnames(temp) <- c("Chrom","Start","End")
temp <- as.data.frame(temp)
temp$Start <- as.numeric(temp$Start)
temp$End <- as.numeric(temp$End)
# chrom <- gsub("chr","",temp$Chrom)
# temp$Chrom <- chrom
# temp$Chrom <- as.numeric(temp$Chrom)
arrange(temp,temp[,1],temp[,2])
link <- paste0(temp$Chrom,"_",temp$Start,"_",temp$End)
myheatmap <- chr_methy[link,]
gdata <- myheatmap[,-1]
gdata <- as.matrix(gdata)
gdata <- t(gdata)
cdata <- colnames(myheatmap[,-1])
cdata <- as.data.frame(cdata)
cdata$group <- c(rep("EAC_T",5),rep("ESCC_T",21),rep("ESCC_N",5),rep("EAC_T",7),rep("EAC_N",7))
cdata$state <- c(rep(1,26),rep(0,5),rep(1,7),rep(0,7))
rownames(cdata) <- cdata$cdata
# cdata <- cdata[,-1]
cdata <- cdata[order(cdata$state),]
gdata <- gdata[cdata$cdata,]
annotation_row <- cdata[,c(2,3)]
temp <- colnames(gdata)
annotation_col <- as.data.frame(temp)
chrom_anno <- c(rep("chr1",length(temp[grep("chr1_",temp)])),
                rep("chr2",length(temp[grep("chr2_",temp)])),
                rep("chr3",length(temp[grep("chr3_",temp)])),
                rep("chr4",length(temp[grep("chr4_",temp)])),
                rep("chr5",length(temp[grep("chr5_",temp)])),
                rep("chr6",length(temp[grep("chr6_",temp)])),
                rep("chr7",length(temp[grep("chr7_",temp)])),
                rep("chr8",length(temp[grep("chr8_",temp)])),
                rep("chr9",length(temp[grep("chr9_",temp)])),
                rep("chr10",length(temp[grep("chr10_",temp)])),
                rep("chr11",length(temp[grep("chr11_",temp)])),
                rep("chr12",length(temp[grep("chr12_",temp)])),
                rep("chr13",length(temp[grep("chr13_",temp)])),
                rep("chr14",length(temp[grep("chr14_",temp)])),
                rep("chr15",length(temp[grep("chr15_",temp)])),
                rep("chr16",length(temp[grep("chr16_",temp)])),
                rep("chr17",length(temp[grep("chr17_",temp)])),
                rep("chr18",length(temp[grep("chr18_",temp)])),
                rep("chr19",length(temp[grep("chr19_",temp)])),
                rep("chr20",length(temp[grep("chr20_",temp)])),
                rep("chr21",length(temp[grep("chr21_",temp)])),
                rep("chr22",length(temp[grep("chr22_",temp)])),
                rep("chrX",length(temp[grep("chrX_",temp)])),
                rep("chrY",length(temp[grep("chrY_",temp)])))
annotation_col$chrom <- chrom_anno
mycolors <- c("
mycolors <- c("
mycolors <- c("
mycolors <- c("
mycolors <- c("
color_bar <- c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A")
pdf("chrY_heatmap_10kb_methylation.pdf")
Heatmap(gdata[,grep("chrY_",colnames(gdata))],show_column_names = F,clustering_method_rows = "complete",
        cluster_rows = F,cluster_columns = F,col = color_bar)
dev.off()
# color_bar <- brewer.pal(11, 'BrBG')
setwd("C:/Users/MyXu/Desktop/result/2")
for (i in 1:24) {
  if(i %in% c(1:22)){
    pdf(paste0("chr",i,"heatmap_10kb_methylation.pdf"))
    Heatmap(gdata[,grep(paste0("chr",i,"_"),colnames(gdata))],show_column_names = F,clustering_method_rows = "complete",
            cluster_rows = F,cluster_columns = F,col = color_bar)
    dev.off()
    }
  else if(i==23){
    pdf("chrX_heatmap.pdf")
    Heatmap(gdata[,grep("chrX_",colnames(gdata))],show_column_names = F,clustering_method_rows = "complete",
            cluster_rows = F,cluster_columns = F,col = color_bar)
    dev.off()
  }
  else{
    pdf("chrY_heatmap.pdf")
    Heatmap(gdata[,grep("chrY_",colnames(gdata))],show_column_names = F,clustering_method_rows = "complete",
            cluster_rows = F,cluster_columns = F,col = color_bar)
    dev.off()
  }
}
library(RColorBrewer)
library(knitr)
brewer.pal(11, "BrBG")
# pdf("heatmap.pdf")
# pheatmap(gdata[,grep("chr13_",colnames(gdata))],
#          # breaks = methyBreaksList,
#          show_colnames = F,
#          show_rownames = T,
#          cluster_rows = F,
#          cluster_cols = F,
#          # annotation_row=annotation_row,
#          # annotation_colors = ann_colors,
#          color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))
#          # (n = length(methyBreaksList))
#          )
# dev.off()
region_select <- data[c(41626:41761),]
rownames(region_select) <- region_select$Region_ID
region_select <- region_select[,-1]
plotdata <- t(region_select)
plotdata <- plotdata[c(1:5,32:38,6:26,27:31,39:45),]
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(9, "Blues"))(50)
pdf("D:/MDD_result/8/genome_paired_heatmap.pdf",width = 20,height = 10)
heatmap(scale(plotdata),Colv = NA, Rowv = NA, scale="column",col = coul)
dev.off()
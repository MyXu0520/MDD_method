# MDD_covered_gene
load("C:/Users/MyXu/Desktop/result/4/cg_TCGA_select.RData")
# hg38_annto <- read.table("D:/EdgeDownload/gencode.v47.annotation (1).gtf/gencode.v47.annotation.gtf",sep = "\t")
library("rtracklayer")
gtf_data = import('D:/EdgeDownload/gencode.v47.annotation (1).gtf/gencode.v47.annotation.gtf')
gtf_data = as.data.frame(gtf_data)
gene_anno <- gtf_data[which(gtf_data$type == "gene"),]
gene_anno <- gene_anno[,c(1,2,3,10)]
write.csv(gene_anno,"gene_anno.csv",col.names = F)
##FWQ,python gene_anno.py /student1/MDD/result/1
EAC_SPEC_gene_anno <- read.csv("D:/MDD_result/5/EAC_SPEC_covered_genes.csv")
ESCC_SPEC_gene_anno <- read.csv("D:/MDD_result/5/ESCC_SPEC_covered_genes.csv")
shared_MDD_gene_anno <- read.csv("D:/MDD_result/5/ESCA_sharedMDDs_covered_genes.csv")
gene <- strsplit(EAC_SPEC_gene_anno$Gene_ID,"\\.")
gene_ID <- c()
for (i in 1:nrow(EAC_SPEC_gene_anno)) {
  gene_ID <- c(gene_ID,gene[[i]][1])
}
EAC_SPEC_gene_anno$gene_ID <- gene_ID
write.csv(EAC_SPEC_gene_anno,"D:/MDD_result/5/EAC_SPEC_covered_genes.csv")
gene <- strsplit(ESCC_SPEC_gene_anno$Gene_ID,"\\.")
gene_ID <- c()
for (i in 1:nrow(ESCC_SPEC_gene_anno)) {
  gene_ID <- c(gene_ID,gene[[i]][1])
}
ESCC_SPEC_gene_anno$gene_ID <- gene_ID
write.csv(ESCC_SPEC_gene_anno,"D:/MDD_result/5/ESCC_SPEC_covered_genes.csv")
gene <- strsplit(shared_MDD_gene_anno$Gene_ID,"\\.")
gene_ID <- c()
for (i in 1:nrow(shared_MDD_gene_anno)) {
  gene_ID <- c(gene_ID,gene[[i]][1])
}
shared_MDD_gene_anno$gene_ID <- gene_ID
write.csv(shared_MDD_gene_anno,"D:/MDD_result/5/ESCA_sharedMDDs_covered_genes.csv")
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(RColorBrewer)
EAC_SPEC_gene <- strsplit(EAC_SPEC_gene_anno$Gene_ID,"\\.")
EAC_gene <- c()
for (i in 1:length(EAC_SPEC_gene)) {
  EAC_gene[i] <- EAC_SPEC_gene[[i]][1]
}
EAC_gene_df <- bitr(EAC_gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
EAC_SPEC_GO_enrich = enrichGO(gene = EAC_gene,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTRIZED",
                     ont = "ALL",
                     pvalueCutoff = 1,qvalueCutoff = 1,
                     readable = T)
EAC_SPEC_GO_enrich  = data.frame(EAC_SPEC_GO_enrich)
write.csv(EAC_SPEC_GO_enrich ,'EAC_SPEC_GO_enrich.csv')
EAC_SPEC_GO_enrich$Description <- factor(EAC_SPEC_GO_enrich$Description,levels = rev(EAC_SPEC_GO_enrich$Description))
# mytheme <- theme(axis.title = element_text(size = 13),
#                  axis.text = element_text(size = 11),
#                  legend.title = element_text(size = 13),
#                  legend.text = element_text(size = 11),
#                  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))
# library(cols4all)
# c4a_gui()
#
# p <- ggplot(data = EAC_SPEC_GO_enrich, aes(x = Count, y = Description, fill = -log10(pvalue))) +
#   labs(x = 'Number of Genes', y = '') +
#   theme_bw() + mytheme
# p
# )
library(tidyverse)
library(ggfun)
library(grid)
EAC_SPEC_result <- read_delim(file = "EAC_SPEC_GO_enrich.csv", col_names = T, delim = ",") %>%
  dplyr::select(2,3,4,7,11) %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()
p <- EAC_SPEC_result %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered = T)) %>%
  ggplot(aes(x = Count, y = Description)) +
  geom_bar(aes(fill = pvalue), stat = "identity", width = 0.75) +
  geom_text(aes(x = Count + 10, label = str_c(Count, str_c("(",sprintf("%.2e", pvalue), ")", sep = ""), sep = " "))) +
  ggtitle(label = "GO enrichment") +
  labs(x = "Count", y = "Description") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "#000000"),
    axis.title = element_text(size = 20, color = "#000000"),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.75),
    legend.background = element_roundrect(color = "#969696"),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  coord_cartesian(clip = "off")
p
p +
  annotation_custom(
    grob = rectGrob(
      x = unit(1, "mm"),
      y = unit(1, "mm"),
      width = unit(45, "mm"),
      height = unit(8, "mm"),
      gp = gpar(fill = "#6baed6", col = NA, alpha = 0.5)
    ),
    xmin = unit(-20, "native"),
    xmax = unit(-20, "native"),
    ymin = unit(0.9, "native"),
    ymax = unit(1.0, "native")
  ) +
  annotation_custom(
    grob = rectGrob(
      x = unit(1, "mm"),
      y = unit(1, "mm"),
      width = unit(80, "mm"),
      height = unit(8, "mm"),
      gp = gpar(fill = "#6baed6", col = NA, alpha = 0.5)
    ),
    xmin = unit(-32, "native"),
    xmax = unit(-32, "native"),
    ymin = unit(4.9, "native"),
    ymax = unit(5.0, "native")
  ) +
  annotation_custom(
    grob = rectGrob(
      x = unit(1, "mm"),
      y = unit(1, "mm"),
      width = unit(80, "mm"),
      height = unit(8, "mm"),
      gp = gpar(fill = "#6baed6", col = NA,alpha = 0.5)
    ),
    xmin = unit(-32, "native"),
    xmax = unit(-32, "native"),
    ymin = unit(12.9, "native"),
    ymax = unit(13.0, "native")
  ) +
  annotation_custom(
    grob = rectGrob(
      x = unit(1, "mm"),
      y = unit(1, "mm"),
      width = unit(110, "mm"),
      height = unit(8, "mm"),
      gp = gpar(fill = "#6baed6", col = NA, alpha = 0.5)
    ),
    xmin = unit(-42, "native"),
    xmax = unit(-42, "native"),
    ymin = unit(14.9, "native"),
    ymax = unit(15.0, "native")
  )
ggsave(filename = "EAC_SPEC_GO_Annotation.pdf",
       height = 10,
       width = 15)
EAC_SPEC_KEGG_enrich <- enrichKEGG(gene = EAC_gene_df$ENTREZID,keyType = "kegg",organism= "hsa", qvalueCutoff = 0.05, pvalueCutoff=0.05)
EAC_SPEC_KEGG_enrich <- as.data.frame(EAC_SPEC_KEGG_enrich)
write.csv(EAC_SPEC_KEGG_enrich ,'EAC_SPEC_KEGG_enrich')
shared_MDD_enrich <- read.table("D:/MDD_result/5/ESCA_SPEC_GO_enrich",sep = ",",header = T)
shared_MDD_enrich <- shared_MDD_enrich[c(1:15),c(4,11)]
library(WeightedTreemaps)
library(RColorBrewer)
oct_coord <- list(  x = sin(seq(0, 2, 2/8)*pi) * 1000 + 1000,  y = cos(seq(0, 2, 2/8)*pi) * 1000 + 1000)
GO<- voronoiTreemap( data = shared_MDD_enrich,
                       levels = c("Count", "Description"),
                       cell_size = "Count",
                       shape = oct_coord,
                       seed = 123)
pdf("GO.PDF",width = 5,height = 5)
drawTreemap(GO, label_size = 10,
            label_color = "black",
            color_palette = c(brewer.pal(5,"Set3"),brewer.pal(5,"Pastel1")) ,
            color_level = 2 )
dev.off()
ESCC_SPEC_enrich <- read.table("D:/MDD_result/5/ESCC_SPEC_GO_enrich.csv",sep = ",",header = T)
ESCC_SPEC_enrich <- ESCC_SPEC_enrich[c(1:15),c(4,11)]
library(WeightedTreemaps)
library(RColorBrewer)
oct_coord <- list(  x = sin(seq(0, 2, 2/8)*pi) * 1000 + 1000,  y = cos(seq(0, 2, 2/8)*pi) * 1000 + 1000)
GO<- voronoiTreemap( data = ESCC_SPEC_enrich,
                     levels = c("Count", "Description"),
                     cell_size = "Count",
                     shape = oct_coord,
                     seed = 123)
pdf("GO.PDF",width = 5,height = 5)
drawTreemap(GO, label_size = 8,
            label_color = "black",
            color_palette = c(brewer.pal(5,"Set3"),brewer.pal(5,"Pastel1")) ,
            color_level = 2 )
dev.off()
EAC_SPEC_enrich <- read.table("D:/MDD_result/5/EAC_SPEC_GO_enrich.csv",sep = ",",header = T)
EAC_SPEC_enrich <- EAC_SPEC_enrich[c(1:15),c(4,11)]
library(WeightedTreemaps)
library(RColorBrewer)
oct_coord <- list(  x = sin(seq(0, 2, 2/8)*pi) * 1000 + 1000,  y = cos(seq(0, 2, 2/8)*pi) * 1000 + 1000)
GO<- voronoiTreemap( data = EAC_SPEC_enrich,
                     levels = c("Count", "Description"),
                     cell_size = "Count",
                     shape = oct_coord,
                     seed = 123)
pdf("GO.PDF",width = 5,height = 5)
drawTreemap(GO, label_size = 8,
            label_color = "black",
            color_palette = c(brewer.pal(5,"Set3"),brewer.pal(5,"Pastel1")) ,
            color_level = 2 )
dev.off()
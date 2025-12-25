# setwd("D:/TCGA_ESCA-RNA/")
# library(easyTCGA)
# ESCA_miRNA<-getmirnaexpr("TCGA-ESCA")
# ESCA_cnv<-getcnv("TCGA-ESCA")
# ESCA_snv<-getsnvmaf("TCGA-ESCA")
#
ESCA_count_matrix <- read.csv("D:/MDD_result/6/COUNT_matrix.csv")
rownames(ESCA_count_matrix) <- ESCA_count_matrix$X
rm()
grep <- grep("^TCGA[.-]([0-9a-zA-Z]{2})[.-]([0-9a-zA-Z]{4})[.-]([0][0-9])",colnames(ESCA_count_matrix))
grep
length(grep)
grep1 <- grep("^TCGA[.-]([0-9a-zA-Z]{2})[.-]([0-9a-zA-Z]{4})[.-]([1][0-9])",colnames(ESCA_count_matrix))
grep1
length(grep1)
tumor <- ESCA_count_matrix[,grep]
normal <- ESCA_count_matrix[,grep1]
library(edgeR)
dim(tumor)
dim(normal)
overall <- cbind(tumor,normal)
dim(overall)
overall[1:5,1:4]
dim(tumor)[2]
dim(normal)[2]
group <- c(rep("Tumor",dim(tumor)[2]),rep("Normal",dim(normal)[2]))
group
o = DGEList(overall,group = group)
o
o <- calcNormFactors(o,method = "TMM")
head(o)
library(limma)
o <- edgeR::cpm(o,log = T)
head(o)
dir()
#log2(CPM+1)
write.csv(o,"ESCA_CPM_ACC.csv")
data <- o
rm(o)
new_name <- c()
temp <- strsplit(rownames(data),"\\.")
for (i in 1:nrow(data)) {
  new_name[i] <- temp[[i]][1]
}
rownames(data) <- new_name
library(readr)
clinical_EAC <- read_tsv("D:/MDD_result/4/EAC/clinical.cart.2024-11-20/clinical.tsv")
clinical_ESCC <- read_tsv("D:/MDD_result/4/ESCC/clinical.cart.2024-11-20/clinical.tsv")
EAC_labels <- clinical_EAC$case_submitter_id
ESCC_labels <- clinical_ESCC$case_submitter_id
colnames(data) <- gsub("\\.","-",colnames(data))
EAC_matrix <- data[,which(substr(colnames(data),1,16) %in% paste0(EAC_labels,"-01A"))]
ESCC_matrix <- data[,which(substr(colnames(data),1,16) %in% paste0(ESCC_labels,"-01A"))]
EAC_SPEC_gene_anno <- read.csv("D:/MDD_result/5/EAC_SPEC_covered_genes.csv")
ESCC_SPEC_gene_anno <- read.csv("D:/MDD_result/5/ESCC_SPEC_covered_genes.csv")
shared_MDD_gene_anno <- read.csv("D:/MDD_result/5/ESCA_sharedMDDs_covered_genes.csv")
# #EAC
# EAC_unique_ID <- paste0(EAC_SPEC_gene_anno$BED_chromosome,"_",EAC_SPEC_gene_anno$BED_start,"_",EAC_SPEC_gene_anno$BED_end)
# EAC_SPEC_gene_anno$unique_ID <- EAC_unique_ID
# EAC_Region_exp <- unique(EAC_SPEC_gene_anno$unique_ID)
# EAC_Region_exp <- as.data.frame(EAC_Region_exp)
#
# temp_exp <- c()
# for (i in 1:length(unique(EAC_SPEC_gene_anno$unique_ID))) {
#   temp <- EAC_SPEC_gene_anno[which(EAC_SPEC_gene_anno$unique_ID == unique(EAC_SPEC_gene_anno$unique_ID)[i]),]
#   gene_select <- c()
#
#   for (j in 1:nrow(temp)) {
#     gene_select[j] <- strsplit(temp$Gene_ID,"\\.")[[j]][1]
#   }
#   EAC_gene_exp <- EAC_matrix[which(rownames(EAC_matrix) %in% gene_select),]
#   temp_exp[i] <- mean(EAC_gene_exp)
# }
#
# EAC_Region_exp$EAC_expression <- temp_exp
#
# #ESCC
# temp_exp <- c()
# for (i in 1:length(unique(EAC_SPEC_gene_anno$unique_ID))) {
#   temp <- EAC_SPEC_gene_anno[which(EAC_SPEC_gene_anno$unique_ID == unique(EAC_SPEC_gene_anno$unique_ID)[i]),]
#   gene_select <- c()
#
#   for (j in 1:nrow(temp)) {
#     gene_select[j] <- strsplit(temp$Gene_ID,"\\.")[[j]][1]
#   }
#   ESCC_gene_exp <- ESCC_matrix[which(rownames(ESCC_matrix) %in% gene_select),]
#   temp_exp[i] <- mean(ESCC_gene_exp)
# }
# EAC_Region_exp$ESCC_expression <- temp_exp
#
#
#
# #ESCC
# #EAC
# ESCC_unique_ID <- paste0(ESCC_SPEC_gene_anno$BED_chromosome,"_",ESCC_SPEC_gene_anno$BED_start,"_",ESCC_SPEC_gene_anno$BED_end)
# ESCC_SPEC_gene_anno$unique_ID <- ESCC_unique_ID
# ESCC_Region_exp <- unique(ESCC_SPEC_gene_anno$unique_ID)
# ESCC_Region_exp <- as.data.frame(ESCC_Region_exp)
#
# temp_exp <- c()
# for (i in 1:length(unique(ESCC_SPEC_gene_anno$unique_ID))) {
#   temp <- ESCC_SPEC_gene_anno[which(ESCC_SPEC_gene_anno$unique_ID == unique(ESCC_SPEC_gene_anno$unique_ID)[i]),]
#   gene_select <- c()
#
#   for (j in 1:nrow(temp)) {
#     gene_select[j] <- strsplit(temp$Gene_ID,"\\.")[[j]][1]
#   }
#   EAC_gene_exp <- EAC_matrix[which(rownames(EAC_matrix) %in% gene_select),]
#   temp_exp[i] <- mean(EAC_gene_exp)
# }
#
# ESCC_Region_exp$EAC_expression <- temp_exp
#
# #ESCC
# temp_exp <- c()
# for (i in 1:length(unique(ESCC_SPEC_gene_anno$unique_ID))) {
#   temp <- ESCC_SPEC_gene_anno[which(ESCC_SPEC_gene_anno$unique_ID == unique(ESCC_SPEC_gene_anno$unique_ID)[i]),]
#   gene_select <- c()
#
#   for (j in 1:nrow(temp)) {
#     gene_select[j] <- strsplit(temp$Gene_ID,"\\.")[[j]][1]
#   }
#   ESCC_gene_exp <- ESCC_matrix[which(rownames(ESCC_matrix) %in% gene_select),]
#   temp_exp[i] <- mean(ESCC_gene_exp)
# }
#
# ESCC_Region_exp$ESCC_expression <- temp_exp
#
#
#
# #Shared
# # EAC
# shared_unique_ID <- paste0(shared_MDD_gene_anno$BED_chromosome,"_",shared_MDD_gene_anno$BED_start,"_",shared_MDD_gene_anno$BED_end)
# shared_MDD_gene_anno$unique_ID <- shared_unique_ID
# shared_Region_exp <- unique(shared_MDD_gene_anno$unique_ID)
# shared_Region_exp <- as.data.frame(shared_Region_exp)
#
# temp_exp <- c()
# for (i in 1:length(unique(shared_MDD_gene_anno$unique_ID))) {
#   temp <- shared_MDD_gene_anno[which(shared_MDD_gene_anno$unique_ID == unique(shared_MDD_gene_anno$unique_ID)[i]),]
#   gene_select <- c()
#
#   for (j in 1:nrow(temp)) {
#     gene_select[j] <- strsplit(temp$Gene_ID,"\\.")[[j]][1]
#   }
#   EAC_gene_exp <- EAC_matrix[which(rownames(EAC_matrix) %in% gene_select),]
#   temp_exp[i] <- mean(EAC_gene_exp)
# }
#
# shared_Region_exp$EAC_expression <- temp_exp
#
# #ESCC
# temp_exp <- c()
# for (i in 1:length(unique(shared_MDD_gene_anno$unique_ID))) {
#   temp <- shared_MDD_gene_anno[which(shared_MDD_gene_anno$unique_ID == unique(shared_MDD_gene_anno$unique_ID)[i]),]
#   gene_select <- c()
#
#   for (j in 1:nrow(temp)) {
#     gene_select[j] <- strsplit(temp$Gene_ID,"\\.")[[j]][1]
#   }
#   ESCC_gene_exp <- ESCC_matrix[which(rownames(ESCC_matrix) %in% gene_select),]
#   temp_exp[i] <- mean(ESCC_gene_exp)
# }
#
# shared_Region_exp$ESCC_expression <- temp_exp
# save(EAC_Region_exp,ESCC_Region_exp,shared_Region_exp,file = "Region_exp.RData")
#
gene_ID <- c()
temp <- strsplit(EAC_SPEC_gene_anno$Gene_ID,"\\.")
for (i in 1:nrow(EAC_SPEC_gene_anno)) {
  gene_ID[i] <- temp[[i]][1]
}
EAC_SPEC_gene_anno$Gene_ID <- gene_ID
gene_ID <- c()
temp <- strsplit(ESCC_SPEC_gene_anno$Gene_ID,"\\.")
for (i in 1:nrow(ESCC_SPEC_gene_anno)) {
  gene_ID[i] <- temp[[i]][1]
}
ESCC_SPEC_gene_anno$Gene_ID <- gene_ID
gene_ID <- c()
temp <- strsplit(shared_MDD_gene_anno$Gene_ID,"\\.")
for (i in 1:nrow(shared_MDD_gene_anno)) {
  gene_ID[i] <- temp[[i]][1]
}
shared_MDD_gene_anno$Gene_ID <- gene_ID
plotData_EAC <- colnames(EAC_matrix)
plotData_EAC <- as.data.frame(plotData_EAC)
colnames(plotData_EAC) <- "Sample"
EAC_SPEC_exp <- c()
ESCC_SPEC_exp <- c()
shared_exp <- c()
for(i in 1:ncol(EAC_matrix)){
  EAC_SPEC_exp[i] <- mean(EAC_matrix[which(rownames(EAC_matrix) %in% EAC_SPEC_gene_anno$Gene_ID),i])
}
plotData_EAC$EXP <- EAC_SPEC_exp
for(i in 1:ncol(EAC_matrix)){
  ESCC_SPEC_exp[i] <- mean(EAC_matrix[which(rownames(EAC_matrix) %in% ESCC_SPEC_gene_anno$Gene_ID),i])
}
temp <- as.data.frame(plotData_EAC$Sample)
colnames(temp) <- "Sample"
temp$EXP <- ESCC_SPEC_exp
plotData_EAC <- rbind(plotData_EAC,temp)
for(i in 1:ncol(EAC_matrix)){
  shared_exp[i] <- mean(EAC_matrix[which(rownames(EAC_matrix) %in% shared_MDD_gene_anno$Gene_ID),i])
}
temp <- as.data.frame(colnames(EAC_matrix))
colnames(temp) <- "Sample"
temp$EXP <- shared_exp
plotData_EAC <- rbind(plotData_EAC,temp)
plotData_EAC$Type <- c(rep("EAC_SPEC",ncol(EAC_matrix)),rep("ESCC_SPEC",ncol(EAC_matrix)),rep("shared_MDD",ncol(EAC_matrix)))
plotData_EAC$Class <- rep("EAC",nrow(plotData_EAC))
plotData_ESCC <- colnames(ESCC_matrix)
plotData_ESCC <- as.data.frame(plotData_ESCC)
colnames(plotData_ESCC) <- "Sample"
EAC_SPEC_exp <- c()
ESCC_SPEC_exp <- c()
shared_exp <- c()
for(i in 1:ncol(ESCC_matrix)){
  EAC_SPEC_exp[i] <- mean(ESCC_matrix[which(rownames(ESCC_matrix) %in% EAC_SPEC_gene_anno$Gene_ID),i])
}
plotData_ESCC$EXP <- EAC_SPEC_exp
for(i in 1:ncol(ESCC_matrix)){
  ESCC_SPEC_exp[i] <- mean(ESCC_matrix[which(rownames(ESCC_matrix) %in% ESCC_SPEC_gene_anno$Gene_ID),i])
}
temp <- as.data.frame(plotData_ESCC$Sample)
colnames(temp) <- "Sample"
temp$EXP <- ESCC_SPEC_exp
plotData_ESCC <- rbind(plotData_ESCC,temp)
for(i in 1:ncol(ESCC_matrix)){
  shared_exp[i] <- mean(ESCC_matrix[which(rownames(ESCC_matrix) %in% shared_MDD_gene_anno$Gene_ID),i])
}
temp <- as.data.frame(colnames(ESCC_matrix))
colnames(temp) <- "Sample"
temp$EXP <- shared_exp
plotData_ESCC <- rbind(plotData_ESCC,temp)
plotData_ESCC$Type <- c(rep("EAC_SPEC",ncol(ESCC_matrix)),rep("ESCC_SPEC",ncol(ESCC_matrix)),rep("shared_MDD",ncol(ESCC_matrix)))
plotData_ESCC$Class <- rep("ESCC",nrow(plotData_ESCC))
plotData <- rbind(plotData_EAC,plotData_ESCC)
save(plotData,file = "boxplotdata_ESCA_exp.RData")
plotData_sorted <- plotData[order(plotData$Sample),]
library(dplyr)
df <- plotData_sorted %>%
  filter(Type %in% c("EAC_SPEC","ESCC_SPEC","shared_MDD")) %>%
  select(Sample,EXP,Class,Type)%>%
  mutate(paired = rep(1:(nrow(plotData_sorted)/3),each=3),Type=factor(Type))
head(df)
library(rstatix)
df %>% group_by(Class, Type) %>% shapiro_test(EXP)
setwd("D:/MDD_result/6/")
library(tidyverse)
library(gaCMInder)
library(ggsci)
library(ggprism)
library(rstatix)
library(ggpubr)
pdf("qplot_all.pdf")
ggqqplot(df, "EXP", ggtheme = theme_bw()) +
  facet_grid(Type ~ Class, labeller = "label_both")
dev.off()
df_p_val1 <- df %>% group_by(Class) %>%
  t_test(EXP  ~Type, paired = TRUE) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj")
df_p_val1
df_p_val1 <- df_p_val1 %>%
  add_xy_position(x = "Type", dodge = 0.8)
pdf("ESCA_region_exp.pdf")
ggplot(df,aes(Type,EXP)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1) +
  geom_boxplot(position=position_dodge(width =0.2),width=0.4) +
  geom_line(aes(group=paired),position = position_dodge(0.2),color="grey80") +
  geom_point(aes(fill= Type,group=paired,size=EXP/2,alpha=EXP),pch=21,
             position = position_dodge(0.2)) +
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=6,hide.ns = T) +
  scale_size_continuous(range=c(1,3)) +
  facet_wrap(.~Class,nrow=1) +
  scale_fill_manual(values = c("#AA2D36","#2F56A3","#64329E"))+
  scale_x_discrete(guide = "prism_bracket") +
  scale_y_continuous(limits = c(-5,2),minor_breaks = seq(-5,2,0.25),guide = "prism_offset_minor") +
  labs(x=NULL,y=NULL) +
  theme_prism(base_line_size =0.5) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(margin = margin(t = -5),color="black",size=10),
        legend.position = "none",
        panel.spacing = unit(0,"lines")) +
  coord_cartesian()
dev.off()
save.image("C:/Users/MyXu/Desktop/result/4/ESCA_paired_boxplot_all.RData")
pdf("ESCA_region_exp2.pdf")
ggplot(df,aes(Type,EXP)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1) +
  geom_boxplot(position=position_dodge(width =0.2),width=0.4) +
  geom_point(aes(fill= Type,size=EXP/2,alpha=EXP),pch=21,
             position = position_dodge(0.2)) +
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=6,hide.ns = T) +
  scale_size_continuous(range=c(1,3)) +
  facet_wrap(.~Class,nrow=1) +
  scale_fill_manual(values = c("#AA2D36","#2F56A3","#64329E"))+
  scale_x_discrete(guide = "prism_bracket") +
  scale_y_continuous(limits = c(-5,2),minor_breaks = seq(-5,2,0.25),guide = "prism_offset_minor") +
  labs(x=NULL,y=NULL) +
  theme_prism(base_line_size =0.5) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(margin = margin(t = -5),color="black",size=10),
        legend.position = "none",
        panel.spacing = unit(0,"lines")) +
  coord_cartesian()
dev.off()
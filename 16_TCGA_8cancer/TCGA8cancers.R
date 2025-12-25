# ESCA/STAD/COAD/READ/LIHC/CHOL/PAAD/HNSC
setwd("E:/33cancer/33cancers/knn_data")
cancer <- c("ESCA","STAD","COAD","READ","LIHC","CHOL","PAAD","HNSC")
library(tsne)
library(pheatmap)
library(dplyr)
library(tidyverse)
# install.packages("Rtsne")
library(Rtsne)
library(data.table)
file_list <- list.files(pattern = "\\.txt$")
temp <- c()
for(i in 1:length(cancer)) {
  temp_list <- file_list[grep(cancer[i],file_list)]
  temp <- c(temp,temp_list)
}
file_list <- temp
all_data <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  # data <- data.table::fread(file,header = F,)
  all_data[[file]] <- data
  print(file)
}
# all_data[[20]] <- read.table(file_list[20], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
save(all_data,file = "E:/33cancer/33cancers/Gastrointestinal_cancer.RData")
load("E:/Gastrointestinal_cancer.RData")
# common_probes <- Reduce(intersect, lapply(all_data, rownames))
common_probes <- rownames(all_data$ESCA.methylation450.knn.txt)
for (i in 1:length(all_data)) {
  common_probes <- intersect(common_probes,rownames(all_data[[i]]))
  if(length(common_probes) == 0){
    print(i)
    break}
}
merged_data <- do.call(cbind, lapply(all_data, function(x) x[common_probes, , drop = FALSE]))
colnames_merge <- c()
for (i in 1:length(all_data)) {
  colnames_merge <- c(colnames_merge,colnames(all_data[[i]]))
}
colnames(merged_data) <- colnames_merge
set.seed(123)
# gsub("\\.","_",colnames(merged_data))
merged_data <- t(merged_data)
tsne_result <- Rtsne(merged_data, dims = 2, perplexity = 30, verbose = FALSE,theta = 0.0,max_iter = 1000)
sample <- c()
group <- c()
for(i in 1:length(all_data)){
  sample <- c(sample,colnames(all_data[[i]]))
  sample_name <- names(all_data[i])
  sample_name <- strsplit(sample_name,"\\.")
  sample_name <- sample_name[[1]][1]
  group <- c(group,rep(sample_name,ncol(all_data[[i]])))
}
sample_df <- as.data.frame(sample)
sample_df$group <- group
save(tsne_result,file = "D:/MDD_result/18/TCGA_Gastrointestinal_cancer_TSNE.RData")
tsne_out = as.data.frame(tsne_result$Y)
colnames(tsne_out) <- c("TSNE_1","TSNE_2")
# ggplot(tsne_out,aes(tSNE1,tSNE2,colour = group))+
#   geom_point()+
install.packages('plotly')
remotes::install_github("plotly/plotly")
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
p <- ggplot(tsne_out, aes(x = TSNE_1, y = TSNE_2)) +
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
pdf("D:/MDD_result/18/TCGA_Gastrointestinal_cancer_TSNE.pdf")
p <- ggplot(tsne_out, aes(x = TSNE_1, y = TSNE_2)) +
  # stat_ellipse(aes(fill = group), level = 0.95, linetype = 1,
  #              show.legend = F, geom = 'polygon', alpha = 0.1) +
  geom_point(aes(colour  =group), size =1.0, alpha = 1) +
  mytheme +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) +
  theme(panel.grid = element_blank())+
# geom_text(data = group, aes(x = TSNE_1, y = TSNE_2, label = group), fontface = "bold", color = 'black', size = 4) +
# guides(color = guide_legend(override.aes = list(size = 5))) +
# scale_fill_manual(values = mycol) +
# scale_fill_brewer(palette = "Spectral")+
  scale_color_brewer(palette = "Spectral")+
  scale_color_manual(values = mycol)
# scale_color_gradientn(values = seq(0,1,0.2),
# colours = c('cyan','blue','green','orange','red'))
p
dev.off()
load("E:/33cancer/33cancers/Gastrointestinal_cancer.RData")
# common_probes <- Reduce(intersect, lapply(all_data, rownames))
common_probes <- rownames(all_data$ESCA.methylation450.knn.txt)
for (i in 1:length(all_data)) {
  common_probes <- intersect(common_probes,rownames(all_data[[i]]))
  if(length(common_probes) == 0){
    print(i)
    break}
}
cg_450K <- read.csv("E:/humanmethylation450.csv",sep = ",",header = T)
cg_450K_locus <- cg_450K[,c(1,12,13,18,22,24,25,26)]
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$Probe_SNPs == ""),]
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$CHR %in% c(as.character(1:22),"X","Y")),]
cg_450K_locus$CHR <- paste0("chr",cg_450K_locus$CHR)
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$IlmnID %in% common_probes),]
#core_MDD
core_MDD <- read.table("D:/MDD_result/12/core_MDD.bed")
colnames(core_MDD) <- c("Chrom","Start","End")
core_MDD_list <- list()
for (i in 1:nrow(core_MDD)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == core_MDD[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= core_MDD[i,2] && temp[j,3] <= core_MDD[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  core_MDD_list[[i]] <- cg_list
}
core_MDD_cg <- c()
for (i in 1:length(core_MDD_list)) {
  temp_cg <- core_MDD_list[[i]]
  core_MDD_cg <- c(core_MDD_cg,temp_cg)
}
sample_methylation_MDD <- list()
for (i in 1:length(all_data)) {
  temp_sample_methylation <- all_data[[i]][core_MDD_cg,]
  sample_methylation_MDD[[i]] <- apply(temp_sample_methylation,2,mean)
}
# sample_methylation_MDD <- list()
# for (i in 1:length(all_data)) {
#   if(names(all_data)[i] == "ESCA.methylation450.knn.txt"){
#     temp_sample_methylation <- all_data[[i]][core_MDD_cg,]
#     sample_methylation_MDD[[i]] <- apply(temp_sample_methylation,2,mean)
#   }
#   else{sample_methylation_MDD[[i]] <- rep(0,ncol(all_data[[i]]))}
#
# }
MDD_anno <- c()
for (i in 1:length(sample_methylation_MDD)) {
  temp_MDD_anno <- sample_methylation_MDD[[i]]
  MDD_anno <- c(MDD_anno,temp_MDD_anno)
}
#ESCA_DMR
hypo_DMR <- dmrs[which(dmrs$diff.Methy <= 0),c(1:3)]
colnames(hypo_DMR) <- c("Chrom","Start","End")
hypo_DMR_list <- list()
for (i in 1:nrow(hypo_DMR)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == hypo_DMR[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= hypo_DMR[i,2] && temp[j,3] <= hypo_DMR[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  hypo_DMR_list[[i]] <- cg_list
}
hypo_DMR_cg <- c()
for (i in 1:length(hypo_DMR_list)) {
  temp_cg <- hypo_DMR_list[[i]]
  hypo_DMR_cg <- c(hypo_DMR_cg,temp_cg)
}
sample_methylation_DMR <- list()
for (i in 1:length(all_data)) {
  temp_sample_methylation <- all_data[[i]][hypo_DMR_cg,]
  sample_methylation_DMR[[i]] <- apply(temp_sample_methylation,2,mean)
}
# sample_methylation_DMR <- list()
# for (i in 1:length(all_data)) {
#   if(names(all_data)[i] == "ESCA.methylation450.knn.txt"){
#     temp_sample_methylation <- all_data[[i]][hypo_DMR_cg,]
#     sample_methylation_DMR[[i]] <- apply(temp_sample_methylation,2,mean)
#   }
#   else{sample_methylation_DMR[[i]] <- rep(0,ncol(all_data[[i]]))}
#
# }
DMR_anno <- c()
for (i in 1:length(sample_methylation_DMR)) {
  temp_DMR_anno <- sample_methylation_DMR[[i]]
  DMR_anno <- c(DMR_anno,temp_DMR_anno)
}
load("D:/MDD_result/18/TCGA_Gastrointestinal_cancer_TSNE.RData")
tsne_out = as.data.frame(tsne_result$Y)
colnames(tsne_out) <- c("TSNE_1","TSNE_2")
mycol <- c4a('10',8)
mytheme <- theme_void() + theme(plot.margin = margin(5.5,15,5.5,5.5))
pdf("D:/MDD_result/18/TCGA_Gastrointestinal_cancer_TSNE_ESCASMDD.pdf")
p <- ggplot(tsne_out, aes(x = TSNE_1, y = TSNE_2)) +
  # stat_ellipse(aes(fill = group), level = 0.95, linetype = 1,
  #              show.legend = F, geom = 'polygon', alpha = 0.1) +
  geom_point(aes(colour=SMDD_anno), size =1.0, alpha = 1) +
  mytheme +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) +
  theme(panel.grid = element_blank())+
  # geom_text(data = group, aes(x = TSNE_1, y = TSNE_2, label = group), fontface = "bold", color = 'black', size = 4) +
  # guides(color = guide_legend(override.aes = list(size = 5))) +
  # scale_fill_manual(values = mycol) +
  # scale_fill_brewer(palette = "Spectral")+
  # scale_color_brewer(palette = "Spectral")+
  # scale_color_manual(values = mycol)
  scale_color_gradientn(values = seq(0,1,0.2),
  # colours = c("gray",'blue','cyan','green','orange','red'))
  colours = c("gray",'blue','yellow','red'))
p
dev.off()
#shared_MDD
shared_MDD <- read.table("D:/MDD_result/1/ESCA_sharedMDDs.bed")
colnames(shared_MDD) <- c("Chrom","Start","End")
shared_MDD_list <- list()
for (i in 1:nrow(shared_MDD)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == shared_MDD[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= shared_MDD[i,2] && temp[j,3] <= shared_MDD[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  shared_MDD_list[[i]] <- cg_list
}
shared_MDD_cg <- c()
for (i in 1:length(shared_MDD_list)) {
  temp_cg <- shared_MDD_list[[i]]
  shared_MDD_cg <- c(shared_MDD_cg,temp_cg)
}
sample_methylation_SMDD <- list()
for (i in 1:length(all_data)) {
  temp_sample_methylation <- all_data[[i]][shared_MDD_cg,]
  sample_methylation_SMDD[[i]] <- apply(temp_sample_methylation,2,mean)
}
# sample_methylation_SMDD <- list()
# for (i in 1:length(all_data)) {
#   if(names(all_data)[i] == "ESCA.methylation450.knn.txt"){
#     temp_sample_methylation <- all_data[[i]][shared_MDD_cg,]
#     sample_methylation_SMDD[[i]] <- apply(temp_sample_methylation,2,mean)
#   }
#   else{sample_methylation_SMDD[[i]] <- rep(0,ncol(all_data[[i]]))}
#
# }
SMDD_anno <- c()
for (i in 1:length(sample_methylation_SMDD)) {
  temp_MDD_anno <- sample_methylation_SMDD[[i]]
  SMDD_anno <- c(SMDD_anno,temp_MDD_anno)
}
#shared_MDD cluster
TCGA8_sample_SMDD <- matrix(nrow = length(shared_MDD_cg))
for (i in 1:length(all_data)) {
  temp_sample <- all_data[[i]][shared_MDD_cg,]
  TCGA8_sample_SMDD <- cbind(TCGA8_sample_SMDD,temp_sample)
}
TCGA8_sample_SMDD <- TCGA8_sample_SMDD[,-1]
set.seed(123)
# gsub("\\.","_",colnames(merged_data))
merged_data_SMDD <- t(TCGA8_sample_SMDD)
tsne_result_SMDD <- Rtsne(merged_data_SMDD, dims = 2, perplexity = 30, verbose = FALSE,theta = 0.0,max_iter = 1000)
tsne_out_SMDD = as.data.frame(tsne_result_SMDD$Y)
colnames(tsne_out_SMDD) <- c("TSNE_1","TSNE_2")
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
p <- ggplot(tsne_out_SMDD, aes(x = TSNE_1, y = TSNE_2)) +
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
#hypo_DMR cluster
TCGA8_sample_DMR <- matrix(nrow = length(hypo_DMR_cg))
for (i in 1:length(all_data)) {
  temp_sample <- all_data[[i]][hypo_DMR_cg,]
  TCGA8_sample_DMR <- cbind(TCGA8_sample_DMR,temp_sample)
}
TCGA8_sample_DMR <- TCGA8_sample_DMR[,-1]
set.seed(123)
# gsub("\\.","_",colnames(merged_data))
merged_data_DMR <- t(TCGA8_sample_DMR)
tsne_result_DMR <- Rtsne(merged_data_DMR, dims = 2, perplexity = 30, verbose = FALSE,theta = 0.0,max_iter = 1000)
tsne_out_DMR = as.data.frame(tsne_result_DMR$Y)
colnames(tsne_out_DMR) <- c("TSNE_1","TSNE_2")
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
p <- ggplot(tsne_out_DMR, aes(x = TSNE_1, y = TSNE_2)) +
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
#core_MDD cluster
TCGA8_sample_CMDD <- matrix(nrow = length(core_MDD_cg))
for (i in 1:length(all_data)) {
  temp_sample <- all_data[[i]][core_MDD_cg,]
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
  geom_point(aes(colour  = type), size =1.0, alpha = 1) +
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
SC <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
AD <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
sample_df
# temp <- grepl("\.","-",sample_df$sample)
new_string <- gsub( "\\.", "-",sample_df$sample)
sample_df$new_sample <- new_string
sample_df$new_sample2 <- substr(new_string,1,12)
type <- c()
for (i in 1:nrow(sample_df)) {
  if(sample_df$new_sample2[i] %in% SC$case_submitter_id){type[i] <- "SC"}
  else if(sample_df$new_sample2[i] %in% AD$case_submitter_id){type[i] <- "AD"}
  else{type[i] <- "OT"}
}
sample_df$type <- type
p <- ggplot(tsne_out_DMR, aes(x = TSNE_1, y = TSNE_2)) +
  # stat_ellipse(aes(fill = group), level = 0.95, linetype = 1,
  #              show.legend = F, geom = 'polygon', alpha = 0.1) +
  geom_point(aes(colour  = type), size =1.0, alpha = 1) +
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
load("E:/33cancer/33cancers/Gastrointestinal_cancer.RData")
common_probes <- rownames(all_data$ESCA.methylation450.knn.txt)
for (i in 1:length(all_data)) {
  common_probes <- intersect(common_probes,rownames(all_data[[i]]))
  if(length(common_probes) == 0){
    print(i)
    break}
}
for (i in 1:length(all_data)) {
  all_data[[i]] <- all_data[[i]][common_probes,]
}
#ESCA_shared_MDD
ESCA_shared_MDD <- read.table("D:/MDD_result/1/ESCA_sharedMDDs.bed")
colnames(ESCA_shared_MDD) <- c("Chrom","Start","End")
#EAC_spec_MDD
EAC_spec_MDD <- read.table("D:/MDD_result/1/EAC_specificMDDs.bed")
colnames(EAC_spec_MDD) <- c("Chrom","Start","End")
#ESCC_spec_MDD
ESCC_spec_MDD <- read.table("D:/MDD_result/1/ESCC_specificMDDs.bed")
colnames(ESCC_spec_MDD) <- c("Chrom","Start","End")
load("D:/MDD_result/4/cg_TCGA_select.RData")
cancer <- c("ESCA","STAD","COAD","READ","LIHC","CHOL","PAAD","HNSC")
all_data_shared <- list()
for (z in 1:length(all_data)) {
  shared_avME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    sample_MDD_ME <- 0
    n <- 0
    for (j in 1:length(shared_list)) {
      if(length(shared_list[[j]]) == 0) {
        temp_MDD_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][shared_list[[j]],i])) == 0){
        temp_MDD_ME <- 0
        n=n+1}
      else{temp_MDD_ME <-  mean(na.omit(all_data[[z]][shared_list[[j]],i]))}
      sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
    }
    shared_avME[i] <- sample_MDD_ME/(length(shared_list)-n)
  }
  shared_avME <- as.data.frame(shared_avME)
  colnames(shared_avME) <- "Meth_level"
  shared_avME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  shared_avME$Type <- rep("shared_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  shared_avME <- cbind(temp,shared_avME)
  all_data_shared[[z]] <- shared_avME
}
save(all_data_shared,file = "D:/MDD_result/18/TCGA_8cancer_MDD.RData")
all_data_EAC_SPEC <- list()
for (z in 1:length(all_data)) {
  EAC_SPEC_avME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    sample_MDD_ME <- 0
    n <- 0
    for (j in 1:length(EAC_spec_list)) {
      if(length(EAC_spec_list[[j]]) == 0) {
        temp_MDD_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][EAC_spec_list[[j]],i])) == 0){
        temp_MDD_ME <- 0
        n=n+1}
      else{temp_MDD_ME <-  mean(na.omit(all_data[[z]][EAC_spec_list[[j]],i]))}
      sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
    }
    EAC_SPEC_avME[i] <- sample_MDD_ME/(length(EAC_spec_list)-n)
  }
  EAC_SPEC_avME <- as.data.frame(EAC_SPEC_avME)
  colnames(EAC_SPEC_avME) <- "Meth_level"
  EAC_SPEC_avME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  EAC_SPEC_avME$Type <- rep("EAC_spec_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  EAC_SPEC_avME <- cbind(temp,EAC_SPEC_avME)
  all_data_EAC_SPEC[[z]] <- EAC_SPEC_avME
}
save(all_data_shared,all_data_EAC_SPEC,file = "D:/MDD_result/18/TCGA_8cancer_MDD.RData")
all_data_ESCC_SPEC <- list()
for (z in 1:length(all_data)) {
  ESCC_SPEC_avME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    sample_MDD_ME <- 0
    n <- 0
    for (j in 1:length(ESCC_spec_list)) {
      if(length(ESCC_spec_list[[j]]) == 0) {
        temp_MDD_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][ESCC_spec_list[[j]],i])) == 0){
        temp_MDD_ME <- 0
        n=n+1}
      else{temp_MDD_ME <-  mean(na.omit(all_data[[z]][ESCC_spec_list[[j]],i]))}
      sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
    }
    ESCC_SPEC_avME[i] <- sample_MDD_ME/(length(ESCC_spec_list)-n)
  }
  ESCC_SPEC_avME <- as.data.frame(ESCC_SPEC_avME)
  colnames(ESCC_SPEC_avME) <- "Meth_level"
  ESCC_SPEC_avME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  ESCC_SPEC_avME$Type <- rep("ESCC_spec_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  ESCC_SPEC_avME <- cbind(temp,ESCC_SPEC_avME)
  all_data_ESCC_SPEC[[z]] <- ESCC_SPEC_avME
}
save(all_data_shared,all_data_EAC_SPEC,all_data_ESCC_SPEC,file = "D:/MDD_result/18/TCGA_8cancer_MDD.RData")
SC <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
AD <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
# sample_df
#all_data_shared
for (i in 1:length(all_data_shared)) {
  temp <- all_data_shared[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_data_shared[[i]] <- temp
}
#all_data_EAC_SPEC
for (i in 1:length(all_data_EAC_SPEC)) {
  temp <- all_data_EAC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_data_EAC_SPEC[[i]] <- temp
}
#all_data_ESCC_SPEC
for (i in 1:length(all_data_ESCC_SPEC)) {
  temp <- all_data_ESCC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_data_ESCC_SPEC[[i]] <- temp
}
#EAC-ESCC delta methylation
EAC_VS_ESCC_delta <- list()
for (i in 1:length(all_data_EAC_SPEC)) {
  delta_meth <- all_data_EAC_SPEC[[i]]$Meth_level-all_data_ESCC_SPEC[[i]]$Meth_level
  EAC_VS_ESCC_delta[[i]] <- all_data_EAC_SPEC[[i]]
  EAC_VS_ESCC_delta[[i]]$Meth_level <- delta_meth
  colnames(EAC_VS_ESCC_delta[[i]]) <- c("Sample","delta_methy","Cancer","Type","Cancer_class")
}
save(EAC_VS_ESCC_delta,file = "D:/MDD_result/21/8cancers_delta_meth.RData")
# sample cancerType      delta        label
# 1       1       ESCC 0.12741917 Pan-Squamous
# 2       1       HNSC 0.08870604       Others
# 3       1       LUSC 0.10726257       Pan-GI
# EAC_VS_ESCC_plot_data <- as.data.frame(matrix(ncol = 5))
# colnames(EAC_VS_ESCC_plot_data) <- colnames(EAC_VS_ESCC_delta[[1]])
# for (i in 1:nrow(EAC_VS_ESCC_delta[[1]])) {
#   temp_plot_data <- as.data.frame(matrix(ncol = 5))
#   colnames(temp_plot_data) <- colnames(EAC_VS_ESCC_delta[[1]])
#   for (j in 1:length(EAC_VS_ESCC_delta)) {
#     temp_plot_data <- rbind(temp_plot_data,EAC_VS_ESCC_delta[[j]][i,])
#   }
#   EAC_VS_ESCC_plot_data <- rbind(EAC_VS_ESCC_plot_data,temp_plot_data)
# }
EAC_VS_ESCC_plot_data <- as.data.frame(matrix(ncol = 4))
colnames(EAC_VS_ESCC_plot_data) <- c("Sample","Cancer","delta_methy","Cancer_class")
for (i in 1:length(EAC_VS_ESCC_delta)) {
  temp_plot_data <- EAC_VS_ESCC_delta[[i]][,c(1,3,2,5)]
  temp_plot_data <- temp_plot_data[which(substr(temp_plot_data$Sample,14,16) == "01A"),]
  temp_plot_data$Sample <- c(1:nrow(temp_plot_data))
  EAC_VS_ESCC_plot_data <- rbind(EAC_VS_ESCC_plot_data,temp_plot_data)
}
EAC_VS_ESCC_plot_data <- EAC_VS_ESCC_plot_data[-1,]
library(ggplot2)
colnames(EAC_VS_ESCC_plot_data) <- c("sample", "cancerType", "delta", "label")
plotTumorATACPointPlot(EAC_VS_ESCC_plot_data,"ESCC only MDDs vs EAC only MDDs", "example_plot2.pdf")
# plotTumorATACPointPlot=function(plotdata2, ylabel, saveFile){
#   plotdata2=plotdata2[order(plotdata2$cancerType, plotdata2$delta),]
#   name<-as.character(unique(plotdata2$cancerType))
#   sum=0
#   for(i in 1:length(name)){
#     tumorRow = nrow(plotdata2[plotdata2$cancerType==name[i],])
#     for(j in 1:tumorRow){
#       plotdata2[sum+j,"Postion"]<- (1000*(i-1)+350+250*j/tumorRow)
#     }
#     sum=sum+tumorRow
#   }
#   for(i in 1:length(name)){
#     plotdata2[plotdata2$cancerType==name[i],"Tick"]<-(1000*(i-1)+350)
#   }
#   meanData=plotdata2[1,]
#   sum2=0
#   for(i in 1:length(name)){
#     tumorRow = nrow(plotdata2[plotdata2$cancerType==name[i],])
#     temp=plotdata2[sum2+ceiling(tumorRow/2),]
#     temp$delta=mean(plotdata2[plotdata2$cancerType==name[i],]$delta)
#     meanData <-rbind(meanData,temp)
#     sum2=sum2+tumorRow
#   }
#   meanData$xmin=meanData$Postion-125
#   meanData$xmax=meanData$Postion+125
#   meanData$ymin=meanData$delta
#   meanData$ymax=meanData$delta
#   meanData=meanData[-1,]
#   write.table(plotdata2, gsub(".pdf", ".data.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
#   p = ggplot(plotdata2, aes(x=Postion, y=delta, color=label))+geom_point(size=0.6,stat="identity")+
#     geom_rect(data=meanData,aes(xmin = xmin, xmax = xmax, ymin = delta, ymax = delta),
#               fill="white",size=0.4,color="orange")
#   p=p+theme_classic()+ylab(ylabel)+ggtitle("TCGA_tumor")
#   p=p+scale_color_manual(name="Cancer groups", values=c("#EA3323","black", "#0000F5"))
#   p=p+scale_x_continuous(name = "", breaks = unique(plotdata2$Tick), labels = unique(plotdata2$cancerType))
#   p=p+geom_hline(yintercept=0, linetype="dashed", color = "grey")
#   p=p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text = element_text(color="black", size=12),
#             plot.title = element_text(hjust = 0.5, size=16, face="bold"), legend.position = "bottom",
#             legend.text = element_text(size=12))
#   pdf(saveFile, width=9, height=4)
#   print(p)
#   dev.off()
# }
###############################################################################################################################
# load("E:/33cancer/33cancers/Gastrointestinal_cancer.RData")
# common_probes <- rownames(all_data$ESCA.methylation450.knn.txt)
# for (i in 1:length(all_data)) {
#   common_probes <- intersect(common_probes,rownames(all_data[[i]]))
#   if(length(common_probes) == 0){
#     print(i)
#     break}
# }
# for (i in 1:length(all_data)) {
#   all_data[[i]] <- all_data[[i]][common_probes,]
# }
#ESCA_shared_DMR
ESCA_shared_DMR <- read.table("D:/MDD_result/7/EAC_ESCC_DMR_bed.bed")
colnames(ESCA_shared_DMR) <- c("Chrom","Start","End")
#EAC_spec_DMR
EAC_spec_DMR <- read.table("D:/MDD_result/7/EAC_SPEC_DMR.bed")
colnames(EAC_spec_DMR) <- c("Chrom","Start","End")
#ESCC_spec_DMR
ESCC_spec_DMR <- read.table("D:/MDD_result/7/ESCC_SPEC_DMR.bed")
colnames(ESCC_spec_DMR) <- c("Chrom","Start","End")
load("D:/MDD_result/18/cg_TCGA_select_DMR.RData")
cancer <- c("ESCA","STAD","COAD","READ","LIHC","CHOL","PAAD","HNSC")
all_DMRdata_shared <- list()
for (z in 1:length(all_data)) {
  shared_DMRavME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    sample_DMR_ME <- 0
    n <- 0
    for (j in 1:length(DMR_shared_list)) {
      if(length(DMR_shared_list[[j]]) == 0) {
        temp_DMR_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][DMR_shared_list[[j]],i])) == 0){
        temp_DMR_ME <- 0
        n=n+1}
      else{temp_DMR_ME <-  mean(na.omit(all_data[[z]][DMR_shared_list[[j]],i]))}
      sample_DMR_ME <- as.numeric(sample_DMR_ME) + as.numeric(temp_DMR_ME)
    }
    shared_DMRavME[i] <- sample_DMR_ME/(length(DMR_shared_list)-n)
  }
  shared_DMRavME <- as.data.frame(shared_DMRavME)
  colnames(shared_DMRavME) <- "Meth_level"
  shared_DMRavME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  shared_DMRavME$Type <- rep("shared_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  shared_DMRavME <- cbind(temp,shared_DMRavME)
  all_DMRdata_shared[[z]] <- shared_DMRavME
}
save(all_DMRdata_shared,file = "D:/MDD_result/21/TCGA_8cancer_DMR.RData")
all_DMRdata_EAC_SPEC <- list()
for (z in 1:length(all_data)) {
  EAC_SPEC_DMRavME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    sample_DMR_ME <- 0
    n <- 0
    for (j in 1:length(DMR_EAC_spec_list)) {
      if(length(DMR_EAC_spec_list[[j]]) == 0) {
        temp_DMR_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][DMR_EAC_spec_list[[j]],i])) == 0){
        temp_DMR_ME <- 0
        n=n+1}
      else{temp_DMR_ME <-  mean(na.omit(all_data[[z]][DMR_EAC_spec_list[[j]],i]))}
      sample_DMR_ME <- as.numeric(sample_DMR_ME) + as.numeric(temp_DMR_ME)
    }
    EAC_SPEC_DMRavME[i] <- sample_DMR_ME/(length(DMR_EAC_spec_list)-n)
  }
  EAC_SPEC_DMRavME <- as.data.frame(EAC_SPEC_DMRavME)
  colnames(EAC_SPEC_DMRavME) <- "Meth_level"
  EAC_SPEC_DMRavME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  EAC_SPEC_DMRavME$Type <- rep("EAC_spec_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  EAC_SPEC_DMRavME <- cbind(temp,EAC_SPEC_DMRavME)
  all_DMRdata_EAC_SPEC[[z]] <- EAC_SPEC_DMRavME
}
save(all_DMRdata_shared,all_DMRdata_EAC_SPEC,file = "D:/MDD_result/21/TCGA_8cancer_DMR.RData")
all_DMRdata_ESCC_SPEC <- list()
for (z in 1:length(all_data)) {
  ESCC_SPEC_DMRavME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    sample_DMR_ME <- 0
    n <- 0
    for (j in 1:length(DMR_ESCC_spec_list)) {
      if(length(DMR_ESCC_spec_list[[j]]) == 0) {
        temp_DMR_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][DMR_ESCC_spec_list[[j]],i])) == 0){
        temp_DMR_ME <- 0
        n=n+1}
      else{temp_DMR_ME <-  mean(na.omit(all_data[[z]][DMR_ESCC_spec_list[[j]],i]))}
      sample_DMR_ME <- as.numeric(sample_DMR_ME) + as.numeric(temp_DMR_ME)
    }
    ESCC_SPEC_DMRavME[i] <- sample_DMR_ME/(length(DMR_ESCC_spec_list)-n)
  }
  ESCC_SPEC_DMRavME <- as.data.frame(ESCC_SPEC_DMRavME)
  colnames(ESCC_SPEC_DMRavME) <- "Meth_level"
  ESCC_SPEC_DMRavME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  ESCC_SPEC_DMRavME$Type <- rep("ESCC_spec_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  ESCC_SPEC_DMRavME <- cbind(temp,ESCC_SPEC_DMRavME)
  all_DMRdata_ESCC_SPEC[[z]] <- ESCC_SPEC_DMRavME
}
save(all_DMRdata_shared,all_DMRdata_EAC_SPEC,all_DMRdata_ESCC_SPEC,file = "D:/MDD_result/21/TCGA_8cancer_DMR.RData")
SC <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
AD <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
# sample_df
#all_DMRdata_shared
for (i in 1:length(all_DMRdata_shared)) {
  temp <- all_DMRdata_shared[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_DMRdata_shared[[i]] <- temp
}
#all_DMRdata_EAC_SPEC
for (i in 1:length(all_DMRdata_EAC_SPEC)) {
  temp <- all_DMRdata_EAC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_DMRdata_EAC_SPEC[[i]] <- temp
}
#all_DMRdata_ESCC_SPEC
for (i in 1:length(all_DMRdata_ESCC_SPEC)) {
  temp <- all_DMRdata_ESCC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_DMRdata_ESCC_SPEC[[i]] <- temp
}
#EAC-ESCC delta methylation
EAC_VS_ESCC_DMRdelta <- list()
for (i in 1:length(all_DMRdata_EAC_SPEC)) {
  delta_meth <- all_DMRdata_EAC_SPEC[[i]]$Meth_level-all_DMRdata_ESCC_SPEC[[i]]$Meth_level
  EAC_VS_ESCC_DMRdelta[[i]] <- all_DMRdata_EAC_SPEC[[i]]
  EAC_VS_ESCC_DMRdelta[[i]]$Meth_level <- delta_meth
  colnames(EAC_VS_ESCC_DMRdelta[[i]]) <- c("Sample","delta_methy","Cancer","Type","Cancer_class")
}
save(EAC_VS_ESCC_DMRdelta,file = "D:/MDD_result/21/8cancers_DMRdelta_meth.RData")
# sample cancerType      delta        label
# 1       1       ESCC 0.12741917 Pan-Squamous
# 2       1       HNSC 0.08870604       Others
# 3       1       LUSC 0.10726257       Pan-GI
# EAC_VS_ESCC_DMRplot_data <- as.data.frame(matrix(ncol = 5))
# colnames(EAC_VS_ESCC_DMRplot_data) <- colnames(EAC_VS_ESCC_DMRdelta[[1]])
# for (i in 1:nrow(EAC_VS_ESCC_DMRdelta[[1]])) {
#   temp_plot_data <- as.data.frame(matrix(ncol = 5))
#   colnames(temp_plot_data) <- colnames(EAC_VS_ESCC_DMRdelta[[1]])
#   for (j in 1:length(EAC_VS_ESCC_DMRdelta)) {
#     temp_plot_data <- rbind(temp_plot_data,EAC_VS_ESCC_DMRdelta[[j]][i,])
#   }
#   EAC_VS_ESCC_DMRplot_data <- rbind(EAC_VS_ESCC_DMRplot_data,temp_plot_data)
# }
EAC_VS_ESCC_DMRplot_data <- as.data.frame(matrix(ncol = 4))
colnames(EAC_VS_ESCC_DMRplot_data) <- c("Sample","Cancer","delta_methy","Cancer_class")
for (i in 1:length(EAC_VS_ESCC_DMRdelta)) {
  temp_plot_data <- EAC_VS_ESCC_DMRdelta[[i]][,c(1,3,2,5)]
  temp_plot_data <- temp_plot_data[which(substr(temp_plot_data$Sample,14,16) == "01A"),]
  temp_plot_data$Sample <- c(1:nrow(temp_plot_data))
  EAC_VS_ESCC_DMRplot_data <- rbind(EAC_VS_ESCC_DMRplot_data,temp_plot_data)
}
EAC_VS_ESCC_DMRplot_data <- EAC_VS_ESCC_DMRplot_data[-1,]
library(ggplot2)
colnames(EAC_VS_ESCC_DMRplot_data) <- c("sample", "cancerType", "delta", "label")
plotTumorATACPointPlot(EAC_VS_ESCC_DMRplot_data,"ESCC only DMRs vs EAC only DMRs", "example_plot2DMR.pdf")
load("D:/MDD_result/4/cg_TCGA_core_select.RData")
cancer <- c("ESCA","STAD","COAD","READ","LIHC","CHOL","PAAD","HNSC")
all_core_MDD <- list()
for (z in 1:length(all_data)) {
  core_MDDavME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    core_MDD_ME <- 0
    n <- 0
    for (j in 1:length(core_MDD_list)) {
      if(length(core_MDD_list[[j]]) == 0) {
        temp_MDD_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][core_MDD_list[[j]],i])) == 0){
        temp_MDD_ME <- 0
        n=n+1}
      else{temp_MDD_ME <-  mean(na.omit(all_data[[z]][core_MDD_list[[j]],i]))}
      core_MDD_ME <- as.numeric(core_MDD_ME) + as.numeric(temp_MDD_ME)
    }
    core_MDDavME[i] <- core_MDD_ME/(length(core_MDD_list)-n)
  }
  core_MDDavME <- as.data.frame(core_MDDavME)
  colnames(core_MDDavME) <- "Meth_level"
  core_MDDavME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  core_MDDavME$Type <- rep("shared_MDD",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  core_MDDavME <- cbind(temp,core_MDDavME)
  all_core_MDD[[z]] <- core_MDDavME
}
save(all_core_MDD,file = "D:/MDD_result/21/TCGA_8cancer_core_MDD.RData")
cancer <- c("ESCA","STAD","COAD","READ","LIHC","CHOL","PAAD","HNSC")
all_core_CMI <- list()
for (z in 1:length(all_data)) {
  core_CMIavME <- c()
  for (i in 1:ncol(all_data[[z]])) {
    core_CMI_ME <- 0
    n <- 0
    for (j in 1:length(core_CMI_list)) {
      if(length(core_CMI_list[[j]]) == 0) {
        temp_CMI_ME <- 0
        n=n+1}
      else if(length(na.omit(all_data[[z]][core_CMI_list[[j]],i])) == 0){
        temp_CMI_ME <- 0
        n=n+1}
      else{temp_CMI_ME <-  mean(na.omit(all_data[[z]][core_CMI_list[[j]],i]))}
      core_CMI_ME <- as.numeric(core_CMI_ME) + as.numeric(temp_CMI_ME)
    }
    core_CMIavME[i] <- core_CMI_ME/(length(core_CMI_list)-n)
  }
  core_CMIavME <- as.data.frame(core_CMIavME)
  colnames(core_CMIavME) <- "Meth_level"
  core_CMIavME$Class <- rep(cancer[z],ncol(all_data[[z]]))
  core_CMIavME$Type <- rep("shared_CMI",ncol(all_data[[z]]))
  temp <- as.data.frame(colnames(all_data[[z]]))
  colnames(temp) <- "Sample"
  core_CMIavME <- cbind(temp,core_CMIavME)
  all_core_CMI[[z]] <- core_CMIavME
}
save(all_core_MDD,all_core_CMI,file = "D:/CMI_result/21/TCGA_8cancer_core_CMI.RData")
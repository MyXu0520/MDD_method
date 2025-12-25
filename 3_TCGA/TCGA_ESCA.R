#TCGA-ESCA
library(readr)
TCGA_ESCA <- read.table("E:/33cancer/33cancers/knn_data/ESCA.methylation450.knn.txt")
clinical_ESCA <- read_tsv("E:/33cancer/33cancers/clinical.cart.2025-12-01/clinical.tsv")
cg_450K <- read.csv("E:/33cancer/33cancers/humanmethylation450.csv",sep = ",",header = T)
cg_450K_locus <- cg_450K[,c(1,12,13,18,22,24,25,26)]
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$Probe_SNPs == ""),]
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$CHR %in% c(as.character(1:22),"X","Y")),]
cg_450K_locus$CHR <- paste0("chr",cg_450K_locus$CHR)
#ESCA_shared_MDD
ESCA_shared_MDD <- read.table("D:/MDD_result/1/ESCA_sharedMDDs.bed")
colnames(ESCA_shared_MDD) <- c("Chrom","Start","End")
#EAC_spec_MDD
EAC_spec_MDD <- read.table("D:/MDD_result/1/EAC_specificMDDs.bed")
colnames(EAC_spec_MDD) <- c("Chrom","Start","End")
#ESCC_spec_MDD
ESCC_spec_MDD <- read.table("D:/MDD_result/1/ESCC_specificMDDs.bed")
colnames(ESCC_spec_MDD) <- c("Chrom","Start","End")
shared_list <- list()
for (i in 1:nrow(ESCA_shared_MDD)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == ESCA_shared_MDD[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= ESCA_shared_MDD[i,2] && temp[j,3] <= ESCA_shared_MDD[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  shared_list[[i]] <- cg_list
}
EAC_spec_list <- list()
for (i in 1:nrow(EAC_spec_MDD)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == EAC_spec_MDD[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= EAC_spec_MDD[i,2] && temp[j,3] <= EAC_spec_MDD[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  EAC_spec_list[[i]] <- cg_list
}
ESCC_spec_list <- list()
for (i in 1:nrow(ESCC_spec_MDD)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == ESCC_spec_MDD[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= ESCC_spec_MDD[i,2] && temp[j,3] <= ESCC_spec_MDD[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  ESCC_spec_list[[i]] <- cg_list
}
save(shared_list,EAC_spec_list,ESCC_spec_list,file = "D:/MDD_result/4/cg_TCGA_select.RData")
clinical_EAC <- read_tsv("D:/MDD_result/4/EAC/clinical.cart.2024-11-20/clinical.tsv")
clinical_ESCC <- read_tsv("D:/MDD_result/4/ESCC/clinical.cart.2024-11-20/clinical.tsv")
EAC_labels <- clinical_EAC$case_submitter_id
ESCC_labels <- clinical_ESCC$case_submitter_id
newcol <- gsub("\\.","-",colnames(TCGA_ESCA))
colnames(TCGA_ESCA) <- newcol
TCGA_EAC_sample <- TCGA_ESCA[,which(substr(colnames(TCGA_ESCA),1,12) %in% EAC_labels)]
TCGA_ESCC_sample <- TCGA_ESCA[,which(substr(colnames(TCGA_ESCA),1,12) %in% ESCC_labels)]
save(TCGA_EAC_sample,TCGA_ESCA,TCGA_ESCC_sample,file = "TCGA_sample_separate.RData")
EAC_shared_avME <- c()
for (i in 1:ncol(TCGA_EAC_sample)) {
  sample_MDD_ME <- 0
  n <- 0
  for (j in 1:length(shared_list)) {
    if(length(shared_list[[j]]) == 0) {
      temp_MDD_ME <- 0
      n=n+1}
    else if(length(na.omit(TCGA_EAC_sample[shared_list[[j]],i])) == 0){
      temp_MDD_ME <- 0
      n=n+1}
    else{temp_MDD_ME <-  mean(na.omit(TCGA_EAC_sample[shared_list[[j]],i]))}
    sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
  }
  EAC_shared_avME[i] <- sample_MDD_ME/(length(shared_list)-n)
}
EAC_shared_avME <- as.data.frame(EAC_shared_avME)
colnames(EAC_shared_avME) <- "Meth_level"
EAC_shared_avME$Class <- rep("EAC",ncol(TCGA_EAC_sample))
EAC_shared_avME$Type <- rep("shared_MDD",ncol(TCGA_EAC_sample))
temp <- as.data.frame(colnames(TCGA_EAC_sample))
colnames(temp) <- "Sample"
EAC_shared_avME <- cbind(temp,EAC_shared_avME)
ESCC_shared_avME <- c()
for (i in 1:ncol(TCGA_ESCC_sample)) {
  sample_MDD_ME <- 0
  n <- 0
  for (j in 1:length(shared_list)) {
    if(length(shared_list[[j]]) == 0) {
      temp_MDD_ME <- 0
      n=n+1}
    else if(length(na.omit(TCGA_ESCC_sample[shared_list[[j]],i])) == 0){
      temp_MDD_ME <- 0
      n=n+1}
    else{temp_MDD_ME <-  mean(na.omit(TCGA_ESCC_sample[shared_list[[j]],i]))}
    sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
  }
  ESCC_shared_avME[i] <- sample_MDD_ME/(length(shared_list)-n)
}
ESCC_shared_avME <- as.data.frame(ESCC_shared_avME)
colnames(ESCC_shared_avME) <- "Meth_level"
ESCC_shared_avME$Class <- rep("ESCC",ncol(TCGA_ESCC_sample))
ESCC_shared_avME$Type <- rep("shared_MDD",ncol(TCGA_ESCC_sample))
temp <- as.data.frame(colnames(TCGA_ESCC_sample))
colnames(temp) <- "Sample"
ESCC_shared_avME <- cbind(temp,ESCC_shared_avME)
# 3)ESCC_ESCC_SPEC
ESCC_ESCC_avME <- c()
for (i in 1:ncol(TCGA_ESCC_sample)) {
  sample_MDD_ME <- 0
  n <- 0
  for (j in 1:length(ESCC_spec_list)) {
    if(length(ESCC_spec_list[[j]]) == 0) {
      temp_MDD_ME <- 0
      n=n+1}
    else if(length(na.omit(TCGA_ESCC_sample[ESCC_spec_list[[j]],i])) == 0){
      temp_MDD_ME <- 0
      n=n+1}
    else{temp_MDD_ME <-  mean(na.omit(TCGA_ESCC_sample[ESCC_spec_list[[j]],i]))}
    sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
  }
  ESCC_ESCC_avME[i] <- sample_MDD_ME/(length(ESCC_spec_list)-n)
}
ESCC_ESCC_avME <- as.data.frame(ESCC_ESCC_avME)
colnames(ESCC_ESCC_avME) <- "Meth_level"
ESCC_ESCC_avME$Class <- rep("ESCC",ncol(TCGA_ESCC_sample))
ESCC_ESCC_avME$Type <- rep("ESCC_SPEC",ncol(TCGA_ESCC_sample))
temp <- as.data.frame(colnames(TCGA_ESCC_sample))
colnames(temp) <- "Sample"
ESCC_ESCC_avME <- cbind(temp,ESCC_ESCC_avME)
# 4)ESCC_EAC_SPEC
ESCC_EAC_avME <- c()
for (i in 1:ncol(TCGA_ESCC_sample)) {
  sample_MDD_ME <- 0
  n <- 0
  for (j in 1:length(EAC_spec_list)) {
    if(length(EAC_spec_list[[j]]) == 0) {
      temp_MDD_ME <- 0
      n=n+1}
    else if(length(na.omit(TCGA_ESCC_sample[EAC_spec_list[[j]],i])) == 0){
      temp_MDD_ME <- 0
      n=n+1}
    else{temp_MDD_ME <-  mean(na.omit(TCGA_ESCC_sample[EAC_spec_list[[j]],i]))}
    sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
  }
  ESCC_EAC_avME[i] <- sample_MDD_ME/(length(EAC_spec_list)-n)
}
ESCC_EAC_avME <- as.data.frame(ESCC_EAC_avME)
colnames(ESCC_EAC_avME) <- "Meth_level"
ESCC_EAC_avME$Class <- rep("ESCC",ncol(TCGA_ESCC_sample))
ESCC_EAC_avME$Type <- rep("EAC_SPEC",ncol(TCGA_ESCC_sample))
temp <- as.data.frame(colnames(TCGA_ESCC_sample))
colnames(temp) <- "Sample"
ESCC_EAC_avME <- cbind(temp,ESCC_EAC_avME)
# 5)EAC_EAC_SPEC
EAC_EAC_avME <- c()
for (i in 1:ncol(TCGA_EAC_sample)) {
  sample_MDD_ME <- 0
  n <- 0
  for (j in 1:length(EAC_spec_list)) {
    if(length(EAC_spec_list[[j]]) == 0) {
      temp_MDD_ME <- 0
      n=n+1}
    else if(length(na.omit(TCGA_EAC_sample[EAC_spec_list[[j]],i])) == 0){
      temp_MDD_ME <- 0
      n=n+1}
    else{temp_MDD_ME <-  mean(na.omit(TCGA_EAC_sample[EAC_spec_list[[j]],i]))}
    sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
  }
  EAC_EAC_avME[i] <- sample_MDD_ME/(length(EAC_spec_list)-n)
}
EAC_EAC_avME <- as.data.frame(EAC_EAC_avME)
colnames(EAC_EAC_avME) <- "Meth_level"
EAC_EAC_avME$Class <- rep("EAC",ncol(TCGA_EAC_sample))
EAC_EAC_avME$Type <- rep("EAC_SPEC",ncol(TCGA_EAC_sample))
temp <- as.data.frame(colnames(TCGA_EAC_sample))
colnames(temp) <- "Sample"
EAC_EAC_avME <- cbind(temp,EAC_EAC_avME)
# 6)EAC_ESCC_SPEC
EAC_ESCC_avME <- c()
for (i in 1:ncol(TCGA_EAC_sample)) {
  sample_MDD_ME <- 0
  n <- 0
  for (j in 1:length(ESCC_spec_list)) {
    if(length(ESCC_spec_list[[j]]) == 0) {
      temp_MDD_ME <- 0
      n=n+1}
    else if(length(na.omit(TCGA_EAC_sample[ESCC_spec_list[[j]],i])) == 0){
      temp_MDD_ME <- 0
      n=n+1}
    else{temp_MDD_ME <-  mean(na.omit(TCGA_EAC_sample[ESCC_spec_list[[j]],i]))}
    sample_MDD_ME <- as.numeric(sample_MDD_ME) + as.numeric(temp_MDD_ME)
  }
  EAC_ESCC_avME[i] <- sample_MDD_ME/(length(ESCC_spec_list)-n)
}
EAC_ESCC_avME <- as.data.frame(EAC_ESCC_avME)
colnames(EAC_ESCC_avME) <- "Meth_level"
EAC_ESCC_avME$Class <- rep("EAC",ncol(TCGA_EAC_sample))
EAC_ESCC_avME$Type <- rep("ESCC_SPEC",ncol(TCGA_EAC_sample))
temp <- as.data.frame(colnames(TCGA_EAC_sample))
colnames(temp) <- "Sample"
EAC_ESCC_avME <- cbind(temp,EAC_ESCC_avME)
EAC_shared_avME <- EAC_shared_avME[which(substr(EAC_shared_avME$Sample,14,16) == "01A"),]
ESCC_shared_avME <- ESCC_shared_avME[which(substr(ESCC_shared_avME$Sample,14,16) == "01A"),]
EAC_EAC_avME <- EAC_EAC_avME[which(substr(EAC_EAC_avME$Sample,14,16) == "01A"),]
ESCC_ESCC_avME <- ESCC_ESCC_avME[which(substr(ESCC_ESCC_avME$Sample,14,16) == "01A"),]
EAC_ESCC_avME <- EAC_ESCC_avME[which(substr(EAC_ESCC_avME$Sample,14,16) == "01A"),]
ESCC_EAC_avME <- ESCC_EAC_avME[which(substr(ESCC_EAC_avME$Sample,14,16) == "01A"),]
#############################
EAC_data_paired_MDD <- rbind(EAC_EAC_avME,EAC_ESCC_avME)
EAC_data_paired_MDD_sorted <- EAC_data_paired_MDD[order(EAC_data_paired_MDD$Sample),]
df <- EAC_data_paired_MDD_sorted %>%
  filter(Type %in% c("EAC_SPEC","ESCC_SPEC")) %>%
  select(Sample,Meth_level,Class,Type)%>%
  mutate(paired = rep(1:(nrow(EAC_data_paired_MDD_sorted)/2),each=2),Type=factor(Type))
head(df)
library(rstatix)
df %>% group_by(Class, Type) %>% shapiro_test(Meth_level)
library(tidyverse)
library(gaCMInder)
library(ggsci)
library(ggprism)
library(rstatix)
library(ggpubr)
pdf("qplot.pdf")
ggqqplot(df, "Meth_level", ggtheme = theme_bw()) +
  facet_grid(Type ~ Class, labeller = "label_both")
dev.off()
df_p_val1 <- df %>% group_by(Class) %>%
  t_test(Meth_level  ~Type, paired = TRUE) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj")
df_p_val1
df_p_val1 <- df_p_val1 %>%
  add_xy_position(x = "Type", dodge = 0.8)
pdf("test_paired.pdf")
ggplot(df,aes( Type,Meth_level)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1) +
  geom_boxplot(position=position_dodge(width =0.2),width=0.4) +
  geom_line(aes(group=paired),position = position_dodge(0.2),color="grey80") +
  geom_point(aes(fill= Type,group=paired,size=Meth_level,alpha=Meth_level),pch=21,
             position = position_dodge(0.2)) +
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=6,hide.ns = T) +
  scale_size_continuous(range=c(1,3)) +
  facet_wrap(.~Class,nrow=1) +
  scale_fill_npg() +
  scale_x_discrete(guide = "prism_bracket") +
  scale_y_continuous(limits = c(0.25,0.75),minor_breaks = seq(0,1,0.1),guide = "prism_offset_minor") +
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
dev.off
save.image("D:/MDD_result/4/ESCA_paired_boxplot.RData")
###################################
data_paired_MDD <- rbind(EAC_EAC_avME,EAC_ESCC_avME,EAC_shared_avME,ESCC_ESCC_avME,ESCC_EAC_avME,ESCC_shared_avME)
data_paired_MDD_sorted <- data_paired_MDD[order(data_paired_MDD$Sample),]
library(dplyr)
df <- data_paired_MDD_sorted %>%
  filter(Type %in% c("EAC_SPEC","ESCC_SPEC","shared_MDD")) %>%
  select(Sample,Meth_level,Class,Type)%>%
  mutate(paired = rep(1:(nrow(data_paired_MDD_sorted)/3),each=3),Type=factor(Type))
head(df)
library(rstatix)
df %>% group_by(Class, Type) %>% shapiro_test(Meth_level)
setwd("D:/MDD_result/4/")
library(tidyverse)
library(gaCMInder)
library(ggsci)
library(ggprism)
library(rstatix)
library(ggpubr)
pdf("qplot_all.pdf")
ggqqplot(df, "Meth_level", ggtheme = theme_bw()) +
  facet_grid(Type ~ Class, labeller = "label_both")
dev.off()
df_p_val1 <- df %>% group_by(Class) %>%
  t_test(Meth_level  ~Type, paired = TRUE) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj")
df_p_val1
df_p_val1 <- df_p_val1 %>%
  add_xy_position(x = "Type", dodge = 0.8)
pdf("test_paired_all.pdf")
ggplot(df,aes( Type,Meth_level)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1) +
  geom_boxplot(position=position_dodge(width =0.2),width=0.4) +
  geom_line(aes(group=paired),position = position_dodge(0.2),color="grey80") +
  geom_point(aes(fill= Type,group=paired,size=Meth_level,alpha=Meth_level),pch=21,
             position = position_dodge(0.2)) +
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=6,hide.ns = T) +
  scale_size_continuous(range=c(1,3)) +
  facet_wrap(.~Class,nrow=1) +
  scale_fill_manual(values = c("#AA2D36","#2F56A3","#64329E"))+
  scale_x_discrete(guide = "prism_bracket") +
  scale_y_continuous(limits = c(0.25,0.75),minor_breaks = seq(0,1,0.1),guide = "prism_offset_minor") +
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
save.image("D:/MDD_result/4/ESCA_paired_boxplot_all.RData")
mycolors <- c("
mycolors <- c("
mycolors <- c("
mycolors <- c("
mycolors <- c("
#EAC
EAC_heatmap <- cbind(EAC_EAC_avME[,c(1,2)])
colnames(EAC_heatmap)[2] <- "EAC_SPEC"
EAC_heatmap <- cbind(EAC_heatmap,EAC_ESCC_avME$Meth_level)
colnames(EAC_heatmap)[3] <- "ESCC_SPEC"
EAC_heatmap <- cbind(EAC_heatmap,EAC_shared_avME$Meth_level)
colnames(EAC_heatmap)[4] <- "shared_MDD"
rownames(EAC_heatmap) <- EAC_heatmap[,1]
EAC_heatmap <- EAC_heatmap[,-1]
library(pheatmap)
p <- pheatmap(as.matrix(EAC_heatmap),cluster_cols = F,cluster_rows = F)
# p1 <- pheatmap(EAC_heatmap, show_rownames = F, cluster_rows = F, cluster_cols = F, border_color = NA,
#             color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(9))
pdf("EAC_heatmap.pdf")
p1 <- pheatmap(EAC_heatmap, show_rownames = F, cluster_rows = F, cluster_cols = F, border_color = NA,
                           color = colorRampPalette(mycolors)(6))
p1
dev.off()
ESCC_heatmap <- cbind(ESCC_ESCC_avME[,c(1,2)])
colnames(ESCC_heatmap)[2] <- "ESCC_SPEC"
ESCC_heatmap <- cbind(ESCC_heatmap,ESCC_ESCC_avME$Meth_level)
colnames(ESCC_heatmap)[3] <- "ESCC_SPEC"
ESCC_heatmap <- cbind(ESCC_heatmap,ESCC_shared_avME$Meth_level)
colnames(ESCC_heatmap)[4] <- "shared_MDD"
rownames(ESCC_heatmap) <- ESCC_heatmap[,1]
ESCC_heatmap <- ESCC_heatmap[,-1]
pdf("ESCC_heatmap.pdf")
p2 <- pheatmap(ESCC_heatmap, show_rownames = F, cluster_rows = F, cluster_cols = F, border_color = NA,
               color = colorRampPalette(mycolors)(6))
p2
dev.off()
library(ggpubr)
library(tidyverse)
library(ggprism)
library(MetBrewer)
library(Rmisc)
EAC_bar <- rbind(EAC_EAC_avME,EAC_ESCC_avME,EAC_shared_avME)
ESCC_bar <- rbind(ESCC_EAC_avME,ESCC_ESCC_avME,ESCC_shared_avME)
bar_plot <- rbind(EAC_bar,ESCC_bar)
bar_plot$Type <- as.factor(bar_plot$Type)
EAC_bar <- summarySE(EAC_bar, measurevar="Meth_level", groupvars=c("Class","Type"))
ESCC_bar <- summarySE(ESCC_bar, measurevar="Meth_level", groupvars=c("Class","Type"))
bar_plot <- summarySE(bar_plot, measurevar="Meth_level", groupvars=c("Class","Type"))
bar_plot$Class <- as.factor(bar_plot$Class)
# Standard error of the mean
# ggplot(bar_plot, aes(x=Type, y=Meth_level,color = Class)) +
#   geom_errorbar(aes(ymin=Meth_level-se, ymax=Meth_level+se), width=.1) +
#   geom_line() +
#   geom_point()+
#   theme_bw()
#
pd <- position_dodge(0.1)
pEAC <- ggplot(EAC_bar, aes(x=Type, y=Meth_level)) +
  geom_errorbar(aes(ymin=Meth_level-se, ymax=Meth_level+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Type") +
  ylab("Methylation_level") +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("OJ", "VC"),
                   labels=c("Orange juice", "Ascorbic acid"),
                   l=40) +                    # Use darker colors, lightness=40
  # ggtitle("") +
  # expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))               # Position legend in bottom right
pESCC <- ggplot(ESCC_bar, aes(x=Type, y=Meth_level)) +
  geom_errorbar(aes(ymin=Meth_level-se, ymax=Meth_level+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Type") +
  ylab("Methylation_level") +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("OJ", "VC"),
                   labels=c("Orange juice", "Ascorbic acid"),
                   l=40) +                    # Use darker colors, lightness=40
  # ggtitle("") +
  # expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))
load("D:/MDD_result/4/TCGA_sample_separate.RData")
# which(substr(colnames(TCGA_ESCA),14,16) == "11A")
TCGA_ESCA <- TCGA_ESCA[]
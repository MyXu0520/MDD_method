library(TCGAbiolinks)
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("minfiData", "sva"))
library(devtools)
install_github("hansenlab/tutorial.450k")
library(minfi)
library(minfiData)
library(sva)
ab <- compartments(TCGA_EAC_sample, chr="chr14", resolution=100*1000)
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("minfiData"))
baseDir <- "/student1/MDD/TCGA/IDAT/EAC/IDAT"
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets = targets)
MSet <- preprocessRaw(RGSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)
# https://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE,
                                     mergeManifest = FALSE, sex = NULL)
ab <- compartments(GRset.quantile, resolution=100*1000, what = "OpenSea",
                       method = c("pearson", "spearman"), keep=TRUE)
load("/student1/MDD/TCGA/IDAT/EAC/EAC_RGSet.RData")
GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE,
                                     mergeManifest = FALSE, sex = NULL)
for (i in 1:22) {
  ab <- compartments(GRset.quantile, resolution=100*1000, what = "OpenSea",
                     method = "pearson",keep=TRUE,chr = paste0("chr",i))
  save(ab,file = paste0("/student1/MDD/TCGA/IDAT/EAC/ab/chr",i,".RData"))
}
load("/student1/MDD/TCGA/IDAT/ESCC/ESCC_RGSet.RData")
GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE,
                                     mergeManifest = FALSE, sex = NULL)
for (i in 1:22) {
  ab <- compartments(GRset.quantile, resolution=100*1000, what = "OpenSea",
                     method = "pearson",keep=TRUE,chr = paste0("chr",i))
  save(ab,file = paste0("/student1/MDD/TCGA/IDAT/ESCC/ab/chr",i,".RData"))
}
load("D:/MDD_result/3/EAC/ab/chr1.RData")
EACab <- ab
load("D:/MDD_result/3/ESCC/ab/chr1.RData")
ESCCab <- ab
rm(ab)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(9, "RdBu"))(50)
coul <- colorRampPalette(brewer.pal(9, "Blues"))(50)
temp1 <- EACab@elementMetadata@listData[["cor.matrix"]]
temp2 <- ESCCab@elementMetadata@listData[["cor.matrix"]]
heatmap(temp1[1:200,1:200],Colv = NA, Rowv = NA, scale="column",col = coul)
heatmap(temp2[1:200,1:200],Colv = NA, Rowv = NA, scale="column",col = coul)
heatmap(temp1[1:100,1:100],Colv = NA, Rowv = NA, scale="column",col = coul)
heatmap(temp2[1:100,1:100],Colv = NA, Rowv = NA, scale="column",col = coul)
library(corrplot)
corrplot(temp1,
         # method = c("circle", "square", "ellipse", "number", "shade", "color", "pie"),
         method = "square",
         # type = c("full", "lower", "upper"),
         type = "full",
         add = FALSE,
         col = NULL, bg = "white", title = "",  is.corr = TRUE,
         diag = TRUE, outline = FALSE, mar = c(0,0,0,0),
         addgrid.col = NULL, addCoef.col = NULL, addCoefasPercent = FALSE,
         order = c("original", "AOE", "FPC", "hclust", "alphabet"),
         # hclust.method = c("complete", "ward", "single", "average",
                           # "mcquitty", "median", "centroid"),
         addrect = NULL, rect.col = "black", rect.lwd = 2,
         tl.pos = NULL, tl.cex = 1,
         tl.col = "red", tl.offset = 0.4, tl.srt = 90,
         cl.pos = NULL, cl.lim = NULL,
         cl.length = NULL, cl.cex = 0.8, cl.ratio = 0.15,
         cl.align.text = "c",cl.offset = 0.5,
         # addshade = c("negative", "positive", "all"),
         addshade = "all",
         shade.lwd = 1, shade.col = "white",
         p.mat = NULL, sig.level = 0.05,
         # insig = c("pch","p-value","blank", "n"),
         pch = 4, pch.col = "black", pch.cex = 3,
         # plotCI = c("n","square", "circle", "rect"),
         # plotCI = "rect"
         ,# lowCI.mat = NULL, uppCI.mat = NULL
         # , ...
         )
# install.packages("gaston")
# library(gaston)
# Load data
# data(AGT)
# x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
# # Compute LD
# ld.x <- LD(x, c(1,ncol(x)))
# # Plot a tiny part of the LD matrix
# LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20] )
#
#
#
#
# EAC_open <- length(which(EACab@elementMetadata@listData[["compartment"]]=="open"))
# EAC_closed <- length(which(EACab@elementMetadata@listData[["compartment"]]=="closed"))
# ESCC_open <- length(which(ESCCab@elementMetadata@listData[["compartment"]]=="open"))
# ESCC_closed <- length(which(ESCCab@elementMetadata@listData[["compartment"]]=="closed"))
#EAC
EAC_B_chr_all <- data.frame()
for (i in c(1:22,"X","Y")) {
  load(paste0("C:/Users/MyXu/Desktop/result/3/EAC/ab/chr",i,".RData"))
  EACab <- ab
  rm(ab)
  EAC_B <- EACab@ranges@start
  EAC_B <- as.data.frame(EAC_B)
  colnames(EAC_B) <- "Start"
  EAC_B$End <- EAC_B$Start+100*1000
  EAC_B$Compartment <- EACab@elementMetadata@listData[["compartment"]]
  EAC_B <- EAC_B[which(EAC_B$Compartment == "closed"),]
  EAC_B$Chr <- paste0("chr",i)
  EAC_B_chr_all <- rbind(EAC_B_chr_all,EAC_B)
}
EAC_B_chr_all <- EAC_B_chr_all[,c(4,1,2)]
# write.table(EAC_B_chr_all,"C:/Users/MyXu/"C:/Users/MyXu/Desktop/result/3/EAC/EAC_B_chr_all.bed"Desktop/result/3/EAC/EAC_B_chr_all.bed",col.names = F,row.names = F,sep = "\t")
write.table(EAC_B_chr_all, file = "C:/Users/MyXu/Desktop/result/3/EAC/EAC_B_chr_all.bed", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
ESCC_B_chr_all <- data.frame()
for (i in c(1:22,"X","Y")) {
  load(paste0("C:/Users/MyXu/Desktop/result/3/ESCC/ab/chr",i,".RData"))
  ESCCab <- ab
  rm(ab)
  ESCC_B <- ESCCab@ranges@start
  ESCC_B <- as.data.frame(ESCC_B)
  colnames(ESCC_B) <- "Start"
  ESCC_B$End <- ESCC_B$Start+100*1000
  ESCC_B$Compartment <- ESCCab@elementMetadata@listData[["compartment"]]
  ESCC_B <- ESCC_B[which(ESCC_B$Compartment == "closed"),]
  ESCC_B$Chr <- paste0("chr",i)
  ESCC_B_chr_all <- rbind(ESCC_B_chr_all,ESCC_B)
}
ESCC_B_chr_all <- ESCC_B_chr_all[,c(4,1,2)]
write.table(ESCC_B_chr_all,"C:/Users/MyXu/Desktop/result/3/ESCC/ESCC_B_chr_all.bed",col.names = F,row.names = F,sep = "\t",quote = FALSE)
#A compartments
# EAC_EAC_overalp <- 0.1638
# EAC_ESCC_overlap <- 0.0293
# EAC_shared_overlap <- 0.3838
#
# ESCC_ESCC_overalp <- 0.0240
# ESCC_EAC_overlap <- 0.1880
# ESCC_shared_overlap <- 0.3898
EAC_EAC_overalp <- 0.1597
EAC_ESCC_overlap <- 0.0201
EAC_shared_overlap <- 0.5467
ESCC_ESCC_overalp <- 0.0254
ESCC_EAC_overlap <- 0.1356
ESCC_shared_overlap <- 0.5407
Group <- c("EAC_SPEC","EAC_SPEC","ESCC_SPEC","ESCC_SPEC","shared_MDD","shared_MDD")
Group <- as.data.frame(Group)
colnames(Group) <- "Group"
SubGroup <- c("EAC","ESCC","EAC","ESCC","EAC","ESCC")
Group$Subgroup <- SubGroup
Index <- c(EAC_EAC_overalp,ESCC_EAC_overlap,EAC_ESCC_overlap,ESCC_ESCC_overalp,EAC_shared_overlap,ESCC_shared_overlap)
Group$Index <- Index
plotdata <- Group
library(plyr)
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col) {
    mean_value <- mean(x[[col]], na.rm = TRUE)
    sd_value <- sd(x[[col]], na.rm = TRUE)
    length_n <- length(x[[col]])
    # standard error= standard deviation/squareroot(n)
    se_value <- sd_value / sqrt(length_n)
    return(c(mean = mean_value, sd = sd_value, se = se_value))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  colnames(data_sum)[which(colnames(data_sum) == "mean")] <- varname
  return(data_sum)
}
plotdata2 <- data_summary(
  data = plotdata,
  varname = "Index",
  groupnames = c("Group", "Subgroup"))
head(plotdata2)
# df_p_val <- plotdata %>%
#   rstatix::group_by(Group) %>%
#   rstatix::t_test(Index ~ Subgroup) %>%
#   rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   rstatix::add_significance(p.col = "p.adj") %>%
#   rstatix::add_xy_position(x = "Group", dodge = 0.8) %>% # important for positioning!
#   dplyr::mutate(p.adj = round(p.adj, 3))
#
# head(df_p_val)
library(ggplot2)
library(ggprism)
ggplot(plotdata, aes(x = Group, y = Index, fill = Subgroup)) +  geom_col (width = 0.5, position = position_dodge(0.7))+
  scale_fill_manual(values = c("#FFCD44","#008F91"))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()
# pl_prism <- ggplot(data = plotdata2, aes(x = Group, y = Index)) +
#   geom_bar(aes(color = Subgroup), stat = "identity",
#            fill = "white",
#            position = position_dodge(width = 0.95)) +
#   # geom_errorbar(aes(color = Subgroup,
#   #                   ymin = Index - sd, ymax = Index + sd),
#   #               width = 0.2, position = position_dodge(width = 0.9),
#   # #               size = 1) +
#   # ggprism::add_pvalue(data = df_p_val,
#   #                     xmin = "xmin",
#   #                     xmax = "xmax",
#   #                     label = "p = {p.adj}",
#   #                     tip.length = 0) +
#   geom_point(data = plotdata,
#              aes(x = Group, y = Index, shape = Subgroup, color = Subgroup),
#              position = position_jitterdodge(jitter.width = 0.5,
#                                              dodge.width = 0.7),
#              size = 2, show.legend = FALSE) +
#   labs(x = "") +
#   scale_y_continuous(guide = "prism_minor", # prism_offset
#                      minor_breaks = seq(0, 1, 0.25),
#                      limits = c(0, 1)) +
#                      # expand = expansion(mult = c(0, 0.1))) +
#   scale_fill_prism() +
#   scale_color_prism() +
#   scale_shape_prism() +
#   guides(shape = "none") +
#   theme_prism() +
#   theme(axis.title = element_text(size = 12, color = "black", face = "bold"),
#         axis.text = element_text(size = 10, color = "black"),
#         text = element_text(size = 9, color = "black"),
#         legend.position = c(0.15, 0.85))
#
# pl_prism
#
#A
EAC_EAC_overalp <- 0.1614
EAC_ESCC_overlap <- 0.2212
EAC_shared_overlap <- 0.3405
ESCC_ESCC_overalp <- 0.2132
ESCC_EAC_overlap <- 0.1527
ESCC_shared_overlap <- 0.3313
Group <- c("EAC_SPEC","EAC_SPEC","ESCC_SPEC","ESCC_SPEC","shared_DMR","shared_DMR")
Group <- as.data.frame(Group)
colnames(Group) <- "Group"
SubGroup <- c("EAC","ESCC","EAC","ESCC","EAC","ESCC")
Group$Subgroup <- SubGroup
Index <- c(EAC_EAC_overalp,ESCC_EAC_overlap,EAC_ESCC_overlap,ESCC_ESCC_overalp,EAC_shared_overlap,ESCC_shared_overlap)
Group$Index <- Index
plotdata <- Group
library(plyr)
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col) {
    mean_value <- mean(x[[col]], na.rm = TRUE)
    sd_value <- sd(x[[col]], na.rm = TRUE)
    length_n <- length(x[[col]])
    # standard error= standard deviation/squareroot(n)
    se_value <- sd_value / sqrt(length_n)
    return(c(mean = mean_value, sd = sd_value, se = se_value))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  colnames(data_sum)[which(colnames(data_sum) == "mean")] <- varname
  return(data_sum)
}
plotdata2 <- data_summary(
  data = plotdata,
  varname = "Index",
  groupnames = c("Group", "Subgroup"))
head(plotdata2)
# df_p_val <- plotdata %>%
#   rstatix::group_by(Group) %>%
#   rstatix::t_test(Index ~ Subgroup) %>%
#   rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   rstatix::add_significance(p.col = "p.adj") %>%
#   rstatix::add_xy_position(x = "Group", dodge = 0.8) %>% # important for positioning!
#   dplyr::mutate(p.adj = round(p.adj, 3))
#
# head(df_p_val)
library(ggplot2)
library(ggprism)
ggplot(plotdata, aes(x = Group, y = Index, fill = Subgroup)) +  geom_col (width = 0.5, position = position_dodge(0.7))+
  scale_fill_manual(values = c("#FFCD44","#008F91"))+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()
#B
EAC_EAC_overalp <- 0.1298
EAC_ESCC_overlap <- 0.2470
EAC_shared_overlap <- 0.2754
ESCC_ESCC_overalp <- 0.2550
ESCC_EAC_overlap <- 0.1384
ESCC_shared_overlap <- 0.2847
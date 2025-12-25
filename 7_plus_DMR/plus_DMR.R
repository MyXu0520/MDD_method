#########################ESCC
#plus_DMR
library(DSS)
library(bsseq)
dir_path <- "/student1/MDD/Tabfile"
files <- list.files(dir_path, pattern = "\\.tab$", full.names = TRUE)
data_list <- list()
sample_names <- c()
cancer_files <- files[grepl("^ESCC_", basename(files))]
for (file in cancer_files) {
  sample_name <- gsub("^ESCC_|\\.tab$", "", basename(file))
  sample_names <- c(sample_names, paste0("C", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
data_listC <- data_list
data_listC <- data_list[c(1:21)]
rm(data_list,sample_names)
data_list <- list()
sample_names <- c()
nonmalignant_files <- files[grepl("^ESCC_Nonmalignant_", basename(files))]
for (file in nonmalignant_files) {
  sample_name <- gsub("^ESCC_Nonmalignant_|\\.tab$", "", basename(file))
  sample_name <- gsub("_", "", sample_name)
  sample_names <- c(sample_names, paste0("N", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
data_listN <- data_list
rm(data_list,sample_names)
data_list <- c(data_listC,data_listN)
BSobj <- makeBSseqData(data_list, names(data_list))
# BSobj <- BSobj[1:1000,]
print(BSobj)
# [1] "C10methyCov"  "C11methyCov"  "C12methyCov"  "C13methyCov"  "C14methyCov"
# [6] "C15methyCov"  "C16methyCov"  "C17methyCov"  "C19methyCov"  "C1methyCov"
# [11] "C20methyCov"  "C21methyCov"  "C22methyCov"  "C2methyCov"   "C3methyCov"
# [16] "C4methyCov"   "C5methyCov"   "C6methyCov"   "C7methyCov"   "C8methyCov"
# [21] "C9methyCov"   "N1methyCov"   "N2methyCov"   "N3methyCov"   "N53FmethyCov"
# [26] "N54MmethyCov"
dmlTest <- DMLtest(BSobj, group1=c("C10methyCov","C11methyCov","C12methyCov","C13methyCov","C14methyCov","C15methyCov","C16methyCov",
                                   "C17methyCov","C19methyCov","C1methyCov","C20methyCov" ,"C21methyCov","C22methyCov","C2methyCov" ,
                                   "C3methyCov","C4methyCov"  ,"C5methyCov" ,"C6methyCov" ,"C7methyCov" ,"C8methyCov","C9methyCov"),
                   group2=c("N1methyCov","N2methyCov","N3methyCov"  ,"N53FmethyCov","N54MmethyCov"))
head(dmlTest)
single = MulticoreParam(workers=1, progressbar=TRUE)
#dmlTest.sm_1 <- DMLtest(BSobj, group1="mdf", group2="mp", equal.disp = T)
dmlTest.sm <- DMLtest(BSobj, group1=c("C10methyCov","C11methyCov","C12methyCov","C13methyCov","C14methyCov","C15methyCov","C16methyCov",
                                        "C17methyCov","C19methyCov","C1methyCov","C20methyCov" ,"C21methyCov","C22methyCov","C2methyCov" ,
                                        "C3methyCov","C4methyCov"  ,"C5methyCov" ,"C6methyCov" ,"C7methyCov" ,"C8methyCov","C9methyCov"),
                        group2=c("N1methyCov","N2methyCov","N3methyCov"  ,"N53FmethyCov","N54MmethyCov"), smoothing=TRUE, BPPARAM=single)
dmlTest.sm <- DMLtest(BSobj, group1=c("C10methyCov","C11methyCov","C12methyCov","C13methyCov","C14methyCov","C15methyCov","C16methyCov",
                                      "C17methyCov","C19methyCov","C1methyCov","C20methyCov" ,"C21methyCov","C22methyCov","C2methyCov" ,
                                      "C3methyCov","C4methyCov"  ,"C5methyCov" ,"C6methyCov" ,"C7methyCov" ,"C8methyCov","C9methyCov"),
                      group2=c("N1methyCov","N2methyCov","N3methyCov"  ,"N53FmethyCov","N54MmethyCov"), smoothing=TRUE,ncores = 2)
head(dmlTest.sm)
save.image(file = "smothing_DMR.RData")
dmls <- callDML(dmlTest, p.threshold=0.001)
# dmls
dmls2 <- callDML(dmlTest, delta=0.1, p.threshold=0.001)
# head(dmls2)
dmls <- callDML(dmlTest.sm,delta=0.1, p.threshold=0.001)
# dmls
dmrs <- callDMR(dmlTest.sm, p.threshold=0.001)
# dmrs
dmrs2 <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
dmrs2
pdf("test.pdf")
showOneDMR(dmrs[1,], BSobj)
dev.off()
load("D:/MDD_result/7/EAC_ESCC_dmrs_dml.RData")
library(qqman)
library(tidyverse)
dmrs <- dmrs[,c(1,2,3,4,8)]
dmrs <- dmrs[order(dmrs$diff.Methy),]
dmrs <- arrange(dmrs,"chr","diff.Methy")
SNP <- paste0("DMR",c(1:nrow(dmrs)))
SNP <- as.data.frame(SNP)
colnames(SNP) <- "SNP"
dmrs <- cbind(SNP,dmrs)
colnames(dmrs) <- c("SNP","CHR","BP","BPcum","total","P")
dmrs <- arrange(dmrs,CHR,P)
# SNP <- paste0("DMR",c(1:nrow(dmrs)))
# SNP <- as.data.frame(SNP)
# colnames(SNP) <- "SNP"
# dmrs$SNP <- SNP
temp <- dmrs$CHR
temp_CHR <- c()
for (i in 1:nrow(dmrs)) {
  if(temp[i] == "chr21"){temp_CHR[i] <- 21}
  else if(temp[i] == "chr1"){temp_CHR[i] <- 1}
  else if(temp[i] == "chr2"){temp_CHR[i] <- 2}
  else if(temp[i] == "chr3"){temp_CHR[i] <- 3}
  else if(temp[i] == "chr4"){temp_CHR[i] <- 4}
  else if(temp[i] == "chr5"){temp_CHR[i] <- 5}
  else if(temp[i] == "chr6"){temp_CHR[i] <- 6}
  else if(temp[i] == "chr7"){temp_CHR[i] <- 7}
  else if(temp[i] == "chr8"){temp_CHR[i] <- 8}
  else if(temp[i] == "chr9"){temp_CHR[i] <- 9}
  else if(temp[i] == "chr10"){temp_CHR[i] <- 10}
  else if(temp[i] == "chr11"){temp_CHR[i] <- 11}
  else if(temp[i] == "chr12"){temp_CHR[i] <- 12}
  else if(temp[i] == "chr13"){temp_CHR[i] <- 13}
  else if(temp[i] == "chr14"){temp_CHR[i] <- 14}
  else if(temp[i] == "chr15"){temp_CHR[i] <- 15}
  else if(temp[i] == "chr16"){temp_CHR[i] <- 16}
  else if(temp[i] == "chr17"){temp_CHR[i] <- 17}
  else if(temp[i] == "chr18"){temp_CHR[i] <- 18}
  else if(temp[i] == "chr19"){temp_CHR[i] <- 19}
  else if(temp[i] == "chr20"){temp_CHR[i] <- 20}
  else{temp_CHR[i] <- 22}
}
dmrs$CHR <- temp_CHR
# DMR$chr <- as.numeric(DMR$chr)
# DMR$start <- as.numeric(DMR$start)
# DMR$end <- as.numeric(DMR$end)
# colnames(DMR)[1] <- "CHR"
# colnames(DMR)[2] <- "BP"
# colnames(DMR)[7] <- "P"
# DMR$SNP <- rownames(DMR)
# DMR <- DMR[,c(11,1:10)]
# DMR <- DMR[,c(1,2,3,8)]
# CpGinterest <- dmrs %>% arrange(P) %>% pull(SNP) %>% head(100)
# pdf("manhatten_DML.pdf")
# # manhattan(dmrs,highlight = CpGinterest)
# manhattan(dmrs)
# dev.off()
dmrs1 <- dmrs[which(dmrs$P >= 0),]
scale_values <- function (x){(x-min(x))/(max(x)-min(x))}
dmrs1$P <- scale_values(dmrs1$P)
dmrs2 <- dmrs[which(dmrs$P < 0),]
dmrs2$P <- abs(dmrs2$P)
dmrs2$P <- scale_values(dmrs2$P)
library(CMplot)
library(openxlsx)
library(ggsci)
# setwd("C:/Users/MyXu/Desktop/result/7")
mycolor <- c(pal_d3("category20",alpha = 0.5)(20),pal_d3("category20",alpha = 0.5)(3))
set.seed(666666)
# dmrs$P <- scale(dmrs$P)
pdf("D:/MDD_result/7/EAC_ESCC_manhatten_DML2.pdf",width = 20,height = 5)
CMplot(dmrs2,
       plot.type="m",
       LOG10=F,
       chr.labels.angle = 45,
       chr.border = F,
       # col= c("#3E0A52", "#423D77","#3F678B",
       #        "#468C8D", "#5FB47F", "#9FD55C","#F9E956"),
       col = mycolor,
       chr.den.col=c("darkgreen", "yellow", "red"),
       highlight = dmrs$SNP,
       highlight.col = NULL,
       highlight.cex = 1,
       highlight.pch = c(15:17),
       # highlight.text = genes,
       # highlight.text.col = "black",
       threshold = 0.05/nrow(dmrs),
       threshold.lty = 2,
       threshold.col = "black",
       threshold.lwd = 3,
       amplify = FALSE,
       file = "pdf",
       # memo = "",
       dpi = 300,
       file.output = F,
       verbose = F)
# width = 28,height = 4)
dev.off()
# library(karyoploteR)
# kp <- plotKaryotype(plot.type=4)
# kp <- kpPlotManhattan(kp, data=dmrs$SNP)
# data <- df_SNP_position %>%
#   mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#   mutate( is_annotate=ifelse(-log10(P)>10, "yes", "no"))
dmrs <- dmrs[order(dmrs$CHR,dmrs$BP),]
dmrs$SNP <-SNP
BP <- c()
for (i in 1:22) {
  BP[i] <- length(which(dmrs$CHR == i))
}
BP_num <- c()
for (i in 1:22) {
  BP_num <- c(BP_num,c(1:BP[i]))
}
dmrs$BP <- BP_num
X_axis <-  dmrs %>% group_by(CHR) %>% summarize(center=( max(BP) +min(BP) ) / 2 )
library(ggforce)
library(ggplot2)
library(tidyverse)
library(ggforce)
library(ggprism)
ggplot(dmrs, aes(x=BP, y=P)) +
  geom_point(aes(color=as.factor(CHR)),alpha=0.8, size=1.5)+
  scale_color_manual(values = rep(c('#30A9DE','#EFDC05','#E53A40','#090707'), 22 ))+
  scale_x_continuous(label = X_axis$CHR, breaks= X_axis$center)+
  scale_y_continuous(expand = c(0, 0) ) +
  geom_hline(yintercept = c(1, 0.05/nrow(dmrs)),
             color = c('green', 'red'),size = 1.2,
             linetype = c("dotted", "twodash"))
# +
  # geom_point(data=subset(dmrs,
                         # is_highlight=="yes"
                         # ), color="green",
             # size=2)
# +
  # facet_zoom(x = BPcum >= 3000 & BPcum <=3500)+
  # theme_prism(palette = "flames",
  #             base_fontface = "plain",
  #             base_family = "serif",
  #             base_size = 16,
  #             base_line_size = 0.8,
  #             axis_text_angle = 45)+
  #       panel.grid = element_blank(),
  #       panel.border = element_blank(),
  #       axis.line.x = element_line(),
  #       axis.line.y = element_line())
# # manhattan(dmls)
# dmls <- dmls[,c(1,2,10)]
# dmls <- dmls[order(dmls$pval),]
# temp <- dmls[which(dmls$pval == 0.000000e+00),]
# dmls <- dmls[-which(rownames(dmls) %in% rownames(temp)),]
#
# dmls <- arrange(dmls,"chr","pval")
# SNP <- paste0("DMR",c(1:nrow(dmls)))
# SNP <- as.data.frame(SNP)
# colnames(SNP) <- "SNP"
# dmls <- cbind(SNP,dmls)
#
# colnames(dmls) <- c("SNP","CHR","BP","P")
# dmls <- arrange(dmls,CHR,P)
#
# SNP <- paste0("DMR",c(1:nrow(dmls)))
# SNP <- as.data.frame(SNP)
# colnames(SNP) <- "SNP"
# dmls$SNP <- SNP
# dmls <- dmls[,-2]
# for (i in 1:22) {
#   for (j in 1:nrow(dmls)) {
#     if(dmls$CHR[j] == paste0("chr",i)){dmls$CHR[j] <-i}
#     else{next}
#   }
#
# }
# temp <- dmls$CHR
# temp_CHR <- c()
# for (i in 1:nrow(dmls)) {
#   if(temp[i] == "chr21"){temp_CHR[i] <- 21}
#   else if(temp[i] == "chr1"){temp_CHR[i] <- 1}
#   else if(temp[i] == "chr2"){temp_CHR[i] <- 2}
#   else if(temp[i] == "chr3"){temp_CHR[i] <- 3}
#   else if(temp[i] == "chr4"){temp_CHR[i] <- 4}
#   else if(temp[i] == "chr5"){temp_CHR[i] <- 5}
#   else if(temp[i] == "chr6"){temp_CHR[i] <- 6}
#   else if(temp[i] == "chr7"){temp_CHR[i] <- 7}
#   else if(temp[i] == "chr8"){temp_CHR[i] <- 8}
#   else if(temp[i] == "chr9"){temp_CHR[i] <- 9}
#   else if(temp[i] == "chr10"){temp_CHR[i] <- 10}
#   else if(temp[i] == "chr11"){temp_CHR[i] <- 11}
#   else if(temp[i] == "chr12"){temp_CHR[i] <- 12}
#   else if(temp[i] == "chr13"){temp_CHR[i] <- 13}
#   else if(temp[i] == "chr14"){temp_CHR[i] <- 14}
#   else if(temp[i] == "chr15"){temp_CHR[i] <- 15}
#   else if(temp[i] == "chr16"){temp_CHR[i] <- 16}
#   else if(temp[i] == "chr17"){temp_CHR[i] <- 17}
#   else if(temp[i] == "chr18"){temp_CHR[i] <- 18}
#   else if(temp[i] == "chr19"){temp_CHR[i] <- 19}
#   else if(temp[i] == "chr20"){temp_CHR[i] <- 20}
#   else{temp_CHR[i] <- 22}
# }
# dmls$CHR <- temp_CHR
#
#
# # DMR$chr <- as.numeric(DMR$chr)
# # DMR$start <- as.numeric(DMR$start)
# # DMR$end <- as.numeric(DMR$end)
# # colnames(DMR)[1] <- "CHR"
# # colnames(DMR)[2] <- "BP"
# # colnames(DMR)[7] <- "P"
# # DMR$SNP <- rownames(DMR)
# # DMR <- DMR[,c(11,1:10)]
# # DMR <- DMR[,c(1,2,3,8)]
# # CpGinterest <- dmls %>% arrange(P) %>% pull(SNP) %>% head(100)
# # pdf("manhatten_DML.pdf")
# # # manhattan(dmls,highlight = CpGinterest)
# # manhattan(dmls)
# # dev.off()
#
#
# library(CMplot)
# library(openxlsx)
# library(ggsci)
# # setwd("C:/Users/MyXu/Desktop/result/7")
# mycolor <- c(pal_d3("category20",alpha = 0.5)(20),pal_d3("category20",alpha = 0.5)(3))
# set.seed(666666)
# pdf("D:/MDD_result/7/EAC_manhatten_DML.pdf",width = 20,height = 5)
# CMplot(dmls,
#        plot.type="m",
#        LOG10=TRUE,
#        chr.labels.angle = 45,
#        chr.border = F,
#        # col= c("#3E0A52", "#423D77","#3F678B",
#        #        "#468C8D", "#5FB47F", "#9FD55C","#F9E956"),
#        col = mycolor,
#        chr.den.col=c("darkgreen", "yellow", "red"),
#        highlight = dmls$SNP,
#        highlight.col = NULL,
#        highlight.cex = 1,
#        highlight.pch = c(15:17),
#        # highlight.text = genes,
#        # highlight.text.col = "black",
#        threshold = 0.05/nrow(dmls),
#        threshold.lty = 2,
#        threshold.col = "black",
#        threshold.lwd = 3,
#        amplify = FALSE,
#        file = "pdf",
#        # memo = "",
#        dpi = 300,
#        file.output = F,
#        verbose = F)
#        # width = 28,height = 4)
# dev.off()
##########################EAC
# plus_DMR
# library(DSS)
# library(bsseq)
#
# dir_path <- "/student1/MDD/Tabfile"
#
# files <- list.files(dir_path, pattern = "\\.tab$", full.names = TRUE)
#
data_list <- list()
sample_names <- c()
cancer_files <- files[grepl("^EAC_", basename(files))]
for (file in cancer_files) {
  sample_name <- gsub("^EAC_|\\.tab$", "", basename(file))
  sample_names <- c(sample_names, paste0("C", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
data_listC2 <- data_list
data_listC <- data_list[c(1:21)]
rm(data_list,sample_names)
#
data_list <- list()
sample_names <- c()
nonmalignant_files <- files[grepl("^GEJ_Nonmalignant_", basename(files))]
for (file in nonmalignant_files) {
  sample_name <- gsub("^GEJ_Nonmalignant_|\\.tab$", "", basename(file))
  sample_name <- gsub("_", "", sample_name)
  sample_names <- c(sample_names, paste0("N", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
data_listN <- data_list
rm(data_list,sample_names)
data_list <- c(data_listC,data_listN)
BSobj <- makeBSseqData(data_list, names(data_list))
# # BSobj <- BSobj[1:1000,]
#
# print(BSobj)
# [1] "C1methyCov" "C2methyCov" "C3methyCov" "C4methyCov" "C6methyCov"
# [6] "C1methyCov" "C2methyCov" "C3methyCov" "C4methyCov" "C5methyCov"
# [11] "C6methyCov" "C7methyCov" "N1methyCov" "N2methyCov" "N3methyCov"
# [16] "N4methyCov" "N5methyCov" "N6methyCov" "N7methyCov"
dmlTest <- DMLtest(BSobj, group1=c("C1methyCov","C2methyCov","C3methyCov","C4methyCov","C5methyCov","C6methyCov",
                                   "C7methyCov","C8methyCov","C9methyCov","C10methyCov","C11methyCov","C12methyCov"),
                   group2=c("N1methyCov","N2methyCov","N3methyCov","N4methyCov","N5methyCov","N6methyCov","N7methyCov"))
head(dmlTest)
# single = MulticoreParam(workers=1, progressbar=TRUE)
# #dmlTest.sm_1 <- DMLtest(BSobj, group1="mdf", group2="mp", equal.disp = T)
# dmlTest.sm <- DMLtest(BSobj, group1=c("C1methyCov","C2methyCov","C3methyCov","C4methyCov","C5methyCov","C6methyCov",
#                                       "C7methyCov","C8methyCov","C9methyCov","C10methyCov","C11methyCov","C12methyCov"),
#                       group2=c("N1methyCov","N2methyCov","N3methyCov","N4methyCov","N5methyCov","N6methyCov","N7methyCov"), smoothing=TRUE, BPPARAM=single)
dmlTest.sm <- DMLtest(BSobj, group1=c("C1methyCov","C2methyCov","C3methyCov","C4methyCov","C5methyCov","C6methyCov",
                                      "C7methyCov","C8methyCov","C9methyCov","C10methyCov","C11methyCov","C12methyCov"),
                      group2=c("N1methyCov","N2methyCov","N3methyCov","N4methyCov","N5methyCov","N6methyCov","N7methyCov"), smoothing=TRUE,ncores = 6)
head(dmlTest.sm)
save.image(file = "smothing_DMR.RData")
dmls <- callDML(dmlTest.sm, p.threshold=0.001)
dmls
dmls <- callDML(dmlTest, p.threshold=0.001)
dmls
dmls <- callDML(dmlTest.sm, delta=0.1, p.threshold=0.001)
# head(dmls2)
dmrs <- callDMR(dmlTest.sm, p.threshold=0.001)
# dmrs
dmrs2 <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
dmrs2
showOneDMR(dmrs[1,], BSobj)
# EAC VS ESCC
# library(DSS)
# library(bsseq)
#
# dir_path <- "/student1/MDD/Tabfile"
#
# files <- list.files(dir_path, pattern = "\\.tab$", full.names = TRUE)
#
data_list <- list()
sample_names <- c()
cancer_files <- files[grepl("^EAC_", basename(files))]
for (file in cancer_files) {
  sample_name <- gsub("^EAC_|\\.tab$", "", basename(file))
  sample_names <- c(sample_names, paste0("C", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
data_listC2 <- data_list
data_listC <- data_list[c(1:21)]
rm(data_list,sample_names)
#
data_list <- list()
sample_names <- c()
nonmalignant_files <- files[grepl("^GEJ_Nonmalignant_", basename(files))]
for (file in nonmalignant_files) {
  sample_name <- gsub("^GEJ_Nonmalignant_|\\.tab$", "", basename(file))
  sample_name <- gsub("_", "", sample_name)
  sample_names <- c(sample_names, paste0("N", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
data_listN <- data_list
rm(data_list,sample_names)
data_list <- c(data_listC,data_listN)
BSobj <- makeBSseqData(data_list, names(data_list))
# # BSobj <- BSobj[1:1000,]
#
# print(BSobj)
# [1] "C1methyCov" "C2methyCov" "C3methyCov" "C4methyCov" "C6methyCov"
# [6] "C1methyCov" "C2methyCov" "C3methyCov" "C4methyCov" "C5methyCov"
# [11] "C6methyCov" "C7methyCov" "N1methyCov" "N2methyCov" "N3methyCov"
# [16] "N4methyCov" "N5methyCov" "N6methyCov" "N7methyCov"
# dmlTest <- DMLtest(BSobj, group1=c("C1methyCov","C2methyCov","C3methyCov","C4methyCov","C5methyCov","C6methyCov",
#                                    "C7methyCov","C8methyCov","C9methyCov","C10methyCov","C11methyCov","C12methyCov"),
#                    group2=c("N1methyCov","N2methyCov","N3methyCov","N4methyCov","N5methyCov","N6methyCov","N7methyCov"))
# head(dmlTest)
# single = MulticoreParam(workers=1, progressbar=TRUE)
# #dmlTest.sm_1 <- DMLtest(BSobj, group1="mdf", group2="mp", equal.disp = T)
# dmlTest.sm <- DMLtest(BSobj, group1=c("C1methyCov","C2methyCov","C3methyCov","C4methyCov","C5methyCov","C6methyCov",
#                                       "C7methyCov","C8methyCov","C9methyCov","C10methyCov","C11methyCov","C12methyCov"),
#                       group2=c("N1methyCov","N2methyCov","N3methyCov","N4methyCov","N5methyCov","N6methyCov","N7methyCov"), smoothing=TRUE, BPPARAM=single)
#
# dmlTest.sm <- DMLtest(BSobj, group1=c("C1methyCov","C2methyCov","C3methyCov","C4methyCov","C5methyCov","C6methyCov",
#                                       "C7methyCov","C8methyCov","C9methyCov","C10methyCov","C11methyCov","C12methyCov"),
#                       group2=c("N1methyCov","N2methyCov","N3methyCov","N4methyCov","N5methyCov","N6methyCov","N7methyCov"), smoothing=TRUE,ncores = 6)
dmlTest.sm <- DMLtest(BSobj, group1=paste0("C",c(1:12),"methyCov"),
                      group2=paste0("C",c(13:33),"methyCov"), smoothing=TRUE,ncores = 4)
head(dmlTest.sm)
save.image(file = "smothing_DMR.RData")
dmls <- callDML(dmlTest.sm, p.threshold=0.001)
dmls
dmls <- callDML(dmlTest, p.threshold=0.001)
dmls
dmls <- callDML(dmlTest.sm, delta=0.1, p.threshold=0.001)
# head(dmls2)
dmrs <- callDMR(dmlTest.sm, p.threshold=0.001)
# dmrs
dmrs2 <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
dmrs2
showOneDMR(dmrs[1,], BSobj)
#ESCC vs EAC
data_list <- c(data_listE,data_listA,data_listB)
BSobj <- makeBSseqData(data_list, names(data_list))
dmlTest.sm <- DMLtest(BSobj, group1=paste0("ESCC",c(1:21),"methyCov"),
                      group2=paste0("EAC",c(1:12),"methyCov"), smoothing=TRUE,ncores = 2)
save.image("DMR_temp.RData")
#core_MDD
# setwd("D:/MDD_result/12/CV/result_topCV/")
dir = "D:/MDD_result/12/CV/result_topCV/SHARED_MDD/result/top_10_percent_cv/"
file = list.files(dir)
core_MDD <- data.frame()
for (i in 1:22) {
  temp_core_MDD <- read.table(paste0(dir,"chr",i,"_top_10_percent_cv.bed"))
  core_MDD <- rbind(core_MDD,temp_core_MDD)
}
core_MDD$V4 <- gsub("CV=","",core_MDD$V4)
core_MDD$V4 <- as.numeric(core_MDD$V4)
core_MDD <- core_MDD[,c(1:3)]
shared_MDD <- read.table("D:/MDD_result/1/ESCA_sharedMDDs.bed")
rownames(core_MDD) <- paste0(core_MDD$V1,"_",core_MDD$V2,"_",core_MDD$V3)
rownames(shared_MDD) <- paste0(shared_MDD$V1,"_",shared_MDD$V2,"_",shared_MDD$V3)
new_core_MDD <- shared_MDD[-which(rownames(shared_MDD) %in% rownames(core_MDD)),]
write.table(new_core_MDD,file = "D:/MDD_result/12/core_MDD.bed",col.names = F,row.names = F,quote = F,sep = "\t")
#
# dir = "D:/MDD_result/12/CV/result_topCV/SHARED_CMI/result/top_10_percent_cv/"
# file = list.files(dir)
# core_CMI <- data.frame()
# for (i in 1:22) {
#   temp_core_CMI <- read.table(paste0(dir,"chr",i,"_top_10_percent_cv.bed"))
#   core_CMI <- rbind(core_CMI,temp_core_CMI)
#
# }
# core_CMI$V4 <- gsub("CV=","",core_CMI$V4)
# core_CMI$V4 <- as.numeric(core_CMI$V4)
# core_CMI <- core_CMI[,c(1:3)]
# shared_CMI <- read.table("D:/MDD_result/1/intersect_shared_CMI_filter.bed")
#
# rownames(core_CMI) <- paste0(core_CMI$V1,"_",core_CMI$V2,"_",core_CMI$V3)
# rownames(shared_CMI) <- paste0(shared_CMI$V1,"_",shared_CMI$V2,"_",shared_CMI$V3)
# new_core_CMI <- shared_CMI[-which(rownames(shared_CMI) %in% rownames(core_CMI)),]
# write.table(new_core_CMI,"core_CMI.bed",col.names = F,row.names = F,quote = F,sep = "\t")
library(xgboost)
set.seed(2021)
core_MDD <- read.table("D:/MDD_result/12/core_MDD.bed")
methylation_file <- read.csv("D:/MDD_result/9/ESCA_MDD_meth/methylation_matrix.csv")
rownames(methylation_file) <- methylation_file$Region
methylation_file <- methylation_file[,-1]
rownames(core_MDD) <- paste0(core_MDD$V1,":",core_MDD$V2,"-",core_MDD$V3)
methylation_file_xgboost <- methylation_file[rownames(core_MDD),]
methylation_file_xgboost <- t(methylation_file_xgboost)
methylation_file_xgboost <- methylation_file_xgboost[order(rownames(methylation_file_xgboost)),]
sample_group <- c(rep(1,26),rep(0,5),rep(1,7),rep(0,7))
methylation_file_xgboost <- cbind(methylation_file_xgboost,sample_group)
methylation_file_xgboost <- as.data.frame(methylation_file_xgboost)
methylation_file_xgboost$sample_group <- as.factor(methylation_file_xgboost$sample_group)
# methylation_file_xgboost[,c(1:ncol(methylation_file_xgboost)-1)] <- as.numeric(methylation_file_xgboost[,c(1:ncol(methylation_file_xgboost)-1)] )
data_train <- methylation_file_xgboost[sample(1:45,0.8*45),]
data_train$sample_group <- as.factor(data_train$sample_group)
data_test <- methylation_file_xgboost[-which(rownames(methylation_file_xgboost) %in% rownames(data_train)),]
data_test$sample_group <- as.factor(data_test$sample_group)
library(caret)
library(future)
plan("multisession",workers=8)
set.seed(1)
rfeControl = rfeControl(functions = rfFuncs,
                        method = "cv",
                        saveDetails = T,
                        number = 10,
                        allowParallel = F
)
# set.seed(1)
# lmProfile <- rfe(methylation_file_xgboost[,c(1:ncol(methylation_file_xgboost)-1)], methylation_file_xgboost$sample_group,
#                  sizes = c(2:446),
#                  rfeControl = rfeControl
# )
set.seed(123)
rfe_results <- rfe(   x = data_train[, -ncol(data_train)],
                                   y = data_train$sample_group,
                                   sizes = c(2:25,30,40,50,60,70,80,90,100,200,300,400),
                                   rfeControl = rfeControl)
# ggplot(data = rfe_results) + theme_bw()
print(rfe_results)
plot(rfe_results, type = c("g", "o"))
selected_features <- predictors(rfe_results)
print(selected_features)
results_df <- rfe_results$results
results_df
library(ggplot2)
library(ggthemes)
ggplot(results_df,
       aes(x = Variables, y = Accuracy, size = AccuracySD, color = Kappa)) +
  geom_point(alpha = 0.6) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size(range = c(3, 15)) +
  labs(title = "Recursive Feature Elimination Results",
       x = "Number of Features",
       y = "Accuracy",
       caption = "Bubble size indicates Accuracy SD; Color indicates Kappa") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
    )
data_train <- data_train[,c(selected_features,"sample_group")]
data_test <- data_test[,c(selected_features,"sample_group")]
library(Matrix)
traindata1 <- as.matrix(data_train[,c(1:ncol(data_train)-1)])
traindata2 <- Matrix(traindata1,sparse = T)
train_y <- as.numeric(data_train[,ncol(data_train)])
traindata <- list(data=traindata2,label=train_y)
# x = model.matrix(sample_group~.,methylation_file_xgboost)[,-1]
# model_martix_train<-model.matrix(sample_group~.,methylation_file_xgboost)[,-1]
# data_train <- xgb.DMatrix(x , label =as.numeric(methylation_file_xgboost$sample_group))
# param <- list(objective = "reg:squarederror")
# xgb_model <- xgb.train(param, data_train, nrounds = 50)
# xgb_model
dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
testset1 <- as.matrix(data_test[,c(1:ncol(data_test)-1)])
testset2 <- Matrix(testset1,sparse=T)
test_y <- as.numeric(data_test[,ncol(data_test)])
testset <- list(data=testset2,label=test_y)
dtest <- xgb.DMatrix(data=testset$data,label=testset$label)
model_xgb <- xgboost(data=dtrain,booster='gbtree',max_depth=6,eta=0.5,objective='multi:softmax',num_class=3,nround=25)
pre <- predict(model_xgb,newdata=dtest)
library(caret)
xgb.cf <-caret::confusionMatrix(as.factor(pre),as.factor(test_y))
xgb.cf
table(testset$label,pre,dnn = c("",""))
require(pROC)
xgboost_roc<-roc(testset$label,as.numeric(pre))
plot(xgboost_roc,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),gird.col=c("green","red"),max.auc.polygon=TRUE,auc.polygon.col="skyblue",print.thres=TRUE,main='xgboost_model_ROC')
save(model_xgb,file = "D:/MDD_result/13/model.RData")
# library("breakDown")
# library(ggplot2)
# nobs <- model_xgb[1L, , drop = FALSE]
#
# explain_2 <- broken(HR_xgb_model, new_observation = nobs,
#                     data = model_martix_train)
# explain_2
#
# vd_xgb <- variable_importance(explainer_xgb, type = "raw")
# plot(vd_xgb)
# load("D:/MDD_result/4/cg_TCGA_select.RData")
shared_MDD <- read.table("D:/MDD_result/1/ESCA_sharedMDDs.bed")
load("D:/MDD_result/4/cg_TCGA_select.RData")
core_MDD <- read.table("D:/MDD_result/12/core_MDD.bed")
names(shared_list) <- paste0(shared_MDD$V1,":",shared_MDD$V2,"-",shared_MDD$V3)
core_list <- shared_list[which(names(shared_list) %in% paste0(core_MDD$V1,":",core_MDD$V2,"-",core_MDD$V3))]
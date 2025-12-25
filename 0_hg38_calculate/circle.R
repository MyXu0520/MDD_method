# MDD_CMI_circle
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
chr.exclude=c("chrX", "chrY")
tracks.inside <- 8
tracks.outside <- 1
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)
RCircos.List.Plot.Parameters()
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
# data(RCircos.Gene.Label.Data)
shared_MDD <- read.table("D:/MDD_result/1/ESCA_sharedMDDs.bed")
colnames(shared_MDD) <- c("Chrom","Start","End")
shared_MDD$group <- rep("shared_MDD",nrow(shared_MDD))
# shared_MDD$End <- shared_MDD$Start+1
EAC_MDD <- read.table("D:/MDD_result/1/EAC_specificMDDs.bed")
colnames(EAC_MDD) <- c("Chrom","Start","End")
EAC_MDD$group <- rep("EAC_SPEC",nrow(EAC_MDD))
# EAC_MDD$End <- EAC_MDD$Start+1
ESCC_MDD <- read.table("D:/MDD_result/1/ESCC_specificMDDs.bed")
colnames(ESCC_MDD) <- c("Chrom","Start","End")
ESCC_MDD$group <- rep("ESCC_SPEC",nrow(ESCC_MDD))
# ESCC_MDD$End <- ESCC_MDD$Start+1
CMI <- read.table("D:/MDD_result/10/intersect_shared_CMI_filtered_unique.bed")
colnames(CMI) <- c("Chrom","Start","End")
CMI$group <- rep("CMI",nrow(CMI))
# CMI$End <_ CMI$Start+1
anno_data <- rbind(shared_MDD,EAC_MDD,ESCC_MDD,CMI)
scale_values <- function (x){(x-min(x))/(max(x)-min(x))}
temp <- sample(1:nrow(anno_data),round(nrow(anno_data)/50),replace = F)
Temp <- temp
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(anno_data[temp,], track.num, side)
name.col <- 4
track.num <- 2
RCircos.Gene.Name.Plot(anno_data[temp,], name.col,track.num, side)
#
#For circle data
setwd("D:/MDD_result/9/")
#ESCC_SPEC
file <- list.files("ESCA_ESCC_MDD_meth/")
file <- file[-length(file)]
temp <- list()
for (i in 1:length(file)) {
  temp[[i]] <- read.table(paste0("ESCA_ESCC_MDD_meth/",file[i]))
}
names(temp) <- file
ESCC_SPEC_sample <- temp[[1]][c(1:3)]
colnames(ESCC_SPEC_sample) <- c("chrom","start","end")
ESCC_SPEC_sample$type <- rep("ESCC_SPEC",nrow(ESCC_SPEC_sample))
methylation_bind <- matrix(nrow = length(temp[[1]]$V1))
methylation_bind <- as.data.frame(methylation_bind)
for (i in 1:length(temp)) {
  methylation_bind <- cbind(methylation_bind,temp[[i]]$V4)
}
methylation_bind <- methylation_bind[,-1]
colnames(methylation_bind) <- file
ESCC_SPEC_sample <- cbind(ESCC_SPEC_sample,methylation_bind)
rownames(ESCC_SPEC_sample) <- paste0(ESCC_SPEC_sample$chrom,"_",ESCC_SPEC_sample$start,"_",ESCC_SPEC_sample$end)
ESCC_SPEC_MDD_mean_methylation <- c()
for (i in 1:nrow(temp[[1]])) {
  temp_mean <- 0
  for (j in 1:length(file)) {
    temp_mean <- temp_mean+temp[[j]]$V4[i]
  }
  temp_mean <- temp_mean/length(file)
  ESCC_SPEC_MDD_mean_methylation[i] <- temp_mean
}
ESCC_SPEC_MDD_region <- temp[[1]][c(1:3)]
ESCC_SPEC_MDD_region$mean <- ESCC_SPEC_MDD_mean_methylation
ESCC_SPEC_sample$mean <- ESCC_SPEC_MDD_region$mean
# Shared_MDD
file <- list.files("ESCA_MDD_meth/")
file <- file[-length(file)]
temp <- list()
for (i in 1:length(file)) {
  temp[[i]] <- read.table(paste0("ESCA_MDD_meth/",file[i]))
}
names(temp) <- file
shared_MDD_sample <- temp[[1]][c(1:3)]
colnames(shared_MDD_sample) <- c("chrom","start","end")
shared_MDD_sample$type <- rep("shared_MDD",nrow(shared_MDD_sample))
methylation_bind <- matrix(nrow = length(temp[[1]]$V1))
methylation_bind <- as.data.frame(methylation_bind)
for (i in 1:length(temp)) {
  methylation_bind <- cbind(methylation_bind,temp[[i]]$V4)
}
methylation_bind <- methylation_bind[,-1]
colnames(methylation_bind) <- file
shared_MDD_sample <- cbind(shared_MDD_sample,methylation_bind)
rownames(shared_MDD_sample) <- paste0(shared_MDD_sample$chrom,"_",shared_MDD_sample$start,"_",shared_MDD_sample$end)
shared_MDD_mean_methylation <- c()
for (i in 1:nrow(temp[[1]])) {
  temp_mean <- 0
  for (j in 1:length(file)) {
    temp_mean <- temp_mean+temp[[j]]$V4[i]
  }
  temp_mean <- temp_mean/length(file)
  shared_MDD_mean_methylation[i] <- temp_mean
}
# shared_MDD_region <- paste0(temp[[1]]$V1,"_",temp[[1]]$V2,"_",temp[[1]]$V3)
# shared_MDD_region <- as.data.frame(shared_MDD_region)
# shared_MDD_region <-
shared_MDD_region <- temp[[1]][c(1:3)]
shared_MDD_region$mean <- shared_MDD_mean_methylation
shared_MDD_sample$mean <- shared_MDD_region$mean
# EAC_SPEC_MDD
#CMI
file <- list.files("ESCA_CMI_meth/")
file <- file[-length(file)]
temp <- list()
for (i in 1:length(file)) {
  temp[[i]] <- read.table(paste0("ESCA_CMI_meth/",file[i]))
}
names(temp) <- file
CMI_sample <- temp[[1]][c(1:3)]
colnames(CMI_sample) <- c("chrom","start","end")
CMI_sample$type <- rep("CMI",nrow(CMI_sample))
methylation_bind <- matrix(nrow = length(temp[[1]]$V1))
methylation_bind <- as.data.frame(methylation_bind)
for (i in 1:length(temp)) {
  methylation_bind <- cbind(methylation_bind,temp[[i]]$V4)
}
methylation_bind <- methylation_bind[,-1]
colnames(methylation_bind) <- file
CMI_sample <- cbind(CMI_sample,methylation_bind)
rownames(CMI_sample) <- paste0(CMI_sample$chrom,"_",CMI_sample$start,"_",CMI_sample$end)
CMI_mean_methylation <- c()
for (i in 1:nrow(temp[[1]])) {
  temp_mean <- 0
  for (j in 1:length(file)) {
    temp_mean <- temp_mean+temp[[j]]$V4[i]
  }
  temp_mean <- temp_mean/length(file)
  CMI_mean_methylation[i] <- temp_mean
}
# shared_MDD_region <- paste0(temp[[1]]$V1,"_",temp[[1]]$V2,"_",temp[[1]]$V3)
# shared_MDD_region <- as.data.frame(shared_MDD_region)
# # shared_MDD_region <-
CMI_region <- temp[[1]][c(1:3)]
CMI_region$mean <- CMI_mean_methylation
CMI_sample$mean <- CMI_region$mean
temp_anno <- anno_data[Temp,]
rownames(temp_anno) <- paste0(temp_anno$Chrom,"_",temp_anno$Start,"_",temp_anno$End)
heatmap_anno <- matrix(ncol = ncol(shared_MDD_sample))
heatmap_anno  <-  as.data.frame(heatmap_anno)
for (i in 1:nrow(temp_anno)) {
  if(temp_anno$group[i] == "shared_MDD"){heatmap_anno[i,] <- shared_MDD_sample[rownames(temp_anno)[i],]}
  # else if(temp_anno$group[i] == "EAC_SPEC"){heatmap_anno[i,] <- EAC_SPEC_sample[rownames(temp_anno)[i],]}
  else if(temp_anno$group[i] == "ESCC_SPEC"){heatmap_anno[i,] <- ESCC_SPEC_sample[rownames(temp_anno)[i],]}
  else{heatmap_anno[i,] <- CMI_sample[rownames(temp_anno)[i],]}
}
heatmap_anno <- na.omit(heatmap_anno)
colnames(heatmap_anno) <- colnames(shared_MDD_sample)
ESCC_mean <- 0
temp_ESCC_mean <- c()
for (i in 1:nrow(heatmap_anno)) {
  ESCC_mean <- ESCC_mean + sum(heatmap_anno[i,c(10:30)])
  ESCC_mean <- ESCC_mean/length(c(10:30))
  temp_ESCC_mean <- c(temp_ESCC_mean,ESCC_mean)
}
heatmap_anno$ESCC_mean <- temp_ESCC_mean
temp_EAC_mean <- c()
EAC_mean <- 0
for (i in 1:nrow(heatmap_anno)) {
  EAC_mean <- EAC_mean + sum(heatmap_anno[i,c(5:9,36:42)])
  EAC_mean <- EAC_mean/length(c(5:9,36:42))
  temp_EAC_mean <- c(temp_EAC_mean,EAC_mean)
}
heatmap_anno$EAC_mean <- temp_EAC_mean
temp_ESCC_N__mean <- c()
ESCC_N__mean <- 0
for (i in 1:nrow(heatmap_anno)) {
  ESCC_N__mean <- ESCC_N__mean + sum(heatmap_anno[i,c(31:35)])
  ESCC_N__mean <- ESCC_N__mean/length(c(31:35))
  temp_ESCC_N__mean <- c(temp_ESCC_N__mean,ESCC_N__mean)
}
heatmap_anno$ESCC_N__mean<- temp_ESCC_N__mean
temp_EAC_N_mean <- c()
EAC_N_mean  <- 0
for (i in 1:nrow(heatmap_anno)) {
  EAC_N_mean  <- EAC_N_mean  + sum(heatmap_anno[i,c(43:49)])
  EAC_N_mean  <- EAC_N_mean /length(c(43:49))
  temp_EAC_N_mean <- c(temp_EAC_N_mean,EAC_N_mean)
}
heatmap_anno$EAC_N_mean <-temp_EAC_N_mean
# scale_values <- function (x){(x-min(x))/(max(x)-min(x))}
# heatmap_anno$mean <- scale_values(heatmap_anno$mean)
# heatmap_anno$ESCC_mean <- scale_values(heatmap_anno$ESCC_mean)
# heatmap_anno$EAC_mean <- scale_values(heatmap_anno$EAC_mean)
# heatmap_anno$ESCC_N__mean <- scale_values(heatmap_anno$ESCC_N__mean)
# heatmap_anno$EAC_N_mean <- scale_values(heatmap_anno$EAC_N_mean)
# data.col <- ncol(heatmap_anno)
# track.num <- 5
# side <- "in"
# RCircos.Heatmap.Plot(heatmap_anno, data.col, track.num, side)
data.col <- 51
track.num <- 5
side <- "in"
RCircos.Heatmap.Plot(heatmap_anno, data.col, track.num, side)
data.col <- 52
track.num <- 6
side <- "in"
RCircos.Heatmap.Plot(heatmap_anno, data.col, track.num, side)
data.col <- 53
track.num <- 8
side <- "in"
RCircos.Heatmap.Plot(heatmap_anno, data.col, track.num, side)
data.col <- 54
track.num <- 9
side <- "in"
RCircos.Heatmap.Plot(heatmap_anno, data.col, track.num, side)
data(RCircos.Link.Data)
RCircos.Link.Data <- RCircos.Link.Data[-which(RCircos.Link.Data$Chromosome %in% c("chrX","chrY")),]
track.num <- 11
RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE)
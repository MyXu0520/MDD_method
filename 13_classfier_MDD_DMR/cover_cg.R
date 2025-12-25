# MDD
#TCGA-ESCA
library(readr)
# TCGA_ESCA <- read.table("E:/33cancer/33cancers/knn_data/ESCA.methylation450.knn.txt")
TCGA_ESCA <- read.table("I:/TCGA/ESCA.methylation450.knn.txt")
clinical_ESCA <- read_tsv("E:/33cancer/33cancers/clinical.cart.2024-11-19/clinical.tsv")
# cg_450K <- read.csv("E:/33cancer/33cancers/humanmethylation450.csv",sep = ",",header = T)
cg_450K <- read.csv("H:/zhi/sun_Meth/data/humanmethylation450.csv",sep = ",",header = T)
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
core_MDD <- read.table("D:/MDD_result/12/core_MDD.bed")
colnames(core_MDD) <- c("Chrom","Start","End")
#EAC_spec_MDD
core_CMI <- read.table("D:/MDD_result/12/core_CMI.bed")
colnames(core_CMI) <- c("Chrom","Start","End")
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
core_CMI_list <- list()
for (i in 1:nrow(core_CMI)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == core_CMI[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= core_CMI[i,2] && temp[j,3] <= core_CMI[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  core_CMI_list[[i]] <- cg_list
}
save(core_MDD_list,core_CMI_list,file = "D:/MDD_result/4/cg_TCGA_core_select.RData")
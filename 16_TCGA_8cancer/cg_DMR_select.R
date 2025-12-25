#TCGA-ESCA
library(readr)
# TCGA_ESCA <- read.table("E:/33cancer/33cancers/knn_data/ESCA.methylation450.knn.txt")
# clinical_ESCA <- read_tsv("E:/33cancer/33cancers/clinical.cart.2024-11-19/clinical.tsv")
cg_450K <- read.csv("E:/33cancer/33cancers/humanmethylation450.csv",sep = ",",header = T)
cg_450K_locus <- cg_450K[,c(1,12,13,18,22,24,25,26)]
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$Probe_SNPs == ""),]
cg_450K_locus <- cg_450K_locus[which(cg_450K_locus$CHR %in% c(as.character(1:22),"X","Y")),]
cg_450K_locus$CHR <- paste0("chr",cg_450K_locus$CHR)
#ESCA_shared_DMR
ESCA_shared_DMR <- read.table("D:/MDD_result/18/EAC_ESCC_DMR_bed.bed")
colnames(ESCA_shared_DMR) <- c("Chrom","Start","End")
#EAC_spec_DMR
EAC_spec_DMR <- read.table("D:/MDD_result/18/EAC_SPEC_DMR.bed")
colnames(EAC_spec_DMR) <- c("Chrom","Start","End")
#ESCC_spec_DMR
ESCC_spec_DMR <- read.table("D:/MDD_result/18/ESCC_SPEC_DMR.bed")
colnames(ESCC_spec_DMR) <- c("Chrom","Start","End")
DMR_shared_list <- list()
for (i in 1:nrow(ESCA_shared_DMR)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == ESCA_shared_DMR[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= ESCA_shared_DMR[i,2] && temp[j,3] <= ESCA_shared_DMR[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  DMR_shared_list[[i]] <- cg_list
}
DMR_EAC_spec_list <- list()
for (i in 1:nrow(EAC_spec_DMR)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == EAC_spec_DMR[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= EAC_spec_DMR[i,2] && temp[j,3] <= EAC_spec_DMR[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  DMR_EAC_spec_list[[i]] <- cg_list
}
DMR_ESCC_spec_list <- list()
for (i in 1:nrow(ESCC_spec_DMR)) {
  temp <- cg_450K_locus[which(cg_450K_locus$CHR == ESCC_spec_DMR[i,1]),]
  cg_list <- c()
  for (j in 1:nrow(temp)) {
    if(temp[j,3] >= ESCC_spec_DMR[i,2] && temp[j,3] <= ESCC_spec_DMR[i,3]){
      cg_list <- c(cg_list,temp$IlmnID[j])
    }
    else{cg_list <- cg_list}
  }
  DMR_ESCC_spec_list[[i]] <- cg_list
}
save(DMR_shared_list,DMR_EAC_spec_list,DMR_ESCC_spec_list,file = "D:/DMR_result/18/cg_TCGA_select_DMR.RData")
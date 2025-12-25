# library(TCGAbiolinks)
#
# query <- GDCquery(
#   project = "TCGA-ESCA",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
#
# GDCdownload(query)
# GDCprepare(query, save = T,save.filename = "D:/MDD_result/16/TCGA-ESCA_SNP.Rdata")
library(maftools)
load("D:/MDD_result/16/TCGA-ESCA_SNP.Rdata")
maf_data <- data[,c(6,7,8)]
write.table(maf_data,"D:/MDD_result/16/maf_data.bed",sep = "\t",quote = F,col.names = F,row.names = F)
# temp <- data
# rownames(temp) <- paste0(temp$Chromosome,"_",temp$Start_Position,"_",temp$End_Position)
shared_MDD <- read.table("D:/MDD_result/16/maf_shared.bed",sep = "\t")
maf_shared <- as.data.frame(matrix(ncol = ncol(data)))
colnames(maf_shared) <- colnames(data)
for (i in 1:length(shared_MDD$V1)) {
  temp <- data[which(data$Chromosome == shared_MDD$V1[i]),]
  temp_point <- temp[which(temp$Start_Position == shared_MDD$V2[i]),]
  maf_shared <- rbind(maf_shared,temp_point)
}
maf.shared <- maf_shared[-1,]
maf1 <- read.maf(maf.shared)
pdf("D:/MDD_result/16/shared_MAF2.pdf")
plotmafSummary(maf = maf1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,color = vc_cols)
dev.off()
pdf("D:/MDD_result/16/shared_MAF.pdf")
oncoplot(maf = maf1, colors = vc_cols, top = 20)
dev.off()
#EAC_SPEC
EAC_SPEC_MDD <- read.table("D:/MDD_result/16/maf_EAC_SPEC.bed",sep = "\t")
maf_EAC_SPEC <- as.data.frame(matrix(ncol = ncol(data)))
colnames(maf_EAC_SPEC) <- colnames(data)
for (i in 1:length(EAC_SPEC_MDD$V1)) {
  temp <- data[which(data$Chromosome == EAC_SPEC_MDD$V1[i]),]
  temp_point <- temp[which(temp$Start_Position == EAC_SPEC_MDD$V2[i]),]
  maf_EAC_SPEC <- rbind(maf_EAC_SPEC,temp_point)
}
maf.EAC_SPEC <- maf_EAC_SPEC[-1,]
maf2 <- read.maf(maf.EAC_SPEC)
plotmafSummary(maf = maf2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,color = vc_cols)
pdf("D:/MDD_result/16/EAC_MAF.pdf")
oncoplot(maf = maf2, colors = vc_cols, top = 20)
dev.off()
#ESCC_SPEC
ESCC_SPEC_MDD <- read.table("D:/MDD_result/16/maf_ESCC_SPEC.bed",sep = "\t")
maf_ESCC_SPEC <- as.data.frame(matrix(ncol = ncol(data)))
colnames(maf_ESCC_SPEC) <- colnames(data)
for (i in 1:length(ESCC_SPEC_MDD$V1)) {
  temp <- data[which(data$Chromosome == ESCC_SPEC_MDD$V1[i]),]
  temp_point <- temp[which(temp$Start_Position == ESCC_SPEC_MDD$V2[i]),]
  maf_ESCC_SPEC <- rbind(maf_ESCC_SPEC,temp_point)
}
maf.ESCC_SPEC <- maf_ESCC_SPEC[-1,]
maf3 <- read.maf(maf.ESCC_SPEC)
plotmafSummary(maf = maf3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,color = vc_cols)
pdf("D:/MDD_result/16/ESCC_MAF.pdf")
oncoplot(maf = maf3, colors = vc_cols, top = 20)
dev.off()
maf.esca <- data
# maf <- read.maf(maf.esca,clinical = "D:/MDD_result/16/clinical.cart.2024-12-23/clinical.tsv")
# esca.pfam <- pfamDomains(maf=maf.esca, AACol="HGVSp_Short", top=10)
# mafSurvival(maf=maf.esca, genes=c("TP53", "TTN", "MUC16"), time="days_to_last_follow_up", Status="days_to_death", isTCGA=TRUE)
# vafclust <- inferHeterogeneity(maf=maf.esca, tsb="TCGA-R6-A8WG")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE)
luad.tnm <- trinucleotideMatrix(maf=luad, ref_genome="BSgenome.Hsapiens.UCSC.hg38")
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
print(vc_cols)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,color = vc_cols)
oncoplot(maf = maf, colors = vc_cols, top = 20)
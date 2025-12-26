###Code for WGBS project######
options(scipen = 20)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(plyr)
library(DelayedMatrixStats)
library(tibble)
library(valr)
library(ggbeeswarm)
library(ggpubr)
library(VennDiagram)
library(dmrseq)
library(bsseq)
library(gridExtra)
library(readr)
library(vipor)
library(dplyr)
library(formattable)
library(stringr)
library(reshape2)
library(BSgenome)
library(MutationalPatterns)
################################################meta_info################################################
EAC_Tumor_sampleList=c(paste0("EAC_", c(1:4,6)), paste0("GEJ_", 1:7))
ESCC_Tumor_sampleList=paste0("ESCC_", c(1:17,19:22))
EAC_Nonmalignant_sampleList=paste0("GEJ_Nonmalignant_", 1:7)
ESCC_Nonmalignant_sampleList=paste0("ESCC_Nonmalignant_", c(1:3, "53F", "54M"))
types=c(rep("ESCC_Tumor",21),rep("ESCC_Nonmalignant", 5), rep("EAC/GEJ_Tumor",12), rep("GEJ_Nonmalignantl",7))

fileList=c(ESCC_Tumor_sampleList, ESCC_Nonmalignant_sampleList, EAC_Tumor_sampleList, EAC_Nonmalignant_sampleList)
chrListTarget=paste("chr",1:22,sep="")

annotation_row=data.frame(Type=c(rep("ESCC_Tumor", 21), rep("ESCC_Nonmalignant", 5), rep("EAC/GEJ_Tumor", 12), rep("GEJ_Nonmalignant", 7)), stringsAsFactors = F)
rownames(annotation_row)=fileList
annotation_row$Type=factor(annotation_row$Type, levels=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor"))
annotation_row=annotation_row[order(annotation_row$Type), ,drop=F]
annotation_col=annotation_row
annotation_col$Sample=rownames(annotation_col)

ann_colors = list(Type = c("#F7CE46", "#75FBFD", "#EA3323", "#0000F5"))
names(ann_colors[[1]])=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor")

my_theme=theme_classic()+theme(axis.text = element_text(color="black", size=10), 
                               axis.title.y=element_text(color="black", size=12),
                               plot.title=element_text(hjust = 0.5, face="bold", size=14))

################################################PMD calling################################################
###ESCA samples
##Call PMDs by MMSeekR
library(MMSeekR.data)
library(MMSeekR)

changeToTabFile=function(sampleIndex){
  data=read_tsv(paste0("Data/ESCA_rmblackList_cov5/", sampleIndex, ".all.cov5.sorted.rmblackList.bed"), col_names = F)
  data=data[,c(1:3,5,4)]
  colnames(data)=c("Chr","Start","End","T","beta")
  data$M=round(data$T*data$beta,0)
  chrListTarget=paste("chr",1:22,sep="")
  data=data[data$Chr%in%chrListTarget,]
  data$Pos=data$Start+1
  data=data[,c(1,7,4,6)]
  write.table(data, file=paste0("Data/MMSeekR_PMDs/TabFile/", sampleIndex, "methyCov.tab"), sep="\t", row.names = F, col.names = F, quote=F)
}
for(fileIndex in fileList){
  print(fileIndex)
  changeToTabFile(fileIndex)
}
methFileList=dir("Data/MMSeekR_PMDs/TabFile/", pattern=".tab", full.names = T)
for(methFile in methFileList){
  fileIndex=gsub(".tab", "", basename(methFile))
  runPMDs("hg38", methFile, paste0("Data/MMSeekR_PMDs/PMDs/", fileIndex))
}


###get ESCA_union_PMDs
cancerTypeCommonPMDs=function(PMDfile, cutoff){
  data=read.table(PMDfile,sep="\t", stringsAsFactors = F,header=T)
  data=data[data$num>=cutoff,]
  outputFile=gsub("_PMDs_multiinter.bed", "_commonPMDs.bed",PMDfile)
  write.table(data[,1:3], file=outputFile, sep="\t", row.names = F, col.names=F, quote=F)
  temp=read_bed(outputFile)
  temp=bed_merge(temp)
  temp=as.data.frame(temp)
  ###merge the region and require the region longer than 2kb
  temp=temp[temp$end-temp$start>2000,]
  print(sum(temp$end-temp$start))
  write.table(temp, file=outputFile, sep="\t", row.names = F, col.names=F, quote=F)
}
cancerTypeCommonPMDs("Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", 8)
cancerTypeCommonPMDs("Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", 14)

###get EAC and ESCC specific PMDs, shared PMD and HMDs
specificPMDs=function(sharedPMDFile, multiinterFile, cutoff, outputFile){
  sharedPMDs=read_bed(sharedPMDFile)
  sharedPMDs=bed_sort(sharedPMDs)
  sharedPMDs=bed_merge(sharedPMDs)
  sharedPMDs_length=sharedPMDs$end-sharedPMDs$start
  sharedPMDs_length_total=sum(sharedPMDs_length)
  print(sharedPMDs_length_total)
  multiinterPMDsAll=read.table(multiinterFile,sep="\t", stringsAsFactors = F,header=T)
  sampleLength=length(colnames(multiinterPMDsAll)[-1:-5])
  multiinterPMDs=multiinterPMDsAll[,1:4]
  multiinterPMDs=multiinterPMDs[multiinterPMDs$num>=(cutoff*sampleLength),]
  multiinterPMDs=tibble(chrom=multiinterPMDs$chrom, start=multiinterPMDs$start, end=multiinterPMDs$end)
  multiinterPMDs=bed_sort(multiinterPMDs)
  multiinterPMDs=bed_merge(multiinterPMDs)
  specificPMDs=bed_subtract(sharedPMDs, multiinterPMDs)
  specificPMDs=bed_sort(specificPMDs)
  specificPMDs=bed_merge(specificPMDs)
  specificPMDs$length=specificPMDs$end-specificPMDs$start
  specificPMDs=specificPMDs[specificPMDs$length>2000,]
  print(sum(specificPMDs$length))
  write.table(specificPMDs[,1:3], file=outputFile, row.names = F, col.names = F, sep="\t", quote=F)
}
#EAC only PMDs
specificPMDs("Data/MMSeekR_PMDs/EAC_commonPMDs.bed", "Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", (1/3), "Data/MMSeekR_PMDs/EAC_specificPMDs.bed")
#ESCC only PMDs
specificPMDs("Data/MMSeekR_PMDs/ESCC_commonPMDs.bed", "Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", (1/3), "Data/MMSeekR_PMDs/ESCC_specificPMDs.bed")
#ESCA sharedPMDs
sharedESCARegions=function(EACPMDFile, ESCCPMDFile, outputFile, type){
  EAC_PMDs=read_bed(EACPMDFile)
  EAC_PMDs_total=as.numeric(sum(EAC_PMDs$end-EAC_PMDs$start))
  ESCC_PMDs=read_bed(ESCCPMDFile)
  ESCC_PMDs_total=as.numeric(sum(ESCC_PMDs$end-ESCC_PMDs$start))
  ESCA_sharedPMDs=bed_intersect(EAC_PMDs, ESCC_PMDs)
  ESCA_sharedPMDs$start=rowMaxs(as.matrix(ESCA_sharedPMDs[,colnames(ESCA_sharedPMDs)%in%c("start.x","start.y")]))
  ESCA_sharedPMDs$end=rowMins(as.matrix(ESCA_sharedPMDs[,colnames(ESCA_sharedPMDs)%in%c("end.x","end.y")]))
  ESCA_sharedPMDs$length=ESCA_sharedPMDs$end-ESCA_sharedPMDs$start
  ESCA_sharedPMDs=ESCA_sharedPMDs[,colnames(ESCA_sharedPMDs)%in%c("chrom", "start", "end")]
  ESCA_sharedPMDs=bed_merge(bed_sort(ESCA_sharedPMDs))
  ESCA_sharedPMDs$length=ESCA_sharedPMDs$end-ESCA_sharedPMDs$start
  ESCA_sharedPMDs=ESCA_sharedPMDs[ESCA_sharedPMDs$length>2000,]
  ESCA_sharedPMDs_total=sum(ESCA_sharedPMDs$length)
  print(EAC_PMDs_total)
  print(ESCC_PMDs_total)
  print(ESCA_sharedPMDs_total)
  write.table(ESCA_sharedPMDs[,1:3], file=outputFile, sep="\t", row.names = F, col.names=F, quote=F)
}
sharedESCARegions("Data/MMSeekR_PMDs/EAC_commonPMDs.bed", "Data/MMSeekR_PMDs/ESCC_commonPMDs.bed","Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", "Tumor")


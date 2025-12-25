# library(TCGAbiolinks)
# ESCAmethy <- GDCquery(
#   project = "TCGA-ESCA",
#   data.category = "DNA Methylation",
#   data.type = "Masked Intensities",
#   platform = "Illumina Human Methylation 450", # Illumina Human Methylation 450
#   legacy = FALSE
# )
# GDCdownload(ESCAmethy)
# GDCprepare(ESCAmethy,save = T,save.filename="ESCA_METHY_idat.Rdata")
#EAC
library(readr)
setwd("D:/TCGA-ESCA-IDAT/ESCC/IDAT")
#SampleSheet
Sample_sheet_Origin <- read_tsv("D:/TCGA-ESCA-IDAT/ESCC/gdc_sample_sheet.2024-11-21.tsv")
temp <- strsplit(Sample_sheet_Origin$`File Name`,"_")
temp_df <- data.frame()
for (i in 1:length(temp)) {
  temp_df[i,c(1:3)] <- temp[[i]][c(1:3)]
}
Sample_Sheet <- paste0(temp_df$V1,"_",temp_df$V2)
Sample_Sheet <- as.data.frame(Sample_Sheet)
colnames(Sample_Sheet) <- "Sample_Name"
Sample_Sheet$Sample_Group <- Sample_sheet_Origin$`Sample Type`
Sample_Sheet$Sentrix_ID <- substr(Sample_sheet_Origin$`Case ID`,1,7)
Sample_Sheet$Sentrix_Position <- substr(Sample_sheet_Origin$`Sample ID`,9,16)
Sample_Sheet$Sample_Name <- paste0(Sample_Sheet$Sentrix_ID,"_",Sample_Sheet$Sentrix_Position)
#
#
# folder <- "D:/TCGA-ESCA-IDAT/EAC/IDAT"
# files <- list.files(folder)
# files <- as.data.frame(files)
# rownames(files) <- files$files
Sample_sheet_filter <- Sample_sheet_Origin[-which(substr(Sample_sheet_Origin$`Sample ID`,14,16) == "01B"),]
Sample_sheet_filter <- Sample_sheet_filter[-which(substr(Sample_sheet_filter$`Sample ID`,14,16) == "11A"),]
Sample_Sheet$select <- paste0(Sample_Sheet$Sentrix_ID,"-",Sample_Sheet$Sentrix_Position)
Sample_Sheet <- Sample_Sheet[which(Sample_Sheet$select %in% Sample_sheet_filter$`Sample ID`),]
#
# new_file_names <- unique(Sample_sheet_filter$`Sample ID`)
# setwd("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/")
# i=1
# while(i<=78){
#   dir.create(new_file_names[i])
#   i=i+1
# }
# setwd("D:/TCGA-ESCA-IDAT/EAC/IDAT/")
# for(i in 1:length(new_file_names)){
#   new_file_folder <- Sample_sheet_filter[which(Sample_sheet_filter$`Sample ID` == new_file_names[i]),]
#   of1 <- file.path(paste0("D:/TCGA-ESCA-IDAT/EAC/IDAT/",new_file_folder$`File ID`[1]),new_file_folder$`File Name`[1])
#   nf1 <- file.path(paste0("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/",new_file_names[i]),new_file_folder$`File Name`[1])
#   file.copy(of1,nf1)
#   of2 <- file.path(paste0("D:/TCGA-ESCA-IDAT/EAC/IDAT/",new_file_folder$`File ID`[2]),new_file_folder$`File Name`[2])
#   nf2 <- file.path(paste0("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/",new_file_names[i]),new_file_folder$`File Name`[2])
#   file.copy(of2,nf2)
#   # file.copy(paste0("D:/TCGA-ESCA-IDAT/EAC/IDAT/",new_file_folder$`File ID`[1],"/",new_file_folder$`File Name`[1]), paste0("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/",new_file_folder$`File Name`[1]))
#   # file.copy(paste0("D:/TCGA-ESCA-IDAT/EAC/IDAT/",new_file_folder$`File ID`[2],"/",new_file_folder$`File Name`[2]), paste0("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/",new_file_folder$`File Name`[2]))
# }
# paste0("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/",new_file_folder$`Sample ID`[2],"_",strsplit(new_file_folder$`File Name`,"_")[[2]][3])
# files <- files[Sample_sheet_filter$`File ID`,]
# for(f in 1:length(files)){
#   newname <- sub(as.character(files[f]),Sample_sheet_filter$`Sample ID`[f],as.character(files[f]))
#   print(newname)
#   file.rename(files[f],newname)
# }
# dir()
# Sample_sheet <- Sample_Sheet$Sample_Name
# folder <- setwd("D:/TCGA-ESCA-IDAT/EAC/IDAT")
folder <- "D:/TCGA-ESCA-IDAT/ESCC/IDAT"
files <- list.files(folder)
files <- as.data.frame(files)
rownames(files) <- files$files
files <- files[Sample_sheet_filter$`File Name`,]
for (i in 1:nrow(Sample_Sheet)) {
  newfile <- files[which(files == Sample_sheet_filter$`File Name`[i])]
  temp <- strsplit(as.character(newfile),"_")
  # paste0(temp[[1]][1],"_",temp[[1]][2])
  newname <- sub(paste0(temp[[1]][1],"_",temp[[1]][2]),Sample_Sheet$Sample_Name[i],as.character(newfile))
  print(newname)
  file.rename(newfile,newname)
}
# for(i in 1:length(files)){
#   newfile <- list.files(paste0("D:/TCGA-ESCA-IDAT/EAC/New_IDAT/",Sample_sheet_filter$`Sample ID`[i],"/"))
#   # newname <- sub(newfile[1],paste0(Sample_Sheet$Sentrix_ID,"_",Sample_Sheet$Sentrix_Position),newfile[1])
#   newname <- paste0(Sample_sheet_filter$`Sample ID`[i],"_",Sample_sheet_Origin$`Sample ID`[i])
#   newfile <- newfile[1]
#   newfile.rename(newfile,newname)
# }
# paste0(Sample_Sheet$Sentrix_ID[f],"_",Sample_Sheet$Sentrix_Position[f])
# Sample_Sheet$Sentrix_ID <- gsub("-","_",Sample_Sheet$Sentrix_ID)
# Sample_Sheet$Sentrix_Position<- gsub("-","_",Sample_Sheet$Sentrix_Position)
# Sample_Sheet$Sample_Name <- strsplit(list.files(),"\\.")
# for (i in 1:nrow(Sample_Sheet)) {
#   Sample_Sheet$Sample_Name[i] <- substr(strsplit(list.files(),"\\.")[[i]][1],1,16)
# }
# for (i in 1:nrow(Sample_Sheet)) {
#   Sample_Sheet$Sample_Name[i] <- paste0(Sample_Sheet$Sentrix_ID[i],"_",Sample_Sheet$Sentrix_Position[i])
# }
Sample_Sheet1 <- Sample_Sheet[,-ncol(Sample_Sheet)]
Sample_Sheet2 <- unique(Sample_Sheet1)
write.csv(Sample_Sheet2,file = "C:/Users/MyXu/Desktop/result/3/Sample_Sheet_ESCC.csv",row.names = F)
# save.image(file = "20211121.RData")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfiData")
library(minfi)
library(minfiData)
library(sva)
library(devtools)
# baseDir <- system.file("D:/TCGA-ESCA-IDAT/EAC/IDAT", package="minfiData")
# targets <- read.450k.sheet(baseDir)
#FWQ
targets <- read.metharray.sheet("/student1/MDD/TCGA/IDAT/EAC/IDAT")
# if (require(minfiData)) {
#   GMset <- mapToGenome(MsetEx)
#   comps <- compartments(GMset)
# }
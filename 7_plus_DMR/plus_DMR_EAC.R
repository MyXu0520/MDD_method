#########################EAC
#plus_DMR
library(DSS)
library(bsseq)
dir_path <- "/student1/MDD/Tabfile"
files <- list.files(dir_path, pattern = "\\.tab$", full.names = TRUE)
data_list <- list()
sample_names <- c()
cancer_files <- files[grepl("^EAC_", basename(files))]
for (file in cancer_files) {
  sample_name <- gsub("^EAC_|\\.tab$", "", basename(file))
  sample_names <- c(sample_names, paste0("C", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
# > names(data_list)
# [1] "C1methyCov" "C2methyCov" "C3methyCov" "C4methyCov" "C6methyCov"
names(data_list) <- sample_names
data_listC1 <- data_list
# data_listC1 <- data_list[c(1:21)]
rm(data_list,sample_names)
library(DSS)
library(bsseq)
dir_path <- "/student1/MDD/Tabfile"
files <- list.files(dir_path, pattern = "\\.tab$", full.names = TRUE)
data_list <- list()
sample_names <- c()
cancer_files <- files[grepl("^GEJ_", basename(files))]
for (file in cancer_files) {
  sample_name <- gsub("^GEJ_|\\.tab$", "", basename(file))
  sample_names <- c(sample_names, paste0("C", sample_name))
  data <- read.table(file, header = FALSE, col.names = c("chr", "pos", "N", "X"))
  data_list <- c(data_list, list(data))
}
names(data_list) <- sample_names
# > names(data_list)
# [1] "C1methyCov"              "C2methyCov"
# [3] "C3methyCov"              "C4methyCov"
# [5] "C5methyCov"              "C6methyCov"
# [7] "C7methyCov"              "CNonmalignant_1methyCov"
# [9] "CNonmalignant_2methyCov" "CNonmalignant_3methyCov"
# [11] "CNonmalignant_4methyCov" "CNonmalignant_5methyCov"
# [13] "CNonmalignant_6methyCov" "CNonmalignant_7methyCov"
# >
data_listC2 <- data_list
data_listC2 <- data_list[c(1:7)]
rm(data_list,sample_names)
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
data_list <- c(data_listC1,data_listC2,data_listN)
BSobj <- makeBSseqData(data_list, names(data_list))
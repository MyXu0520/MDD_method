require(multiROC)
library(multiROC)
library(ggplot2)
library(pROC)
library(rpart)
library(rpart.plot)
library(caret)
load("D:/MDD_result/21/TCGA_8cancer_DMR.RData")
load("D:/MDD_result/18/TCGA_8cancer_MDD.RData")
SC <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
AD <- read.table("D:/MDD_result/15/clinical.cart.2024-12-25/clinical.tsv",header = T,sep = "\t",fill = T)
#all_data_shared
for (i in 1:length(all_data_shared)) {
  temp <- all_data_shared[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_data_shared[[i]] <- temp
}
#all_data_EAC_SPEC
for (i in 1:length(all_data_EAC_SPEC)) {
  temp <- all_data_EAC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_data_EAC_SPEC[[i]] <- temp
}
#all_data_ESCC_SPEC
for (i in 1:length(all_data_ESCC_SPEC)) {
  temp <- all_data_ESCC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_data_ESCC_SPEC[[i]] <- temp
}
for (i in 1:length(all_DMRdata_shared)) {
  temp <- all_DMRdata_shared[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_DMRdata_shared[[i]] <- temp
}
#all_DMRdata_EAC_SPEC
for (i in 1:length(all_DMRdata_EAC_SPEC)) {
  temp <- all_DMRdata_EAC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_DMRdata_EAC_SPEC[[i]] <- temp
}
#all_DMRdata_ESCC_SPEC
for (i in 1:length(all_DMRdata_ESCC_SPEC)) {
  temp <- all_DMRdata_ESCC_SPEC[[i]]
  new_string <- gsub( "\\.", "-",temp$Sample)
  temp$Sample <- new_string
  cancer_class <- c()
  for (j in 1:nrow(temp)) {
    if(substr(temp$Sample[j],1,12) %in% SC$case_submitter_id){cancer_class[j] <- "SC"}
    else if(substr(temp$Sample[j],1,12) %in% AD$case_submitter_id){cancer_class[j] <- "AD"}
    else{cancer_class[j] <- "OT"}
  }
  temp$Cancer_class <- cancer_class
  all_DMRdata_ESCC_SPEC[[i]] <- temp
}
atlas_MDD_EAC_SPEC <- as.data.frame(matrix(ncol = 5))
colnames(atlas_MDD_EAC_SPEC) <- colnames(all_data_EAC_SPEC[[1]])
for (i in 1:length(all_data_EAC_SPEC)) {
  temp_atlas <- all_data_EAC_SPEC[[i]]
  atlas_MDD_EAC_SPEC <- rbind(atlas_MDD_EAC_SPEC,temp_atlas)
}
atlas_MDD_EAC_SPEC <- atlas_MDD_EAC_SPEC[-1,]
colnames(atlas_MDD_EAC_SPEC)[2] <- "methy_level_EAC"
atlas_MDD_ESCC_SPEC <- as.data.frame(matrix(ncol = 5))
colnames(atlas_MDD_ESCC_SPEC) <- colnames(all_data_ESCC_SPEC[[1]])
for (i in 1:length(all_data_ESCC_SPEC)) {
  temp_atlas <- all_data_ESCC_SPEC[[i]]
  atlas_MDD_ESCC_SPEC <- rbind(atlas_MDD_ESCC_SPEC,temp_atlas)
}
atlas_MDD_ESCC_SPEC <- atlas_MDD_ESCC_SPEC[-1,]
colnames(atlas_MDD_ESCC_SPEC)[2] <- "methy_level_ESCC"
atlas_MDD_SPEC <- merge(atlas_MDD_EAC_SPEC,atlas_MDD_ESCC_SPEC,by = "Sample")
atlas_MDD_SPEC <- atlas_MDD_SPEC[,c(1,2,6,5)]
colnames(atlas_MDD_SPEC)[4] <- "Cancer_class"
atlas_MDD_SPEC <- atlas_MDD_SPEC[,-1]
atlas_MDD_SPEC$Cancer_class <- as.factor(atlas_MDD_SPEC$Cancer_class)
set.seed(123456)
total_number <- nrow(atlas_MDD_SPEC)
train_idx <- sample(total_number, round(total_number*0.7))
train_df <- atlas_MDD_SPEC[train_idx, ]
test_df <- atlas_MDD_SPEC[-train_idx, ]
# train_df$Cancer_class <- as.factor(train_df$Cancer_class)
rf_res <- randomForest::randomForest(Cancer_class ~ ., data = train_df[,-1], ntree = 100)
rf_res
rf_pred <- predict(rf_res, test_df, type = 'prob')
rf_pred <- data.frame(rf_pred)
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")
mn_res <- nnet::multinom(Cancer_class ~., data = train_df)
mn_pred <- predict(mn_res, test_df, type = 'prob')
mn_pred <- data.frame(mn_pred)
colnames(mn_pred) <- paste(colnames(mn_pred), "_pred_MN")
true_label <- dummies::dummy(test_df$Cancer_class, sep = ".")
## Warning in model.matrix.default(~x - 1, model.frame(~x - 1), contrasts =
## FALSE): non-list contrasts argument ignored
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, rf_pred, mn_pred)
head(final_df)
roc_res <- multi_roc(final_df, force_diag=T)
plot_roc_df <- plot_roc_data(roc_res)
head(plot_roc_df)
require(ggplot2)
## Loading required package: ggplot2
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               colour='grey', linetype = 'dotdash') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(),
        legend.background = element_rect(fill=NULL, size=0.5,
                                         linetype="solid", colour ="black"))
auc_list <- roc_res$AUC
auc_df <- do.call(rbind, lapply(names(auc_list), function(method) {
  method_auc <- auc_list[[method]]
  data.frame(
    Method = method,
    Group = names(method_auc),
    AUC = as.numeric(method_auc),
    row.names = NULL
  )
}))
# auc_df <- subset(auc_df, !(Group %in% c("macro", "micro")))
auc_df$Group <- trimws(auc_df$Group)
print(auc_df)
atlas_DMR_EAC_SPEC <- as.data.frame(matrix(ncol = 5))
colnames(atlas_DMR_EAC_SPEC) <- colnames(all_DMRdata_EAC_SPEC[[1]])
for (i in 1:length(all_DMRdata_EAC_SPEC)) {
  temp_atlas <- all_DMRdata_EAC_SPEC[[i]]
  atlas_DMR_EAC_SPEC <- rbind(atlas_DMR_EAC_SPEC,temp_atlas)
}
atlas_DMR_EAC_SPEC <- atlas_DMR_EAC_SPEC[-1,]
colnames(atlas_DMR_EAC_SPEC)[2] <- "methy_level_EAC"
atlas_DMR_ESCC_SPEC <- as.data.frame(matrix(ncol = 5))
colnames(atlas_DMR_ESCC_SPEC) <- colnames(all_DMRdata_ESCC_SPEC[[1]])
for (i in 1:length(all_DMRdata_ESCC_SPEC)) {
  temp_atlas <- all_DMRdata_ESCC_SPEC[[i]]
  atlas_DMR_ESCC_SPEC <- rbind(atlas_DMR_ESCC_SPEC,temp_atlas)
}
atlas_DMR_ESCC_SPEC <- atlas_DMR_ESCC_SPEC[-1,]
colnames(atlas_DMR_ESCC_SPEC)[2] <- "methy_level_ESCC"
atlas_DMR_SPEC <- merge(atlas_DMR_EAC_SPEC,atlas_DMR_ESCC_SPEC,by = "Sample")
atlas_DMR_SPEC <- atlas_DMR_SPEC[,c(1,2,6,5)]
colnames(atlas_DMR_SPEC)[4] <- "Cancer_class"
atlas_DMR_SPEC <- atlas_DMR_SPEC[,-1]
atlas_DMR_SPEC$Cancer_class <- as.factor(atlas_DMR_SPEC$Cancer_class)
set.seed(123456)
total_number <- nrow(atlas_DMR_SPEC)
train_idx <- sample(total_number, round(total_number*0.7))
train_df <- atlas_DMR_SPEC[train_idx, ]
test_df <- atlas_DMR_SPEC[-train_idx, ]
train_df$Cancer_class <- as.factor(train_df$Cancer_class)
rf_res <- randomForest::randomForest(Cancer_class ~ ., data = train_df[,-1], ntree = 100)
rf_res
rf_pred <- predict(rf_res, test_df, type = 'prob')
rf_pred <- data.frame(rf_pred)
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")
mn_res <- nnet::multinom(Cancer_class ~., data = train_df)
mn_pred <- predict(mn_res, test_df, type = 'prob')
mn_pred <- data.frame(mn_pred)
colnames(mn_pred) <- paste(colnames(mn_pred), "_pred_MN")
true_label <- dummies::dummy(test_df$Cancer_class, sep = ".")
## Warning in model.matrix.default(~x - 1, model.frame(~x - 1), contrasts =
## FALSE): non-list contrasts argument ignored
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, rf_pred, mn_pred)
head(final_df)
roc_res <- multi_roc(final_df, force_diag=T)
plot_roc_df <- plot_roc_data(roc_res)
head(plot_roc_df)
require(ggplot2)
## Loading required package: ggplot2
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               colour='grey', linetype = 'dotdash') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(),
        legend.background = element_rect(fill=NULL, size=0.5,
                                         linetype="solid", colour ="black"))
auc_list <- roc_res$AUC
auc_df <- do.call(rbind, lapply(names(auc_list), function(method) {
  method_auc <- auc_list[[method]]
  data.frame(
    Method = method,
    Group = names(method_auc),
    AUC = as.numeric(method_auc),
    row.names = NULL
  )
}))
# auc_df <- subset(auc_df, !(Group %in% c("macro", "micro")))
auc_df$Group <- trimws(auc_df$Group)
print(auc_df)
rownames(atlas_DMR_SPEC) <- rownames(atlas_DMR_EAC_SPEC)
rownames(atlas_MDD_SPEC) <- rownames(atlas_MDD_EAC_SPEC)
colnames(atlas_DMR_SPEC) <- c("methy_level_EAC_DMR","methy_level_ESCC_DMR","Cancer_class")
colnames(atlas_MDD_SPEC) <- c("methy_level_EAC_MDD","methy_level_ESCC_MDD","Cancer_class")
all_atlas <- cbind(atlas_DMR_SPEC,atlas_MDD_SPEC[,-ncol(atlas_MDD_SPEC)])
all_atlas <- all_atlas[,c(1,2,4,5,3)]
set.seed(123456)
total_number <- nrow(all_atlas)
train_idx <- sample(total_number, round(total_number*0.7))
train_df <- all_atlas[train_idx, ]
test_df <- all_atlas[-train_idx, ]
# train_df$Cancer_class <- as.factor(train_df$Cancer_class)
rf_res <- randomForest::randomForest(Cancer_class ~ ., data = train_df[,-1], ntree = 100)
rf_res
rf_pred <- predict(rf_res, test_df, type = 'prob')
rf_pred <- data.frame(rf_pred)
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")
mn_res <- nnet::multinom(Cancer_class ~., data = train_df)
mn_pred <- predict(mn_res, test_df, type = 'prob')
mn_pred <- data.frame(mn_pred)
colnames(mn_pred) <- paste(colnames(mn_pred), "_pred_MN")
true_label <- dummies::dummy(test_df$Cancer_class, sep = ".")
## Warning in model.matrix.default(~x - 1, model.frame(~x - 1), contrasts =
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, rf_pred, mn_pred)
head(final_df)
colnames(true_label) <- gsub("Cancer_class\\.", "", colnames(true_label))
colnames(true_label) <- trimws(true_label)
colnames(true_label) <- paste0(colnames(true_label), "_true")
colnames(rf_pred) <- trimws(colnames(rf_pred))
colnames(rf_pred) <- paste0(gsub(" _pred_RF", "", colnames(rf_pred)), "_pred_RF")
colnames(mn_pred) <- trimws(colnames(mn_pred))
colnames(mn_pred) <- paste0(gsub(" _pred_MN", "", colnames(mn_pred)), "_pred_MN")
final_df <- cbind(true_label, rf_pred, mn_pred)
roc_res <- multi_roc(final_df, force_diag = TRUE)
# roc_res <- multi_roc(final_df, force_diag=T)
plot_roc_df <- plot_roc_data(roc_res)
head(plot_roc_df)
require(ggplot2)
## Loading required package: ggplot2
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               colour='grey', linetype = 'dotdash') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(),
        legend.background = element_rect(fill=NULL, size=0.5,
                                         linetype="solid", colour ="black"))
auc_table <- roc_res[roc_res$Metric == "AUC", c("Group", "Method", "Value")]
print(auc_table)
# library(ggplot2)
#
# ggplot(plot_roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
#   geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
#   geom_text(aes(label = paste("AUC = ", round(AUC, 2))),
#             size = 3, vjust = -0.5, hjust = 1.1,
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.justification = c(1, 0),
#         legend.position = c(.95, .05),
#         legend.title = element_blank(),
#         legend.background = element_rect(fill = NULL, size = 0.5,
#                                          linetype = "solid", colour = "black")) +
#   labs(title = "ROC Curve with AUC values")
#
auc_list <- roc_res$AUC
auc_df <- do.call(rbind, lapply(names(auc_list), function(method) {
  method_auc <- auc_list[[method]]
  data.frame(
    Method = method,
    Group = names(method_auc),
    AUC = as.numeric(method_auc),
    row.names = NULL
  )
}))
# auc_df <- subset(auc_df, !(Group %in% c("macro", "micro")))
auc_df$Group <- trimws(auc_df$Group)
print(auc_df)
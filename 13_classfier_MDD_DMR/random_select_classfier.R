# D:/MDD_result/4
load("D:/MDD_result/4/TCGA_sample_separate.RData")
load("D:/MDD_result/4/cg_TCGA_select.RData")
load("D:/MDD_result/4/shared_ESCA.RData")
load("D:/MDD_result/4/shared_Normal.RData")
Normal_sample <- TCGA_ESCA[,which(substr(colnames(TCGA_ESCA),14,15) == "11" )]
TCGA_EAC_sample <- TCGA_EAC_sample[,which(substr(colnames(TCGA_EAC_sample),14,15) == "01" )]
TCGA_ESCC_sample <- TCGA_ESCC_sample[,which(substr(colnames(TCGA_ESCC_sample),14,15) == "01" )]
cancer_sample <- cbind(TCGA_EAC_sample,TCGA_ESCC_sample)
# shared_ESCA <- matrix(nrow = length(shared_list),ncol = ncol(cancer_sample))
# for (i in 1:ncol(cancer_sample)) {
#   for (j in 1:length(shared_list)) {
#     temp_cg <- shared_list[[j]]
#     shared_ESCA[j,i] <- mean(na.omit(cancer_sample[temp_cg,i]))
#   }
# }
# colnames(shared_ESCA) <- colnames(cancer_sample)
# rownames(shared_ESCA) <- paste0("MDD",c(1:length(shared_list)))
# shared_Normal <- matrix(nrow = length(shared_list),ncol = ncol(Normal_sample))
# for (i in 1:ncol(Normal_sample)) {
#   for (j in 1:length(shared_list)) {
#     temp_cg <- shared_list[[j]]
#     shared_Normal[j,i] <- mean(na.omit(Normal_sample[temp_cg,i]))
#   }
# }
# colnames(shared_Normal) <- colnames(Normal_sample)
# rownames(shared_Normal) <- paste0("MDD",c(1:length(shared_list)))
shared_MDD_sample <- cbind(shared_Normal,shared_ESCA)
shared_MDD_sample <- t(shared_MDD_sample)
sample_group <- c(rep(0,ncol(shared_Normal)),rep(1,ncol(shared_ESCA)))
shared_MDD_sample <- as.data.frame(shared_MDD_sample)
shared_MDD_sample$group <- sample_group
set.seed(2021)
shared_MDD_train <- shared_MDD_sample[sample(1:nrow(shared_MDD_sample),0.8*nrow(shared_MDD_sample),replace = F),]
shared_MDD_test <- shared_MDD_train[-which(rownames(shared_MDD_sample) %in% rownames(shared_MDD_train)),]
shared_MDD_train <- shared_MDD_train[,colSums(is.na(shared_MDD_train)) == 0]
shared_MDD_test <- shared_MDD_test[,colSums(is.na(shared_MDD_test)) == 0]
# if (!requireNamespace("randomForest", quietly = TRUE)) {
#   install.packages("randomForest")
# }
library(randomForest)
set.seed(2021)
rf_model <- randomForest(x = as.matrix(shared_MDD_train[, -ncol(shared_MDD_train)]),
                         y = shared_MDD_train$group,
                         ntree = 1000,
                         mtry = floor(sqrt(ncol(shared_MDD_train) - 1)))
importance_scores <- importance(rf_model)
num_features_to_select <- ceiling(0.05 * length(importance_scores))
top_features <- order(importance_scores, decreasing = TRUE)[1:num_features_to_select]
shared_MDD_train_selected <- shared_MDD_train[, c(top_features, ncol(shared_MDD_train))]
shared_MDD_test_selected <- shared_MDD_test[, top_features]
# shared_MDD_train_selected <- shared_MDD_train[,c(top_features,group)]
dtrain_selected <- xgb.DMatrix(data = as.matrix(shared_MDD_train_selected[, -ncol(shared_MDD_train_selected)]),
                               label = shared_MDD_train_selected$group)
dtest_selected <- xgb.DMatrix(data = as.matrix(shared_MDD_test_selected))
param <- list(objective = "binary:logistic")
bst_selected <- xgb.train(param, dtrain_selected, nrounds = 50)
pred_selected <- predict(bst_selected, newdata = dtest_selected)
pred_selected_binary <- round(pred_selected)
table(shared_MDD_test$group, pred_selected_binary, dnn = c("", ""))
require(pROC)
xgboost_roc_selected <- roc(shared_MDD_test$group, as.numeric(pred_selected_binary))
pdf("xgboost_MDD_selected_features.pdf")
plot(xgboost_roc_selected, print.auc = TRUE, auc.polygon = TRUE, grid = c(0.1, 0.2), grid.col = c("green", "red"),
     max.auc.polygon = TRUE, auc.polygon.col = "skyblue", print.thres = TRUE, main = 'XGBoostROC（）')
dev.off()
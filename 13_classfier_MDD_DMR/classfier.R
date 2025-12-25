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
                      sizes = c(2:25,30,40,50,60,70,80,90,100,150,200,250,300,350,400),
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
# data_test <- data_test[,c(selected_features)]
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
pre <- predict(model_xgb,newdata = dtest)
library(caret)
xgb.cf <-caret::confusionMatrix(as.factor(pre),as.factor(test_y))
xgb.cf
table(testset$label,pre,dnn = c("",""))
require(pROC)
xgboost_roc<-roc(testset$label,as.numeric(pre))
plot(xgboost_roc,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),gird.col=c("green","red"),max.auc.polygon=TRUE,auc.polygon.col="skyblue",print.thres=TRUE,main='xgboost_model_ROC')
save.image("D:/MDD_result/15/MDD_ESCA_classfier.RData")
#MDD_core
#DMR
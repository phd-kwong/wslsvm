library(pROC)
library(e1071)
source('extract_intersect_external_breast.R')

# WSLSVM
# Data  = read.table("Feature_GSE10780_WSLSVM.txt ", header = T, check.names = FALSE)
Data  = read.table("Feature_GSE21422_WSLSVM.txt ", header = T, check.names = FALSE)
# Data  = read.table("Feature_GSE38959_WSLSVM.txt ", header = T, check.names = FALSE)
# Data  = read.table("Feature_GSE42568_WSLSVM.txt ", header = T, check.names = FALSE)

Data2 = read.table("Feature_breast_WSLSVM.txt", header = T, check.names = FALSE)



#####################################  svm model -- built model makes prediction #########
x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
y.train[y.train==0] <- -1
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]
y.test[y.test==0] <- -1

tuned <- tune.svm(x.train,y.train, gamma = 10^(-6:-1), cost = 10^(1:2)) # tune
summary(tuned) # to select best gamma and cost

model <- svm(x.train, y.train, kernel =  "linear", cost = tuned$best.parameters$cost, gamma=tuned$best.parameters$gamma,  scale = FALSE)

summary(model)




# Predict -----------------------------------------------------------------
p_test <- predict(model, x.test)
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
# setwd('D:\\E\\...\\Data\\RTCGA\\result\\A_test')
# write.csv(A_test,"A_test_20437.csv", row.names=F)
names(A_test) <- c("p", "outcome")
# jpeg(file = "pAUC_train.jpg")
plot.roc(A_test$outcome, A_test$p, print.auc=T)
ans <- roc(A_test$outcome, A_test$p, print.auc=T)
# dev.off()
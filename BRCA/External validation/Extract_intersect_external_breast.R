library(openxlsx)
library(caret)
library(tidyverse)
library(corrplot)
library(readr)
library(openxlsx)
library(data.table)
library(MASS)
library(logADMMmodify1)
library(logADMMspline)
X_node <- read.xlsx('all_gene_2026120.xlsx')
# XX_all <- read.table("GSE10780_scale_breast.txt", header = T, check.names = FALSE)
XX_all <- read.table("GSE21422_scale_breast.txt", header = T, check.names = FALSE)
# XX_all <- read.table("GSE38959_scale_breast.txt", header = T, check.names = FALSE)
# XX_all <- read.table("GSE42568_scale_breast.txt", header = T, check.names = FALSE)

XX_all <- t(XX_all)
common_columns <- intersect(X_node$WSLSVM, colnames(XX_all))
colnames <- common_columns
colnames
sort(colnames)
length(colnames)

# XX1 <- read.table("breast_scale_net_adjp_net.txt", header = T, check.names = FALSE)
# XX1 <- data.frame(t(XX1))
# XX1 <- read.table("E:/SVM/Realdata/RD1_SVM/Data_train/1.txt", header = T, check.names = FALSE)
XX1 <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene.txt", header = T, check.names = FALSE)
XX1End <- data.frame(t(XX1))
# XX1End <- data.frame(XX1)
XXX1 <- as.matrix(XX1End[,-1])
data_X <- XXX1
data_Y <- XX1End[,1]
XX <- cbind(data_Y,data_X)

selected_data <- XX[ , colnames]
XXX1 <- as.matrix(selected_data)
data_X <- XXX1

data_Y <- XX1End[,1]
XX_end <- as.data.frame(cbind(data_Y,data_X))

write.table(t(XX_end), file = "Feature_breast_WSLSVM.txt", sep = "\t", col.names = TRUE, row.names = TRUE)

XX_true <- XX_all[ , colnames]
# Y10780 <- read.xlsx('GSE10780.xlsx')
Y10780 <- read.xlsx('GSE21422.xlsx')
# Y10780 <- read.xlsx('GSE38959.xlsx')
# Y10780 <- read.xlsx('GSE42568.xlsx')

test_data <- as.data.frame(cbind(Y10780[,2],XX_true))
# write.table(t(test_data), file = "Feature_GSE10780_WSLSVM.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(t(test_data), file = "Feature_GSE21422_WSLSVM.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(t(test_data), file = "Feature_GSE38959_WSLSVM.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(t(test_data), file = "Feature_GSE42568_WSLSVM.txt", sep = "\t", col.names = TRUE, row.names = TRUE)


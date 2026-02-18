library(openxlsx)
library(caret)
library(tidyverse)
library(corrplot)
library(readr)
library(openxlsx)
library(data.table)
library(MASS)
library(Rcpp)
library(readr)
library(MASS)
library(glmnet)
library(newSVM)
library(penalizedSVM)
library(newSVMspline)
library(newSVMlasso)
library(sparseSVM)
source("svmfun.R")


XX1 <- read.xlsx('LSVT_sorted_matrix.xlsx')
XXX1 <- as.matrix(XX1[,-1])
data_X <- XXX1
data_Y <- XX1[,1]
data_Y[which(data_Y==0)] <- -1
XX <- cbind(data_Y,data_X)
### proprocessing ###

Xx0 <- XX[which(XX[,1]==-1),]
Xx1 <- XX[which(XX[,1]==1),]
n0=round(dim(Xx0)[1]*0.6)
n1=round(dim(Xx1)[1]*0.6)
nv0=round(dim(Xx0)[1]*0.2)
nv1=round(dim(Xx1)[1]*0.2)

nt0=dim(Xx0)[1]-n0-nv0
nt1=dim(Xx1)[1]-n1-nv1
n=n0+n1
nv=nv0+nv1
nt=nt0+nt1
p <- dim(data_X)[2]
### end ###


G_or <- c(rep(1,p))
G1 <- diag(1,p,p)
diag(G1[-1,-p]) = -1
G1 <- t(G1)
G <- G1[-p,]
G2 <- t(G)%*%G
G_2 <- G2[c(-1,-p),]
G <- G_2
# G <- diag(p)
m=dim(G)[1]

ga1 <- 0.95
theta <- 0.05


lm1 <- qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*0.1
alpha.v=c(seq(0.3, 0.01, length.out = 20))
lg=length(alpha.v)
bMatrix=matrix(0,p+4,100)
bMatrixm=matrix(0,p+3,100)
sr_1 <- matrix(0,1,100)
for (M in 1:100) {
  X= read.csv(file = paste0('E:/SVM/training data/X/X_',M,'.csv'),encoding = 'UTF-8')[,-1]     
  Y= read.csv(file = paste0('E:/SVM/training data/Y/Y_',M,'.csv'),encoding = 'UTF-8')[,-1]        
  
  Xva= read.csv(file = paste0('E:/SVM/validation data/X/Xva_',M,'.csv'),encoding = 'UTF-8')[,-1]     
  Yva= read.csv(file = paste0('E:/SVM/validation data/Y/Yva_',M,'.csv'),encoding = 'UTF-8')[,-1]        
  
  Xtest= read.csv(file = paste0('E:/SVM/testing data/X/Xtest_',M,'.csv'),encoding = 'UTF-8')[,-1]          
  Ytest= read.csv(file = paste0('E:/SVM/testing data/Y/Ytest_',M,'.csv'),encoding = 'UTF-8')[,-1] 
  
  bmat=matrix(0,p+4,lg)
  for (j in 1:lg) {
    pen=c(rep(1,2*p))
    lmnew=lm1
    ANS1=newSVMADMMspline(Y,X,c(alpha.v[j]*lmnew*pen[1:p],(1-alpha.v[j])*lmnew*pen[p+1:m]),
                          G,n,p,ItM=10*10^5,eabsn=10^(-7),erel=10^(-7))
    bet1=ANS1[[1]][-1]
    b0=ANS1[[1]][1]
    cvnew=targetf(Yva,Xva,bet1,(1-alpha.v[j])*lmnew*G,c(alpha.v[j]*lmnew*pen[1:p],pen[p+1:m]),b0)
    bmat[,j] <- c(round(c(b0,bet1,alpha.v[j],cvnew[2]),5), ANS1[[2]])##+loss
    print(paste0(M,'-',j,'-',ANS1[[2]]))
  }
  ind=which.min(bmat[p+3,1:lg]) ##loss----alpha
  b0 = bmat[1,ind]
  bM=bmat[2:(p+1),ind]
  cv1 <- bmat[p+3,ind]
  alp=alpha.v[ind]
  bMatrix[,M]=c(b0,bM,cv1,alp,lmnew)
  
  Y_test1G <- Xtest%*%bMatrix[2:(p+1),M] + bMatrix[1,M]
  result_1G <- ifelse(round(Y_test1G,2) >= 0, 1, -1)
  sr_1[M] <- length(Ytest) - length(which(Ytest==result_1G))
  
}

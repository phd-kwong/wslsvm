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
library(sparseSVM)
source("svmfun.R")

XX1 <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene.txt", header = T, check.names = FALSE)
XX1 <- data.frame(t(XX1))
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

G_weight <- read.csv('structural_G.csv')[,-1]
G <- as.matrix(G_weight)
m=dim(G)[1]

ga1 <- 0.95
theta <- 0.05

lm1 <- qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*0.05
lg=length(alpha.v)
Gnew=G
Gnew1=t(Gnew)%*%Gnew
bMatrix=matrix(0,p+4,100)
sr_1 <- matrix(0,1,100)


for (i in 1:100) {
  X= read.csv(file = paste0('E:/SVM/DATA/Data_train/X/M',i,'X.csv'),encoding = 'UTF-8')[,-1]       
  Y= read.csv(file = paste0('E:/SVM/DATA/Data_train/Y/M',i,'Y.csv'),encoding = 'UTF-8')[,-1]        
  
  Xva= read.csv(file = paste0('E:/SVM/DATA/Data_val/X/M',i,'X.csv'),encoding = 'UTF-8')[,-1]       
  Yva= read.csv(file = paste0('E:/SVM/DATA/Data_val/Y/M',i,'Y.csv'),encoding = 'UTF-8')[,-1]        
  
  Xtest= read.csv(file = paste0('E:/SVM/DATA/Data_test/X/M',i,'X.csv'),encoding = 'UTF-8')[,-1]         
  Ytest= read.csv(file = paste0('E:/SVM/DATA/Data_test/Y/M',i,'Y.csv'),encoding = 'UTF-8')[,-1]         
  
  bmat=matrix(0,p+4,lg)
  for (j in 1:lg) {
    Xwav=X%*%solve(alpha.v[j]*diag(p)+(1-alpha.v[j])*Gnew1)%*%cbind(diag(p),t(Gnew))
    pen=round(sqrt(colMeans(Xwav^2)),4)
    lmnew=lm1
    ANS1=newSVMADMM(Y,X,c(alpha.v[j]*lmnew*pen[1:p],lmnew*pen[p+1:m]*n^(-4/5)),
                    (1-alpha.v[j])*G*n^(4/5),n,p,ItM=3*10^5,eabsn=10^(-5),erel=10^(-5))
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

library(MASS)
library(lda)
library(Rcpp)
library(readr)
library(MASS)
library(glmnet)
library(newSVM)
library(penalizedSVM)
library(newSVMspline)

source("svmfun.R")
n=200
p=500
nv=200
nt=100
G_or <- c(rep(1,p))
G1 <- diag(1,p,p)
diag(G1[-1,-p]) = -1
G1 <- t(G1)
G <- G1[-p,]
G2 <- t(G)%*%G
G_2 <- G2[c(-1,-p),]
# G <- G_2
ratio_minor <- 0.05
n_minor <- max(1, as.integer(n * ratio_minor))  # 少数类样本数（至少1个）
n_major <- n - n_minor 

nv_minor <- max(1, as.integer(nv * ratio_minor))  # 少数类样本数（至少1个）
nv_major <- nv - nv_minor 

nt_minor <- max(1, as.integer(nt * ratio_minor))  # 少数类样本数（至少1个）
nt_major <- nt - nt_minor 

out1=23
out2=23
out3=100
beta=c(rep(0,out3),0.333,rep(0.167,out1),0.333,rep(0,out3),0.333,rep(0.167,out2),0.333,rep(0,p-250))


m=dim(G)[1]
EG=eigen(t(G)%*%G)[[2]][,-c(1:m)]
AN=rbind(G,t(EG))
ga1 <- 0.95
theta <- 0.05
alpha.v=c(exp(seq(log(0.1),log(0.35),length.out= 15)))
lg=length(alpha.v)
lm1 <- qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*0.5
Gnew=rbind(G,t(eigen(t(G)%*%G)[[2]][,-(1:m)]))
bMatrix1=matrix(0,p+4,100)
cv_test<- matrix(0,1,100)

Y <- read.csv('Y005.csv')[,-1]
Yva <- read.csv('Yva005.csv')[,-1]
Ytest <- read.csv('Ytest005.csv')[,-1]

for (M in 1:100) {
 X = as.matrix(read.csv(file = paste0("E:/SVM/LDA_005/X/M",M,'X.csv'),encoding = 'UTF-8')[,-1])
  Xva = as.matrix(read.csv(file = paste0("E:/SVM/LDA_005/Xva/M",M,'Xva.csv'),encoding = 'UTF-8')[,-1])
  Xtest = as.matrix(read.csv(file = paste0("E:/SVM/LDA_005/Xtest/M",M,'Xtest.csv'),encoding = 'UTF-8')[,-1])

  bmat=matrix(0,p+3,lg)
  for (j in 1:lg) {
    pen=c(rep(1,2*p))
    lmnew=lm1
    ANS1=newSVMADMMspline(Y,X,c(alpha.v[j]*lmnew*pen[1:p],(1-alpha.v[j])*lmnew*pen[p+1:m]),
                    G,n,p,ItM=3*10^5,eabsn=10^(-7),erel=10^(-7))
    bet1=ANS1[[1]][-1]
    b0=ANS1[[1]][1]
    cvnew=targetf(Yva,Xva,bet1,(1-alpha.v[j])*lmnew*G,c(alpha.v[j]*lmnew*pen[1:p],pen[p+1:m]),b0)
    bmat[,j] <- c(round(c(b0,bet1),5), cvnew[2], ANS1[[2]])
    print(paste0(M,'-',j,'-',ANS1[[2]]))
  }
  ind=which.min(bmat[p+2,1:lg])
  b0 = bmat[1,ind]
  bM=bmat[2:(p+1),ind]
  cv1 <- bmat[p+2,ind]
  alp=alpha.v[ind]
  bMatrix1[,M]=c(b0,bM,cv1,alp,lmnew)

 Y_test1G111 <- Xtest%*%bM + b0
 Y_test1G <- ifelse(Y_test1G111 >=0, 1, -1)
 cv_test[M]= 1-length(which(Y_test1G==Ytest))/length(Ytest)



}

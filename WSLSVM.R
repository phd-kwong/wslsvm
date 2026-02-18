library(MASS)
library(lda)
library(Rcpp)
library(readr)
library(MASS)
library(glmnet)
library(newSVM)
library(penalizedSVM)
library(newSVMspline)
library(mvtnorm)

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
beta=c(rep(0,out3),0.3333333,rep(0.1666667,out1),0.3333333,rep(0,out3),
       0.3333333,rep(0.1666667,out2),0.3333333,rep(0,p-250))

m=dim(G)[1]
EG=eigen(t(G)%*%G)[[2]][,-c(1:m)]
AN=rbind(G,t(EG))
ga1 <- 0.95
theta <- 0.05



alpha.v=c(exp(seq(log(0.0015),log(0.01),length.out= 5))) ##1G 601700
lg=length(alpha.v)
lm1 <- qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*(0.125)

Gnew=rbind(G,t(eigen(t(G)%*%G)[[2]][,-(1:m)]))
Gnew1=t(Gnew)%*%Gnew
bMatrix=matrix(0,p+4,100)
bMatrixm=bMatrix

Y <- read.csv('Y005.csv')[,-1]
Yva <- read.csv('Yva005.csv')[,-1]
Ytest <- read.csv('Ytest005.csv')[,-1]
cv_test<- matrix(0,1,100)
cv_test1<- matrix(0,1,100)
for (M in 3:100) { #2 
  # M <- IID[i]
  X = as.matrix(read.csv(file = paste0("E:/SVM/LDA_005/X/M",M,'X.csv'),encoding = 'UTF-8')[,-1])
  Xva = as.matrix(read.csv(file = paste0("E:/SVM/LDA_005/Xva/M",M,'Xva.csv'),encoding = 'UTF-8')[,-1])
  Xtest = as.matrix(read.csv(file = paste0("E:/SVM/LDA_005/Xtest/M",M,'Xtest.csv'),encoding = 'UTF-8')[,-1])
  
   bmat=matrix(0,p+3,lg)
  for (j in 1:lg) {
    Xwav=X%*%solve(alpha.v[j]*diag(p)+(1-alpha.v[j])*Gnew1)%*%cbind(diag(p),t(Gnew))
    pen=round(sqrt(colMeans(Xwav^2)),4)
    lmnew=lm1
    ANS1=newSVMADMM(Y,X,c(alpha.v[j]*lmnew*pen[1:p],lmnew*pen[p+1:m]*n^(-4/5)),
                    (1-alpha.v[j])*G*n^(4/5),n,p,ItM=10*10^5,eabsn=10^(-6),erel=10^(-6))
    bet1=ANS1[[1]][-1]
    b0=ANS1[[1]][1]
    cvnew=targetf(Yva,Xva,bet1,(1-alpha.v[j])*lmnew*G,c(alpha.v[j]*lmnew*pen[1:p],pen[p+1:m]),b0)
    bmat[,j] <- c(round(c(b0,bet1),5), cvnew[2],ANS1[[2]])
    print(paste0(M,'-',j,'-',ANS1[[2]]))
  }
  ind=which.min(bmat[p+2,1:lg])
  b0 = bmat[1,ind]
  bM=bmat[2:(p+1),ind]
  cv1 <- bmat[p+2,ind]
  alp=alpha.v[ind]
  bMatrix[,M]=c(b0,bM,cv1,alp,lmnew)
  
  Y_test1G111 <- Xtest%*%bM + b0
  Y_test1G <- ifelse(Y_test1G111 >=0, 1, -1)
  cv_test[M]= 1-length(which(Y_test1G==Ytest))/length(Ytest)
  
  
  Xwav1=X%*%solve(alp*diag(p)+(1-alp)*t(Gnew)%*%Gnew)%*%cbind(diag(p),t(Gnew))
  pen1=round(sqrt(colMeans(Xwav1^2)),4)
  lmnew1=lm1
  Gnonzeroind=which(abs(G%*%bmat[2:(p+1),ind])>=0.0001)
  if(length(Gnonzeroind)==0){
    bet2=bmat[2:(p+1),ind]
    b02=bmat[1,ind]
  }else{
    m1=dim(G[-Gnonzeroind,])[1]
    # lmmax=max(norm(t(X)%*%Y,'i'),norm(t(X%*%solve(AN))%*%Y,'i'))
    # lmmax=sqrt(log(2*(2*p+1)/theta)/(2*n))/ga1
    lmmax=10*qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*0.5
    # ANS2=newSVMADMM(Y,X,c(lmnew1*alp*pen1[1:p]*n,rep(50,m1)*n),10*lmmax*G[-Gnonzeroind,],n,p,ItM=3*1e5,eabsn=10^(-8),erel=10^(-8))
    ANS2=newSVMADMM(Y,X,c(lmnew1*alp*pen1[1:p],rep(lmmax,m1)*n),G[-Gnonzeroind,],n,p,ItM=3*1e5,eabsn=10^(-6),erel=10^(-5))
    bet2=ANS2[[1]][-1]
    b02=ANS2[[1]][1]
  }
  cvnew2=targetf1(Yva,Xva,bet2,b02)
  bMatrixm[,M]=c(b02,bet2,cvnew2,alp,lmnew1)
  
  Y_test2<- Xtest%*%bMatrixm[2:(p+1),M] + bMatrixm[1,M]
  Y_test2 <- ifelse(Y_test2 >=0, 1, -1)
  cv_test1[M]= 1-length(which(Y_test2==Ytest))/length(Ytest)
  
}

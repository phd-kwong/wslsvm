
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
# G=G_2

out1=25
t1 <- seq(0, pi, length.out= out1)
beta11 <- 0.1*sin(t1)+ 0.1
out2=25
out3=100
beta=c(rep(0,out3), beta11 ,rep(0,out3),rep(-0.25,out2),rep(0,p-(out1+2*out3+out2)))


# X=as.matrix(read.csv('X_Imd05_gev.csv')[,-1])
# Xva=as.matrix(read.csv('Xva_Imd05_gev.csv')[,-1])
# Xtest=as.matrix(read.csv('Xtest_Imd05_gev.csv')[,-1])
# 
# Y1=as.matrix(read.csv('Y_Imd05_gev.csv')[,-1])
# Yva1=as.matrix(read.csv('Yva_Imd05_gev.csv')[,-1])
# Ytest=as.matrix(read.csv('Ytest_Imd05_gev.csv')[,-1])


# X=as.matrix(read.csv('X_Imd07_gev.csv')[,-1])
# Xva=as.matrix(read.csv('Xva_Imd07_gev.csv')[,-1])
# Xtest=as.matrix(read.csv('Xtest_Imd07_gev.csv')[,-1])
# 
# Y1=as.matrix(read.csv('Y_Imd07_gev.csv')[,-1])
# Yva1=as.matrix(read.csv('Yva_Imd07_gev.csv')[,-1])
# Ytest=as.matrix(read.csv('Ytest_Imd07_gev.csv')[,-1])



m=dim(G)[1]
EG=eigen(t(G)%*%G)[[2]][,-c(1:m)]
AN=rbind(G,t(EG))
ga1 <- 0.95
theta <- 0.05


alpha.v=c(exp(seq(log(0.005),log(0.05),length.out= 5))) #1G
lg=length(alpha.v)
lm1 <- qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*(0.5)

Gnew=rbind(G,t(eigen(t(G)%*%G)[[2]][,-(1:m)]))
Gnew1=t(Gnew)%*%Gnew
bMatrix=matrix(0,p+4,1500)
bMatrixm=bMatrix
cv_test<- matrix(0,1,1500)
cv_test1<- matrix(0,1,1500)
for (M in 1:100) {
  Y <- Y1[,M]
  Yva <- Yva1[,M]
  bmat=matrix(0,p+3,lg)
  for (j in 1:lg) {
    Xwav=X%*%solve(alpha.v[j]*diag(p)+(1-alpha.v[j])*Gnew1)%*%cbind(diag(p),t(Gnew))
    pen=round(sqrt(colMeans(Xwav^2)),4)
    lmnew=lm1
    ANS1=newSVMADMM(Y,X,c(alpha.v[j]*lmnew*pen[1:p],lmnew*pen[p+1:m]*n^(-4/5)),
                    (1-alpha.v[j])*G*n^(4/5),n,p,ItM=3*10^5,eabsn=10^(-6),erel=10^(-6))
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
  cv_test[M]= 1-length(which(Y_test1G==Ytest[,M]))/length(Ytest[,M])
  
  
  Xwav1=X%*%solve(alp*diag(p)+(1-alp)*Gnew1)%*%cbind(diag(p),t(Gnew))
  pen1=round(sqrt(colMeans(Xwav1^2)),4)
  lmnew1=lm1
  
  Gnonzeroind=which(abs(G%*%bmat[2:(p+1),ind])>=0.0001)
  if(length(Gnonzeroind)==0){
    bet2=bmat[2:(p+1),ind]
    b02=bmat[1,ind]
  }else{
    m1=dim(G[-Gnonzeroind,])[1]
    lmmax=10*qnorm(1-theta/(2*(2*p+1)))/(ga1*sqrt(n))*0.5
     ANS2=newSVMADMM(Y,X,c(lmnew1*alp*pen1[1:p],rep(lmmax,m1)*n),G[-Gnonzeroind,],n,p,ItM=3*1e5,eabsn=10^(-6),erel=10^(-5))
    bet2=ANS2[[1]][-1]
    b02=ANS2[[1]][1]
  }
  cvnew2=targetf1(Yva,Xva,bet2,b02)
  bMatrixm[,M]=c(b02,bet2,cvnew2,alp,lmnew1)
  
  Y_test21<- Xtest%*%bMatrixm[2:(p+1),M] + bMatrixm[1,M]
  Y_test2 <- ifelse(Y_test21 >=0, 1, -1)
  cv_test1[M]= 1-length(which(Y_test2==Ytest[,M]))/length(Ytest[,M])

}

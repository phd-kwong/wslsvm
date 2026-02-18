

newSVMADMM<-function(Y,X,lm,G,n,p,numax=1000,numin=0.0001,nuCMIT=30000,ItM=10^6,eabsn=10^(-6),erel=10^(-8)){
  It<-0;ST<-0;
  m=dim(G)[1]
  Xwave <- cbind(rep(1,n),X)
  Gwav=cbind(rep(0,m),sqrt(n)*diag(sqrt(lm[-(1:p)]))%*%G)
  b<-rep(0,p+1);
  # nu<-sqrt(p);eta=2; mu=10
  # nu<-20; eta=0.1; mu=5
  nu<- 1; eta=2; mu=5
  lmG=c(sqrt(lm[-(1:p)]))/nu
  lmb=c(0,sqrt(n)*lm[1:p])/nu

  # nu<-sqrt(p);eta=2; mu=10#2G
  # lmG=c(sqrt(lm[-(1:p)]))/(sqrt(p)*nu)
  # lmb=c(0,sqrt(n)*lm[1:p])/(sqrt(p)*nu)
  
  
  out=ADMMnewsvmC(Xwave, diag(Y), n, diag(p+1), Gwav, as.matrix(c(rep(1,length(Y)))),
                  lmG, lmb, b, b+2, numax,numin,eta, mu, nu, ItM, nuCMIT, eabsn, erel)
  out[[1]]=round(out[[1]],6)
  return(out)
}


newSVMADMMspline<-function(Y,X,lm,G,n,p,numax=1000,numin=0.0001,nuCMIT=30000,ItM=10^6,eabsn=10^(-6),erel=10^(-8)){
  It<-0;ST<-0;
  m=dim(G)[1]
  Xwave <- cbind(rep(1,n),X)
  Gwav=cbind(rep(0,m),sqrt(n)*diag(sqrt(lm[-(1:p)]))%*%G)
  b<-rep(0,p+1);
  # nu<-sqrt(p);eta=2; mu=10
  # nu<- 1.1; eta=1.6; mu=10
  nu<- 1; eta=2; mu=5
  lmG=c(sqrt(lm[-(1:p)]))/nu
  lmb=c(0,sqrt(n)*lm[1:p])/nu
  # lmb=c(lm[1]+lm[p+1],lm[1:p])/nu
  out=ADMMnewsvmsplineC(Xwave, diag(Y), n, diag(p+1), Gwav, as.matrix(c(rep(1,length(Y)))),
                  lmG, lmb, b, b+2, numax,numin,eta, mu, nu, ItM, nuCMIT, eabsn, erel)
  out[[1]]=round(out[[1]],6)
  return(out)
}

newSVMADMMEnet<-function(Y,X,lm,G,n,p,numax=1000,numin=0.0001,nuCMIT=30000,ItM=10^6,eabsn=10^(-6),erel=10^(-8)){
  It<-0;ST<-0;
  m=dim(G)[1]
  Xwave <- cbind(rep(1,n),X)
  Gwav=cbind(rep(0,m),sqrt(n)*diag(sqrt(lm[-(1:p)]))%*%G)
  b<-rep(0,p+1);
  # nu<-sqrt(p);eta=2; mu=10
  # nu<- 1.1; eta=1.6; mu=10
  nu<- 1; eta=2; mu=5
  lmG=c(sqrt(lm[-(1:p)]))/nu
  lmb=c(0,sqrt(n)*lm[1:p])/nu
  # lmb=c(lm[1]+lm[p+1],lm[1:p])/nu
  out=ADMMnewsvmsplineC(Xwave, diag(Y), n, diag(p+1), Gwav, as.matrix(c(rep(1,length(Y)))),
                        lmG, lmb, b, b+2, numax,numin,eta, mu, nu, ItM, nuCMIT, eabsn, erel)
  out[[1]]=round(out[[1]],6)
  return(out)
}


targetf<-function(Y,X,beta,G,lm,b0){
  loss= c(rep(1,length(Y)))-diag(Y)%*%(c(rep(1,length(Y)))*b0+X%*%beta)
  loss=ifelse(loss>0,loss,0)
  loss=loss/length(Y)
  penalty=sum(lm[1:p]*abs(beta))+sum(lm[-c(1:p)]*abs(G%*%beta))
  return(c(loss+penalty,loss,penalty))
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
#include<iostream>
#include<RcppEigen.h>
#include<math.h>
#include<Eigen/Core>
#include<Eigen/Dense>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]

SEXP  ADMMnewsvmC(Eigen::MatrixXd X,Eigen::MatrixXd Y, int n,
                Eigen::MatrixXd Ep,Eigen::MatrixXd G,
                Eigen::VectorXd En,
                Eigen::VectorXd lmG,Eigen::VectorXd lmb,
                Eigen::VectorXd b,Eigen::VectorXd s, double numax,
                double numin, double eta, double mu, double nu,
                int ItM, int nuCMIT  ,double eabsn, double erel) {
  Eigen::MatrixXd XGEI=(X.transpose()*Y.transpose()*Y*X+G.transpose()*G+Ep.transpose()*Ep).inverse();
  Eigen::MatrixXd X1=XGEI*X.transpose()*Y.transpose();
  Eigen::MatrixXd X2=XGEI*Ep;
  Eigen::MatrixXd X3=XGEI*G.transpose();
  Eigen::VectorXd z;
  Eigen::VectorXd gam;
  Eigen::VectorXd u;
  Eigen::VectorXd t;
  Eigen::VectorXd zold;
  Eigen::VectorXd bold;
  Eigen::VectorXd gamold;
  Eigen::VectorXd Znew;
  Eigen::VectorXd betanew;
  Eigen::VectorXd Gamnew;
  Rcpp::NumericVector bV;
  Rcpp::NumericVector lmbV;
  Rcpp::NumericVector gamV;
  Rcpp::NumericVector lmGV;
  Rcpp::NumericVector xminuspenaltyb;
  Rcpp::NumericVector xminuspenaltygam;
  Rcpp::NumericVector signalb;
  Rcpp::NumericVector signalgam;
  Eigen::VectorXd r1;
  Eigen::VectorXd r2;
  Eigen::VectorXd r3;
  
  NumericVector ZupN;
  NumericVector aND;
  Eigen::MatrixXd Ynew;
  
  double ep1;
  double epnew;
  double ep;
  double epm;
  double ed;
  double edm;
  
  Ynew=Y*X;
  z=En-Ynew*b;
  gam=G*b;
  u=z;
  t=gam;
  double it=0;
  for (int It=0;It<ItM;++It) {
    it +=1;
    zold=z;
    bold=b;
    gamold=gam;
    
    aND=wrap(z-u);
    ZupN=ifelse(aND>1/(sqrt(n)*nu),aND-1/(sqrt(n)*nu),ifelse(aND>0,0,aND));
    Znew =as<Eigen::Map<Eigen::VectorXd> >(ZupN);
    
    
    bV=wrap(b-s);
    lmbV=wrap(lmb);
    xminuspenaltyb =ifelse(abs(bV)>lmbV,abs(bV)-lmbV,0);
    signalb=(wrap(sign(bV)));
    bV=xminuspenaltyb*signalb;
    betanew=as<Eigen::Map<Eigen::VectorXd> >(bV);
    
    gamV=wrap(gam-t);
    lmGV=wrap(lmG);
    xminuspenaltygam =ifelse(abs(gamV)>lmGV,abs(gamV)-lmGV,0);
    signalgam=(wrap(sign(gamV)));
    gamV=xminuspenaltygam*signalgam;
    Gamnew=as<Eigen::Map<Eigen::VectorXd> >(gamV);
    
    r1=eta*Znew+(1-eta)*zold+u;
    r2=eta*betanew+(1-eta)*bold+s;
    r3=eta*Gamnew+(1-eta)*gamold+t;
    
    b=X1*(En-r1) +X2*(r2) +X3*(r3);
    z=En-Ynew*b;
    gam=G*b;
    
    
    
    ep=sqrt((z-Znew).dot(z-Znew)+(b-betanew).dot(b-betanew)+
      (gam-Gamnew).dot(gam-Gamnew));
    ep1= sqrt(z.dot(z)+b.dot(b)+gam.dot(gam));
    epnew= sqrt(Znew.dot(Znew)+betanew.dot(betanew)+Gamnew.dot(Gamnew));
    if(ep1>epnew){
      epm=eabsn+erel*ep1;
    }else{
      epm=eabsn+erel*epnew;
    }
    
    ed=nu*sqrt((zold-z).dot(zold-z)+
      (bold-b).dot(bold-b)+
      (gamold-gam).dot(gamold-gam));
    edm=eabsn+erel*nu*sqrt(u.dot(u)+s.dot(s)+t.dot(t));
    
    u=r1-z;
    s=r2-b;
    t=r3-gam;
    
    if(ep<epm&ed<edm){break;}
    if (It % 5000 == 0){Rcpp::checkUserInterrupt();}
    if((ep/epm)>(mu*ed/edm) & It<nuCMIT){
      if(nu<numax){
        eta=0.8;
      }
    }else if(ed/edm>mu*ep/epm & It<nuCMIT){
      if(nu>numin){
        eta=1.2;
      }
    }
    else{eta=1;}
  }
  List output=Rcpp::List::create(Named("b")=wrap(b),
                                 Named("it")=it);
  return (  output);
}





// [[Rcpp::export]]

Eigen::MatrixXd SOVE(Eigen::MatrixXd x1){
  Eigen::MatrixXd res=x1.inverse();
  return res;
}

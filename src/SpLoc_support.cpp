
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void set_seed(unsigned int seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SpLocMeanC(arma::mat& ymat, arma::sp_mat& NNmatrix, int nperm, int s, int side){
  int q=NNmatrix.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  double sd;
  double mn;
  
  arma::mat permU(q, nperm); 
  
  U=NNmatrix*ymat*onevec;  

  arma::mat flipmat(n,nperm); 
  set_seed(s);
  flipmat.randn(); 
  flipmat=sign(flipmat);
  permU=NNmatrix*ymat*flipmat;

  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }
  
  if (side==2){
    permU=permU%permU;
    U=U%U;  
  }

  arma::vec permMax(nperm);
    
  if (side==-1){
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).min();
    }
  } else{
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).max();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SpLocDiffC(arma::mat& ymat, arma::sp_mat& NNmatrix, arma::vec group, int nperm, int s, int side){
  int q=NNmatrix.n_rows;
  int n=group.size();
  arma::vec U(q);
  double sd;
  double mn;

  arma::mat permU(q, nperm);

  U=NNmatrix*ymat*group;  
  
  arma::mat permmat(n,nperm);
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permmat.col(i)=shuffle(group);
  }

  permU=NNmatrix*ymat*permmat;
    
  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }

  if (side==2){
    permU=permU%permU;
    U=U%U;  
  }

  arma::vec permMax(nperm);

  if (side==-1){
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).min();
    }
  } else{
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).max();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List MassiveMeanC(arma::mat ymat, int nperm,  int s, int side){
  int q=ymat.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  double sd;
  double mn;  

  arma::mat permU(q, nperm); 
  
  U=ymat*onevec;  

  arma::mat flipmat(n,nperm); 
  set_seed(s);
  flipmat.randn(); 
  flipmat=sign(flipmat);
  permU=ymat*flipmat;
  
  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }
  
  if (side==2){
    permU=permU%permU;
    U=U%U;  
  }
  
  arma::vec permMax(nperm);

  if (side==-1){
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).min();
    }
  } else{
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).max();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List MassiveDiffC(arma::mat ymat, arma::vec group, int nperm,  int s, int side){
  int q=ymat.n_rows;
  int n=group.size();
  arma::vec U(q);
  double sd;
  double mn;
  arma::mat permU(q, nperm);

  U=ymat*group;  

  arma::mat permmat(n,nperm);
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permmat.col(i)=shuffle(group);
  }
  permU=ymat*permmat;
  

  for (int k=0; k<q; ++k){
    mn=mean(permU.row(k));
    sd=stddev(permU.row(k));
    permU.row(k)=(permU.row(k)-mn)/sd;
    U(k)=(U(k)-mn)/sd;
  }

  if (side==2){
    permU=permU%permU;
    U=U%U;  
  }

  arma::vec permMax(nperm);

  if (side==-1){
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).min();
    }
  } else{
    for (int i=0; i<nperm; ++i){
      permMax(i)=permU.col(i).max();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}

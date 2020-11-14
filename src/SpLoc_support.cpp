
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;


arma::vec avg_rank(arma::vec x) {
  arma::uvec w = arma::stable_sort_index(x, "descend");
  R_xlen_t sz = x.size();
  arma::vec r(sz);
  
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      // r[w[i + k]] = i + (n + 1) / 2.;
      r[w[i + k]] = i + (n ) ;
    }
  }
  return r;
}

// [[Rcpp::export]]
double quantileC(arma::vec Tstatvec, double alpha){
  int n=Tstatvec.size();
  arma::vec avgrank=avg_rank(Tstatvec);
  int thres=floor(n*alpha);
  int index=0;
  for (int i=0; i<n; ++i){
    if (avgrank(i)==thres){
      index=i;
      break;
    }
  }
  double out=Tstatvec(index);
  return out;
}

// [[Rcpp::export]]
void set_seed(unsigned int seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SpLocMeanC(arma::mat& ymat, arma::sp_mat& NNmatrix, int nperm, double alpha, int s){
  int q=NNmatrix.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  double sd;
  
  arma::mat permU(q, nperm); 
  
  U=NNmatrix*ymat*onevec;  

  arma::mat flipmat(n,nperm); 
  set_seed(s);
  flipmat.randn(); 
  flipmat=sign(flipmat);
  permU=NNmatrix*ymat*flipmat;
  
  for (int k=0; k<q; ++k){
    sd=stddev(permU.row(k));
    permU.row(k)=permU.row(k)/sd;
    U(k)=U(k)/sd;
  }
  
  permU=permU%permU;
  U=U%U;  
  
  arma::vec permMax(nperm);
  for (int i=0; i<nperm; ++i){
    permMax(i)=permU.col(i).max();
  }
  
  double qt=quantileC(permMax, alpha);  
  
  return Rcpp::List::create(Rcpp::Named("threshold")=qt,
                            Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List SpLocDiffC(arma::mat& ymat, arma::sp_mat& NNmatrix, arma::vec group, int nperm, double alpha, int s){
  int q=NNmatrix.n_rows;
  int n=group.size();
  arma::vec U(q);
  double sd;
  arma::mat permU(q, nperm);

  U=NNmatrix*ymat*group;  
  
  arma::mat permmat(n,nperm);
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permmat.col(i)=shuffle(group);
  }

  permU=NNmatrix*ymat*permmat;
  
  for (int k=0; k<q; ++k){
    sd=stddev(permU.row(k));
    permU.row(k)=permU.row(k)/sd;
    U(k)=U(k)/sd;
  }

  permU=permU%permU;
  U=U%U;

  arma::vec permMax(nperm);
  for (int i=0; i<nperm; ++i){
    permMax(i)=permU.col(i).max();
  }
  
  double qt=quantileC(permMax, alpha);
  
  return Rcpp::List::create(Rcpp::Named("threshold")=qt,
                            Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List MassiveMeanC(arma::mat ymat, int nperm, double alpha, int s){
  int q=ymat.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);
  arma::vec U(q);
  double sd;
  
  arma::mat permU(q, nperm); 
  
  U=ymat*onevec;  

  arma::mat flipmat(n,nperm); 
  set_seed(s);
  flipmat.randn(); 
  flipmat=sign(flipmat);
  permU=ymat*flipmat;
  
  for (int k=0; k<q; ++k){
    sd=stddev(permU.row(k));
    permU.row(k)=permU.row(k)/sd;
    U(k)=U(k)/sd;
  }
  
  permU=permU%permU;
  U=U%U;  
  
  arma::vec permMax(nperm);
  for (int i=0; i<nperm; ++i){
    permMax(i)=permU.col(i).max();
  }
  
  double qt=quantileC(permMax, alpha);  
  
  return Rcpp::List::create(Rcpp::Named("threshold")=qt,
                            Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List MassiveDiffC(arma::mat ymat, arma::vec group, int nperm, double alpha, int s){
  int q=ymat.n_rows;
  int n=group.size();
  arma::vec U(q);
  double sd;
  arma::mat permU(q, nperm);

  U=ymat*group;  
  
  arma::mat permmat(n,nperm);
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permmat.col(i)=shuffle(group);
  }

  permU=ymat*permmat;
  
  for (int k=0; k<q; ++k){
    sd=stddev(permU.row(k));
    permU.row(k)=permU.row(k)/sd;
    U(k)=U(k)/sd;
  }

  permU=permU%permU;
  U=U%U;

  arma::vec permMax(nperm);
  for (int i=0; i<nperm; ++i){
    permMax(i)=permU.col(i).max();
  }
  
  double qt=quantileC(permMax, alpha);
  
  return Rcpp::List::create(Rcpp::Named("threshold")=qt,
                            Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax,
                            Rcpp::Named("nperm")=nperm);
}


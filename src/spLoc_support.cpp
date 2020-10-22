#define ARMA_64BIT_WORD 1


#include <armadillo>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
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
Rcpp::List SpLocC(arma::sp_mat NN, arma::mat ymat, int nperm, double alpha){
  int q=NN.n_rows;
  int p=ymat.n_rows;
  int n=ymat.n_cols;
  arma::mat permU(q,nperm);
  arma::vec permy(p);
  arma::vec rand(n); 
  arma::vec sdvec(q); sdvec.fill(0);
  arma::vec y(p);y.fill(0);
  arma::vec U(q);

  for (int subj=0; subj<n; ++subj){
    y=y+ymat.col(subj);
  }
  U=NN*y;  
  
  for (int i=0; i<nperm; ++i){
    permy.fill(0);
    rand.randn();
    rand=rand/abs(rand);
    for (int subj=0; subj<n; ++subj){
      permy=permy+rand(subj)*ymat.col(subj);
    }
    permU.col(i)=NN*permy;
  }
  
  for (int k=0; k<q; ++k){
    double sd=stddev(permU.row(k));
    permU.row(k)=permU.row(k)/sd;
    U(k)=U(k)/sd;
  }
  permU=permU%permU;
  
  arma::vec permMax(nperm);
  for (int i=0; i<nperm; ++i){
    permMax(i)=permU.col(i).max();
  }
  
  double qt=quantileC(permMax, alpha);
  
  U=U%U;

  return Rcpp::List::create(Rcpp::Named("threshold")=qt,
                            Rcpp::Named("Tstat")=U,
                            Rcpp::Named("permMax")=permMax);
}




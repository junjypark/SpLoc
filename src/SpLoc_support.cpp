
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

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
Rcpp::List SpLocMeanC(arma::sp_mat& NNmatrix, arma::mat& ymat, int nperm, double alpha, int s, SEXP pU){
  int q=NNmatrix.n_rows;
  int n=ymat.n_cols;
  arma::vec onevec(n); onevec.fill(1);

  arma::vec U(q);
  double sd;
  
  XPtr<BigMatrix> xpMat(pU);
  arma::mat permU = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
  
  U=NNmatrix*ymat*onevec;  

  arma::mat randmat(n,nperm); 
  set_seed(s);
  randmat.randn(); 
  randmat=sign(randmat);
  permU=NNmatrix*ymat*randmat;
  
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
Rcpp::List SpLocMeanC2(arma::sp_mat& NNmatrix, arma::mat& ymat, int nperm, double alpha, int s, SEXP pU, SEXP pY){
  int q=NNmatrix.n_rows;
  int n=ymat.n_cols;
  arma::vec rand(n); rand.fill(1);
  arma::vec U(q);
  double sd;

  XPtr<BigMatrix> xpMat(pU);
  arma::mat permU = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);

  XPtr<BigMatrix> xp2Mat(pY);
  arma::mat permY = arma::Mat<double> ( (double *)xp2Mat->matrix(), xp2Mat->nrow(), xp2Mat->ncol(), false);

  U=NNmatrix*ymat*rand;  
  
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    rand.randn();
    rand=sign(rand);
    permY.col(i)=ymat*rand;
  }

  permU=NNmatrix*permY;
  
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
Rcpp::List SpLocDiffC(arma::sp_mat& NNmatrix, arma::mat& ymat, arma::vec group, int nperm, double alpha, int s, SEXP pU){
  int q=NNmatrix.n_rows;
  int n=group.size();
  arma::mat permgroup(n,nperm);
  arma::vec U(q);
  double sd;

  XPtr<BigMatrix> xpMat(pU);
  arma::mat permU = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);

  U=NNmatrix*ymat*group;  
  
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permgroup.col(i)=shuffle(group);
  }

  permU=NNmatrix*ymat*permgroup;
  
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
Rcpp::List SpLocDiffC2(arma::sp_mat& NNmatrix, arma::mat& ymat, arma::vec group, int nperm, double alpha, int s, SEXP pU, SEXP pY){
  int q=NNmatrix.n_rows;
  int p=group.size();
  arma::vec permgroup(p);
  arma::vec U(q);
  double sd;

  XPtr<BigMatrix> xpMat(pU);
  arma::mat permU = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);

  XPtr<BigMatrix> xp2Mat(pY);
  arma::mat permY = arma::Mat<double> ( (double *)xp2Mat->matrix(), xp2Mat->nrow(), xp2Mat->ncol(), false);

  U=NNmatrix*ymat*group;  
  
  set_seed(s);
  for (int i=0; i<nperm; ++i){
    permgroup=shuffle(group);
    permY.col(i)=ymat*permgroup;
  }

  permU=NNmatrix*permY;
  
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





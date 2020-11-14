// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// quantileC
double quantileC(arma::vec Tstatvec, double alpha);
RcppExport SEXP _SpLoc_quantileC(SEXP TstatvecSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Tstatvec(TstatvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(quantileC(Tstatvec, alpha));
    return rcpp_result_gen;
END_RCPP
}
// set_seed
void set_seed(unsigned int seed);
RcppExport SEXP _SpLoc_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// SpLocMeanC
Rcpp::List SpLocMeanC(arma::sp_mat& NNmatrix, arma::mat& ymat, int nperm, double alpha, int s);
RcppExport SEXP _SpLoc_SpLocMeanC(SEXP NNmatrixSEXP, SEXP ymatSEXP, SEXP npermSEXP, SEXP alphaSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocMeanC(NNmatrix, ymat, nperm, alpha, s));
    return rcpp_result_gen;
END_RCPP
}
// SpLocDiffC
Rcpp::List SpLocDiffC(arma::sp_mat& NNmatrix, arma::mat& ymat, arma::vec group, int nperm, double alpha, int s);
RcppExport SEXP _SpLoc_SpLocDiffC(SEXP NNmatrixSEXP, SEXP ymatSEXP, SEXP groupSEXP, SEXP npermSEXP, SEXP alphaSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocDiffC(NNmatrix, ymat, group, nperm, alpha, s));
    return rcpp_result_gen;
END_RCPP
}
// SpLocMeanC2
Rcpp::List SpLocMeanC2(arma::sp_mat& NNmatrix, arma::mat& ymat, int nperm, double alpha, int s, SEXP pU);
RcppExport SEXP _SpLoc_SpLocMeanC2(SEXP NNmatrixSEXP, SEXP ymatSEXP, SEXP npermSEXP, SEXP alphaSEXP, SEXP sSEXP, SEXP pUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pU(pUSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocMeanC2(NNmatrix, ymat, nperm, alpha, s, pU));
    return rcpp_result_gen;
END_RCPP
}
// SpLocDiffC2
Rcpp::List SpLocDiffC2(arma::sp_mat& NNmatrix, arma::mat& ymat, arma::vec group, int nperm, double alpha, int s, SEXP pU);
RcppExport SEXP _SpLoc_SpLocDiffC2(SEXP NNmatrixSEXP, SEXP ymatSEXP, SEXP groupSEXP, SEXP npermSEXP, SEXP alphaSEXP, SEXP sSEXP, SEXP pUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pU(pUSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocDiffC2(NNmatrix, ymat, group, nperm, alpha, s, pU));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpLoc_quantileC", (DL_FUNC) &_SpLoc_quantileC, 2},
    {"_SpLoc_set_seed", (DL_FUNC) &_SpLoc_set_seed, 1},
    {"_SpLoc_SpLocMeanC", (DL_FUNC) &_SpLoc_SpLocMeanC, 5},
    {"_SpLoc_SpLocDiffC", (DL_FUNC) &_SpLoc_SpLocDiffC, 6},
    {"_SpLoc_SpLocMeanC2", (DL_FUNC) &_SpLoc_SpLocMeanC2, 6},
    {"_SpLoc_SpLocDiffC2", (DL_FUNC) &_SpLoc_SpLocDiffC2, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpLoc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

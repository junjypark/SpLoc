// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
Rcpp::List SpLocMeanC(arma::mat& ymat, arma::sp_mat& NNmatrix, int nperm, int s, int side);
RcppExport SEXP _SpLoc_SpLocMeanC(SEXP ymatSEXP, SEXP NNmatrixSEXP, SEXP npermSEXP, SEXP sSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocMeanC(ymat, NNmatrix, nperm, s, side));
    return rcpp_result_gen;
END_RCPP
}
// SpLocMeanC2
Rcpp::List SpLocMeanC2(arma::mat& ymat, arma::sp_mat& NNmatrix, int nperm, int s, int side);
RcppExport SEXP _SpLoc_SpLocMeanC2(SEXP ymatSEXP, SEXP NNmatrixSEXP, SEXP npermSEXP, SEXP sSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocMeanC2(ymat, NNmatrix, nperm, s, side));
    return rcpp_result_gen;
END_RCPP
}
// SpLocDiffC
Rcpp::List SpLocDiffC(arma::mat& ymat, arma::sp_mat& NNmatrix, arma::vec group, int nperm, int s, int side);
RcppExport SEXP _SpLoc_SpLocDiffC(SEXP ymatSEXP, SEXP NNmatrixSEXP, SEXP groupSEXP, SEXP npermSEXP, SEXP sSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type NNmatrix(NNmatrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(SpLocDiffC(ymat, NNmatrix, group, nperm, s, side));
    return rcpp_result_gen;
END_RCPP
}
// MassiveMeanC
Rcpp::List MassiveMeanC(arma::mat ymat, int nperm, int s, int side);
RcppExport SEXP _SpLoc_MassiveMeanC(SEXP ymatSEXP, SEXP npermSEXP, SEXP sSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(MassiveMeanC(ymat, nperm, s, side));
    return rcpp_result_gen;
END_RCPP
}
// MassiveDiffC
Rcpp::List MassiveDiffC(arma::mat ymat, arma::vec group, int nperm, int s, int side);
RcppExport SEXP _SpLoc_MassiveDiffC(SEXP ymatSEXP, SEXP groupSEXP, SEXP npermSEXP, SEXP sSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(MassiveDiffC(ymat, group, nperm, s, side));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpLoc_set_seed", (DL_FUNC) &_SpLoc_set_seed, 1},
    {"_SpLoc_SpLocMeanC", (DL_FUNC) &_SpLoc_SpLocMeanC, 5},
    {"_SpLoc_SpLocMeanC2", (DL_FUNC) &_SpLoc_SpLocMeanC2, 5},
    {"_SpLoc_SpLocDiffC", (DL_FUNC) &_SpLoc_SpLocDiffC, 6},
    {"_SpLoc_MassiveMeanC", (DL_FUNC) &_SpLoc_MassiveMeanC, 4},
    {"_SpLoc_MassiveDiffC", (DL_FUNC) &_SpLoc_MassiveDiffC, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpLoc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

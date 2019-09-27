// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// soft_thr
arma::colvec soft_thr(arma::colvec input, double lambda, int pos, arma::rowvec sp_locs);
RcppExport SEXP _ssJIVE_soft_thr(SEXP inputSEXP, SEXP lambdaSEXP, SEXP posSEXP, SEXP sp_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type input(inputSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type pos(posSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type sp_locs(sp_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_thr(input, lambda, pos, sp_locs));
    return rcpp_result_gen;
END_RCPP
}
// sfpca_bic
List sfpca_bic(arma::mat x, int K, arma::mat Omegv, arma::rowvec sp_locsv, arma::rowvec lamvs, arma::rowvec alphavs);
RcppExport SEXP _ssJIVE_sfpca_bic(SEXP xSEXP, SEXP KSEXP, SEXP OmegvSEXP, SEXP sp_locsvSEXP, SEXP lamvsSEXP, SEXP alphavsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omegv(OmegvSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type sp_locsv(sp_locsvSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type lamvs(lamvsSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type alphavs(alphavsSEXP);
    rcpp_result_gen = Rcpp::wrap(sfpca_bic(x, K, Omegv, sp_locsv, lamvs, alphavs));
    return rcpp_result_gen;
END_RCPP
}
// soft_thr
arma::colvec soft_thr(arma::colvec input, double lambda, int pos, arma::rowvec sp_locs);
RcppExport SEXP _ssJIVE_soft_thr(SEXP inputSEXP, SEXP lambdaSEXP, SEXP posSEXP, SEXP sp_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type input(inputSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type pos(posSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type sp_locs(sp_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_thr(input, lambda, pos, sp_locs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ssJIVE_soft_thr", (DL_FUNC) &_ssJIVE_soft_thr, 4},
    {"_ssJIVE_sfpca_bic", (DL_FUNC) &_ssJIVE_sfpca_bic, 6},
    {"_ssJIVE_soft_thr", (DL_FUNC) &_ssJIVE_soft_thr, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ssJIVE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
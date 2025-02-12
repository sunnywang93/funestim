// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// epaKernelSmoothingCurve
arma::vec epaKernelSmoothingCurve(const arma::vec& U, const arma::vec& T, const arma::vec& Y, const arma::vec& b, const double& n_obs_min);
RcppExport SEXP _funestim_epaKernelSmoothingCurve(SEXP USEXP, SEXP TSEXP, SEXP YSEXP, SEXP bSEXP, SEXP n_obs_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type n_obs_min(n_obs_minSEXP);
    rcpp_result_gen = Rcpp::wrap(epaKernelSmoothingCurve(U, T, Y, b, n_obs_min));
    return rcpp_result_gen;
END_RCPP
}
// uniKernelSmoothingCurve
arma::vec uniKernelSmoothingCurve(const arma::vec& U, const arma::vec& T, const arma::vec& Y, const arma::vec& b, const double& n_obs_min);
RcppExport SEXP _funestim_uniKernelSmoothingCurve(SEXP USEXP, SEXP TSEXP, SEXP YSEXP, SEXP bSEXP, SEXP n_obs_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type n_obs_min(n_obs_minSEXP);
    rcpp_result_gen = Rcpp::wrap(uniKernelSmoothingCurve(U, T, Y, b, n_obs_min));
    return rcpp_result_gen;
END_RCPP
}
// biweightKernelSmoothingCurve
arma::vec biweightKernelSmoothingCurve(const arma::vec& U, const arma::vec& T, const arma::vec& Y, const arma::vec& b, const double& n_obs_min);
RcppExport SEXP _funestim_biweightKernelSmoothingCurve(SEXP USEXP, SEXP TSEXP, SEXP YSEXP, SEXP bSEXP, SEXP n_obs_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type n_obs_min(n_obs_minSEXP);
    rcpp_result_gen = Rcpp::wrap(biweightKernelSmoothingCurve(U, T, Y, b, n_obs_min));
    return rcpp_result_gen;
END_RCPP
}
// estimateSigma
double estimateSigma(const List& curves);
RcppExport SEXP _funestim_estimateSigma(SEXP curvesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type curves(curvesSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateSigma(curves));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_funestim_epaKernelSmoothingCurve", (DL_FUNC) &_funestim_epaKernelSmoothingCurve, 5},
    {"_funestim_uniKernelSmoothingCurve", (DL_FUNC) &_funestim_uniKernelSmoothingCurve, 5},
    {"_funestim_biweightKernelSmoothingCurve", (DL_FUNC) &_funestim_biweightKernelSmoothingCurve, 5},
    {"_funestim_estimateSigma", (DL_FUNC) &_funestim_estimateSigma, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_funestim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

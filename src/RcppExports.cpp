// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getMLD
std::map<int,double> getMLD(int& L, int& K1, double& rho1, double& rho2, IntegerVector& times, NumericVector& minus_log_runif_L, int& nT, double& dt);
RcppExport SEXP _SimulateHAT_getMLD(SEXP LSEXP, SEXP K1SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP timesSEXP, SEXP minus_log_runif_LSEXP, SEXP nTSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< int& >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< double& >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< double& >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type times(timesSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type minus_log_runif_L(minus_log_runif_LSEXP);
    Rcpp::traits::input_parameter< int& >::type nT(nTSEXP);
    Rcpp::traits::input_parameter< double& >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(getMLD(L, K1, rho1, rho2, times, minus_log_runif_L, nT, dt));
    return rcpp_result_gen;
END_RCPP
}
// SimulateOnce
IntegerVector SimulateOnce(long int& L, long int& T, long int& K, long int seed);
RcppExport SEXP _SimulateHAT_SimulateOnce(SEXP LSEXP, SEXP TSEXP, SEXP KSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< long int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< long int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< long int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(SimulateOnce(L, T, K, seed));
    return rcpp_result_gen;
END_RCPP
}
// SimulateIndirect
IntegerVector SimulateIndirect(long int& L, long int& T, long int& K, long int seed);
RcppExport SEXP _SimulateHAT_SimulateIndirect(SEXP LSEXP, SEXP TSEXP, SEXP KSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< long int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< long int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< long int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(SimulateIndirect(L, T, K, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SimulateHAT_getMLD", (DL_FUNC) &_SimulateHAT_getMLD, 8},
    {"_SimulateHAT_SimulateOnce", (DL_FUNC) &_SimulateHAT_SimulateOnce, 4},
    {"_SimulateHAT_SimulateIndirect", (DL_FUNC) &_SimulateHAT_SimulateIndirect, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SimulateHAT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

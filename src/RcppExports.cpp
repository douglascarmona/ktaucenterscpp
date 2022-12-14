// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// distance_to_centers
List distance_to_centers(NumericMatrix data, NumericMatrix centers);
RcppExport SEXP _ktaucenterscpp_distance_to_centers(SEXP dataSEXP, SEXP centersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    rcpp_result_gen = Rcpp::wrap(distance_to_centers(data, centers));
    return rcpp_result_gen;
END_RCPP
}
// normal_consistency_constants
double normal_consistency_constants(int p);
RcppExport SEXP _ktaucenterscpp_normal_consistency_constants(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_consistency_constants(p));
    return rcpp_result_gen;
END_RCPP
}
// c1
double c1();
RcppExport SEXP _ktaucenterscpp_c1() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(c1());
    return rcpp_result_gen;
END_RCPP
}
// c2
double c2(int p);
RcppExport SEXP _ktaucenterscpp_c2(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(c2(p));
    return rcpp_result_gen;
END_RCPP
}
// median_cpp
double median_cpp(NumericVector x);
RcppExport SEXP _ktaucenterscpp_median_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(median_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// mscale
double mscale(NumericVector u, double c, double b);
RcppExport SEXP _ktaucenterscpp_mscale(SEXP uSEXP, SEXP cSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mscale(u, c, b));
    return rcpp_result_gen;
END_RCPP
}
// rho_opt
NumericVector rho_opt(NumericVector t, double c);
RcppExport SEXP _ktaucenterscpp_rho_opt(SEXP tSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_opt(t, c));
    return rcpp_result_gen;
END_RCPP
}
// psi_opt
NumericVector psi_opt(NumericVector t, double c);
RcppExport SEXP _ktaucenterscpp_psi_opt(SEXP tSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(psi_opt(t, c));
    return rcpp_result_gen;
END_RCPP
}
// derpsi_opt
NumericVector derpsi_opt(NumericVector t, double c);
RcppExport SEXP _ktaucenterscpp_derpsi_opt(SEXP tSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(derpsi_opt(t, c));
    return rcpp_result_gen;
END_RCPP
}
// tau_scale
double tau_scale(NumericVector u, double c, double s, double b);
RcppExport SEXP _ktaucenterscpp_tau_scale(SEXP uSEXP, SEXP cSEXP, SEXP sSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(tau_scale(u, c, s, b));
    return rcpp_result_gen;
END_RCPP
}
// wni
NumericVector wni(NumericVector u, double c1, double c2, double s, double b1, double b2);
RcppExport SEXP _ktaucenterscpp_wni(SEXP uSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP sSEXP, SEXP b1SEXP, SEXP b2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< double >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    rcpp_result_gen = Rcpp::wrap(wni(u, c1, c2, s, b1, b2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ktaucenterscpp_distance_to_centers", (DL_FUNC) &_ktaucenterscpp_distance_to_centers, 2},
    {"_ktaucenterscpp_normal_consistency_constants", (DL_FUNC) &_ktaucenterscpp_normal_consistency_constants, 1},
    {"_ktaucenterscpp_c1", (DL_FUNC) &_ktaucenterscpp_c1, 0},
    {"_ktaucenterscpp_c2", (DL_FUNC) &_ktaucenterscpp_c2, 1},
    {"_ktaucenterscpp_median_cpp", (DL_FUNC) &_ktaucenterscpp_median_cpp, 1},
    {"_ktaucenterscpp_mscale", (DL_FUNC) &_ktaucenterscpp_mscale, 3},
    {"_ktaucenterscpp_rho_opt", (DL_FUNC) &_ktaucenterscpp_rho_opt, 2},
    {"_ktaucenterscpp_psi_opt", (DL_FUNC) &_ktaucenterscpp_psi_opt, 2},
    {"_ktaucenterscpp_derpsi_opt", (DL_FUNC) &_ktaucenterscpp_derpsi_opt, 2},
    {"_ktaucenterscpp_tau_scale", (DL_FUNC) &_ktaucenterscpp_tau_scale, 4},
    {"_ktaucenterscpp_wni", (DL_FUNC) &_ktaucenterscpp_wni, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ktaucenterscpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

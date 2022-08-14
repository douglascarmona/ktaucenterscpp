#ifndef RHO_OPT_H
#define RHO_OPT_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector rho_opt(NumericVector, double);
NumericVector psi_opt(NumericVector, double);
NumericVector derpsi_opt(NumericVector, double);

#endif
#ifndef RHO_OPT_H
#define RHO_OPT_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector rho_opt();
NumericVector psi_opt();
NumericVector derpsi_opt();

#endif
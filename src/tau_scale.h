#ifndef TAU_SCALE_H
#define TAU_SCALE_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector tau_scale(NumericVector, double, double, double);
NumericVector wni(NumericVector, double, double, double, double, double);

#endif
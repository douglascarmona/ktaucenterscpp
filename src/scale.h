#ifndef SCALE_H
#define SCALE_H

#include <Rcpp.h>
using namespace Rcpp;

double normal_consistency_constants(int);
double const_c2(std::size_t);
double median_cpp(NumericVector);
double mscale(NumericVector, const double, const double);
double tau_scale(NumericVector, const double, const double);
NumericVector wni(NumericVector, double, double, double);
NumericVector get_weights(NumericVector, IntegerVector);
NumericMatrix get_new_centers(NumericMatrix, NumericVector, IntegerVector,
                              NumericVector);

#endif
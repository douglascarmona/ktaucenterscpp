#ifndef SCALE_H
#define SCALE_H

#include <Rcpp.h>
using namespace Rcpp;

double normal_consistency_constants(int);
double const_c1();
double const_c2(std::size_t);
double median_cpp(NumericVector);
double mscale(NumericVector, double, double);
double tau_scale(NumericVector, double, double);
NumericVector wni(NumericVector, double, double, double);
NumericVector get_weights(NumericVector, IntegerVector);
NumericMatrix get_new_centers(NumericMatrix, NumericVector, IntegerVector,
                              const int, NumericVector);

#endif
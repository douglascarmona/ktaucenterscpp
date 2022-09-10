#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

// IntegerVector tabulatecpp(const IntegerVector, const unsigned);
// NumericMatrix matrix_mult(NumericMatrix, NumericVector);
// NumericMatrix row_sum_cpp(NumericMatrix, IntegerVector);
IntegerVector top_index(NumericVector, int);
double median_cpp(NumericVector);
IntegerVector tabulatecpp(IntegerVector, const unsigned);
double max_tolerance(NumericMatrix, NumericMatrix);
#endif
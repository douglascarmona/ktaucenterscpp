#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

IntegerVector top_index(NumericVector, int);
double median_cpp(NumericVector);
IntegerVector tabulatecpp(IntegerVector, const std::size_t);
double max_tolerance(NumericMatrix, NumericMatrix);
#endif
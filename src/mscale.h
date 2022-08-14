#ifndef MSCALE_H
#define MSCALE_H

#include <Rcpp.h>
using namespace Rcpp;

double normal_consistency_constants(int);
double c1();
double c2(int);
double median_cpp(NumericVector);
double mscale(NumericVector, double, double);

#endif
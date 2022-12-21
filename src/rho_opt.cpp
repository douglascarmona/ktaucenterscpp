#include "rho_opt.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Quasi optimal rho function following reference [1]
//'
//' @param x numeric vector.
//' @param c tunning constant.
//'
//' @return
//' Numeric vector with with quasi optimal rho computation for each element of
//' x.
//'
//'@references
//' [1] Salibian-Barrera, M., Willems, G., & Zamar, R. (2008). The fast-tau
//' estimator for regression. Journal of Computational and GraphicalStatistics,
//' 17(3), 659-682.
//'
// [[Rcpp::export]]
NumericVector rho_opt(NumericVector x, double c) {
  NumericVector out(no_init(x.size()));

  std::size_t idx = 0;
  for (const auto &x_it : x) {
    if (fabs(x_it) <= 2 * c) {
      out[idx] = 0.5 * pow(x_it, 2) / (3.25 * pow(c, 2));
    } else if (fabs(x_it) <= 3 * c) {
      out[idx] =
          (1.792 - 0.972 * pow(x_it, 2) / pow(c, 2) +
           0.432 * pow(x_it, 4) / pow(c, 4) - 0.052 * pow(x_it, 6) / pow(c, 6) +
           0.002 * pow(x_it, 8) / pow(c, 8)) /
          3.25;
    } else {
      out[idx] = 1.0;
    }
    ++idx;
  }

  return out;
}

//' Implementation of the derivative of quasi optimal rho function
//'
//' @param x a numeric vector
//' @param c a tunning constant.
//'
//' @return
//' Numeric vector with with the derivative of the quasi optimal rho computation
//' for each element of x.
//'
// [[Rcpp::export]]
NumericVector psi_opt(NumericVector x, double c) {
  NumericVector out(no_init(x.size()));

  std::size_t idx = 0;
  for (const auto &x_it : x) {
    if (fabs(x_it) <= 2 * c) {
      out[idx] = x_it / (3.25 * pow(c, 2));
    } else if (fabs(x_it) <= 3 * c) {
      out[idx] = (-1.944 * x_it / pow(c, 2) + 1.728 * pow(x_it, 3) / pow(c, 4) -
                  0.312 * pow(x_it, 5) / pow(c, 6) +
                  0.016 * pow(x_it, 7) / pow(c, 8)) /
                 3.25;
    } else {
      out[idx] = 0.0;
    }
    ++idx;
  }

  return out;
}

//' Implementation of the second derivative of the rho function
//'
//' @param x a numeric vector
//' @param c a tunning constant.
//'
//' @return
//' Numeric vector with with the second derivative of the quasi optimal rho
//' computation for each element of x.
//'
// [[Rcpp::export]]
NumericVector derpsi_opt(NumericVector x, double c) {
  NumericVector out(no_init(x.size()));

  std::size_t idx = 0;
  for (const auto &x_it : x) {
    if (fabs(x_it) <= 2 * c) {
      out[idx] = 1.0 / (3.25 * pow(c, 2));
    } else if (fabs(x_it) <= 3 * c) {
      out[idx] =
          (-1.944 / pow(c, 2) + 5.184 * pow(x_it, 2) / pow(c, 4) -
           1.56 * pow(x_it, 4) / pow(c, 6) + 0.112 * pow(x_it, 6) / pow(c, 8)) /
          3.25;
    } else {
      out[idx] = 0.0;
    }
    ++idx;
  }

  return out;
}

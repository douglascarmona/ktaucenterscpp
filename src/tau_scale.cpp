#include "mscale.h"
#include "rho_opt.h"
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

// TODO: Add docs
//'rho_opt function
//'@export
// [[Rcpp::export]]
double tau_scale(NumericVector u, double c, double s, double b) {
  return s * sqrt(mean(rho_opt(u / s, c))) / sqrt(b);
}

// TODO: Add docs
//'Wni function
//'@export
// [[Rcpp::export]]
NumericVector wni(NumericVector u, double c1, double c2, double s, double b1,
                  double b2) {

  NumericVector dnor = u / s;
  double A = mean(2 * rho_opt(dnor, c2) - psi_opt(dnor, c2) * dnor);
  double B = mean(psi_opt(dnor, c1) * dnor);
  return ifelse(u == 0.0, A * derpsi_opt(0.0, c1) + B * derpsi_opt(0.0, c2),
                (A * psi_opt(dnor, c1) + B * psi_opt(dnor, c2)) / dnor);
}
#include "rho_opt.h"
#include <Rcpp.h>
using namespace Rcpp;

//'rho_opt function
//'
//' An implementation of quasi optimal rho functions following reference [1]
//'
//'@param t numeric vector.
//'@param c tunning constant.
//'@return rho(\code{t}/\code{c})
//'@examples val <- rho_opt(t = 0.5, c = 1.0)
//'@references [1] Salibian-Barrera, M., Willems, G., & Zamar,
//'   R. (2008). The fast-tau estimator for regression. Journal of
//'   Computational and Graphical Statistics, 17(3), 659-682.
//'
//'@export
// [[Rcpp::export]]
NumericVector rho_opt(NumericVector t, double c) {
  NumericVector out(Rf_allocVector(REALSXP, t.size()));

  R_xlen_t idx = 0;
  for (const auto &t_it : t) {
    if (fabs(t_it) <= 2 * c) {
      out[idx] = 0.5 * pow(t_it, 2) / (3.25 * pow(c, 2));
    } else if (fabs(t_it) <= 3 * c) {
      out[idx] =
          (1.792 - 0.972 * pow(t_it, 2) / pow(c, 2) +
           0.432 * pow(t_it, 4) / pow(c, 4) - 0.052 * pow(t_it, 6) / pow(c, 6) +
           0.002 * pow(t_it, 8) / pow(c, 8)) /
          3.25;
    } else {
      out[idx] = 1.0;
    }
    ++idx;
  }

  return out;
}

//'psi_opt function
//'
//' An implementation of the derivative of quasi optimal rho function
//'
//'@param t a numeric vector
//'@param c a tunning constant.
//'@return psi(\code{t}/ \code{c})
//'@examples val <- psi_opt(t = 0.5, c = 1)
//'
//'@export
// [[Rcpp::export]]
NumericVector psi_opt(NumericVector t, double c) {
  NumericVector out(Rf_allocVector(REALSXP, t.size()));

  R_xlen_t idx = 0;
  for (const auto &t_it : t) {
    if (fabs(t_it) <= 2 * c) {
      out[idx] = t_it / (3.25 * pow(c, 2));
    } else if (fabs(t_it) <= 3 * c) {
      out[idx] = (-1.944 * t_it / pow(c, 2) + 1.728 * pow(t_it, 3) / pow(c, 4) -
                  0.312 * pow(t_it, 5) / pow(c, 6) +
                  0.016 * pow(t_it, 7) / pow(c, 8)) /
                 3.25;
    } else {
      out[idx] = 0.0;
    }
    ++idx;
  }

  return out;
}

//' derpsiOpt
//' the derivative of the psi function
//'
//' @param t a numeric vector
//' @param c a tunning constant.
//' @return rho'(\code{x}/ \code{c})
//' @examples val <- derpsi_opt(t = 0.5, c = 1.0)
//'
//' @export
// [[Rcpp::export]]
NumericVector derpsi_opt(NumericVector t, double c) {
  NumericVector out(Rf_allocVector(REALSXP, t.size()));

  R_xlen_t idx = 0;
  for (const auto &t_it : t) {
    if (fabs(t_it) <= 2 * c) {
      out[idx] = 1.0 / (3.25 * pow(c, 2));
    } else if (fabs(t_it) <= 3 * c) {
      out[idx] =
          (-1.944 / pow(c, 2) + 5.184 * pow(t_it, 2) / pow(c, 4) -
           1.56 * pow(t_it, 4) / pow(c, 6) + 0.112 * pow(t_it, 6) / pow(c, 8)) /
          3.25;
    } else {
      out[idx] = 0.0;
    }
    ++idx;
  }

  return out;
}

#include <Rcpp.h>
using namespace Rcpp;

//'rho_aux function
//'
//'Auxiliary function to compute rho
//'@param t numeric vector.
//'@param c tunning constant.
// [[Rcpp::export]]
NumericVector rho_aux(NumericVector t, const double c) {
  return (1.792 - 0.972 * pow(t, 2) / pow(c, 2) +
          0.432 * pow(t, 4) / pow(c, 4) - 0.052 * pow(t, 6) / pow(c, 6) +
          0.002 * pow(t, 8) / pow(c, 8)) /
         3.25;
}

//'rho_opt function
//'
//' An implementation of quasi optimal rho functions following reference [1]
//'
//'@param t numeric vector.
//'@param c tunning constant.
//'@return rho(\code{t}/\code{c}\code{c})
//'@examples val <- rho_opt(t = 0.5, c = 1.0)
//'@references [1] Salibian-Barrera, M., Willems, G., & Zamar,
//'   R. (2008). The fast-tau estimator for regression. Journal of
//'   Computational and Graphical Statistics, 17(3), 659-682.
//'
//'@export
// [[Rcpp::export]]
NumericVector rho_opt(NumericVector t, const double c) {
  NumericVector out(t.size());

  out = 0.5 * pow(t, 2) / (3.25 * pow(c, 2));
  out[abs(t) >= 2 * c] = rho_aux(t[abs(t) >= 2 * c], c);
  out[abs(t) > 3 * c] = 1.0;

  return out;
}
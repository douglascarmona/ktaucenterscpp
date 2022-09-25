#include "distance.h"
#include <Rcpp.h>
using namespace Rcpp;

// TODO: Add docs
// distance function
//'
//'@export
// [[Rcpp::export]]
List distance_to_centers(NumericMatrix x, NumericMatrix centers) {
  /* x is a matrix with observations(rows) and variables,
   centers is a matrix with cluster centers coordinates (rows)
   */

  const std::size_t k = centers.rows();
  const std::size_t n = x.rows();
  const std::size_t p = x.cols();
  IntegerVector membership(
      n); // membership vector TODO: Change initialization to Ralloc
  NumericVector min_distance(n); // distance vector to closest center TODO:
                                 // Change initialization to Ralloc

  for (std::size_t n_iter = 0; n_iter < n; ++n_iter) {
    double min_dist_aux = R_PosInf;
    for (std::size_t k_iter = 0; k_iter < k; ++k_iter) {
      double dist = 0.0;
      for (std::size_t p_iter = 0; p_iter < p; ++p_iter) {
        dist += pow(x(n_iter, p_iter) - centers(k_iter, p_iter), 2);
      }
      dist = sqrt(dist);
      if (dist < min_dist_aux) {
        min_distance(n_iter) = dist;
        membership(n_iter) = k_iter + 1;
        min_dist_aux = dist;
      }
    }
  }
  return (List::create(_["membership"] = membership,
                       _["min_distance"] = min_distance));
}

#include "distance.h"
#include <Rcpp.h>
using namespace Rcpp;

// TODO: Add docs
// distance function
// [[Rcpp::export]]
List distance_to_centers(NumericMatrix x, NumericMatrix centers) {
  /* x is a matrix with observations(rows) and variables,
   centers is a matrix with cluster's centers coordinates (rows)
   */

  const std::size_t k = centers.rows();
  const std::size_t n = x.rows();
  const std::size_t p = x.cols();
  IntegerVector clusters(no_init(n));
  NumericVector distances_min(no_init(n));

  for (std::size_t n_iter = 0; n_iter < n; ++n_iter) {
    double min_dist_aux = R_PosInf;
    for (std::size_t k_iter = 0; k_iter < k; ++k_iter) {
      double dist = 0.0;
      for (std::size_t p_iter = 0; p_iter < p; ++p_iter) {
        dist += pow(x(n_iter, p_iter) - centers(k_iter, p_iter), 2);
      }
      dist = sqrt(dist);
      if (dist < min_dist_aux) {
        distances_min(n_iter) = dist;
        clusters(n_iter) = k_iter + 1;
        min_dist_aux = dist;
      }
    }
  }
  return (List::create(_["clusters"] = clusters,
                       _["distances_min"] = distances_min));
}

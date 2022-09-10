#include "distance.h"
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

// TODO: Add docs
// distance function
//'
//'@export
// [[Rcpp::export]]
List distance_to_centers(NumericMatrix x, NumericMatrix centers) {
  /* x is a matrix with observations(rows) and variables,
   centers is a matrix with cluster centers coordinates (rows)
   */

  const int k = centers.rows(); // number of centers
  const int n =
      x.rows(); // number of observations, TODO: Change type to R_xlen_t
  const int p = x.cols(); // number of variables, TODO: Change type to R_xlen_t
  IntegerVector membership(
      n); // membership vector TODO: Change initialization to Ralloc
  NumericVector min_distance(n); // distance vector to closest center TODO:
                                 // Change initialization to Ralloc
  NumericMatrix distance_matrix(
      n, k); // return value TODO: Change initialization to Ralloc

  for (int n_iter = 0; n_iter < n; ++n_iter) {
    double min_dist_aux = R_PosInf;
    for (int k_iter = 0; k_iter < k; ++k_iter) {
      double dist = 0.0;
      for (int p_iter = 0; p_iter < p; ++p_iter) {
        dist += pow(x(n_iter, p_iter) - centers(k_iter, p_iter), 2);
      }
      dist = sqrt(dist);
      distance_matrix(n_iter, k_iter) = dist;
      if (dist < min_dist_aux) {
        min_distance(n_iter) = dist;
        membership(n_iter) = k_iter + 1;
        min_dist_aux = dist;
      }
    }
  }
  return (List::create(_["distance_matrix"] = distance_matrix,
                       _["membership"] = membership,
                       _["min_distance"] = min_distance));
}

#include "distance.h"
#include "scale.h"
#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// TODO: Add docs
// ktaucenters_run function
//'
//'@export
// [[Rcpp::export]]
List ktaucenters_run(NumericMatrix x, NumericMatrix centers, double tolerance,
                     int max_iter) {
  int n_clusters = centers.rows();

  // TODO: Change type to R_xlen_t
  int n = x.rows();
  int p = x.cols();
  double c1 = const_c1();
  double c2 = const_c2(p);
  // TODO: Replace and use b1 and b2 as function parameters
  double b1 = 0.5;
  double b2 = 1.0;

  int iter = 0;
  double tol = tolerance + 1.0;
  // TODO: Replace initializacion to avoid copying and insert values by
  // iteration instead of 0s
  NumericVector tau_path(max_iter, 0.0);
  NumericVector weights(n);
  IntegerVector clusters(n);
  while (iter < max_iter && tol > tolerance) {
    // Step 1: (re)compute labels
    List dists = distance_to_centers(x, centers);

    // TODO: Think a better way to replace List in c++ internal functions
    // only use Rcpp classes for exported functions
    NumericVector distances_min = dists["min_distance"];
    clusters = dists["membership"];
    double s = mscale(distances_min, c1, b1);
    tau_path[iter] = tau_scale(distances_min, c2, s);

    // Step 2: (re)compute centers
    NumericMatrix old_centers = centers;
    NumericVector Wni = wni(distances_min, c1, c2, s);
    weights = get_weights(Wni, clusters);

    centers = get_new_centers(x, weights, clusters, n_clusters, distances_min);

    tol = max_tolerance(old_centers, centers);
    iter += 1;
  }

  return (List::create(_["tau_path"] = tau_path, _["iter"] = iter,
                       _["centers"] = centers, _["clusters"] = clusters,
                       _["tol"] = tol, _["weights"] = weights));
}
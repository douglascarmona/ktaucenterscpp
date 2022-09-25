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
List ktaucenters_run(NumericMatrix x, NumericMatrix centers,
                     const double tolerance, const int iter_max,
                     const std::string method) {

  const int n_clusters = centers.rows();
  const int n = x.rows();
  const int p = x.cols();
  const double c1 = const_c1();
  const double c2 = const_c2(p);
  // TODO: Replace and use b1 and b2 as function parameters
  const double b1 = 0.5;
  const double b2 = 1.0;

  int iter = 0;
  double tol = tolerance + 1.0;
  NumericVector weights(n);
  NumericVector distances_min(n);
  IntegerVector clusters(n);
  double tau;
  while (iter < iter_max && tol > tolerance) {
    // Step 1: (re)compute labels
    List dists = distance_to_centers(x, centers);

    // TODO: Think a better way to replace List in c++ internal functions
    // only use Rcpp classes for exported functions
    distances_min = dists["min_distance"];
    clusters = dists["membership"];
    double s = mscale(distances_min, c1, b1);
    tau = tau_scale(distances_min, c2, s);
    // Step 2: (re)compute centers
    NumericMatrix old_centers = centers;
    NumericVector Wni = wni(distances_min, c1, c2, s);
    weights = get_weights(Wni, clusters);
    centers = get_new_centers(x, weights, clusters, n_clusters, distances_min);

    tol = max_tolerance(old_centers, centers);
    iter += 1;
  }

  return (List::create(_["tau"] = tau, _["iter"] = iter,
                       _["di"] = distances_min, _["centers"] = centers,
                       _["clusters"] = clusters, _["p"] = p,
                       _["weights"] = weights));
}
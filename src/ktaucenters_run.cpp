#include "cluster.h"
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
                     const double tolerance, const unsigned int iter_max,
                     const std::string method) {

  const std::size_t n_clusters = centers.rows();
  const std::size_t n = x.rows();
  const std::size_t p = x.cols();
  const double c1 = 1.0;
  const double c2 = const_c2(p);
  // TODO: Replace and use b1 and b2 as function parameters
  const double b1 = 0.5;
  const double b2 = 1.0;

  int iter = 0;
  double tol = tolerance + 1.0;
  NumericVector weights(n);
  NumericVector distance_min(n);
  IntegerVector clusters(n);
  double tau;
  while (iter < iter_max && tol > tolerance) {
    // Step 1: (re)compute labels
    List cluster_loc = cluster_location(x, centers);
    distance_min = cluster_loc["distance"];
    clusters = cluster_loc["clusters"];

    double s = mscale(distance_min, c1, b1);
    tau = tau_scale(distance_min, c2, s);

    // Step 2: (re)compute centers
    NumericMatrix old_centers = centers;
    NumericVector Wni = wni(distance_min, c1, c2, s);
    weights = get_weights(Wni, clusters);
    centers = get_new_centers(x, weights, clusters, distance_min);

    tol = max_tolerance(old_centers, centers);
    iter += 1;
  }

  return (List::create(_["tau"] = tau, _["iter"] = iter, _["di"] = distance_min,
                       _["centers"] = centers, _["clusters"] = clusters,
                       _["p"] = p, _["weights"] = weights));
}
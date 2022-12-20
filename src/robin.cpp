#include "knn.h"
#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

//' Estimates the local points density.
//'
//' @param D a distance matrix, which contains the distances between the rows of
//' a matrix.
//' @param k number of neighbors to calculate local point density.
//'
//' @return
//' A vector containing the density values for each point.
//'
// [[Rcpp::export]]
NumericVector point_density(NumericMatrix D, const std::size_t k) {
  const std::size_t n = D.nrow();

  List knn = dist_to_kNN(D, k);
  IntegerMatrix id = knn["id"];
  NumericMatrix distances = knn["dist"];

  NumericVector out(no_init(n));

  for (std::size_t i = 0; i < n; ++i) {
    NumericVector max_distance(k);
    for (std::size_t j = 0; j < k; ++j) {
      max_distance[j] = std::max(distances(id(i, j), k - 1), distances(i, j));
    }
    out[i] = k / sum(max_distance);
  }
  return out;
}

//' Utility function to estimate robinden center
//'
//' @param idp a vector with containing the inverse density each point.
//' @param indexes vector with sorted indexes.
//' @param crit_robin critical robin value.
//'
//' @return
//' Index of the cluster center
//'
// [[Rcpp::export]]
std::size_t robin_center(NumericVector idp, IntegerVector indexes,
                         const double crit_robin) {

  NumericVector idp_sorted_points = idp[indexes];

  // Sometimes all idp_sorted_points are greater than the
  // crit_robin value, then we take the nearest point to crit_robin
  LogicalVector comp = idp_sorted_points <= crit_robin;

  if (is_true(any(comp))) {
    IntegerVector tmp = indexes[comp];
    return tmp[0];
  } else {
    return indexes[which_min(idp_sorted_points - crit_robin)];
  }
}

//' Robust Initialization based on Inverse Density estimator (ROBINDEN)
//'
//' Searches for k initial cluster seeds for k-means based clustering methods.
//'
//' @param D a distance matrix, which contains the distances between the rows of
//' a matrix.
//' @param k number of cluster centers to find.
//' @param mp number of nearest neighbors to find dense regions by LOF
//'
//' @return A list with the following components:
//' \item{centers }{A numeric vector with initial cluster centers indexes
//' centers}
//' \item{idpoints }{A real vector containing the inverse of point density
//' estimation}
//'
//' @details
//' The centers are the observations located in the most dense region
//' and far away from each other at the same time.
//' In order to find the observations in the highly dense region, this function
//' uses point density estimation (instead of Local Outlier Factor, Breunig et
//' al (2000)), see more details.
//'
//' @note This is a slightly modified version of ROBIN algorithm
//' implementation done by Sarka Brodinova <sarka.brodinova@tuwien.ac.at>.
//' @author Juan Domingo Gonzalez <juanrst@hotmail.com>
//'
//' @references Hasan AM, et al. Robust partitional clustering by
//' outlier and density insensitive seeding. Pattern Recognition Letters,
//' 30(11), 994-1002, 2009.
//'
//'@export
// [[Rcpp::export]]
List robinden(NumericMatrix D, const std::size_t n_clusters,
              const std::size_t mp) {

  const std::size_t n = D.nrow();

  // Compute the inverse density points.
  NumericVector idp = 1.0 / point_density(D, mp);

  // Outliers have a high idp value. In unbalanced cases and when the number of
  // clusters increases, all the observations from a group might be above the
  // crit_robin. So we need to increase the crit_robin in order to
  // avoid two initials centers from the same group.

  // Minus 1 to get 0 based index
  const std::size_t position =
      trunc(std::max(0.5, 0.96 * (1 - (1.5 / n_clusters))) * n) - 1;
  NumericVector sorted_idp = clone(idp).sort(false);
  const double crit_robin = sorted_idp[position];

  // Start with a point with maximum density
  std::size_t r = which_min(idp);
  IntegerVector sorted_idxs = top_index(D.column(r), n, true);

  IntegerVector centers(n_clusters);
  centers[0] = robin_center(idp, sorted_idxs, crit_robin);

  for (std ::size_t iter = 1; iter < n_clusters; ++iter) {
    NumericVector minimum_values(no_init(n));

    for (std::size_t column = 0; column < n; ++column) {
      NumericVector c = D.column(column);
      NumericVector tmp = c[centers[Range(0, iter - 1)]];
      minimum_values[column] = min(tmp);
    }

    sorted_idxs = top_index(minimum_values, n, true);
    centers[iter] = robin_center(idp, sorted_idxs, crit_robin);
  }
  return List::create(_["centers"] = centers, _["idpoints"] = idp);
}

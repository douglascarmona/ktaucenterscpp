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
  const std::size_t size = idp.size();
  NumericVector idp_sort_points = idp[indexes];
  bool flag = false;
  std::size_t id;
  NumericVector diff(no_init(size));

  for (std::size_t i = 0; i < size; ++i) {
    diff[i] = idp_sort_points[i] - crit_robin;

    if (idp_sort_points[i] <= crit_robin) {
      id = i;
      flag = true;
      break;
    }
  }

  if (flag == false) {
    // Sometimes all idp_sort_points are greater than the
    // crit_robin value, then we take the nearest point to crit_robin
    id = which_min(diff);
  }

  return indexes[id];
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
List robinden(NumericMatrix D, const std::size_t k, const std::size_t mp) {

  const std::size_t n = D.nrow();

  NumericVector idp = 1 / point_density(D, mp);

  // Outliers have a high idp value. In unbalanced cases and when k increases,
  // all the observations from a group might be above the crit_robin.
  // So we need to increase the crit_robin in order to avoid two initials
  // centers from the same group.
  NumericVector cloned_idp = clone(idp);
  const int nth_element = trunc(std::max(0.5, 0.96 * (1 - (1.5 / k))) * n);
  std::nth_element(cloned_idp.begin(), cloned_idp.begin() + nth_element - 1,
                   cloned_idp.end());
  const double crit_robin = cloned_idp[nth_element - 1];

  // Start with a point with maximum density
  std::size_t min_idx = which_min(idp);
  IntegerVector sorted_idxs = top_index(D.column(min_idx), n, true);

  IntegerVector centers(no_init(k));
  centers[0] = robin_center(idp, sorted_idxs, crit_robin);

  for (std ::size_t m = 1; m < k; ++m) {
    NumericVector minimum_values(no_init(n));

    for (std::size_t column = 0; column < n; ++column) {
      NumericVector c = D.column(column);
      NumericVector tmp = c[centers[Range(0, m - 1)]];
      minimum_values[column] = min(tmp);
    }

    sorted_idxs = top_index(minimum_values, n, true);
    centers[m] = robin_center(idp, sorted_idxs, crit_robin);
  }
  return List::create(_["centers"] = centers, _["idpoints"] = idp);
}

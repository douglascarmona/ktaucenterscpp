#include <Rcpp.h>
#include <algorithm>
#include <queue>
using namespace Rcpp;

// // [[Rcpp::export]]
// NumericMatrix row_sum_cpp(NumericMatrix X, IntegerVector G) {

//   IntegerVector gr = sort_unique(G);
//   int gr_n = gr.size();
//   int nrow = X.nrow(), ncol = X.ncol();

//   // This constructor zero-initializes memory (kind of like
//   // making a copy). You should use:
//   //
//   //     NumericMatrix out = no_init(gr_n, ncol)
//   //
//   // to ensure the memory is allocated, but not zeroed.
//   //
//   // EDIT: We don't have no_init for matrices right now, but you can hack
//   // around that with:
//   //
//   //     NumericMatrix out(Rf_allocMatrix(REALSXP, gr_n, ncol));
//   NumericMatrix out(Rf_allocMatrix(REALSXP, gr_n, ncol));

//   for (int g = 0; g < gr_n; g++) {

//     // subsetting with operator[] is cheaper, so use gr[g] when
//     // you can be sure bounds checks are not necessary
//     int g_id = gr[g];

//     for (int j = 0; j < ncol; j++) {
//       double total = 0;
//       for (int i = 0; i < nrow; i++) {

//         // similarily here
//         if (G[i] != gr[g])
//           continue; // not sure how else to do this
//         total += X(i, j);
//       }
//       // IIUC, you are filling the matrice row-wise. This is slower as
//       // R matrices are stored in column-major format, and so filling
//       // matrices column-wise will be faster.
//       out[g, j] = total;
//     }
//   }
//   return out;
// }

// TODO: Add docs
// tabulateC function
//'
//'@export
// [[Rcpp::export]]
IntegerVector tabulatecpp(IntegerVector x, const std::size_t max) {
  IntegerVector counts(max);
  for (const auto &now : x) {
    if (now > 0 && now <= max)
      counts[now - 1] += 1;
  }
  return counts;
}

// TODO: Add cdocs
// distance function
//'
//'@export
// [[Rcpp::export]]
double median_cpp(NumericVector x) {
  // TODO: Move function to utils. Replace function with partial_sort
  std::size_t size = x.size();
  std::sort(x.begin(), x.end());
  if (size % 2 == 0)
    return (x[size / 2 - 1] + x[size / 2]) / 2.0;
  return x[size / 2];
}

class IndexComparator {
public:
  IndexComparator(const NumericVector &data_) : data(data_.begin()) {}

  inline bool operator()(int i, int j) const {
    return data[i] > data[j] || (data[i] == data[j] && j > i);
  }

private:
  const Vector<REALSXP>::const_iterator data;
};

class IndexQueue {
public:
  typedef std::priority_queue<int, std::vector<int>, IndexComparator> Queue;

  IndexQueue(const NumericVector &data_)
      : comparator(data_), q(comparator), data(data_) {}

  inline operator IntegerVector() {
    int n = q.size();
    IntegerVector res(n);
    for (int i = n - 1; i >= 0; --i) {
      // Add +1 for 1-based R indexing
      res[i] = q.top();
      q.pop();
    }
    return res;
  }
  inline void input(int i) {
    // if( data[ q.top() ] < data[i] ) {
    if (comparator(i, q.top())) {
      q.pop();
      q.push(i);
    }
  }
  inline void pop() { q.pop(); }
  inline void push(int i) { q.push(i); }

private:
  IndexComparator comparator;
  Queue q;
  const Vector<REALSXP> &data;
};

// TODO: Add docs
// top_index function
//'
//'@export
// [[Rcpp::export]]
IntegerVector top_index(NumericVector v, int n) {
  int size = v.size();

  // not interesting case. Less data than n
  if (size < n) {
    return seq(0, n - 1);
  }

  IndexQueue q(v);
  for (int i = 0; i < n; ++i)
    q.push(i);
  for (int i = n; i < size; ++i)
    q.input(i);
  return q;
}

// TODO: Add cdocs
// tolerance function
//'
//'@export
// [[Rcpp::export]]
double max_tolerance(NumericMatrix old_centers, NumericMatrix centers) {

  const int p = centers.cols();
  const int k = centers.rows();
  NumericMatrix diff(k, p);

  std::transform(old_centers.begin(), old_centers.end(), centers.begin(),
                 diff.begin(), [](double old_center, double new_center) {
                   return pow(old_center - new_center, 2);
                 });

  NumericVector tolerances(k);

  for (int column = 0; column != p; ++column) {
    for (int row = 0; row != k; ++row) {
      tolerances[row] += diff(row, column);
    }
  }

  std::transform(tolerances.begin(), tolerances.end(), tolerances.begin(),
                 ::sqrt);

  return max(tolerances);
}

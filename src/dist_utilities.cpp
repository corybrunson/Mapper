// dist_utilities.cpp
#include "dist_utilities.h"

// Given a 'dist' object (vector) and an integer vector of indices (1-based), create a new dist object
// represented by the subset of 'dist' indexed by 'idx'
// [[Rcpp::export]]
NumericVector dist_subset(const NumericVector& dist, IntegerVector idx){
  const int n = dist.attr("Size");
  const int cl_n = idx.length();
  if (n == cl_n) { return dist; } // if the full index set was specified. Assumes idx is unique.
  NumericVector new_dist = NumericVector((cl_n * (cl_n - 1))/2);
  int ii = 0;
  for (IntegerVector::iterator i = idx.begin(); i != idx.end(); ++i){
    for (IntegerVector::iterator j = i; j != idx.end(); ++j){
      if (*i == *j) { continue; }
      const int ij_idx = INDEX_TF(n, (*i < *j ? *i : *j) - 1, (*i < *j ? *j : *i) - 1);
      new_dist[ii++] = dist[ij_idx];
    }
  }
  new_dist.attr("Size") = cl_n;
  new_dist.attr("class") = "dist";
  return(new_dist);
}

// Given a 'dist' object (vector) and two (n, m)-length integer vectors of indices (1-based), 
// create a (n x m) matrix of the pairwise distance between the indices of points represented 
// in the 'dist' object. Does not check dimensions or indices. 
// [[Rcpp::export]]
NumericMatrix dist_subset_pw(const NumericVector& dist, IntegerVector idx1, IntegerVector idx2){
  const int dist_n = dist.attr("Size");
  const size_t n = idx1.size();
  const size_t m = idx2.size();
  NumericMatrix row_dist = NumericMatrix(n, m);
  size_t ii, jj; 
  for (size_t i = 0; i < n; ++i){
    for (size_t j = 0; j < m; ++j){
      ii = idx1[i] - 1;
      jj = idx2[j] - 1;
      const size_t ij_idx = INDEX_TF(dist_n, (ii < jj ? ii : jj), (ii < jj ? jj : ii));
      row_dist(i, j) = dist[ij_idx];
    }
  }
  CharacterVector rn = CharacterVector(idx1.begin(), idx1.end());
  CharacterVector cn = CharacterVector(idx2.begin(), idx2.end());
  rownames(row_dist) = rn; 
  colnames(row_dist) = cn; 
  return(row_dist);
}


NumericVector combine(const NumericVector& t1, const NumericVector& t2){
  std::size_t n = t1.size() + t2.size();
  NumericVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}

// Given a distance vector 'dist' and a new set of pairwise distances 'new_dists' from 'n_new'
// points, concatenate two dist objects together.
// [[Rcpp::export]]
NumericVector concatDist(const NumericVector& dist, const int n_new, const NumericVector& new_dists){
  const int n = dist.attr("Size");
  NumericVector::const_iterator new_dist = new_dists.begin();
  NumericVector::const_iterator old_dist = dist.begin();
  const std::size_t new_size = dist.size() + new_dists.size();
  const int n_pts = n + n_new;
  NumericVector output = Rcpp::no_init(new_size);
  for (int i = 0; i < (int) new_size; ++i){
    const int to = INDEX_TO(i, n_pts);
    const int from = INDEX_FROM(i, n_pts, to);
    // Rcout << "i: " << i << ", from: " << from << ", to: " << to << std::endl;
    if (to >= n || from >= n){ output[i] = *new_dist; new_dist++; }
    else { output[i] = *old_dist; old_dist++; }
  }
  output.attr("Size") = n_pts;
  output.attr("class") = "dist";
  return(output);
}

// For each point in Q, computes the pairwise euclidean distance to every point in R.
// Assumes the two are disjoint.
// [[Rcpp::export]]
NumericVector dist_from_to(const NumericMatrix& X_query, const NumericMatrix& X_ref){
  const int q_n = X_query.nrow();
  const int r_n = X_ref.nrow();
  NumericVector output = Rcpp::no_init(q_n*r_n + q_n*(q_n - 1)/2);
  // if (add_self){ output = ; }
  // else { output = Rcpp::no_init(q_n*r_n); }

  // Compute all pairwise distances
  int counter = 0;
  for (int i = 0; i < r_n; ++i){
    for (int j = 0; j < q_n; ++j){
      NumericMatrix::ConstRow query_pt = X_query(j, _);
      NumericMatrix::ConstRow ref_pt = X_ref(i, _);
      output[counter++] = euc_dist(query_pt, ref_pt);
    }
  }

  // If add_self is true, then continuously add distances between query points as well
  for (int i = 0; i < q_n; ++i){
    for (int j = i+1; j < q_n; ++j){
    NumericMatrix::ConstRow query_pt1 = X_query(i, _);
    NumericMatrix::ConstRow query_pt2 = X_query(j, _);
    output[counter++] = euc_dist(query_pt1, query_pt2);
    }
  }
  return(output);
}


/*** R
# Rcpp::sourceCpp('src/dist_utilities.cpp', embeddedR = FALSE)

## Testing
# iris_mat <- as.matrix(iris[, 1:4])
# x <- dist(iris_mat)
# ref_x <- as.matrix(iris_mat[1:49, 1:4])
# query_x <- as.matrix(iris_mat[50:150, 1:4])
# new_dists <- dist_from_to(X_query = query_x, X_ref = ref_x, add_self = 1)
# ref_dist <- dist(ref_x)
# test_dist <- concatDist(dist = ref_dist, n_new = nrow(query_x), new_dists = new_dists)
# all(x == test_dist)

## Benchmarks
# x <- cbind(rnorm(10000), rnorm(10000), rnorm(10000), rnorm(10000))
# distx <- dist(x)
# ref_x <- x[1:9900,]
# query_x <- x[9901:10000,]
# new_dists <- dist_from_to(X_query = query_x, X_ref = ref_x)
# ref_dist <- dist(ref_x)
# microbenchmark::microbenchmark(dist(x), times = 5L)
# # Unit: milliseconds
# # expr      min       lq     mean   median       uq      max neval
# #   dist(x) 573.0578 582.6419 602.2179 597.9135 610.0308 647.4452     5
# microbenchmark::microbenchmark(concatDist(dist = ref_dist, n_new = nrow(query_x), new_dists = dist_from_to(X_query = query_x, X_ref = ref_x)), times = 5L)
# # Unit: milliseconds
# # expr      min       lq     mean   median
# # concatDist(dist = ref_dist, n_new = nrow(query_x), new_dists = dist_from_to(X_query = query_x,      X_ref = ref_x)) 939.5041 943.7769 963.5275 958.7722
# # uq      max neval
# # 978.9124 996.6717     5
# microbenchmark::microbenchmark(dist_subset(distx, idx = 9901:10000), times = 5L)
# # Unit: microseconds
# # expr    min    lq     mean median     uq     max neval
# #   dist_subset(distx, idx = 9901:10000) 28.243 28.25 218.3678 28.863 58.526 947.957     5
*/

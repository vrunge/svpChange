#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Smallest Valid Partitioning with a User-Defined Validity Test
//'
//' @title Smallest Valid Partitioning with a User-Defined Validity Test
//' @description This function uses dynamic programming to find a partition of
//' a univariate signal whose segments all pass a user-defined validity test.
//' It first minimizes the number of segments and, among partitions with the
//' same number of segments, minimizes the within-segment sum of squares.
//'
//' @param data Numeric vector containing the univariate signal to segment.
//' @param gamma Numeric threshold passed to `test`.
//' @param test Function of the form `function(segment, gamma)` returning one
//'   non-missing logical value: `TRUE` if the segment is valid. The function
//'   is not called for singleton segments; they are always valid.
//' @param subtests Character scalar controlling validity-based pruning. The
//'   choices are `"both"` (default), `"right"`, and `"none"`.
//'   `"right"` removes a candidate start when its current segment is invalid;
//'   `"both"` additionally removes that start and all older starts; and
//'   `"none"` applies neither rule. Invalid candidates are never used for the
//'   current optimum.
//' @param PELT_pruning Logical; whether to apply the additional cost-based
//'   candidate pruning rule.
//'
//' @details
//' A candidate boundary `s` at endpoint `t` represents the R segment
//' `data[(s + 1):t]`; `s` is zero-based and `t` is one-based. The quadratic
//' cost is the residual sum of squares around the segment mean. Singleton
//' segments are always valid, regardless of the result of `test`.
//'
//' @examples
//' range_test <- function(segment, gamma) {
//'   diff(range(segment)) <= gamma
//' }
//' set.seed(1)
//' data <- ts_generator(
//'   chpts = c(30, 60), parameters = c(0, 4),
//'   sd_noise = 0.25, type = "gauss"
//' )
//' fit <- svp0(data, gamma = 2, test = range_test,
//'            subtests = "both")
//' fit$changepoints
//'
//' @return A list with the following components:
//' \describe{
//'   \item{changepoints}{Numeric vector of increasing segment-ending
//'     positions, including `length(data)`. Internal values are estimated
//'     change-point positions.}
//'   \item{lastIndexSet}{Numeric vector of candidate boundaries remaining at
//'     termination. These are zero-based boundaries and are returned in
//'     decreasing order.}
//'   \item{nb}{Numeric vector in time order. Element `t` is the number of
//'     candidate entries examined at endpoint `t`, before pruning.}
//'   \item{costQ}{Always `NULL`; cumulative costs are stored in the first
//'     column of `R`.}
//'   \item{R}{Numeric matrix with `length(data)` rows and three columns. Row
//'     `t` describes the optimum ending at observation `t`:
//'     \describe{
//'       \item{Q}{cumulative cost}
//'       \item{K}{number of segments in Q}
//'       \item{s}{zero-based previous boundary; the final segment is
//'         `data[(s + 1):t]`}
//'     }
//'   }
//' }
//'
//' @export
// [[Rcpp::export]]
List svp0(std::vector<double> data,
          double gamma,
          Function test,
          std::string subtests = "both",
          bool PELT_pruning = false)
{
  size_t n = data.size();

  //
  // STOPS
  //
  if (n == 0)
  {
    stop("'data' must contain at least one observation.");
  }
  if (subtests != "both" && subtests != "right" && subtests != "none")
  {
    stop("Invalid value for 'subtests'. Must be 'both', 'right', or 'none'.");
  }

  // Initialization of elements
  NumericMatrix R(n + 1, 3); // Q, K, s
  R(0, 0) = 0.0; // Default value
  R(0, 1) = 0.0; // Default value
  R(0, 2) = 0.0; // Default value
  // nb of candidates examined at each t
  std::vector<size_t> nb(n);

  //
  // PREPROCESSING
  //
  // Cumulative sum for optimized calculations
  std::vector<double> S1(n + 1, 0);
  std::vector<double> S2(n + 1, 0);
  for (size_t i = 0; i < n; i++)
  {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }

  bool valid; // for validity test

  double best_Q;
  size_t best_K;
  size_t best_s = 0;
  size_t s;
  std::vector<double> seg; // for the validity
  double candidate_Q; // for the lex. comparison
  size_t candidate_K; // for the lex. comparison

  std::vector<size_t> INDEX = {0};
  std::vector<size_t> valid_INDEX; //indices that pass the validity test
  std::vector<size_t> non_pruned_INDEX; // candidates retained by pruning

  //
  // MAIN LOOP
  //
  for (size_t t = 1; t < n + 1; t++)
  {
    nb[t - 1] = INDEX.size(); // number of candidates examined at each previous time
    // initialization
    best_Q = std::numeric_limits<double>::infinity(); // Inf
    best_K = std::numeric_limits<size_t>::max(); // Inf
    best_s = 0;

    valid_INDEX.clear(); // Reuse storage for valid candidates.
    size_t max_non_valid_index = 0; // the largest index of a candidate that fails the validity test

    for (size_t k = 0; k < INDEX.size(); k++)
    {
      s = INDEX[k];

      // Singleton segments are always valid. Other segments are tested by R.
      if (t - s == 1)
      {
        valid = true; // singleton segments are always valid
      }
      else
      {
        seg.assign(data.begin() + s, data.begin() + t);
        valid = as<bool>(test(seg, gamma)); // test against the threshold gamma
      }

      // IF VALID, do the comparisons and save s in valid_INDEX
      if (valid == true)
      {
        valid_INDEX.push_back(s);

        // compute for s, the overall cost and number of segments
        candidate_Q = R(s, 0) + (S2[t] - S2[s]) -
          (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        candidate_K = R(s, 1) + 1;

        // lexicographic order, if the new candidate is better
        // update the best_Q, best_K and best_s
        if (candidate_K < best_K ||
          (candidate_K == best_K && candidate_Q < best_Q))
        {
          best_Q = candidate_Q;
          best_K = candidate_K;
          best_s = s;
        }
      }
      else if (s > max_non_valid_index)
      {
        max_non_valid_index = s;
      }
    }

    // write the best answer in R
    R(t, 0) = best_Q;
    R(t, 1) = best_K;
    R(t, 2) = best_s;
    // Right subtests retain only candidates whose segment from s to t is valid.
    if (subtests == "both" || subtests == "right")
    {
      INDEX.swap(valid_INDEX);
    }

    // In "both", the inclusive left rule removes candidates older than the
    // largest invalid candidate (the invalid candidate was removed above).
    if (subtests == "both")
    {
      non_pruned_INDEX.clear();
      for (size_t k = 0; k < INDEX.size(); k++)
      {
        s = INDEX[k];
        if (s >= max_non_valid_index)
        {
          non_pruned_INDEX.push_back(s);
        }
      }
      INDEX.swap(non_pruned_INDEX);
    }


    //
    //  PRUNING PELT
    //
    if (PELT_pruning == true)
    {
      non_pruned_INDEX.clear(); // Reuse storage for retained candidates.
      for (size_t k = 0; k < INDEX.size(); k++)
      {
        s = INDEX[k];
        candidate_Q = R(s, 0) + (S2[t] - S2[s]) -
          (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        candidate_K = R(s, 1);

        if (!((candidate_Q > best_Q) && (candidate_K == best_K)))
        {
          non_pruned_INDEX.push_back(s);
        }

      }
      non_pruned_INDEX.push_back(t);
      INDEX.swap(non_pruned_INDEX); // Keep candidates retained by PELT pruning.
    }
    else
    {
      INDEX.push_back(t);
    }
  }

  //
  // BACKTRACKING
  //
  // Change points reconstruction
  std::vector<size_t> changepoints;
  size_t i = n;
  while (i > 0)
  {
    changepoints.push_back(i);
    i = R(i, 2);
  }
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(INDEX.begin(), INDEX.end());

  return List::create(
    Named("changepoints") = changepoints,
    Named("lastIndexSet") = INDEX,
    Named("nb") = nb,
    Named("costQ") = R_NilValue,
    Named("R") = R(Range(1, R.nrow() - 1), Range(0, R.ncol() - 1))
  );
}

#include "tests.h"
#include <Rcpp.h>
#include <functional>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Smallest Valid Partitioning with FOCUS validity test
//'
//' @title Smallest Valid Partitioning with Validation and Pruning
//' @description This function implements a dynamic programming approach to segment a univariate signal into the smallest number of valid segments, according to a user-defined validation function. Each segment must pass a validity test (e.g., based on variance, range, etc.). The algorithm minimizes a quadratic cost subject to this constraint.
//'
//' @param data A numeric vector representing the univariate signal to be segmented.
//' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
//' @param test A function of the form `function(data, gamma)` returning TRUE if the segment is valid. Default is `valid_OP`.
//' @param prune_if_unvalid Logical. If TRUE (default), the algorithm applies *segment-wise validation*:
//' at each time step, it tests whether the candidate segment \code{data[(s+1):t]} is valid using the
//' user-defined function \code{test}. If the segment fails the test, the candidate \code{s} is removed
//' (pruned) from the set of possible changepoints. This accelerates computation by avoiding invalid
//' segment extensions. If FALSE, the algorithm skips this validation and considers all candidate
//' segments without checking their validity (which can be faster but may return invalid segments).
//'
//' @return A list with the following components :
//' \describe{
//'   \item{changepoints}{Integer vector indicating the ending index of each segment (i.e., positions of changepoints).}
//'  \item{lastIndexSet}{Integer vector containing the non-pruned indices at the end of the algorithm execution}
//'   \item{nb}{Integer vector of length \code{length(data)}. At each position \code{t}, it records the number of candidates tested.}
//'   \item{costQ}{Numeric vector of length \code{length(data)}. Quadratic cost value at each time step. Set to NULL as it is recorded into matrix R}
//'   \item{R}{A matrix of dimension \code{(length(data)+1) x 3} containing, for each time step :
//'     \describe{
//'       \item{Q}{cumulative cost}
//'       \item{K}{number of segments in Q}
//'       \item{s}{previous changepoint}
//'     }
//'   }
//' }
//'
//'
//' @export
// [[Rcpp::export]]
List SVP(std::vector<double> data,
         double gamma,
         std::string test,
         bool prune_if_unvalid = true)
{
  size_t n = data.size();

  // Initialization of elements in the return
  NumericMatrix R(n + 1, 3); // Q, K, s
  R(0, 0) = 0.0;
  R(0, 1) = 0.0;
  R(0, 2) = 0.0;

  std::vector<size_t> nb(n); // nb of candidates examined at each t

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

  //
  // Define a generic initializer for test objects
  //
  std::function<std::unique_ptr<TestBase>()> newTest;
  if (test == "gaussian_mean")
  {
    newTest = []() { return std::make_unique<GaussianMean>();};
  }
  // Add more cases here for other tests, e.g.:
  // else if (test == "bernoulli_mean") {
  //   newTest = []() { return std::make_unique<BernoulliMean>(); };
  // }
  else
  {
    stop("Unknown test type");
  }

  //
  //
  //
  std::vector<size_t> INDEX = {0};
  std::vector<size_t> valid_INDEX;  // indices that pass the validity test

  std::vector<std::unique_ptr<TestBase>> TESTS;
  TESTS.push_back(newTest());

  double best_Q;
  size_t best_K;
  size_t best_s = 0;
  size_t s;
  bool valid;         // for validity test
  double candidate_Q; // for the lex. comparison
  size_t candidate_K; // for the lex. comparison

  for (size_t t = 1; t < n + 1; ++t)
  {
    nb[t - 1] = INDEX.size();

    // Initialization
    best_Q = std::numeric_limits<double>::infinity();
    best_K = std::numeric_limits<size_t>::max();

    ///
    /// the elements to be saved
    ///
    valid_INDEX.clear(); // set to length 0 this vector, fill it with valid indices
    std::vector<std::unique_ptr<TestBase>> valid_TESTS; // to IMPROVE

    for (size_t k = 0; k < INDEX.size(); ++k)
    {
      s = INDEX[k];
      auto& test_instance = TESTS[k];

      test_instance->update(data[t - 1]);
      valid = test_instance->statistic() < gamma;

      if (valid == true)
      {
        ////
        //// save the s and the instance
        ////
        valid_INDEX.push_back(s);
        valid_TESTS.push_back(std::move(test_instance));

        candidate_Q = R(s, 0) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        candidate_K = R(s, 1) + 1;

        if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
        {
          best_Q = candidate_Q;
          best_K = candidate_K;
          best_s = s;
        }
      }
    }
    R(t, 0) = best_Q;
    R(t, 1) = best_K;
    R(t, 2) = best_s;

    //
    //  PRUNING if prune_if_unvalid == true
    //
    if (prune_if_unvalid == true)
    {
      INDEX.swap(valid_INDEX); // index now contains the valid_INDEX only
      TESTS = std::move(valid_TESTS);
    }
    TESTS.push_back(newTest());
    INDEX.push_back(t);
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
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = INDEX,
    _["nb"] = nb,
    _["costQ"] = NULL,
    _["R"] = R(Range(1, R.nrow() - 1), Range(0, R.ncol() - 1))
  );
}
